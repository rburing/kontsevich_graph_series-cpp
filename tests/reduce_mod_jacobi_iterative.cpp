#include "../kontsevich_graph_series.hpp"
#include "../leibniz_graph.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include "../util/continued_fraction.hpp"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
using namespace std;
using namespace GiNaC;

typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double> SparseMatrix;

double threshold = 1e-5;

int main(int argc, char* argv[])
{
    if (argc == 1 || (argc == 2 && string(argv[1]) == "--help"))
    {
        cout << "Usage: " << argv[0] << " <graph-series-filename> [--max-jac-indegree=k] [--skew-leibniz] [--solve]\n\n"
             << "Accepts only homogeneous power series: graphs with n internal vertices at order n.\n\n"
             << "--max-iterations=i     perform at most i iterations.\n"
             << "--max-jac-indegree=k   restricts the number of arrows falling on Jacobiators to be <= k.\n"
             << "--skew-leibniz         skew-symmetrizes each Leibniz graph before subtracting it with an undetermined coefficient.\n"
             << "--leibniz-in=filename  input graph series already contains Leibniz graphs, with encodings in filename.\n"
             << "--leibniz-out=filename store Leibniz graph encodings in filename.\n"
             << "--coeff-prefix=c       let the coefficients of leibniz graphs be c_n.\n"
             << "--solve                the undetermined variables in the input are added to the linear system to-be-solved.\n"
             << "--solve-numerically    try to solve the linear system using Eigen.\n"
             << "--interactive          ask whether to continue to the next iteration, and whether to solve numerically.\n";
        return 1;
    }

    bool interactive = false;
    bool solve = false;
    bool solve_numerically = false;
    size_t max_jac_indegree = numeric_limits<size_t>::max();
    size_t max_iterations = numeric_limits<size_t>::max();
    bool skew_leibniz = false;
    string leibniz_in_filename = "";
    string leibniz_out_filename = "";
    string coefficient_prefix = "c";

    // Process arguments
    for (int idx = 2; idx < argc; ++idx)
    {
        string argument = argv[idx];
        size_t equals_pos = argument.find('=');
        if (equals_pos != string::npos)
        {
            string key = argument.substr(0, equals_pos);
            string value = argument.substr(equals_pos+1);
            if (key == "--max-iterations")
                max_iterations = stoi(value);
            else if (key == "--max-jac-indegree")
                max_jac_indegree = stoi(value);
            else if (key == "--coeff-prefix")
                coefficient_prefix = value;
            else if (key == "--leibniz-in")
                leibniz_in_filename = value;
            else if (key == "--leibniz-out")
                leibniz_out_filename = value;
            else {
                cout << "Unrecognized option: " << argument << "\n";
                return 1;
            }
        }
        else
        {
            if (argument == "--skew-leibniz")
                skew_leibniz = true;
            else if (argument == "--solve")
                solve = true;
            else if (argument == "--solve-numerically")
                solve_numerically = true;
            else if (argument == "--interactive")
                interactive = true;
            else {
                cout << "Unrecognized option: " << argument << "\n";
                return 1;
            }
        }
    }

    cout << "Options: "
         << "max-iterations = ";
    if (max_iterations == numeric_limits<size_t>::max())
        cout << "none";
    else
        cout << max_iterations;
    cout << ", max-jac-indegree = ";
    if (max_jac_indegree == numeric_limits<size_t>::max())
        cout << "none";
    else
        cout << max_jac_indegree;
    cout << ", solve = " << (solve ? "yes" : "no")
         << ", solve-numerically = " << (solve_numerically ? "yes" : "no")
         << ", skew-leibniz = " << (skew_leibniz ? "yes" : "no")
         << ", leibniz-in = " << (leibniz_in_filename == "" ? "none" : leibniz_in_filename)
         << ", coeff-prefix = " << coefficient_prefix
         << ", interactive = " << (interactive ? "yes" : "no") << "\n";

    // Reading in Leibniz graphs
    parser coefficient_reader;
    map< LeibnizGraph<ex>, ex> leibniz_graphs;

    if (leibniz_in_filename != "")
    {
        ifstream leibniz_in_file(leibniz_in_filename);
        leibniz_graphs = LeibnizGraph<ex>::map_from_istream(leibniz_in_file,
                                                                [&coefficient_reader](string s) -> ex
                                                                {
                                                                    return coefficient_reader(s);
                                                                });
    }
    std::vector<ex> leibniz_coeffs;
    for (auto const& pair : leibniz_graphs)
        leibniz_coeffs.push_back(pair.second);

    // Reading in graph series
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    bool homogeneous = true;
    map<size_t, set< vector<size_t> > > in_degrees;
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file,
        [&coefficient_reader](string s) -> ex { return coefficient_reader(s); },
        [&homogeneous, &in_degrees](KontsevichGraph graph, size_t order) -> bool
                                   {
                                       in_degrees[order].insert(graph.in_degrees());
                                       return homogeneous &= graph.internal() == order;
                                   }
    );
    size_t order = graph_series.precision();
    if (!homogeneous)
    {
        cerr << "Only accepting homogeneous power series: graphs with n internal vertices at order n.\n";
        return 1;
    }

    graph_series.reduce_mod_skew();

    vector<symbol> unknowns_list;
    for (auto const& namevar : coefficient_reader.get_syms())
        if (find(leibniz_coeffs.begin(), leibniz_coeffs.end(), namevar.second) == leibniz_coeffs.end())
            unknowns_list.push_back(ex_to<symbol>(namevar.second));

    if (unknowns_list.size() != 0 && !solve)
    {
        cerr << "Input contains unknowns; remove them or specify --solve to solve for them.\n";
        return 1;
    }

    vector<symbol> coefficient_list = unknowns_list;
    coefficient_list.reserve(leibniz_coeffs.size());
    for (auto& leibniz_coeff : leibniz_coeffs)
        coefficient_list.push_back(ex_to<symbol>(leibniz_coeff));

    KontsevichGraphSeries<ex> leibniz_graph_series = graph_series;

    set<KontsevichGraph> processed_graphs;

    size_t counter = 0;
    bool converged = false;
    size_t step = 0;
    while (!converged && ++step <= max_iterations)
    {
        size_t old_counter = counter;
        for (size_t n = 0; n <= order; ++n)
        {
            size_t termcounter = 0;
            for (auto& term : graph_series[n])
            {
                cerr << "\rProcessing term " << ++termcounter << " / " << graph_series[n].size() << ". ";
                if (skew_leibniz)
                    cerr << "Skew-";
                cerr << "Leibniz graphs so far: " << counter;

                KontsevichGraph& graph = term.second;

                // Check if this graph has already been processed
                // TODO: optimization: take minimum of skew symmetrization, if skew_leibniz
                if (find(processed_graphs.begin(), processed_graphs.end(), graph) != processed_graphs.end())
                    continue;

                // Use that the graph series is reduced (graphs are in normal form with sign +1)
                processed_graphs.insert(graph);

                // Subtract Leibniz graphs that yield this Kontsevich graph
                for (LeibnizGraph<ex> leibniz_graph : LeibnizGraph<ex>::those_yielding_kontsevich_graph(graph, skew_leibniz))
                {
                    leibniz_graph.normalize();
                    if (leibniz_graphs.find(leibniz_graph) != leibniz_graphs.end())
                        continue;
                    symbol coefficient(coefficient_prefix + "_" + to_string(counter));
                    KontsevichGraphSum<ex> graph_sum = leibniz_graph.expansion(coefficient);
                    if (graph_sum.size() != 0)
                    {
                        coefficient_list.push_back(coefficient);
                        leibniz_graphs[leibniz_graph] = coefficient;
                        ++counter;
                        leibniz_graph_series[n] -= graph_sum;
                    }
                }
            }
        }

        cout << "\nNumber of Leibniz graphs: " << leibniz_graphs.size() << "\n";
        cout << "\nNumber of terms (before reducing): " << leibniz_graph_series[order].size() << "\n";

        cout.flush();

        cerr << "\nReducing...\n";
        leibniz_graph_series.reduce_mod_skew();

        cout << "\nNumber of terms (after reducing): " << leibniz_graph_series[order].size() << "\n";

        lst equations;

        for (size_t n = 0; n <= order; ++n)
            for (auto& term : leibniz_graph_series[n])
            {
                cout << term.second.encoding() << "    " << term.first << "==0\n";
                equations.append(term.first);
            }

        size_t rows = equations.nops();
        size_t cols = coefficient_list.size();

        cout << "Got linear system of size " << rows << " x " << cols << ".\n";

        char solve_input = 'N';
        if (interactive)
        {
            cerr << "Solve? (Y/N) ";
            cin >> solve_input;
            solve_numerically = solve_input == 'Y';
        }

        if (solve_numerically)
        {
            // Set up sparse matrix linear system

            cerr << "Setting up linear system for numerical solution...\n";
            Eigen::VectorXd b(rows);
            SparseMatrix matrix(rows,cols);

            vector<Triplet> tripletList;
            size_t idx = 0;
            for (ex equation : equations)
            {
                if (!is_a<add>(equation))
                    equation = lst({equation});
                for (ex term : equation)
                {
                    if (!is_a<mul>(term))
                        term = lst({term});
                    double prefactor = 1;
                    symbol coefficient("one");
                    for (ex factor : term)
                    {
                        if (is_a<numeric>(factor))
                            prefactor *= ex_to<numeric>(factor).to_double();
                        else if (is_a<symbol>(factor))
                            coefficient = ex_to<symbol>(factor);
                    }
                    if (coefficient.get_name() == "one") // constant term
                        b(idx) = -prefactor;
                    else
                        tripletList.push_back(Triplet(idx,find(coefficient_list.begin(), coefficient_list.end(), coefficient) - coefficient_list.begin(), prefactor));
                        // NB: Eigen uses zero-based indices (contrast MATLAB, Mathematica)
                }
                ++idx;
            }

            matrix.setFromTriplets(tripletList.begin(), tripletList.end());
            
            cerr << "Solving linear system (" << rows << " x " << cols << ") numerically...\n";

            Eigen::SparseQR< SparseMatrix, Eigen::COLAMDOrdering<int> > qr(matrix);
            Eigen::VectorXd x = qr.solve(b);
            
            cerr << "Residual norm = " << (matrix * x - b).squaredNorm() << "\n";

            cerr << "Rounding...\n";
            x = x.unaryExpr([](double elem) { return fabs(elem) < threshold ? 0.0 : elem; });

            double residual_norm = (matrix * x - b).squaredNorm();
            cerr << "Still a solution? Residual norm = " << residual_norm << "\n";

            if (!isnan(residual_norm) && residual_norm < 1.0)
            {
                cerr << "Approximating numerical solution by rational solution...\n";

                lst zero_substitution;
                lst solution_substitution;
                for (int i = 0; i != x.size(); i++)
                {
                    ex result = best_rational_approximation(x.coeff(i), threshold);
                    if (result == 0)
                        zero_substitution.append(coefficient_list[i] == 0);
                    else
                        solution_substitution.append(coefficient_list[i] == result);
                }

                cerr << "Substituting zeros...\n";

                KontsevichGraphSeries<ex> leibniz_graph_series_copy = leibniz_graph_series; // TODO: inefficient?

                for (auto& order: leibniz_graph_series_copy)
                    for (auto& term : leibniz_graph_series_copy[order.first])
                        term.first = term.first.subs(zero_substitution);

                cerr << "Reducing zeros...\n";

                leibniz_graph_series_copy.reduce_mod_skew();

                for (size_t n = 0; n <= leibniz_graph_series_copy.precision(); ++n)
                {
                    cout << "h^" << n << ":\n";
                    for (auto& term : leibniz_graph_series_copy[n])
                    {
                        cout << term.second.encoding() << "    " << term.first << "\n";
                    }
                }

                cerr << "Verifying solution...\n";

                for (auto& order: leibniz_graph_series_copy)
                    for (auto& term : leibniz_graph_series_copy[order.first])
                        term.first = term.first.subs(solution_substitution);

                leibniz_graph_series_copy.reduce_mod_skew();

                cout << "Do we really have a solution? " << (leibniz_graph_series_copy == 0 ? "Yes" : "No") << "\n";

                if (leibniz_graph_series_copy == 0)
                {
                    for (auto& pair : leibniz_graphs)
                        if (pair.second != ex(pair.second).subs(solution_substitution))
                            cout << pair.first.encoding() << "    " << pair.second << "==" << ex(pair.second).subs(solution_substitution) << "\n";
                    for (auto var : unknowns_list)
                        cout << var << "==" << ex(var).subs(solution_substitution) << "\n";
                }
            }
        }

        char iterate = 'Y';
        if (interactive)
        {
            cerr << "Next iteration? (Y/N) ";
            cin >> iterate;
        }
        if (iterate != 'Y')
            break;

        graph_series = leibniz_graph_series;

        converged = counter == old_counter;
    }
    if (converged)
        cout << "\nConverged in " << step << " steps.\n\n";
    if (skew_leibniz)
        cout << "Skew-";
    cout << "Leibniz graphs:";
    ostream* leibniz_out_stream = &cout;
    ofstream leibniz_out_fstream;
    if (leibniz_out_filename != "")
    {
        cout << " writing to " << leibniz_out_filename << "\n";
        leibniz_out_fstream.open(leibniz_out_filename);
        leibniz_out_stream = &leibniz_out_fstream;
    }
    else
        cout << "\n";
    for (auto& pair : leibniz_graphs)
    {
        (*leibniz_out_stream) << pair.first.encoding() << "    " << pair.second << "\n";
    }
    leibniz_out_fstream.close();
}
