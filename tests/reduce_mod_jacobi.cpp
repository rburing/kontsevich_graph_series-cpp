#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
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

double threshold = 1e-13;

int main(int argc, char* argv[])
{
    if (argc != 2 && argc != 3)
    {
        cout << "Usage: " << argv[0] << " <graph-series-filename> [max-jacobiators]\n\n"
             << "Accepts only homogeneous power series: graphs with n internal vertices at order n.\n"
             << "Optional argument max-jacobiators restricts the number of Jacobiators in each generated differential consequence.\n";
        return 1;
    }

    size_t max_jacobiators = numeric_limits<size_t>::max();
    if (argc == 3)
        max_jacobiators = stoi(argv[2]);

    // Reading in graph series
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    parser coefficient_reader;
    bool homogeneous = true;
    map<size_t, set< vector<size_t> > > in_degrees;
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file,
        [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); },
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

    graph_series.reduce();

    size_t counter = 0;
    std::vector<symbol> coefficient_list;

    for (size_t n = 2; n <= order; ++n) // need at least 2 internal vertices for Jacobi
    {
        if (graph_series[n].size() == 0)
            continue;

        cout << "h^" << n << ":\n";

        // First we choose the target vertices i,j,k of the Jacobiators (which contain 2 bivectors), in increasing order (without loss of generality)
        
        // Jacobi must have three distinct arguments, and not act on itself, but can act on other Jacobi

        for (size_t k = 1; k <= min(n/2, max_jacobiators); ++k)
        {
            std::vector<size_t> jacobi_vertices(3*k, n + 3);
            CartesianProduct jacobi_indices(jacobi_vertices);
            for (auto jacobi_index = jacobi_indices.begin(); jacobi_index != jacobi_indices.end(); ++jacobi_index)
            {
                bool accept = true;
                for (size_t i = 0; i != k; ++i)
                {
                    if ((*jacobi_index)[i*3] >= (*jacobi_index)[i*3 + 1] || (*jacobi_index)[i*3 + 1] >= (*jacobi_index)[i*3 + 2]) // not strictly increasing
                    {
                        accept = false;
                        break;
                    }
                }
                if (!accept)
                    continue;

                // Then we choose the target vertices of the remaining n - 2*k bivectors, stored in a multi-index of length 2*(n-2*k)
                // Here we have k fewer possible targets: out of the last 2*k internal vertices, the first k act as placeholders for the respective Jacobiators,
                // to be replaced by the Leibniz rule later on
                std::vector<size_t> remaining_edges(2*(n-2*k), n + 3 - k);
                CartesianProduct indices(remaining_edges);
                for (auto multi_index = indices.begin(); multi_index != indices.end(); ++multi_index)
                {
                    bool accept = true;
                    for (size_t idx = 0; idx != n - 2*k; ++idx)
                    {
                        if ((*multi_index)[2*idx] >= (*multi_index)[2*idx+1])
                        {
                            accept = false; // accept only strictly increasing indices
                            break;
                        }
                        // TODO: filter out tadpoles, maybe?
                    }
                    if (!accept)
                        continue;

                    // We build the list of targets for the graph, as described above (using i,j,k and the multi-index)
                    std::vector<KontsevichGraph::VertexPair> targets(n);
                    // first part:
                    for (size_t idx = 0; idx != n - 2*k; ++idx)
                        targets[idx] = {(*multi_index)[2*idx], (*multi_index)[2*idx+1]};
                    // second part:
                    for (size_t i = 0; i != k; ++i)
                    {
                        targets[n - 2*k + 2*i].first = KontsevichGraph::Vertex((*jacobi_index)[3*i]);
                        targets[n - 2*k + 2*i].second = KontsevichGraph::Vertex((*jacobi_index)[3*i + 1]);
                        targets[n - 2*k + 2*i + 1].first = KontsevichGraph::Vertex(n + 3 - 2*k + 2*i);
                        targets[n - 2*k + 2*i + 1].second = KontsevichGraph::Vertex((*jacobi_index)[3*i + 2]);
                    }

                    // Make vector of references to bad targets: those in first part with target >= (n + 3 - 2*k), the placeholders for the Jacobiators:
                    std::map<KontsevichGraph::Vertex*, int> bad_targets;
                    for (size_t idx = 0; idx != n - 2*k; ++idx) // look for bad targets in first part
                    {
                        if ((int)targets[idx].first >= (int)n + 3 - 2*(int)k)
                            bad_targets[&targets[idx].first] = (int)targets[idx].first - (n + 3 - 2*k);
                        if ((int)targets[idx].second >= (int)n + 3 - 2*(int)k)
                            bad_targets[&targets[idx].second] = (int)targets[idx].second - (n + 3 - 2*k);
                    }

                    KontsevichGraphSum<ex> graph_sum;

                    // Replace bad targets by Leibniz rule:
                    map< vector<size_t>, symbol > coefficients;
                    std::vector<size_t> leibniz_sizes(bad_targets.size(), 2);
                    CartesianProduct leibniz_indices(leibniz_sizes);
                    for (auto leibniz_index = leibniz_indices.begin(); leibniz_index != leibniz_indices.end(); ++leibniz_index)
                    {
                        size_t idx = 0;
                        for (auto& bad_target : bad_targets)
                            *(bad_target.first) = KontsevichGraph::Vertex(3 + n - 2*k + 2*(bad_target.second) + (*leibniz_index)[idx++]);

                        for (size_t i = 0; i != k; ++i)
                        {
                            for (auto jacobi_targets_choice : std::vector< std::vector<KontsevichGraph::Vertex> >({ { targets[n-2*k+2*i].first, targets[n-2*k+2*i].second, targets[n-2*k+2*i+1].second },
                                                                                                                    { targets[n-2*k+2*i].second, targets[n-2*k+2*i+1].second, targets[n-2*k+2*i].first },
                                                                                                                    { targets[n-2*k+2*i+1].second, targets[n-2*k+2*i].first, targets[n-2*k+2*i].second } }))
                            {
                                // Set Jacobiator targets to one of the three permutatations
                                targets[n-2*k+2*i].first = jacobi_targets_choice[0];
                                targets[n-2*k+2*i].second = jacobi_targets_choice[1];
                                targets[n-2*k+2*i+1].second = jacobi_targets_choice[2];

                                KontsevichGraph graph(n, 3, targets);

                                vector<size_t> indegrees = graph.in_degrees();

                                if (in_degrees[n].find(indegrees) == in_degrees[n].end()) // skip terms
                                    continue;
                                if (coefficients.find(indegrees) == coefficients.end())
                                {
                                    symbol coefficient("c_" + to_string(k) + "_" + to_string(counter) + "_" + to_string(indegrees[0]) + to_string(indegrees[1]) + to_string(indegrees[2]));
                                    coefficients[indegrees] = coefficient;
                                }
                                graph_sum += KontsevichGraphSum<ex>({ { coefficients[indegrees], graph } });
                            }
                        }
                    }

                    graph_sum.reduce();
                    if (graph_sum.size() != 0)
                    {
                        cerr << "\r" << ++counter;
                        for (auto& pair : coefficients)
                        {
                            coefficient_list.push_back(pair.second);
                        }
                    }
                    graph_series[n] -= graph_sum;
                }
            }
        }
    }

    cout << "\nNumber of coefficients: " << coefficient_list.size() << "\n";
    cout << "\nNumber of terms: " << graph_series[order].size() << "\n";
    cout << "\nNumber of terms per coefficient: " << (float)graph_series[order].size()/coefficient_list.size() << "\n";

    cout.flush();

    cerr << "\nReducing...\n";
    graph_series.reduce();

    lst equations;

    for (size_t n = 0; n <= order; ++n)
        for (auto& term : graph_series[n])
        {
            cerr << term.first << "==0\n";
            equations.append(term.first);
        }

    // Set up sparse matrix linear system

    cerr << "Setting up linear system...\n";
    size_t rows = equations.nops();
    size_t cols = coefficient_list.size();

    Eigen::VectorXd b(rows);
    SparseMatrix matrix(rows,cols);

    std::vector<Triplet> tripletList;
    size_t idx = 0;
    for (ex equation : equations)
    {
        if (!is_a<add>(equation))
            equation = lst(equation);
        for (ex term : equation)
        {
            if (!is_a<mul>(term))
                term = lst(term);
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
    
    cerr << "Solving linear system numerically...\n";

    Eigen::SparseQR< SparseMatrix, Eigen::COLAMDOrdering<int> > qr(matrix);
    Eigen::VectorXd x = qr.solve(b);
    
    cerr << "Residual norm = " << (matrix * x - b).squaredNorm() << "\n";

    cerr << "Rounding...\n";
    x = x.unaryExpr([](double elem) { return fabs(elem) < threshold ? 0.0 : elem; });

    cerr << "Still a solution? Residual norm = " << (matrix * x - b).squaredNorm() << "\n";

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

    for (auto& order: graph_series)
        for (auto& term : graph_series[order.first])
            term.first = term.first.subs(zero_substitution);

    cerr << "Reducing zeros...\n";

    graph_series.reduce();

    for (size_t n = 0; n <= graph_series.precision(); ++n)
    {
        cout << "h^" << n << ":\n";
        for (auto& term : graph_series[n])
        {
            cout << term.second.encoding() << "    " << term.first << "\n";
        }
    }

    cerr << "Verifying solution...\n";

    for (auto& order: graph_series)
        for (auto& term : graph_series[order.first])
            term.first = term.first.subs(solution_substitution);

    graph_series.reduce();

    cout << "Do we really have a solution? " << (graph_series == 0 ? "Yes" : "No") << "\n";

    for (ex subs : solution_substitution)
        cout << subs << "\n";
}
