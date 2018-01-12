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
             << "--max-jac-indegree=k   restricts the number of arrows falling on Jacobiators to be <= k.\n"
             << "--skew-leibniz         skew-symmetrizes each Leibniz graph before subtracting it with an undetermined coefficient.\n"
             << "--leibniz-in=filename  input graph series already contains Leibniz graphs, with encodings in filename.\n"
             << "--coeff-prefix=c       let the coefficients of leibniz graphs be c_n.\n"
             << "--solve                the undetermined variables in the input are added to the linear system to-be-solved.\n"
             << "--interactive          ask whether to continue to the next iteration, and whether to solve numerically.\n";
        return 1;
    }

    bool interactive = false;
    bool solve = false;
    size_t max_jac_indegree = numeric_limits<size_t>::max();
    bool skew_leibniz = false;
    string leibniz_in_filename = "";
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
            if (key == "--max-jac-indegree")
                max_jac_indegree = stoi(value);
            else if (key == "--coeff-prefix")
                coefficient_prefix = value;
            else if (key == "--leibniz-in")
                leibniz_in_filename = value;
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
            else if (argument == "--interactive")
                interactive = true;
            else {
                cout << "Unrecognized option: " << argument << "\n";
                return 1;
            }
        }
    }

    cout << "Options: "
         << "max-jac-indegree = ";
    if (max_jac_indegree == numeric_limits<size_t>::max())
        cout << "none";
    else
        cout << max_jac_indegree;
    cout << ", solve = " << (solve ? "yes" : "no")
         << ", skew-leibniz = " << (skew_leibniz ? "yes" : "no")
         << ", leibniz-in = " << (leibniz_in_filename == "" ? "none" : leibniz_in_filename)
         << ", coeff-prefix = " << coefficient_prefix
         << ", interactive = " << (interactive ? "yes" : "no") << "\n";

    // Reading in Leibniz graphs
    parser coefficient_reader;
    map<LeibnizGraph, symbol> leibniz_graphs;

    if (leibniz_in_filename != "")
    {
        ifstream leibniz_in_file(leibniz_in_filename);
        leibniz_graphs = LeibnizGraph::map_from_istream<symbol>(leibniz_in_file,
                                                                [&coefficient_reader](string s) -> symbol
                                                                {
                                                                    return ex_to<symbol>(coefficient_reader(s));
                                                                });
    }

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

    // Find number of external vertices
    // TODO: move to method in KontsevichGraphSum / KontsevichGraphSeries?
    size_t external = 0;
    for (size_t n = 0; n <= order; ++n)
    {
        if (graph_series[n].size() != 0)
        {
            external = graph_series[n].front().second.external();
            break;
        }
    }
    graph_series.reduce_mod_skew();

    vector<symbol> unknowns_list;
    for (auto namevar : coefficient_reader.get_syms())
        unknowns_list.push_back(ex_to<symbol>(namevar.second));

    if (unknowns_list.size() != 0 && !solve)
    {
        cerr << "Input contains unknowns; remove them or specify --solve to solve for them.\n";
        return 1;
    }

    vector<symbol> coefficient_list = unknowns_list;

    KontsevichGraphSeries<ex> leibniz_graph_series = graph_series;

    set<KontsevichGraph> processed_graphs;

    size_t counter = 0;
    bool converged = false;
    size_t step = 0;
    while (!converged)
    {
        ++step;
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

                for (KontsevichGraph::Vertex v : graph.internal_vertices())
                {
                    for (KontsevichGraph::Vertex w : graph.neighbors_in(v))
                    {
                        vector<KontsevichGraph::VertexPair> targets_template = graph.targets();
                        KontsevichGraph::VertexPair& target_pair_v = targets_template[(size_t)v - external];
                        // Check that there is no loop between v and w
                        if (target_pair_v.first == w || target_pair_v.second == w)
                            continue;
                        KontsevichGraph::VertexPair& target_pair_w = targets_template[(size_t)w - external];
                        // Check that the "Jacobiator" consisting of v and w falls on 3 distinct targets
                        KontsevichGraph::Vertex& a = target_pair_v.first;
                        KontsevichGraph::Vertex& b = target_pair_v.second;
                        KontsevichGraph::Vertex& c = (target_pair_w.first == v) ? target_pair_w.second : target_pair_w.first;
                        if (c == a || c == b)
                            continue;
                        
                        // Build the list of references to bad targets (incoming edges on v or w, except the edge w -> v)
                        set<KontsevichGraph::Vertex*> bad_targets;
                        for (KontsevichGraph::Vertex u : graph.internal_vertices())
                        {
                            if (u == w)
                                continue;
                            KontsevichGraph::VertexPair& target_pair = targets_template[(size_t)u - external];
                            if (target_pair.first == v || target_pair.first == w)
                                bad_targets.insert(&target_pair.first);
                            if (target_pair.second == v || target_pair.second == w)
                                bad_targets.insert(&target_pair.second);
                        }
                        if (bad_targets.size() > max_jac_indegree)
                            continue;

                        symbol coefficient(coefficient_prefix + "_" + to_string(counter));

                        // Normal form of Leibniz graph (three permutations of Jacobi, take minimal encoding, remember where Jacobiator is)
                        for (auto& bad_target : bad_targets)
                            *bad_target = v;
                        vector< vector<KontsevichGraph::Vertex> > jacobi_targets_choices({ { a, b, c },
                                                                                           { b, c, a },
                                                                                           { c, a, b } });
                        vector<LeibnizGraph> my_leibniz_graphs;

                        vector<KontsevichGraph::Vertex> ground_vertices { 0, 1, 2};
                        do
                        {
                            for (auto jacobi_targets_choice : jacobi_targets_choices)
                            {
                                // Set Jacobiator targets to one of the three permutatations
                                a = jacobi_targets_choice[0];
                                b = jacobi_targets_choice[1];
                                c = jacobi_targets_choice[2];

                                vector<KontsevichGraph::VertexPair> d_targets = targets_template;

                                size_t idx = 0;
                                for (KontsevichGraph::VertexPair& target_pair : d_targets)
                                {
                                    if ((size_t)target_pair.first < graph.external())
                                        target_pair.first = ground_vertices[idx++];
                                    if ((size_t)target_pair.second < graph.external())
                                        target_pair.second = ground_vertices[idx++];
                                }

                                // Find permutation of vertex labels such that the list of targets is minimal with respect to the defined ordering
                                std::vector<KontsevichGraph::VertexPair> global_minimum = d_targets;
                                sort_pairs(global_minimum.begin(), global_minimum.end());

                                size_t d_internal = graph.internal();
                                size_t d_external = graph.external();

                                KontsevichGraph::VertexPair new_vw = {v, w};
                                std::vector<KontsevichGraph::Vertex> vertices(d_external + d_internal);
                                std::iota(vertices.begin(), vertices.end(), 0);
                                while (std::next_permutation(vertices.begin() + d_external, vertices.end()))
                                {
                                    std::vector<KontsevichGraph::VertexPair> local_minimum = d_targets;
                                    apply_permutation(d_internal, d_external, local_minimum, vertices);
                                    if (local_minimum < global_minimum)
                                    {
                                        global_minimum = local_minimum;
                                        // Find where Jacobiator is
                                        new_vw = { vertices[(size_t)v], vertices[(size_t)w] };
                                    }
                                }
                                d_targets = global_minimum;

                                KontsevichGraph leibniz_graph(d_targets.size(), d_external, d_targets, 1, true);
                                my_leibniz_graphs.push_back({ leibniz_graph, { new_vw } });
                            }
                        } while (skew_leibniz && std::next_permutation(ground_vertices.begin(), ground_vertices.end()));
                        LeibnizGraph leibniz_normal_form = *min_element(my_leibniz_graphs.begin(), my_leibniz_graphs.end());

                        if (leibniz_graphs.find(leibniz_normal_form) != leibniz_graphs.end())
                            continue;

                        KontsevichGraphSum<ex> graph_sum;

                        // Replace bad targets by Leibniz rule:
                        vector<size_t> leibniz_sizes(bad_targets.size(), 2);
                        CartesianProduct leibniz_indices(leibniz_sizes);
                        for (auto leibniz_index = leibniz_indices.begin(); leibniz_index != leibniz_indices.end(); ++leibniz_index)
                        {
                            size_t idx = 0;
                            for (auto& bad_target : bad_targets)
                                *bad_target = (*leibniz_index)[idx++] == 0 ? v : w;
                            for (auto jacobi_targets_choice : jacobi_targets_choices)
                            {
                                // Set Jacobiator targets to one of the three permutatations
                                a = jacobi_targets_choice[0];
                                b = jacobi_targets_choice[1];
                                c = jacobi_targets_choice[2];

                                KontsevichGraph new_graph(targets_template.size(), external, targets_template, graph.sign());

                                graph_sum += KontsevichGraphSum<ex>({ { coefficient, new_graph } });
                            }
                        }

                        if (skew_leibniz)
                            graph_sum = graph_sum.skew_symmetrization();

                        graph_sum.reduce_mod_skew();
                        if (graph_sum.size() != 0)
                        {
                            ++counter;
                            coefficient_list.push_back(coefficient);
                            leibniz_graphs[leibniz_normal_form] = coefficient;
                        }
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

        char solve = 'N';
        if (interactive)
        {
            cerr << "Solve? (Y/N) ";
            cin >> solve;
        }

        if (solve == 'Y')
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
                    for (auto pair : leibniz_graphs)
                        if (pair.second != ex(pair.second).subs(solution_substitution))
                            cout << pair.first.first.encoding() << "    " << pair.second << "==" << ex(pair.second).subs(solution_substitution) << "\n";
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
    cout << "\nConverged in " << step << " steps.\n\n";
    if (skew_leibniz)
        cout << "Skew-";
    cout << "Leibniz graphs:\n";
    for (auto pair : leibniz_graphs)
    {
        cout << pair.first.encoding() << "    " << pair.second << "\n";
    }
}
