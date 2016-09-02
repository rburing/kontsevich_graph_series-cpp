#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
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
    if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " <graph-series-filename>\n\n"
             << "Accepts only homogeneous power series: graphs with n internal vertices at order n.\n";
        return 1;
    }

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

        // First we choose the target vertices i,j,k of the Jacobiator (which contains 2 bivectors), in increasing order (without loss of generality)
        for (size_t i = 0; i != n + 3 - 4; ++i)
        {
            for (size_t j = i + 1; j != n + 3 - 3; ++j)
            {
                for (size_t k = j + 1; k != n + 3 - 2; ++k)
                {
                    // Then we choose the target vertices of the remaining n - 2 bivectors, stored in a multi-index of length 2*(n-2)
                    // Here we have one fewer possible target: the internal vertex (n + 3 - 2) acts as a placeholder for the Jacobiator, to be replaced by the Leibniz rule later on
                    std::vector<size_t> remaining_edges(2*(n-2), n + 3 - 1);
                    CartesianProduct indices(remaining_edges);
                    for (auto multi_index = indices.begin(); multi_index != indices.end(); ++multi_index)
                    {
                        bool accept = true;
                        for (size_t idx = 0; idx != n - 2; ++idx)
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
                        for (size_t idx = 0; idx != n - 2; ++idx)
                            targets[idx] = {(*multi_index)[2*idx], (*multi_index)[2*idx+1]};
                        // second part:
                        targets[n-1].first = static_cast<KontsevichGraph::Vertex>(n + 3 - 2);
                        targets[n-1].second = k;
                        targets[n-2].first = i;
                        targets[n-2].second = j;

                        std::vector<KontsevichGraph::Vertex*> jacobi_targets { &targets[n-2].first, &targets[n-2].second, &targets[n-1].second };
                        std::vector<KontsevichGraph::Vertex> jacobi_vertices { // to be used for edges incoming on the Jacobiator, applying the Leibniz rule
                            KontsevichGraph::Vertex(i),
                            KontsevichGraph::Vertex(j),
                            KontsevichGraph::Vertex(k),
                            KontsevichGraph::Vertex(n + 3 - 2),
                            KontsevichGraph::Vertex(n + 3 - 1)
                        };

                        // Make vector of references to bad targets: those in first part with target equal to (n + 3 - 2), the placeholder for the Jacobiator:
                        std::vector<KontsevichGraph::Vertex*> bad_targets;
                        for (size_t idx = 0; idx != n - 2; ++idx) // look for bad targets in first part
                        {
                            if (targets[idx].first == KontsevichGraph::Vertex(n + 3 - 2))
                                bad_targets.push_back(&targets[idx].first);
                            if (targets[idx].second == KontsevichGraph::Vertex(n + 3 - 2))
                                bad_targets.push_back(&targets[idx].second);
                        }

                        KontsevichGraphSum<ex> graph_sum;

                        map< std::pair< vector<size_t>, vector<size_t> >, symbol > coefficients;
                        for (auto jacobi_targets_choice : std::vector< std::vector<KontsevichGraph::Vertex> >({ { targets[n-2].first, targets[n-2].second, targets[n-1].second },
                                                                                                                { targets[n-2].second, targets[n-1].second, targets[n-2].first },
                                                                                                                { targets[n-1].second, targets[n-2].first, targets[n-2].second } }))
                        {
                            // Set Jacobiator targets to one of the three permutatations
                            targets[n-2].first = jacobi_targets_choice[0];
                            targets[n-2].second = jacobi_targets_choice[1];
                            targets[n-1].second = jacobi_targets_choice[2];
                            // Replace bad targets by Leibniz rule:
                            std::vector<size_t> leibniz_sizes(bad_targets.size(), jacobi_vertices.size());
                            CartesianProduct leibniz_indices(leibniz_sizes);
                            for (auto leibniz_index = leibniz_indices.begin(); leibniz_index != leibniz_indices.end(); ++leibniz_index)
                            {
                                for (size_t idx = 0; idx != bad_targets.size(); ++idx)
                                {
                                    *bad_targets[idx] = jacobi_vertices[(*leibniz_index)[idx]];
                                }
                                KontsevichGraph graph(n, 3, targets);
                                vector<size_t> indegrees = graph.in_degrees();
                                vector<size_t> jacobi_indegrees({ 1, 1, 1 });
                                for (size_t idx = 0; idx != bad_targets.size(); ++idx)
                                {
                                    if ((*leibniz_index)[idx] < 3)
                                    {
                                        ++jacobi_indegrees[(*leibniz_index)[idx]]; // this is correct, because the targets of Jacobi are permuted in place
                                    }
                                }
                                if (jacobi_indegrees != vector<size_t>({ 1, 1, 1})) // it suffices to restrict the internal Jacobi to in-degree 1,1,1
                                    continue;
                                if (in_degrees[n].find(indegrees) == in_degrees[n].end()) // skip terms
                                    continue;
                                if (coefficients.find({ indegrees, jacobi_indegrees }) == coefficients.end())
                                {
                                    symbol coefficient("c_" + to_string(counter) + "_" + to_string(indegrees[0]) + to_string(indegrees[1]) + to_string(indegrees[2]) + "_" + to_string(jacobi_indegrees[0]) + to_string(jacobi_indegrees[1]) + to_string(jacobi_indegrees[2]));
                                    coefficients[{ indegrees, jacobi_indegrees}] = coefficient;
                                }
                                graph_sum += KontsevichGraphSum<ex>({ { coefficients[{ indegrees, jacobi_indegrees } ], graph } });
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
            equations.append(term.first);

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
