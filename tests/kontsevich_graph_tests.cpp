#include "../kontsevich_graph_series.hpp"
#include <iostream>
#include <sstream>
#include <cmath>
using namespace std;

int main()
{
    KontsevichGraph g(2, 2, { {0, 1}, {0, 2} });
    cout << g << "\n";
    cout << "Edges: ";
    for (size_t v : g.internal_vertices())
    {
        KontsevichGraph::VertexPair targets = g.targets(v);
        cout << "(" << v << ", " << (size_t)targets.first << ", 'L') ";
        cout << "(" << v << ", " << (size_t)targets.second<< ", 'R') ";
    }
    cout << "\n";
    cout << "In-degrees of external vertices: ";
    for (size_t d : g.in_degrees())
    {
        cout << d << " ";
    }
    cout << "\n";
    cout << "Incoming neighbors of external vertices:\n";
    for (size_t i = 0; i < g.external(); ++i)
    {
        cout << i << ": ";
        for (size_t n : g.neighbors_in(i))
        {
            cout << n << " ";
        }
        cout << "\n";
    }
    cout << "Sign: " << g.sign() << "\n";
    KontsevichGraph g2(2, 2, { {1, 0}, {1, 0} });
    cout << "Equality testing " << (g == g2 ? "fails" : "works") << ".\n";

    KontsevichGraph g3(2, 2, { {1, 0}, {0, 2} });

    KontsevichGraphSum<int> gs({ {3, g}, {2, g2}, {1, g2}, {-1, g}, {5, g3} });
    cout << gs << "\n";
    gs.reduce_mod_skew();
    cout << gs << "\n";

    cout << "Series: ";
    KontsevichGraphSeries<int> star({{0, gs}, {1, gs}});
    star.precision(1);
    star.reduce_mod_skew();
    cout << star << "\n";
    cout << "Series precision: " << star.precision() << "\n";

    cout << "Composition:\n";
    KontsevichGraph p(1, 2, { {0, 1} });
    KontsevichGraphSum<int> sum({ { 1, p } });
    KontsevichGraphSum<int> composition = sum({ sum, sum });
    composition.reduce_mod_skew();
    cout << composition << "\n";
    cout << composition.size() << "\n";
    for (auto& term : composition)
    {
        cout << term.first * term.second.sign() << "\t";
        for (size_t v : term.second.internal_vertices())
        {
            KontsevichGraph::VertexPair targets = term.second.targets(v);
            cout << "(" << v << ", " << targets.first << ", 'L'), ";
            cout << "(" << v << ", " << targets.second<< ", 'R'), ";
        }
        cout << "\n";
    }
    stringstream ss;
    ss << "2 2 1 0 1 0 2" << std::endl;
    KontsevichGraph g_read;
    ss >> g_read;
    cout << g_read << "\n";
    cout << "Reading in graphs " << (g_read == g ? "works" : "fails") << "\n";

    // Compare with composition from SAGE.
    KontsevichGraphSum<int> total;
    ss << "1 	4 3 1 	5 6	0 1	2 3 "
       << "1 	4 3 1 	5 2	0 1	2 3 "
       << "1 	4 3 1 	5 3	0 1	2 3 "
       << "1 	4 3 1 	0 6	0 1	2 3 "
       << "1 	4 3 1 	0 2	0 1	2 3 "
       << "1 	4 3 1 	0 3	0 1	2 3 "
       << "1 	4 3 1 	1 6	0 1	2 3 "
       << "1 	4 3 1 	1 2	0 1	2 3 "
       << "1 	4 3 1 	1 3	0 1	2 3 " << std::endl;
    ss >> total;
    cout << "Total: " << total << "\n";
    cout << "Total size: " << total.size() << "\n";
    cout << "Do we agree with SAGE? " << (total == composition ? "Yes" : "No") << "\n";
    KontsevichGraphSum<int> difference = composition - total;
    difference.reduce_mod_skew();
    cout << "Difference: " << difference << "\n";
    cout << "Difference size: " << difference.size() << "\n";

    cout << "Series composition: ";
    KontsevichGraph one(0, 2, {});
    KontsevichGraphSum<int> onesum({ {1, one} });
    KontsevichGraphSeries<int> oneseries({ {0, onesum} });
    KontsevichGraphSeries<int> pseries({ {0, sum} });
    KontsevichGraphSeries<int> result = pseries ({ oneseries, oneseries });
    result.reduce_mod_skew();
    cout << result;
    cout << "\n";

    cout << "Generating graphs:\n";
    size_t n = 4;
    size_t num_graphs = 0;
    std::set<KontsevichGraph> graphs = KontsevichGraph::graphs(n, 2, true, false, nullptr, nullptr);
    for (auto& g : graphs)
    {
        std::vector<KontsevichGraph::VertexPair> targets = g.abs().second;
        cout << "graph: ";
        for (size_t i = 0; i != targets.size(); ++i)
            cout << targets[i].first << " " << targets[i].second << "\t";
        cout << " (sign " << g.sign() << ", ";
        cout << (g.is_prime() ? "" : "not ") << "prime";
        cout << ", multiplicity " << g.multiplicity() << ")"; 
        cout << "\n";
        num_graphs += g.multiplicity();
    }
    cout << "Number of graphs: " << num_graphs << "\n";
    cout << "(n*(n+1))^n = " << pow(n*(n+1), n) << "\n";
}
