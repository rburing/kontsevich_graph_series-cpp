#include "../kontsevich_graph_series.hpp"
#include <iostream>
using namespace std;

int main()
{
    KontsevichGraph g(2, 2, { {0, 1}, {0, 2} });
    cout << g << "\n";
    cout << "Edges: ";
    for (size_t v : g.internal_vertices())
    {
        std::pair<size_t, size_t> targets = g.targets(v);
        cout << "(" << v << ", " << targets.first << ", 'L') ";
        cout << "(" << v << ", " << targets.second<< ", 'R') ";
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
    gs.reduce();
    cout << gs << "\n";

    KontsevichGraphSeries<int> star({{0, gs}, {1, gs}});
    cout << star << "\n";

    cout << "Composition:\n";
    KontsevichGraph p(1, 2, { {0, 1} });
    KontsevichGraphSum<int> sum({ { 1, p } });
    KontsevichGraphSum<int> composition = sum({ sum, sum });
    composition.reduce();
    cout << composition << "\n";
    cout << composition.size() << "\n";
    for (auto& term : composition)
    {
        cout << term.first * term.second.sign() << "\t";
        for (size_t v : term.second.internal_vertices())
        {
            std::pair<size_t, size_t> targets = term.second.targets(v);
            cout << "(" << v << ", " << targets.first << ", 'L'), ";
            cout << "(" << v << ", " << targets.second<< ", 'R'), ";
        }
        cout << "\n";
    }
}
