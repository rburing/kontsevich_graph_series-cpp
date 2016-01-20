#include <kontsevich_graph.hpp>
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
    cout << "Sign: " << g.sign() << "\n";
    KontsevichGraph g2(2, 2, { {1, 0}, {2, 0} });
    cout << "Equality testing " << (g == g2 ? "works" : "fails") << ".\n";
}
