#include <kontsevich_graph.hpp>
#include <iostream>
using namespace std;

int main()
{
    KontsevichGraph g(2, 2, { {0, 1}, {0, 2} });
    cout << g << "\n";
    cout << "Targets: \n";
    for (size_t v : g.internal_vertices())
    {
        std::pair<size_t, size_t> targets = g.targets(v);
        cout << v << ": " << std::get<0>(targets) << " " << std::get<1>(targets) << "\n";
    }
}
