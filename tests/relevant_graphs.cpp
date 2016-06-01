#include "../kontsevich_graph_series.hpp"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " <order>\n";
        return 1;
    }
    size_t order = stoi(argv[1]);
    set<KontsevichGraph> relevants = KontsevichGraph::graphs(order, 2, true, true,
                    [](KontsevichGraph g) -> bool
                    {
                        return g.positive_differential_order() && g.is_prime();
                    });
    size_t counter = 0;
    for (KontsevichGraph const& g : relevants)
    {
        cout << g.encoding() << "    ";
        if (g.is_zero())
            cout << "0\n";
        else
            cout << "w_" << order << "_" << (++counter) << "\n";
    }
}
