#include "../kontsevich_graph_series.hpp"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 2 && argc != 3)
    {
        cout << "Usage: " << argv[0] << " <internal> [external]\n";
        return 1;
    }
    size_t internal = stoi(argv[1]);
    size_t external = 2;
    if (argc == 3)
        external = stoi(argv[2]);

    size_t counter = 0;
    KontsevichGraph::graphs(internal, external, true, true,
        [internal, &counter](KontsevichGraph g)
        {
            cout << g.encoding() << "    ";
            if (g.is_zero())
                cout << "0\n";
            else
                cout << "w_" << internal << "_" << (++counter) << "\n";
        },
        [](KontsevichGraph g) -> bool
        {
            return g.positive_differential_order() && g.is_prime();
        }
    );
}
