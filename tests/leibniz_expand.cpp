#include "../leibniz_graph.hpp"
#include "../kontsevich_graph_sum.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " <leibniz-graph-series-filename>\n";
        return 1;
    }

    string leibniz_in_filename(argv[1]);

    // Reading in Leibniz graphs
    parser coefficient_reader;
    map< LeibnizGraph<ex>, ex> leibniz_graphs;

    ifstream leibniz_in_file(leibniz_in_filename);
    leibniz_graphs = LeibnizGraph<ex>::map_from_istream(leibniz_in_file,
                                                        [&coefficient_reader](string s) -> ex
                                                        {
                                                            return coefficient_reader(s);
                                                        });

    KontsevichGraphSum<ex> graph_sum;
    for (auto& pair : leibniz_graphs)
    {
        LeibnizGraph<ex> L(pair.first);
        L.skew(true); // XXX: temporary, assume all Leibniz graphs are skew-Leibniz graphs
        KontsevichGraphSum<ex> leibniz_expansion = L.expansion(pair.second);
        for (auto& term : leibniz_expansion)
        {
            cout << term.second.encoding() << "    " << term.first << "\n";
        }
    }
}
