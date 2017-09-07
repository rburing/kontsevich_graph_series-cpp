#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits>
#include <sstream>
using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 2 && argc != 3)
    {
        cout << "Usage: " << argv[0] << " <graph-series-filename> [--print-differential-orders]\n";
        return 1;
    }

    bool print_differential_orders = false;
    if (argc == 3 && string(argv[2]) == "--print-differential-orders")
        print_differential_orders = true;
    
    // Reading in graph series:
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });
    graph_series.reduce_mod_skew();

    for (size_t n = 0; n <= graph_series.precision(); ++n)
    {
        if (graph_series[n] != 0 || n == graph_series.precision())
            cout << "h^" << n << ":\n";
        for (auto& indegree : graph_series[n].in_degrees(true))
        {
            if (print_differential_orders)
            {
                cout << "# ";
                for (size_t in : indegree)
                    cout << in << " ";
                cout << "\n";
            }
            for (auto& term : graph_series[n][indegree])
            {
                cout << term.second.encoding() << "    " << term.first << "\n";
            }
        }
    }
}
