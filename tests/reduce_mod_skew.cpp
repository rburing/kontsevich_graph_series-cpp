#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 2 && argc != 3 && argc != 4 && argc != 5)
    {
        cout << "Usage: " << argv[0] << " <graph-series-filename> [--print-differential-orders] [--modulo-reversion] [--print-variables]\n";
        return 1;
    }

    bool print_differential_orders = false;
    bool modulo_reversion = false;
    bool print_variables = false;
    if (argc >= 3)
        for (int j = 2; j < argc; ++j)
        {
            string argument(argv[j]);
            if (argument == "--print-differential-orders")
                print_differential_orders = true;
            else if (argument == "--modulo-reversion")
                modulo_reversion = true;
            else if (argument == "--print-variables")
                print_variables = true;
        }
    
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
            if (modulo_reversion)
            {
                auto reversed_indegree = indegree;
                reverse(reversed_indegree.begin(), reversed_indegree.end());
                if (reversed_indegree < indegree)
                    continue;
            }
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
    if (print_variables)
    {
        cerr << "Number of variables: " << coefficient_reader.get_syms().size() << "\n";
        cerr << "Variables: {";
        size_t cnt = 0;
        for (auto pair: coefficient_reader.get_syms())
        {
            cerr << pair.first;
            if (++cnt != coefficient_reader.get_syms().size())
                cerr << ", ";
        }
        cerr << "}\n";
    }
}
