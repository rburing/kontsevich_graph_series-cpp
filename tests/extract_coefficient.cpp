#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cout << "Usage: " << argv[0] << " <graph-series-filename> <variable-name>\n";
        return 1;
    }

    // Reading in graph series:
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });
    graph_series.reduce();

    ex variable = coefficient_reader(string(argv[2]));

    for (size_t n = 0; n <= graph_series.precision(); ++n)
    {
        if (graph_series[n] != 0 || n == graph_series.precision())
            cout << "h^" << n << ":\n";
        for (auto& term : graph_series[n])
        {
            ex coefficient = term.first.coeff(variable);
            if (coefficient != 0)
                cout << term.second.encoding() << "    " << coefficient << "\n";
        }
    }
}
