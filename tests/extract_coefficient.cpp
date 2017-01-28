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
        cout << "Usage: " << argv[0] << " <graph-series-filename> <expression>\n";
        return 1;
    }

    // Reading in graph series:
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });

    ex expression = coefficient_reader(string(argv[2]));

    for (size_t n = 0; n <= graph_series.precision(); ++n)
        for (auto& term : graph_series[n])
            term.first = term.first.coeff(expression);

    graph_series.reduce();

    for (size_t n = 0; n <= graph_series.precision(); ++n)
    {
        if (graph_series[n] != 0 || n == graph_series.precision())
            cout << "h^" << n << ":\n";
        for (auto& term : graph_series[n])
            cout << term.second.encoding() << "    " << term.first << "\n";
    }
}
