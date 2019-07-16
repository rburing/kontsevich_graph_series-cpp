#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " <graph-series-filename1> <graph-series-filename2>\n";
        return 1;
    }

    // Reading in graph series:
    parser coefficient_reader;
    string graph_series_filename1(argv[1]);
    ifstream graph_series_file1(graph_series_filename1);
    KontsevichGraphSeries<ex> graph_series1 = KontsevichGraphSeries<ex>::from_istream(graph_series_file1, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });
    string graph_series_filename2(argv[2]);
    ifstream graph_series_file2(graph_series_filename2);
    KontsevichGraphSeries<ex> graph_series2 = KontsevichGraphSeries<ex>::from_istream(graph_series_file2, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });

    size_t order = min(graph_series1.precision(), graph_series2.precision());

    KontsevichGraphSeries<ex> gerstenhaber_bracket_series = gerstenhaber_bracket(graph_series1, graph_series2);

    for (size_t n = 0; n <= order; ++n)
    {
        if (gerstenhaber_bracket_series[n] != 0 || n == order)
            cout << "h^" << n << ":\n";
        for (auto& term : gerstenhaber_bracket_series[n])
        {
            cout << term.second.encoding() << "    " << term.first << "\n";
        }
    }
}
