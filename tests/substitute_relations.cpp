#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
#include <limits>

using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cout << "Usage: " << argv[0] << " <graphs-series-filename> <relations-filename>\n";
        return 1;
    }
    string graph_series_filename(argv[1]), relations_filename(argv[2]);

    // Reading in coefficient relations:
    ifstream relations_file(relations_filename);
    lst relations;
    parser coefficient_reader;
    for (string lhs, rhs; getline(relations_file, lhs, '=') && relations_file.ignore(1) && getline(relations_file, rhs); )
    {
        relations.append(coefficient_reader(lhs) == coefficient_reader(rhs));
    }

    // Reading in graphs and their (possibly symbolic) coefficients:
    ifstream graph_series_file(graph_series_filename);
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });

    for (size_t n = 0; n <= graph_series.precision(); ++n)
    {
        cout << "h^" << n << ":\n";
        for (auto& term : graph_series[n])
        {
            cout << term.second.encoding() << "    " << term.first.subs(relations) << "\n";
        }
    }
}
