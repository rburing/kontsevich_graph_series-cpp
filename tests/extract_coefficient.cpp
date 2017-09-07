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
        cout << "Usage: " << argv[0] << " <graph-series-filename> <expression>\n\n"
             << "If <expression>=1, the constant part is taken.\n";
        return 1;
    }

    // Reading in graph series:
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });

    ex expression = coefficient_reader(string(argv[2]));

    lst allzero_substitution;
    if (expression == 1)
        for (auto named_symbol : coefficient_reader.get_syms())
            allzero_substitution.append(named_symbol.second==0);

    for (size_t n = 0; n <= graph_series.precision(); ++n)
    {
        if (graph_series[n] != 0 || n == graph_series.precision())
            cout << "h^" << n << ":\n";
        for (auto& term : graph_series[n])
        {
            if (expression == 1)
                term.first = term.first.subs(allzero_substitution);
            else
                term.first = term.first.coeff(expression);
            if (term.first != 0)
                cout << term.second.encoding() << "    " << term.first << "\n";
        }
    }
}
