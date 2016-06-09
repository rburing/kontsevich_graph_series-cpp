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
    if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " <graph-series-filename>\n";
        return 1;
    }
    
    // Reading in graph series:
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> graph_series;
    KontsevichGraphSum<ex> term;
    size_t order = 0;
    for (string line; getline(graph_series_file, line); )
    {
        if (line.length() == 0)
            continue;
        if (line[0] == 'h')
        {
            graph_series[order] = term;
            term = KontsevichGraphSum<ex>({ });
            order = stoi(line.substr(2));
        }
        else
        {
            KontsevichGraph graph;
            stringstream ss(line);
            ss >> graph;
            string coefficient_str;
            ss >> coefficient_str;
            ex coefficient = coefficient_reader(coefficient_str);
            term += KontsevichGraphSum<ex>({ { coefficient, graph } });
        }
    }
    graph_series[order] = term; // the last one

    graph_series.reduce();

    for (size_t n = 0; n <= order; ++n)
    {
        cout << "h^" << n << ":\n";
        for (auto& term : graph_series[n])
        {
            cout << term.second.encoding() << "    " << term.first << "\n";
        }
    }
}
