#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
using namespace std;
using namespace GiNaC;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "Usage: " << argv[0] << " <graph-series-filename>\n";
        return 1;
    }

    // Reading in graph series:
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });
    size_t order = graph_series.precision();

    for (size_t n = 0; n <= order; ++n)
    {
        if (graph_series[n] != 0 || n == order)
        {
            if (n > 0)
                cout << "+\\hbar";
            if (n >= 2)
                cout << "^{" << n << "}";
            if (graph_series[n].size() > 1)
                cout << "\\big(\n";
        }
        size_t term_number = 0;
        for (auto& term : graph_series[n])
        {
            KontsevichGraph& graph = term.second;
            ex coefficient = graph.sign() * term.first;
            if (graph.sign() == 0 || coefficient == 0)
                continue;
            graph.sign(1);
            if (coefficient != 1)
            {
                stringstream coefficient_stream;
                coefficient_stream << latex << coefficient;
                string coefficient_string = coefficient_stream.str();
                if (coefficient_string[0] != '-' && term_number != 0)
                    coefficient_string = "+" + coefficient_string;
                cout << coefficient_string << " ";
            }
            vector< vector<string> > indices = { { "i", "j" }, { "k", "\\ell" }, { "m", "n" }, { "p", "q" }, { "r", "s" }, { "t", "v"} };
            vector<string> functions = { "f", "g", "h"};
            vector<string> factors(graph.vertices());
            for (size_t idx = 0; idx != graph.external(); ++idx)
                factors[idx] = functions[idx];
            for (size_t idx = 0; idx != graph.internal(); ++idx)
                factors[graph.external() + idx] = "\\cP^{" + indices[idx][0] + indices[idx][1] + "}";
            for (KontsevichGraph::Vertex v : graph.internal_vertices())
            {
                size_t idx = (size_t)v - graph.external();
                KontsevichGraph::VertexPair targets = graph.targets(v);
                factors[targets.first]  = "\\partial_{" + indices[idx][0] + "} " + factors[targets.first];
                factors[targets.second] = "\\partial_{" + indices[idx][1] + "} " + factors[targets.second];
            }
            for (size_t idx = 0; idx != graph.internal(); ++idx)
                cout << factors[graph.external() + idx] << " ";
            for (size_t idx = 0; idx != graph.external(); ++idx)
                cout << factors[idx] << " ";
            cout << "\n";
            ++term_number;
        }
        if (graph_series[n] != 0 || n == order)
        {
            if (graph_series[n].size() > 1)
                cout << "\\big)\n";
        }
    }
}
