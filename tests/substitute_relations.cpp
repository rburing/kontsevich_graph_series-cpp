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
        cout << "Usage: " << argv[0] << " <graphs-and-weights-filename> <weight-relations-filename>\n";
        return 1;
    }
    string graphs_weights_filename(argv[1]), weight_relations_filename(argv[2]);

    // Reading in known weight relations:
    ifstream weight_relations_file(weight_relations_filename);
    lst weight_relations;
    parser weights_reader;
    for (string lhs, rhs; getline(weight_relations_file, lhs, '=') && weight_relations_file.ignore(1) && getline(weight_relations_file, rhs); )
    {
        weight_relations.append(weights_reader(lhs) == weights_reader(rhs));
    }

    // Reading in known graphs and their (possibly symbolic) weights:
    ifstream weights_file(graphs_weights_filename);
    weights_file.ignore(numeric_limits<streamsize>::max(), '\n'); // ignore first line with column headers
    KontsevichGraph graph;
    string weight_str;
    cout << "Graph\t\t\tWeight\n";
    while (weights_file >> graph >> weight_str)
    {
        ex weight = weights_reader(weight_str);
        cout << graph.encoding() << "    " << weight.subs(weight_relations) << "\n";
    }
}
