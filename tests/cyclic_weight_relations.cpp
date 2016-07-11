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
        cout << "Usage: " << argv[0] << " <star-product-filename>\n";
        return 1;
    }

    // Reading in star product:
    string star_product_filename(argv[1]);
    ifstream star_product_file(star_product_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> star_product = KontsevichGraphSeries<ex>::from_istream(star_product_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });
    size_t order = star_product.precision();

    // Make a list of weight variables
    lst weight_vars;
    for (std::pair<string, ex> pair : coefficient_reader.get_syms())
        weight_vars.append(pair.second);

    // Computing cyclic linear weight relations:
    lst weight_relations;
    for (auto& term : star_product[order])
    {
        KontsevichGraph graph = term.second;
        // coefficient = multiplicity * weight / n! => weight = n! * coefficient / multiplicity
        ex linear_combination = (graph.internal() % 2 ? 1 : -1) * star_product[order][graph] / graph.multiplicity();
        // iterate over subsets of set of edges: 2^(2n) elements
        vector<size_t> selector(2*graph.internal(), 2);
        CartesianProduct edge_selector(selector);
        for (auto edge_selection = edge_selector.begin(); edge_selection != edge_selector.end(); ++edge_selection)
        {
            std::vector<KontsevichGraph::VertexPair> targets = graph.targets();

            // skip those where an edge is already connected to the first ground vertex
            bool skip = false;
            for (size_t choice = 0; choice != 2*graph.internal(); ++choice)
                if ((*edge_selection)[choice] && ((choice % 2 == 0 && targets[choice/2].first == 0) || (choice % 2 == 1 && targets[choice/2].second == 0)))
                    skip = true;
            if (skip)
                continue;

            // replace edges by edges to the first ground vertex:
            for (size_t choice = 0; choice != 2*graph.internal(); ++choice)
            {
                if ((*edge_selection)[choice]) // replace edge
                {
                    if (choice % 2 == 0)
                        targets[choice/2].first = 0;
                    else
                        targets[choice/2].second = 0;
                }
            }
            KontsevichGraph modified_graph(graph.internal(), graph.external(), targets, graph.sign());
            // build linear relation
            linear_combination -= (modified_graph.in_degree(0) % 2 ? 1 : -1) * star_product[order][modified_graph] / modified_graph.multiplicity();
        }
        if (linear_combination != 0)
        {
            for (std::pair<string, ex> pair : coefficient_reader.get_syms())
            {
                if (linear_combination.coeff(pair.second) != 0) // first weight variable occurring in expression
                {
                    linear_combination /= linear_combination.coeff(pair.second); // divide by its coefficient
                    break;
                }
            }
            weight_relations.append(linear_combination == 0);
        }
    }
    for (ex const& weight_relation : lsolve(weight_relations, weight_vars))
    {
        cout << weight_relation << "\n";
    }
}
