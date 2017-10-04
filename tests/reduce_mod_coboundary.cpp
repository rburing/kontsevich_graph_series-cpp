#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
#include <string>
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
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });

    KontsevichGraphSeries<ex> graph_series_copy = graph_series;

    size_t counter = 0;

    KontsevichGraphSeries<ex> potential_primitives;
    potential_primitives.precision(graph_series.precision());
    for (size_t n = 0; n <= graph_series.precision(); ++n)
    {
        set<KontsevichGraph> normal_forms;
        for (auto& term : graph_series[n])
        {
            KontsevichGraph& graph = term.second;
            // See whether this graph could come from a coboundary
            KontsevichGraphSum<ex> potential_primitive; // the "primitive" is X in Q = d_P(X) = [[P, X]]
            // First possibility: the existence of a "top"
            vector<KontsevichGraph::VertexPair> tops;
            for (KontsevichGraph::Vertex v : graph.internal_vertices())
            {
                if (graph.in_degree(v) != 0) // TODO: optimize
                    continue;
                pair<KontsevichGraph::Vertex, KontsevichGraph::Vertex> targets = graph.targets(v);
                if ((size_t)targets.first < graph.external() && (size_t)targets.second >= graph.external())
                    tops.push_back({ v, targets.first });
                else if ((size_t)targets.second < graph.external() && (size_t)targets.first >= graph.external())
                    tops.push_back({ v, targets.second });
            }
            for (KontsevichGraph::VertexPair top : tops)
            {
                // swap the top and the last vertex
                KontsevichGraph::Vertex last = graph.internal() + graph.external() - 1;
                vector<KontsevichGraph::VertexPair> new_targets = graph.targets();
                for (KontsevichGraph::VertexPair& pair : new_targets)
                {
                    if (pair.first == top.first) pair.first = last;
                    else if (pair.first == last) pair.first = top.first;
                    if (pair.second == top.first) pair.second = last;
                    else if (pair.second == last) pair.second = top.first;
                }
                new_targets[top.first - graph.external()] = new_targets[last - graph.external()];
                // remove the last vertex
                new_targets.pop_back();
                // relabel the ground vertices
                for (KontsevichGraph::VertexPair& pair : new_targets)
                {
                    if (pair.first >= top.second) pair.first--;
                    if (pair.second >= top.second) pair.second--;
                }
                KontsevichGraph new_graph(graph.internal() - 1, graph.external() - 1, new_targets);
                new_graph.sign(1);
                KontsevichGraphSum<ex> kgs { {1, new_graph } };
                kgs = kgs.skew_symmetrization();
                kgs.reduce_mod_skew();
                pair< ex, KontsevichGraph> normal_form = *min_element(kgs.begin(), kgs.end());
                if (normal_forms.find(normal_form.second) == normal_forms.end())
                {
                    normal_forms.insert(normal_form.second);
                    symbol coefficient("t_" + to_string(counter++));
                    for (auto& term : kgs)
                        term.first *= coefficient;
                    potential_primitive += kgs;
                }
                /*
                cout << "OLD: " << graph.encoding() << " (top at " << top.first << ")\n";
                cout << "NEW: " << new_graph.encoding() << "\n";
                */
            }
            // Second possibility: the existence of a "wedge" on ground vertices
            vector<KontsevichGraph::Vertex> wedges;
            for (KontsevichGraph::Vertex v : graph.internal_vertices())
            {
                if (graph.in_degree(v) != 1) // TODO: optimize
                    continue;
                pair<KontsevichGraph::Vertex, KontsevichGraph::Vertex> targets = graph.targets(v);
                if ((size_t)targets.first < graph.external() && (size_t)targets.second < graph.external())
                    wedges.push_back(v);
            }
            for (KontsevichGraph::Vertex wedge : wedges)
            {
                KontsevichGraph::VertexPair wedge_targets = graph.targets(wedge);
                // swap the wedge with the last vertex
                KontsevichGraph::Vertex last = graph.internal() + graph.external() - 1;
                vector<KontsevichGraph::VertexPair> new_targets = graph.targets();
                for (KontsevichGraph::VertexPair& pair : new_targets)
                {
                    if (pair.first == wedge) pair.first = last;
                    else if (pair.first == last) pair.first = wedge;
                    if (pair.second == wedge) pair.second = last;
                    else if (pair.second == last) pair.second = wedge;
                }
                new_targets[wedge - graph.external()] = new_targets[last - graph.external()];
                // remove the last vertex
                new_targets.pop_back();
                // replace edges to the last vertex by edges to the first target of the wedge
                for (KontsevichGraph::VertexPair& pair : new_targets)
                {
                    if (pair.first == last) pair.first = wedge_targets.first;
                    if (pair.second == last) pair.second = wedge_targets.first;
                }
                // relabel the ground vertices
                for (KontsevichGraph::VertexPair& pair : new_targets)
                {
                    if (pair.first >= wedge_targets.second) pair.first--;
                    if (pair.second >= wedge_targets.second) pair.second--;
                }
                KontsevichGraph new_graph(graph.internal() - 1, graph.external() - 1, new_targets);
                new_graph.sign(1);
                /*
                cout << "OLD: " << graph.encoding() << " (wedge at " << wedge << ")\n";
                cout << "NEW: " << new_graph.encoding() << "\n";
                */
                KontsevichGraphSum<ex> kgs { {1, new_graph } };
                kgs = kgs.skew_symmetrization();
                kgs.reduce_mod_skew();
                pair< ex, KontsevichGraph> normal_form = *min_element(kgs.begin(), kgs.end());
                if (normal_forms.find(normal_form.second) == normal_forms.end())
                {
                    normal_forms.insert(normal_form.second);
                    symbol coefficient("w_" + to_string(counter++));
                    for (auto& term : kgs)
                        term.first *= coefficient;
                    potential_primitive += kgs;
                }
            }
            KontsevichGraphSum<ex> coboundary = schouten_bracket(KontsevichGraphSum<ex> { {1, KontsevichGraph(1, 2, { {0, 1} }) } }, potential_primitive);
            coboundary.reduce_mod_skew();
            graph_series_copy[n] -= coboundary;
            potential_primitives[n] += potential_primitive;
        }
    }
    cout << "Reducing...\n";
    cout.flush();
    graph_series_copy.reduce_mod_skew();
    for (size_t n = 0; n <= graph_series_copy.precision(); ++n)
    {
        if (graph_series_copy[n] != 0 || n == graph_series.precision())
            cout << "h^" << n << ":\n";
        for (auto& term : graph_series_copy[n])
        {
            cout << term.second.encoding() << "    " << term.first << "\n";
        }
    }
    cout << "Potential primitive:\n";
    for (size_t n = 0; n <= potential_primitives.precision(); ++n)
    {
        if (potential_primitives[n] != 0 || n == potential_primitives.precision())
            cout << "h^" << n << ":\n";
        for (auto& term : potential_primitives[n])
        {
            cout << term.second.encoding() << "    " << term.first << "\n";
        }
    }
}
