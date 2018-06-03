#include "leibniz_graph.hpp"
#include <sstream>

LeibnizGraph::LeibnizGraph(KontsevichGraph graph, std::vector<KontsevichGraph::VertexPair> jacobiators, bool skew_leibniz)
: std::pair<KontsevichGraph, std::vector<KontsevichGraph::VertexPair> >::pair(graph, jacobiators), skew(skew_leibniz)
{
}

std::string LeibnizGraph::encoding() const
{
    std::stringstream ss;
    ss << this->second.size();
    ss << "   ";
    ss << this->first.encoding();
    ss << "   ";
    for (KontsevichGraph::VertexPair jacobiator : this->second)
        ss << jacobiator.first << " " << jacobiator.second;
    return ss.str();
}

std::istream& operator>>(std::istream& is, LeibnizGraph& g)
{
    size_t jacobiators;
    is >> jacobiators;
    is >> g.first;
    g.second.clear();
    std::pair<size_t, size_t> jacobiator;
    size_t jacobiator_count = 0;
    while (jacobiator_count++ < jacobiators && is >> jacobiator.first >> jacobiator.second)
        g.second.push_back(jacobiator);
    return is;
}

template<class T>
std::map<LeibnizGraph, T> LeibnizGraph::map_from_istream(std::istream& is, std::function<T(std::string)> const& parser)
{
    std::map<LeibnizGraph, T> result;
    if (parser == nullptr)
        return result;
    for (std::string line; getline(is, line);)
    {
        if (line.length() == 0 || line[0] == '#') // also skip comments
            continue;
        std::stringstream ss(line);
        LeibnizGraph g;
        ss >> g;
        std::string coefficient_str;
        ss >> coefficient_str;
        T coefficient = parser(coefficient_str);
        result[g] = coefficient;
    }
    return result;
}

std::set<LeibnizGraph> LeibnizGraph::those_yielding_kontsevich_graph(KontsevichGraph& graph, bool skew_leibniz)
{
    // TODO: multiple Jacobiators
    std::set<LeibnizGraph> leibniz_graphs;
    size_t external = graph.external();
    for (KontsevichGraph::Vertex v : graph.internal_vertices())
    {
        for (KontsevichGraph::Vertex w : graph.neighbors_in(v))
        {
            std::vector<KontsevichGraph::VertexPair> targets = graph.targets();
            KontsevichGraph::VertexPair& target_pair_v = targets[(size_t)v - external];
            // Check that there is no loop between v and w
            if (target_pair_v.first == w || target_pair_v.second == w)
                continue;
            KontsevichGraph::VertexPair& target_pair_w = targets[(size_t)w - external];
            // Check that the "Jacobiator" consisting of v and w falls on 3 distinct targets
            KontsevichGraph::Vertex& a = target_pair_v.first;
            KontsevichGraph::Vertex& b = target_pair_v.second;
            KontsevichGraph::Vertex& c = (target_pair_w.first == v) ? target_pair_w.second : target_pair_w.first;
            if (c == a || c == b)
                continue;
            leibniz_graphs.insert(LeibnizGraph(graph, { { v, w } }, skew_leibniz));
        }
    }
    return leibniz_graphs;
}
