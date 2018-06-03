#include "leibniz_graph.hpp"
#include <sstream>
#include <tuple>

LeibnizGraph::LeibnizGraph(KontsevichGraph graph, std::vector<KontsevichGraph::VertexPair> jacobiators, bool skew)
: KontsevichGraph(graph), d_jacobiators(jacobiators), d_skew(skew)
{
    std::map<KontsevichGraph::Vertex, size_t> which_jacobiator;
    for (size_t j = 0; j != d_jacobiators.size(); ++j)
    {
        which_jacobiator[d_jacobiators[j].first]  = j;
        which_jacobiator[d_jacobiators[j].second] = j;
    }

    // Start building the sets of references to Leibniz targets (incoming edges on Jacobiator vertices)
    d_leibniz_targets.resize(jacobiators.size());
    for (KontsevichGraph::VertexPair& target_pair : d_targets)
    {
        auto it1 = which_jacobiator.find(target_pair.first);
        if (it1 != which_jacobiator.end())
            d_leibniz_targets[it1->second].insert(&target_pair.first);
        auto it2 = which_jacobiator.find(target_pair.second);
        if (it2 != which_jacobiator.end())
            d_leibniz_targets[it2->second].insert(&target_pair.second);
    }

    // Build the sets of three Jacobiator targets each
    d_jacobiator_targets.resize(jacobiators.size());
    d_max_jac_indegree = 0;
    size_t j = 0;
    for (KontsevichGraph::VertexPair& jacobiator : d_jacobiators)
    {
        KontsevichGraph::Vertex v = jacobiator.first;
        KontsevichGraph::Vertex w = jacobiator.second;
        KontsevichGraph::VertexPair& target_pair_v = d_targets[(size_t)v - d_external];
        KontsevichGraph::VertexPair& target_pair_w = d_targets[(size_t)w - d_external];
        KontsevichGraph::Vertex* a = &target_pair_v.first;
        KontsevichGraph::Vertex* b = &target_pair_v.second;
        KontsevichGraph::Vertex* c = (target_pair_w.first == v) ? &target_pair_w.second : &target_pair_w.first;
        // Remove internal Jacobiator edge from Leibniz targets
        d_leibniz_targets[j].erase(target_pair_w.first == v ? &target_pair_w.first : &target_pair_w.second);
        // Update maximum Jacobiator indegree
        size_t jac_indegree = d_leibniz_targets[j].size();
        if (jac_indegree > d_max_jac_indegree)
            d_max_jac_indegree = jac_indegree;
        // Set Jacobiator targets
        d_jacobiator_targets[j++] = { a, b, c };
    }
}

bool LeibnizGraph::operator<(const LeibnizGraph& rhs) const
{
    return std::tie(this->d_skew, this->d_external, this->d_internal, this->d_targets, this->d_jacobiators, this->d_sign) < \
           std::tie(rhs.d_skew, rhs.d_external, rhs.d_internal, rhs.d_targets, rhs.d_jacobiators, rhs.d_sign);
}

size_t LeibnizGraph::max_jac_indegree() const
{
    return d_max_jac_indegree;
}

std::string LeibnizGraph::encoding() const
{
    std::stringstream ss;
    ss << this->d_jacobiators.size();
    ss << "   ";
    ss << KontsevichGraph::encoding();
    ss << "   ";
    for (KontsevichGraph::VertexPair jacobiator : this->d_jacobiators)
        ss << jacobiator.first << " " << jacobiator.second;
    return ss.str();
}

std::istream& operator>>(std::istream& is, LeibnizGraph& g)
{
    size_t jacobiators;
    is >> jacobiators;
    is >> (KontsevichGraph&)g;
    g.d_jacobiators.clear();
    std::pair<size_t, size_t> jacobiator;
    size_t jacobiator_count = 0;
    while (jacobiator_count++ < jacobiators && is >> jacobiator.first >> jacobiator.second)
        g.d_jacobiators.push_back(jacobiator);
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
    std::vector<KontsevichGraph::VertexPair> targets = graph.targets();
    for (KontsevichGraph::Vertex v : graph.internal_vertices())
    {
        for (KontsevichGraph::Vertex w : graph.neighbors_in(v))
        {
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
