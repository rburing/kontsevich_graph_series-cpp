#include "leibniz_graph.hpp"
#include "util/cartesian_product.hpp"
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
    std::map<size_t, size_t> jac_indegree;
    for (KontsevichGraph::VertexPair& target_pair : d_targets)
    {
        auto it1 = which_jacobiator.find(target_pair.first);
        if (it1 != which_jacobiator.end())
        {
            d_leibniz_targets[&target_pair.first] = it1->second;
            ++jac_indegree[it1->second];
        }
        auto it2 = which_jacobiator.find(target_pair.second);
        if (it2 != which_jacobiator.end())
        {
            d_leibniz_targets[&target_pair.second] = it2->second;
            ++jac_indegree[it2->second];
        }
    }
    d_max_jac_indegree = 0;
    for (auto indegree : jac_indegree)
        if (indegree.second - 1 > d_max_jac_indegree)
            d_max_jac_indegree = indegree.second - 1;

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
        d_leibniz_targets.erase(target_pair_w.first == v ? &target_pair_w.first : &target_pair_w.second);
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

void LeibnizGraph::normalize()
{
    // Normal form of Leibniz graph (three permutations of each Jacobiator, take minimal encoding, remember where Jacobiators are)

    // TODO: remember sign?
    // TODO: save partial expansion?
    // TODO: construct LeibnizGraph (call the constructor) only once; let the intermediates be tuples or something

    // Set Leibniz targets to "bottom" vertex in Jacobiator, i.e. v in { v, w } (the Jacobiator edge is v <-- w)
    for (auto& leibniz_target : d_leibniz_targets)
        *leibniz_target.first = d_jacobiators[leibniz_target.second].first;

    // Fix some ordering of Jacobiator arguments (as a vector, instead of a set)
    std::vector< std::vector<KontsevichGraph::Vertex> > jacobiator_arguments(d_jacobiators.size());
    for (size_t j = 0; j != d_jacobiators.size(); ++j)
    {
        jacobiator_arguments[j].resize(3);
        size_t k = 0;
        for (auto it = d_jacobiator_targets[j].begin(); it != d_jacobiator_targets[j].end(); ++it)
            jacobiator_arguments[j][k++] = **it;
    }

    std::vector<LeibnizGraph> leibniz_graphs;

    std::vector<KontsevichGraph::Vertex> ground_vertices(d_external);
    std::iota(ground_vertices.begin(), ground_vertices.end(), 0);
    do
    {
        // Choose shifts (by 0, 1, or 2) in Jacobiator arguments
        std::vector<size_t> shifts_max(d_jacobiators.size(), 3);
        CartesianProduct all_shifts(shifts_max);

        for (CartesianProduct shifts = all_shifts.begin(); shifts != all_shifts.end(); ++shifts)
        {
            // Set Jacobiator arguments
            for (size_t j = 0; j != d_jacobiators.size(); ++j)
            {
                size_t k = 0;
                for (auto it = d_jacobiator_targets[j].begin(); it != d_jacobiator_targets[j].end(); ++it)
                    **it = jacobiator_arguments[j][(k++ + (*shifts)[j]) % 3];
            }

            // Permute ground vertices (if skew)
            for (KontsevichGraph::VertexPair& target_pair : d_targets)
            {
                if ((size_t)target_pair.first < d_external)
                    target_pair.first = ground_vertices[target_pair.first];
                if ((size_t)target_pair.second < d_external)
                    target_pair.second = ground_vertices[target_pair.second];
            }

            // Find permutation of vertex labels such that the list of targets is minimal with respect to the defined ordering
            std::vector<KontsevichGraph::VertexPair> global_minimum = d_targets;

            sort_pairs(global_minimum.begin(), global_minimum.end());

            std::vector<KontsevichGraph::VertexPair> new_jacobiators = d_jacobiators;

            std::vector<KontsevichGraph::Vertex> vertices(d_external + d_internal);
            std::iota(vertices.begin(), vertices.end(), 0);

            while (std::next_permutation(vertices.begin() + d_external, vertices.end()))
            {
                std::vector<KontsevichGraph::VertexPair> local_minimum = d_targets;
                apply_permutation(d_internal, d_external, local_minimum, vertices);
                if (local_minimum < global_minimum)
                {
                    global_minimum = local_minimum;
                    // Find where Jacobiators are
                    for (KontsevichGraph::VertexPair& new_jacobiator : new_jacobiators)
                        new_jacobiator = { vertices[(size_t)new_jacobiator.first], vertices[(size_t)new_jacobiator.second] };
                }
            }

            KontsevichGraph leibniz_graph(d_internal, d_external, global_minimum, 1, true);
            leibniz_graphs.push_back(LeibnizGraph(leibniz_graph, new_jacobiators, d_skew));
        }
    } while (d_skew && std::next_permutation(ground_vertices.begin(), ground_vertices.end()));
    LeibnizGraph leibniz_normal_form = *min_element(leibniz_graphs.begin(), leibniz_graphs.end());
    std::swap(leibniz_normal_form, *this);
}
