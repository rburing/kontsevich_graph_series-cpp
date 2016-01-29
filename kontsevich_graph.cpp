#include "kontsevich_graph.hpp"
#include "util/sort_pairs.hpp"
#include <algorithm>

KontsevichGraph::KontsevichGraph(size_t internal, size_t external, std::vector< std::pair<size_t, size_t> > targets, int sign, bool normalized)
: d_internal(internal), d_external(external), d_targets(targets)
{
    if (!normalized)
    {
        size_t exchanges = sort_pairs(d_targets);
        sort(d_targets.begin(), d_targets.end());
        d_sign = sign * (exchanges % 2 == 0) ? 1 : -1;
    }
    else
        d_sign = sign;
}

std::vector<size_t> KontsevichGraph::internal_vertices() const
{
    std::vector<size_t> vertices(d_internal);
    std::iota(vertices.begin(), vertices.end(), d_external);
    return vertices;
}

std::vector< std::pair<size_t, size_t> > KontsevichGraph::targets() const
{
    return d_targets;
}

std::pair<size_t, size_t> KontsevichGraph::targets(size_t internal_vertex) const
{
    return d_targets[internal_vertex - d_external];
}

int KontsevichGraph::sign() const
{
    return d_sign;
}

int KontsevichGraph::sign(int new_sign)
{
    return d_sign = new_sign;
}

std::pair< size_t, std::vector< std::pair<size_t, size_t> > > KontsevichGraph::abs() const
{
    return { d_external, d_targets };
}

size_t KontsevichGraph::internal() const
{
    return d_internal;
}

size_t KontsevichGraph::external() const
{
    return d_external;
}

std::vector<size_t> KontsevichGraph::in_degrees() const
{
    std::vector<size_t> indegrees(d_external);
    for (auto& target_pair : d_targets)
    {
        if (target_pair.first < d_external)
            indegrees[target_pair.first]++;
        if (target_pair.second < d_external)
            indegrees[target_pair.second]++;
    }
    return indegrees;
}

std::vector<size_t> KontsevichGraph::neighbors_in(size_t vertex) const
{
    std::vector<size_t> neighbors;
    for (size_t idx = 0; idx < d_internal; ++idx)
    {
        if (d_targets[idx].first == vertex || d_targets[idx].second == vertex)
            neighbors.push_back(d_external + idx);
    }
    return neighbors;
}

bool operator==(const KontsevichGraph &lhs, const KontsevichGraph &rhs)
{
    return (lhs.d_external == rhs.d_external) && (lhs.d_sign == rhs.d_sign) && (lhs.d_targets == rhs.d_targets);
}

bool operator!=(const KontsevichGraph &lhs, const KontsevichGraph &rhs)
{
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream &os, const KontsevichGraph& g)
{
    return os << "Kontsevich graph with " << g.d_internal << " vertices on " << g.d_external << " ground vertices";
}

