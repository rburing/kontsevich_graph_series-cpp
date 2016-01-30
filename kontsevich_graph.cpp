#include "kontsevich_graph.hpp"
#include "util/sort_pairs.hpp"
#include <algorithm>

KontsevichGraph::KontsevichGraph(size_t internal, size_t external, std::vector< std::pair<size_t, size_t> > targets, int sign, bool normalized)
: d_internal(internal), d_external(external), d_targets(targets), d_sign(sign)
{
    if (!normalized)
    {
        std::vector< std::pair<size_t, size_t> > global_minimum = d_targets;
        std::vector<size_t> vertices(d_external + d_internal);
        std::iota(vertices.begin(), vertices.end(), 0);
        size_t exchanges = 0;
        while (std::next_permutation(vertices.begin() + d_external, vertices.end()))
        {
            std::vector< std::pair<size_t, size_t> > local_minimum = d_targets;
            // Relabel elements of target pairs
            for (size_t i = 0; i != d_internal; ++i) {
                local_minimum[i].first = vertices[local_minimum[i].first];
                local_minimum[i].second = vertices[local_minimum[i].second];
            }
            // Reorder target pairs
            for (size_t i = 0; i != d_internal; ++i) {
                std::swap(local_minimum[i], local_minimum[vertices[d_external + i] - d_external]);
            }
            size_t local_exchanges = sort_pairs(local_minimum.begin(), local_minimum.end());
            if (local_minimum < global_minimum) {
                global_minimum = local_minimum;
                exchanges = local_exchanges;
            }
        }
        d_targets = global_minimum;
        d_sign *= (exchanges % 2 == 0) ? 1 : -1;
    }
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

size_t KontsevichGraph::vertices() const
{
    return d_internal + d_external;
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

