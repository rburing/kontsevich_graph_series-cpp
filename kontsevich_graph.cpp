#include "kontsevich_graph.hpp"
#include <algorithm>

inline std::pair<size_t, size_t> exchange_pair(std::pair<size_t, size_t> p)
{
    return { p.second, p.first };
}

KontsevichGraph::KontsevichGraph(size_t internal, size_t external, std::vector< std::pair<size_t, size_t> > targets)
: d_internal(internal), d_external(external), d_targets(targets), d_targets_normalized(targets)
{
    size_t exchanges = 0;
    for (auto target = d_targets_normalized.begin(); target != d_targets_normalized.end(); target++)
    {
        std::pair<size_t, size_t> exchanged = exchange_pair(*target);
        if (exchanged < *target)
        {
            *target = exchanged;
            exchanges++;
        }
    }
    sort(d_targets_normalized.begin(), d_targets_normalized.end());
    d_sign = (exchanges % 2 == 0) ? 1 : -1;
}

std::vector<size_t> KontsevichGraph::internal_vertices() const
{
    std::vector<size_t> vertices(d_internal);
    std::iota(vertices.begin(), vertices.end(), d_external);
    return vertices;
}

std::pair<size_t, size_t> KontsevichGraph::targets(size_t internal_vertex) const
{
    return d_targets[internal_vertex - d_external];
}

int KontsevichGraph::sign() const
{
    return d_sign;
}

bool operator==(const KontsevichGraph &lhs, const KontsevichGraph &rhs)
{
    return (lhs.d_sign == rhs.d_sign) && (lhs.d_targets_normalized == rhs.d_targets_normalized);
}

bool operator!=(const KontsevichGraph &lhs, const KontsevichGraph &rhs)
{
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream &os, const KontsevichGraph& g)
{
    return os << "Kontsevich graph with " << g.d_internal << " vertices on " << g.d_external << " ground vertices";
}

