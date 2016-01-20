#include "kontsevich_graph.hpp"
#include <algorithm>

KontsevichGraph::KontsevichGraph(size_t internal, size_t external, std::vector< std::pair<size_t, size_t> > targets)
: d_internal(internal), d_external(external), d_targets(targets)
{
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

std::ostream& operator<<(std::ostream &os, const KontsevichGraph& g)
{
    return os << "Kontsevich graph with " << g.d_internal << " vertices on " << g.d_external << " ground vertices";
}
