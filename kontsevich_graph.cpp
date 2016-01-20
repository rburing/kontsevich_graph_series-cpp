#include "kontsevich_graph.hpp"

KontsevichGraph::KontsevichGraph(size_t internal, size_t external, std::vector< std::pair<size_t, size_t> > targets)
: d_internal(internal), d_external(external), d_targets(targets)
{
}

std::ostream& operator<<(std::ostream &os, const KontsevichGraph& g)
{
    return os << "Kontsevich graph with " << g.d_internal << " vertices on " << g.d_external << " ground vertices";
}
