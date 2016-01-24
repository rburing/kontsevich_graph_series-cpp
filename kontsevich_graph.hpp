#ifndef INCLUDED_KONTSEVICH_GRAPH_H_
#define INCLUDED_KONTSEVICH_GRAPH_H_

#include <vector>
#include <utility>
#include <cstddef>
#include <ostream>

class KontsevichGraph
{
    size_t d_internal;
    size_t d_external;
    std::vector< std::pair<size_t, size_t> > d_targets;
    int d_sign;

    public:
    KontsevichGraph(size_t internal, size_t external, std::vector< std::pair<size_t, size_t> > targets);
    std::vector<size_t> internal_vertices() const;
    std::pair<size_t, size_t> targets(size_t internal_vertex) const;
    int sign() const;
    int sign(int new_sign);
    std::pair< size_t, std::vector< std::pair<size_t, size_t> > > abs() const;

    private:
    friend std::ostream& operator<<(std::ostream &os, const KontsevichGraph& g);
    friend bool operator==(const KontsevichGraph &lhs, const KontsevichGraph& rhs);
    friend bool operator!=(const KontsevichGraph &lhs, const KontsevichGraph& rhs);
};

#endif
