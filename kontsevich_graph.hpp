#ifndef INCLUDED_KONTSEVICH_GRAPH_H_
#define INCLUDED_KONTSEVICH_GRAPH_H_

#include <vector>
#include <utility>
#include <cstddef>
#include <iostream>
#include <set>
#include <functional>

class KontsevichGraph
{
    size_t d_internal = 0;
    size_t d_external = 0;
    std::vector< std::pair<char, char> > d_targets;
    int d_sign = 1;

    public:
    typedef char Vertex;
    typedef std::pair<Vertex, Vertex> VertexPair;

    KontsevichGraph() = default;
    KontsevichGraph(size_t internal, size_t external, std::vector<VertexPair> targets, int sign = 1, bool normalized = false);
    void normalize();
    std::vector<Vertex> internal_vertices() const;
    std::vector<VertexPair> targets() const;
    VertexPair targets(Vertex internal_vertex) const;
    int sign() const;
    int sign(int new_sign);
    std::pair< size_t, std::vector<VertexPair> > abs() const;
    size_t internal() const;
    size_t external() const;
    size_t vertices() const;
    size_t multiplicity() const;
    bool is_zero() const;
    std::vector<size_t> in_degrees() const;
    std::vector<Vertex> neighbors_in(Vertex vertex) const;
    bool operator<(const KontsevichGraph& rhs) const;
    KontsevichGraph& operator*=(const KontsevichGraph& rhs);
    bool is_prime() const;
    KontsevichGraph mirror_image() const;
    bool positive_differential_order() const;

    static std::set<KontsevichGraph> graphs(size_t internal, size_t external = 2, bool modulo_signs = false, bool modulo_mirror_images = false, std::function<bool(KontsevichGraph)> const& filter = nullptr);

    private:
    friend std::ostream& operator<<(std::ostream &os, const KontsevichGraph& g);
    friend std::istream& operator>>(std::istream& is, KontsevichGraph& g);
    friend bool operator==(const KontsevichGraph &lhs, const KontsevichGraph& rhs);
    friend bool operator!=(const KontsevichGraph &lhs, const KontsevichGraph& rhs);
};

KontsevichGraph operator*(KontsevichGraph lhs, const KontsevichGraph& rhs);

std::ostream& operator<<(std::ostream &os, const KontsevichGraph::Vertex v);

#endif
