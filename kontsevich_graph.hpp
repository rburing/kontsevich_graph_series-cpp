#ifndef INCLUDED_KONTSEVICH_GRAPH_H_
#define INCLUDED_KONTSEVICH_GRAPH_H_

#include <vector>
#include <utility>
#include <cstddef>
#include <iostream>
#include <set>
#include <functional>
#include "util/sort_pairs.hpp"

class KontsevichGraph
{
    protected:
    size_t d_internal = 0;
    size_t d_external = 0;
    std::vector< std::pair<char, char> > d_targets;
    int d_sign = 1;

    public:
    typedef char Vertex;
    typedef std::pair<Vertex, Vertex> VertexPair;

    KontsevichGraph() = default;
    KontsevichGraph(size_t internal, size_t external, std::vector<VertexPair> targets, int sign = 1, bool normalized = false);
    std::vector<VertexPair> targets() const;
    VertexPair targets(Vertex internal_vertex) const;
    int sign() const;
    int sign(int new_sign);
    size_t internal() const;
    size_t external() const;
    size_t vertices() const;
    std::vector<Vertex> internal_vertices() const;
    std::pair< size_t, std::vector<VertexPair> > abs() const;
    size_t multiplicity() const;
    size_t in_degree(KontsevichGraph::Vertex vertex) const;
    std::vector<size_t> in_degrees() const;
    std::vector<Vertex> neighbors_in(Vertex vertex) const;
    KontsevichGraph mirror_image() const;
    std::string as_sage_expression() const;
    std::string encoding() const;
    std::vector< std::tuple<KontsevichGraph, int, int> > permutations() const;
    void normalize();
    KontsevichGraph& operator*=(const KontsevichGraph& rhs);
    bool operator<(const KontsevichGraph& rhs) const;
    bool is_zero() const;
    bool is_prime() const;
    bool positive_differential_order() const;
    bool has_cycles() const;
    bool has_tadpoles() const;
    bool has_multiple_edges() const;
    bool has_max_internal_indegree(size_t max_indegree) const;

    static std::set<KontsevichGraph> graphs(size_t internal, size_t external = 2, bool modulo_signs = false, bool modulo_mirror_images = false, std::function<void(KontsevichGraph&)> const& callback = nullptr, std::function<bool(KontsevichGraph&)> const& filter = nullptr);

    private:
    friend std::ostream& operator<<(std::ostream &os, const KontsevichGraph& g);
    friend std::istream& operator>>(std::istream& is, KontsevichGraph& g);
    friend bool operator==(const KontsevichGraph &lhs, const KontsevichGraph& rhs);
    friend bool operator!=(const KontsevichGraph &lhs, const KontsevichGraph& rhs);
};

inline size_t apply_permutation(size_t internal, size_t external, std::vector<KontsevichGraph::VertexPair>& targets, std::vector<KontsevichGraph::Vertex>& permutation)
{
    // Relabel elements of target pairs
    for (size_t i = 0; i != internal; ++i) {
        targets[i].first = permutation[targets[i].first];
        targets[i].second = permutation[targets[i].second];
    }
    // Apply permutation to list of target pairs
    std::vector<KontsevichGraph::VertexPair> permuted(targets.size());
    for (size_t i = 0; i != internal; ++i)
    {
        permuted[permutation[external + i] - external] = targets[i];
    }
    targets.swap(permuted);
    // Sort elements of target pairs
    return sort_pairs(targets.begin(), targets.end());
}

KontsevichGraph operator*(KontsevichGraph lhs, const KontsevichGraph& rhs);

std::ostream& operator<<(std::ostream &os, const KontsevichGraph::Vertex v);

#endif
