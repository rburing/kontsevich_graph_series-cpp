#include "kontsevich_graph.hpp"
#include "util/sort_pairs.hpp"
#include "util/cartesian_product.hpp"
#include "util/factorial.hpp"
#include <algorithm>
#include <tuple>
#include <stack>
#include <sstream>
#include <cmath>

KontsevichGraph::KontsevichGraph(size_t internal, size_t external, std::vector<KontsevichGraph::VertexPair> targets, int sign, bool normalized)
: d_internal(internal), d_external(external), d_targets(targets), d_sign(sign)
{
    if (!normalized)
        normalize();
}

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

void KontsevichGraph::normalize()
{
    std::vector<KontsevichGraph::VertexPair> global_minimum = d_targets;
    size_t exchanges = sort_pairs(global_minimum.begin(), global_minimum.end());
    std::vector<KontsevichGraph::Vertex> vertices(d_external + d_internal);
    std::iota(vertices.begin(), vertices.end(), 0);
    while (std::next_permutation(vertices.begin() + d_external, vertices.end()))
    {
        std::vector<KontsevichGraph::VertexPair> local_minimum = d_targets;
        size_t local_exchanges = apply_permutation(d_internal, d_external, local_minimum, vertices);
        if (local_exchanges % 2 == 1 && local_minimum == d_targets)
            d_sign = 0;
        if (local_minimum < global_minimum) {
            global_minimum = local_minimum;
            exchanges = local_exchanges;
        }
    }
    d_targets = global_minimum;
    d_sign *= (exchanges % 2 == 0) ? 1 : -1;
}

std::vector<KontsevichGraph::Vertex> KontsevichGraph::internal_vertices() const
{
    std::vector<KontsevichGraph::Vertex> vertices(d_internal);
    std::iota(vertices.begin(), vertices.end(), d_external);
    return vertices;
}

std::vector<KontsevichGraph::VertexPair> KontsevichGraph::targets() const
{
    return d_targets;
}

KontsevichGraph::VertexPair KontsevichGraph::targets(KontsevichGraph::Vertex internal_vertex) const
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

std::pair< size_t, std::vector<KontsevichGraph::VertexPair> > KontsevichGraph::abs() const
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

size_t KontsevichGraph::multiplicity() const
{
    std::set< std::vector<KontsevichGraph::VertexPair> > seen;
    std::vector<KontsevichGraph::Vertex> vertices(d_external + d_internal);
    std::iota(vertices.begin(), vertices.end(), 0);
    do
    {
        std::vector<KontsevichGraph::VertexPair> permuted = d_targets;
        apply_permutation(d_internal, d_external, permuted, vertices); // apply internal vertex permutation
        sort_pairs(permuted.begin(), permuted.end()); // ignore edge labeling
        if (seen.find(permuted) == seen.end())
            seen.insert(permuted);
    }
    while (std::next_permutation(vertices.begin() + d_external, vertices.end()));
    return pow(2, d_internal) * seen.size();
}

bool KontsevichGraph::is_zero() const
{
    std::vector<KontsevichGraph::Vertex> vertices(d_external + d_internal);
    std::iota(vertices.begin(), vertices.end(), 0);
    while (std::next_permutation(vertices.begin() + d_external, vertices.end()))
    {
        std::vector<KontsevichGraph::VertexPair> permuted = d_targets;
        size_t exchanges = apply_permutation(d_internal, d_external, permuted, vertices);
        if (permuted == d_targets && exchanges % 2 == 1)
            return true;
    }
    return false;
}

size_t KontsevichGraph::in_degree(KontsevichGraph::Vertex vertex) const
{
    size_t count = 0;
    for (auto& target_pair : d_targets)
    {
        if (target_pair.first == vertex)
            ++count;
        if (target_pair.second == vertex)
            ++count;
    }
    return count;
}

std::vector<size_t> KontsevichGraph::in_degrees() const
{
    std::vector<size_t> indegrees(d_external);
    for (auto& target_pair : d_targets)
    {
        if ((size_t)target_pair.first < d_external)
            indegrees[target_pair.first]++;
        if ((size_t)target_pair.second < d_external)
            indegrees[target_pair.second]++;
    }
    return indegrees;
}

std::vector<KontsevichGraph::Vertex> KontsevichGraph::neighbors_in(KontsevichGraph::Vertex vertex) const
{
    std::vector<KontsevichGraph::Vertex> neighbors;
    for (size_t idx = 0; idx < d_internal; ++idx)
    {
        if (d_targets[idx].first == vertex || d_targets[idx].second == vertex)
            neighbors.push_back(d_external + idx);
    }
    return neighbors;
}

bool KontsevichGraph::operator<(const KontsevichGraph& rhs) const
{
    return std::tie(this->d_external, this->d_internal, this->d_targets, this->d_sign) < std::tie(rhs.d_external, rhs.d_internal, rhs.d_targets, rhs.d_sign);
}

KontsevichGraph& KontsevichGraph::operator*=(const KontsevichGraph& rhs)
{
    // TODO: maybe check if d_external == rhs.d_external
    d_targets.reserve(d_targets.size() + rhs.d_targets.size());
    // Concatenate lists of targets
    d_targets.insert(d_targets.end(), rhs.d_targets.begin(), rhs.d_targets.end());
    // Add offsets to RHS' internal targets
    for (size_t i = 0; i != rhs.d_internal; ++i)
    {
        if ((size_t)d_targets[d_internal + i].first >= d_external)
            d_targets[d_internal + i].first += d_internal;
        if ((size_t)d_targets[d_internal + i].second >= d_external)
            d_targets[d_internal + i].second += d_internal;
    }
    d_internal += rhs.d_internal;
    d_sign *= rhs.d_sign;
    normalize();
    return *this;
}

KontsevichGraph operator*(KontsevichGraph lhs, const KontsevichGraph& rhs)
{
    lhs *= rhs;
    return lhs;
}

bool KontsevichGraph::is_prime() const
{
    size_t vertex_count = 0;
    std::set<KontsevichGraph::Vertex> seen;
    // Choose a vertex connected to the ground, if possible
    KontsevichGraph::Vertex start = 0;
    for (size_t i = 0; i != d_internal; ++i)
    {
        if ((size_t)d_targets[i].first < d_external || (size_t)d_targets[i].second < d_external)
        {
            start = d_external + i;
            break;
        }
    }
    // Otherwise, choose the first internal vertex
    if (start == 0)
        start = d_external;
    // Breadth-first search
    std::stack< std::pair<KontsevichGraph::Vertex, int> > vertex_stack;
    vertex_stack.push({ start, 0 });
    while (vertex_stack.size() > 0)
    {
        std::pair<KontsevichGraph::Vertex, int> vertex = vertex_stack.top();
        vertex_stack.pop();
        if (seen.find(vertex.first) == seen.end())
        {
            ++vertex_count;
            seen.insert(vertex.first);
            // Outgoing neighbors
            KontsevichGraph::VertexPair target_pair = d_targets[vertex.first - d_external];
            if ((size_t)target_pair.first >= d_external)
                if (seen.find(target_pair.first) == seen.end())
                    vertex_stack.push({ target_pair.first, vertex.second + 1});
            if ((size_t)target_pair.second >= d_external)
                if (seen.find(target_pair.second) == seen.end())
                    vertex_stack.push({ target_pair.second, vertex.second + 1});
            // Incoming neighbors
            for (auto neighbor : this->neighbors_in(vertex.first))
                if (seen.find(neighbor) == seen.end())
                    vertex_stack.push({ neighbor, vertex.second + 1});
        }
    }
    return vertex_count == d_internal;
}

KontsevichGraph KontsevichGraph::mirror_image() const
{
    std::vector<KontsevichGraph::VertexPair> targets = d_targets;
    // Reverse the ground vertices
    for (auto& target_pair : targets)
    {
        if ((size_t)target_pair.first < d_external)
            target_pair.first = d_external - 1 - target_pair.first;
        if ((size_t)target_pair.second < d_external)
            target_pair.second = d_external - 1 - target_pair.second;
    }
    return KontsevichGraph(d_internal, d_external, targets, d_sign);
}

bool KontsevichGraph::positive_differential_order() const
{
    std::set<KontsevichGraph::Vertex> seen;
    for (auto& target_pair : d_targets)
    {
        if ((size_t)target_pair.first < d_external)
            seen.insert(target_pair.first);
        if ((size_t)target_pair.second < d_external)
            seen.insert(target_pair.second);
    }
    return seen.size() == d_external;
}

std::string KontsevichGraph::as_sage_expression() const
{
    std::stringstream ss;
    ss << "KontsevichGraph([";
    for (size_t v : this->internal_vertices())
    {
        KontsevichGraph::VertexPair targets = this->targets(v);
        ss << "(" << v << ", " << (size_t)targets.first << ", 'L'), ";
        ss << "(" << v << ", " << (size_t)targets.second<< ", 'R')";
        if (v != this->vertices() - 1)
            ss << ", ";
    }
    ss << "], ground_vertices=(";
    for (size_t v = 0; v < d_external; ++v)
    {
        ss << v;
        if (v != d_external-1)
            ss << ", ";
    }
    ss << "), immutable=True).normalize_vertex_labels(inplace=False)";
    return ss.str();
}

std::set<KontsevichGraph> KontsevichGraph::graphs(size_t internal, size_t external, bool modulo_signs, bool modulo_mirror_images, std::function<bool(KontsevichGraph)> const& filter)
{
    std::set<KontsevichGraph> result;
    std::vector<size_t> ends(2*internal);
    for (size_t i = 0; i != 2*internal; ++i)
    {
        ends[i] = internal + external;
    }
    CartesianProduct graph_encodings(ends);
    std::vector<KontsevichGraph::VertexPair> targets(internal);
    for (auto graph_encoding = graph_encodings.begin(); graph_encoding != graph_encodings.end(); ++graph_encoding)
    {
        bool skip = false;
        for (size_t i = 0; i != internal; ++i)
        {
            KontsevichGraph::VertexPair target_pair = { (*graph_encoding)[2*i], (*graph_encoding)[2*i + 1] };
            // Avoid double edges and tadpoles:
            if (target_pair.first == target_pair.second || (size_t)target_pair.first == external + i || (size_t)target_pair.second == external + i)
            {
                skip = true;
                break;
            }
            targets[i] = target_pair;
        }
        if (!skip)
        {
            KontsevichGraph graph(internal, external, targets);
            if (filter && !filter(graph))
                continue;
            if (modulo_mirror_images)
            {
                KontsevichGraph mirror = graph.mirror_image();
                if (mirror < graph)
                    graph = mirror;
            }
            if (modulo_signs)
                graph.sign(1);
            result.insert(graph);
        }
    }
    return result;
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

std::istream& operator>>(std::istream& is, KontsevichGraph& g)
{
    is >> g.d_external;
    is >> g.d_internal;
    is >> g.d_sign;
    g.d_targets.clear();
    std::string line;
    std::pair<size_t, size_t> target_pair;
    size_t pair_count = 0;
    while (pair_count++ < g.d_internal && is >> target_pair.first >> target_pair.second)
        g.d_targets.push_back(target_pair);
    g.d_internal = g.d_targets.size();
    g.normalize();
    return is;
}

std::string KontsevichGraph::encoding() const
{
    std::stringstream ss;
    ss << d_external << " " << d_internal << " " << d_sign << "   ";
    for (auto& target_pair : d_targets)
    {
        ss << target_pair.first << " " << target_pair.second;
        if (d_targets.size() > 0 && &target_pair != &d_targets[d_targets.size()-1])
            ss << " ";
    }
    return ss.str();
}

bool KontsevichGraph::has_cycles() const
{
    for (KontsevichGraph::Vertex v : this->internal_vertices())
    {
        VertexPair targets = this->targets(v);
        for (KontsevichGraph::Vertex w : { targets.first, targets.second })
        {
            KontsevichGraph::VertexPair targets = this->targets(w);
            if (targets.first == v || targets.second == v)
                return true;
        }
    }
    return false;
}

bool KontsevichGraph::has_tadpoles() const
{
    for (KontsevichGraph::Vertex v : this->internal_vertices())
    {
        VertexPair targets = this->targets(v);
        if (targets.first == v || targets.second == v)
            return true;
    }
    return false;
}

bool KontsevichGraph::has_multiple_edges() const
{
    for (auto const& target_pair : d_targets)
        if (target_pair.first == target_pair.second)
            return true;
    return false;
}

std::ostream& operator<<(std::ostream &os, const KontsevichGraph::Vertex v)
{
    return os << (size_t)v;
}
