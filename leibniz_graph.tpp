#include "leibniz_graph.hpp"
#include "util/cartesian_product.hpp"
#include <sstream>
#include <tuple>

template<class T>
LeibnizGraph<T>::LeibnizGraph(KontsevichGraph graph, std::vector<KontsevichGraph::VertexPair> jacobiators, bool skew)
: KontsevichGraph(graph), d_jacobiators(jacobiators), d_skew(skew)
{
    set_jacobiator_and_leibniz_targets();
}

template<class T>
LeibnizGraph<T>::LeibnizGraph(const LeibnizGraph<T>& other)
: KontsevichGraph(other), d_jacobiators(other.d_jacobiators), d_skew(other.d_skew)
{
    set_jacobiator_and_leibniz_targets();
}

template<class T>
LeibnizGraph<T>& LeibnizGraph<T>::operator=(const LeibnizGraph<T>& other)
{
    KontsevichGraph::operator=(other);
    d_jacobiators = other.d_jacobiators;
    d_skew = other.d_skew;
    set_jacobiator_and_leibniz_targets();
    return *this;
}

template<class T>
void LeibnizGraph<T>::set_jacobiator_and_leibniz_targets()
{
    std::map<KontsevichGraph::Vertex, size_t> which_jacobiator;
    for (size_t j = 0; j != d_jacobiators.size(); ++j)
    {
        which_jacobiator[d_jacobiators[j].first]  = j;
        which_jacobiator[d_jacobiators[j].second] = j;
    }

    // Start building the sets of references to Leibniz targets (incoming edges on Jacobiator vertices)
    d_leibniz_targets.clear();
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
    d_jacobiator_targets.clear();
    d_jacobiator_targets.resize(d_jacobiators.size());
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

template<class T>
bool LeibnizGraph<T>::skew() const
{
    return d_skew;
}

template<class T>
bool LeibnizGraph<T>::skew(bool new_skew)
{
    d_skew = new_skew;
    return d_skew;
}

template<class T>
bool LeibnizGraph<T>::operator<(const LeibnizGraph<T>& rhs) const
{
    return std::tie(this->d_skew, this->d_external, this->d_internal, this->d_targets, this->d_jacobiators, this->d_sign) < \
           std::tie(rhs.d_skew, rhs.d_external, rhs.d_internal, rhs.d_targets, rhs.d_jacobiators, rhs.d_sign);
}

template<class T>
size_t LeibnizGraph<T>::max_jac_indegree() const
{
    return d_max_jac_indegree;
}

template<class T>
std::string LeibnizGraph<T>::encoding() const
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

template<class T>
std::istream& operator>>(std::istream& is, LeibnizGraph<T>& g)
{
    size_t jacobiators;
    is >> jacobiators;
    is >> (KontsevichGraph&)g;
    g.d_jacobiators.clear();
    std::pair<size_t, size_t> jacobiator;
    size_t jacobiator_count = 0;
    while (jacobiator_count++ < jacobiators && is >> jacobiator.first >> jacobiator.second)
        g.d_jacobiators.push_back(jacobiator);
    g.set_jacobiator_and_leibniz_targets();
    return is;
}

template<class T>
std::map< LeibnizGraph<T>, T> LeibnizGraph<T>::map_from_istream(std::istream& is, std::function<T(std::string)> const& parser)
{
    std::map< LeibnizGraph<T>, T> result;
    if (parser == nullptr)
        return result;
    for (std::string line; getline(is, line);)
    {
        if (line.length() == 0 || line[0] == '#') // also skip comments
            continue;
        std::stringstream ss(line);
        LeibnizGraph<T> g;
        ss >> g;
        std::string coefficient_str;
        ss >> coefficient_str;
        T coefficient = parser(coefficient_str);
        result[g] = coefficient;
    }
    return result;
}

template<class T>
std::set< LeibnizGraph<T> > LeibnizGraph<T>::those_yielding_kontsevich_graph(KontsevichGraph& graph, bool skew_leibniz)
{
    // TODO: multiple Jacobiators
    std::set< LeibnizGraph<T> > leibniz_graphs;
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
            leibniz_graphs.insert(LeibnizGraph<T>(graph, { { v, w } }, skew_leibniz));
        }
    }
    return leibniz_graphs;
}

template<class T>
void LeibnizGraph<T>::normalize()
{
    // Normal form of Leibniz graph (three permutations of each Jacobiator, take minimal encoding, remember where Jacobiators are)

    // TODO: remember sign?
    // TODO: save partial expansion?

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

    std::vector< std::pair< std::vector<KontsevichGraph::VertexPair>, std::vector<KontsevichGraph::VertexPair> > > leibniz_graphs;

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

            std::vector<KontsevichGraph::VertexPair> d_targets_template = d_targets;

            // Permute ground vertices (if skew)
            for (KontsevichGraph::VertexPair& target_pair : d_targets_template)
            {
                if ((size_t)target_pair.first < d_external)
                    target_pair.first = ground_vertices[target_pair.first];
                if ((size_t)target_pair.second < d_external)
                    target_pair.second = ground_vertices[target_pair.second];
            }

            // Find permutation of vertex labels such that the list of targets is minimal with respect to the defined ordering
            std::vector<KontsevichGraph::VertexPair> global_minimum = d_targets_template;

            sort_pairs(global_minimum.begin(), global_minimum.end());

            std::vector<KontsevichGraph::VertexPair> new_jacobiators = d_jacobiators;

            std::vector<KontsevichGraph::Vertex> vertices(d_external + d_internal);
            std::iota(vertices.begin(), vertices.end(), 0);

            while (std::next_permutation(vertices.begin() + d_external, vertices.end()))
            {
                std::vector<KontsevichGraph::VertexPair> local_minimum = d_targets_template;
                apply_permutation(d_internal, d_external, local_minimum, vertices);
                if (local_minimum < global_minimum)
                {
                    global_minimum = local_minimum;
                    // Find where Jacobiators are
                    new_jacobiators = d_jacobiators;
                    for (KontsevichGraph::VertexPair& new_jacobiator : new_jacobiators)
                        new_jacobiator = { vertices[(size_t)new_jacobiator.first], vertices[(size_t)new_jacobiator.second] };
                }
            }

            leibniz_graphs.push_back({ global_minimum, new_jacobiators });
        }
    } while (d_skew && std::next_permutation(ground_vertices.begin(), ground_vertices.end()));
    auto& leibniz_normal_form = *min_element(leibniz_graphs.begin(), leibniz_graphs.end());
    d_targets = leibniz_normal_form.first;
    d_jacobiators = leibniz_normal_form.second;
    set_jacobiator_and_leibniz_targets();
}

template<class T>
KontsevichGraphSum<T> LeibnizGraph<T>::expansion(T prefactor)
{
    // Back-up d_targets:
    std::vector<KontsevichGraph::VertexPair> d_targets_values = d_targets;

    KontsevichGraphSum<T> graph_sum;

    // Fix some ordering of Jacobiator arguments (as a vector, instead of a set)
    std::vector< std::vector<KontsevichGraph::Vertex> > jacobiator_arguments(d_jacobiators.size());
    for (size_t j = 0; j != d_jacobiators.size(); ++j)
    {
        jacobiator_arguments[j].resize(3);
        size_t k = 0;
        for (auto it = d_jacobiator_targets[j].begin(); it != d_jacobiator_targets[j].end(); ++it)
            jacobiator_arguments[j][k++] = **it;
    }

    std::vector<size_t> leibniz_sizes(d_leibniz_targets.size(), 2);
    CartesianProduct leibniz_indices(leibniz_sizes);
    for (CartesianProduct leibniz_index = leibniz_indices.begin(); leibniz_index != leibniz_indices.end(); ++leibniz_index)
    {
        // Leibniz rule
        size_t idx = 0;
        for (auto& leibniz_target : d_leibniz_targets)
            *(leibniz_target.first) = (*leibniz_index)[idx++] == 0 ? d_jacobiators[leibniz_target.second].first : d_jacobiators[leibniz_target.second].second;

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

            KontsevichGraph new_graph(d_internal, d_external, d_targets, d_sign);
            graph_sum += KontsevichGraphSum<T>({ { prefactor, new_graph } });
        }
    }

    if (d_skew)
        graph_sum = graph_sum.skew_symmetrization();

    graph_sum.reduce_mod_skew();

    // Restore back-up of d_targets:
    for (size_t j = 0; j != d_internal; ++j)
        d_targets[j] = d_targets_values[j];

    return graph_sum;
}
