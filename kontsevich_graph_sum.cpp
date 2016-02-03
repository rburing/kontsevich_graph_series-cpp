#include "kontsevich_graph_sum.hpp"
#include "util/cartesian_product.hpp"
#include "util/sort_pairs.hpp"
#include <algorithm>

template <class T>
void KontsevichGraphSum<T>::reduce()
{
    auto current_term = this->begin();
    while (current_term < this->end())
    {
        current_term->first *= current_term->second.sign();
        current_term->second.sign(1);
        auto subsequent_term = current_term + 1;
        while (subsequent_term < this->end())
        {
            if (subsequent_term->second.abs() == current_term->second.abs())
            {
                current_term->first += subsequent_term->first * subsequent_term->second.sign();
                subsequent_term = this->erase(subsequent_term);
            }
            else
                subsequent_term++;
        }
        if (current_term->first == 0)
            current_term = this->erase(current_term);
        else
            current_term++;
    }
}

template <class T>
std::ostream& operator<<(std::ostream& os, const std::pair<T, KontsevichGraph>& term)
{
    os << term.first * term.second.sign() << "*(" << term.second << ")";
    return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const KontsevichGraphSum<T>& gs)
{
    if (gs.size() == 0)
        return os << "0";
    for (auto& term : gs)
    {
        os << term;
        if (&term != &gs.back())
            os << " + ";
    }
    return os;
}

template <class T>
std::istream& operator>>(std::istream& is, KontsevichGraphSum<T>& sum)
{
    typename KontsevichGraphSum<T>::Term term;
    while (is >> term.first >> term.second) {
        sum.push_back(term);
    }
    return is;
}

template <class T>
bool KontsevichGraphSum<T>::operator==(const KontsevichGraphSum<T> &other)
{
    KontsevichGraphSum<T> difference = *this - other;
    difference.reduce();
    return difference.size() == 0;
}

template <class T>
bool KontsevichGraphSum<T>::operator!=(const KontsevichGraphSum<T> &other)
{
    return !(*this == other);
}

template <class T>
KontsevichGraphSum<T>& KontsevichGraphSum<T>::operator+=(const KontsevichGraphSum<T>& rhs)
{
    this->reserve(this->size() + rhs.size());
    this->insert(this->end(), rhs.begin(), rhs.end());
    return *this;
}

template <class T>
KontsevichGraphSum<T> operator+(KontsevichGraphSum<T> lhs, const KontsevichGraphSum<T> &rhs)
{
    lhs += rhs;
    return lhs;
}

template <class T>
KontsevichGraphSum<T>& KontsevichGraphSum<T>::operator-=(const KontsevichGraphSum<T>& rhs)
{
    size_t my_size = this->size();
    // Add the lists of terms
    *this += rhs;
    // Flip the appropriate signs
    for (auto it = this->begin() + my_size; it != this->end(); ++it)
        it->first = -it->first;
    return *this;
}

template <class T>
KontsevichGraphSum<T> operator-(KontsevichGraphSum<T> lhs, const KontsevichGraphSum<T>& rhs)
{
    lhs -= rhs;
    return lhs;
}

template <class T>
KontsevichGraphSum<T> KontsevichGraphSum<T>::operator()(std::vector< KontsevichGraphSum<T> > arguments)
{
    KontsevichGraphSum<T> total; // TODO: pre-compute size?

    for (auto& main_term : *this) // Linearity
    {
        // TODO: check arguments.size() equals main_term.d_external?

        std::vector<size_t> argument_sizes(arguments.size());
        for (size_t i = 0; i != arguments.size(); ++i)
            argument_sizes[i] = arguments[i].size();
        CartesianProduct multilinearity_indices(argument_sizes);

        for (auto arg_indices = multilinearity_indices.begin(); arg_indices != multilinearity_indices.end(); ++arg_indices) // Multi-linearity
        {
            // Now main_term acts on arguments indicated by arg_indices
            T coeff = main_term.first;
            size_t internal = main_term.second.internal(), external = 0;
            int sign = main_term.second.sign();
            for (size_t i = 0; i != arguments.size(); ++i)
            {
                coeff *= arguments[i][(*arg_indices)[i]].first;
                KontsevichGraph* graph = &arguments[i][(*arg_indices)[i]].second;
                internal += graph->internal();
                external += graph->external();
                sign *= graph->sign();
            }
            // Prepare to concatenate targets
            std::vector< std::pair<size_t, size_t> > new_targets(internal);
            auto new_targets_it = new_targets.begin();
            // Relabel argument targets and copy them
            size_t start_internal = external;
            size_t start_external = 0;
            std::vector<size_t> start_internal_vec(arguments.size());
            std::vector<size_t> start_external_vec(arguments.size());
            for (size_t i = 0; i != arguments.size(); ++i)
            {
                KontsevichGraphSum<T>::Term arg_term = arguments[i][(*arg_indices)[i]];
                std::vector< std::pair<size_t, size_t> > arg_targets = arg_term.second.targets();
                std::vector<size_t> offsets(arg_term.second.vertices());
                std::fill(offsets.begin(), offsets.begin() + arg_term.second.external(), start_external); 
                std::fill(offsets.begin() + arg_term.second.external(), offsets.begin() + arg_term.second.vertices(), start_internal - arg_term.second.external());
                for (auto& target_pair : arg_targets)
                {
                    target_pair.first += offsets[target_pair.first];
                    target_pair.second += offsets[target_pair.second];
                }
                std::copy(arg_targets.begin(), arg_targets.end(), new_targets_it);
                start_external_vec[i] = start_external;
                start_external += arg_term.second.external();
                start_internal_vec[i] = start_internal;
                start_internal += arg_term.second.internal();
                new_targets_it += arg_targets.size();
            }
            // Relabel main_term internal targets, but not the external ones
            std::vector< std::pair<size_t, size_t> > main_targets = main_term.second.targets();
            for (auto& target_pair : main_targets)
            {
                for (size_t* target : {&target_pair.first, &target_pair.second})
                    if (*target >= main_term.second.external()) // internal
                        *target += (start_internal - main_term.second.external());
            }
            std::copy(main_targets.begin(), main_targets.end(), new_targets_it);
            // Prepare to apply the Leibniz rule
            std::vector<size_t> indegrees = main_term.second.in_degrees();
            size_t incoming_edges = 0;
            for (size_t i = 0; i != indegrees.size(); ++i)
                incoming_edges += indegrees[i];
            std::vector<size_t> leibniz_factors(incoming_edges);
            // Build list of pointers to targets that should be replaced to apply the Leibniz rule
            std::vector< size_t* > targets(incoming_edges);
            size_t target_idx = 0;
            size_t main_idx = new_targets.size() - main_targets.size();
            for (size_t i = 0; i != arguments.size(); ++i)
            {
                for (size_t n : main_term.second.neighbors_in(i))
                {
                    KontsevichGraph* graph = &arguments[i][(*arg_indices)[i]].second;
                    leibniz_factors[target_idx] = graph->vertices();
                    auto& source = new_targets[main_idx + n - main_term.second.external()];
                    // Use the fact that main_term's ground vertices were not relabeled:
                    targets[target_idx] = (source.first == i) ? &source.first : &source.second;
                    target_idx++;
                }
            }
            auto leibniz_indices = CartesianProduct(leibniz_factors);
            for (auto target_indices = leibniz_indices.begin(); target_indices != leibniz_indices.end(); ++target_indices) // Leibniz rule
            {
                // Now assign targets according to the Leibniz rule
                size_t target_idx = 0;
                for (size_t i = 0; i != arguments.size(); ++i)
                {
                    KontsevichGraph* graph = &arguments[i][(*arg_indices)[i]].second;
                    for (size_t j = 0; j != indegrees[i]; ++j)
                    {
                        if ((*target_indices)[target_idx] >= graph->external()) // internal
                            *targets[target_idx] = start_internal_vec[i] - graph->external() + (*target_indices)[target_idx];
                        else // external
                            *targets[target_idx] = start_external_vec[i] + (*target_indices)[target_idx];
                        target_idx++;
                    }
                }
                // TODO: can we do something here to normalize the result more efficiently?
                KontsevichGraph graph(internal, external, new_targets, sign);
                total.push_back({ coeff, graph });
            }
        }
    }
    return total;
}

template class KontsevichGraphSum<int>;
template std::ostream& operator<<(std::ostream& os, const KontsevichGraphSum<int>& gs);
template std::istream& operator>>(std::istream& is, KontsevichGraphSum<int>& sum);
template KontsevichGraphSum<int> operator-(KontsevichGraphSum<int> lhs, const KontsevichGraphSum<int>& rhs);
template KontsevichGraphSum<int> operator+(KontsevichGraphSum<int> lhs, const KontsevichGraphSum<int>& rhs);
