#include "kontsevich_graph_series.hpp"
#include "util/cartesian_product.hpp"
#include <sstream>

template <class T>
size_t KontsevichGraphSeries<T>::precision() const
{
    return d_precision;
}

template <class T>
size_t KontsevichGraphSeries<T>::precision(size_t new_precision)
{
    return d_precision = new_precision;
}

template <class T>
void KontsevichGraphSeries<T>::reduce()
{
    auto current_term = this->begin();
    while (current_term != this->end())
    {
        current_term->second.reduce();
        if (current_term->second.size() == 0)
            current_term = this->erase(current_term);
        else
            ++current_term;
    }
}

template <class T>
void KontsevichGraphSeries<T>::reduce_mod_permutations()
{
    auto current_term = this->begin();
    while (current_term != this->end())
    {
        current_term->second.reduce_mod_permutations();
        if (current_term->second.size() == 0)
            current_term = this->erase(current_term);
        else
            ++current_term;
    }
}

template <class T>
KontsevichGraphSeries<T> KontsevichGraphSeries<T>::operator()(std::vector< KontsevichGraphSeries<T> > arguments) const
{
    KontsevichGraphSeries<T> result;
    // Theoretical precision of the result (may be the theoretical maximum):
    size_t new_precision = precision();
    for (auto& argument : arguments)
        new_precision = std::min(new_precision, argument.precision());
    result.precision(new_precision);
    // Return zero if the series itself or any of its arguments are zero:
    if (this->empty())
        return result;
    for (auto& argument : arguments)
        if (argument.empty())
            return result;
    // Practical precision (actually considering the data available):
    size_t practical_precision = std::min(this->rbegin()->first, new_precision);
    std::vector<size_t> argument_sizes(arguments.size());
    for (size_t i = 0; i != arguments.size(); ++i)
    {
        argument_sizes[i] = arguments[i].rbegin()->first + 1;
    }
    // Actual composition:
    for (size_t n = 0; n <= practical_precision; ++n)
    {
        auto entry = this->find(n);
        if (entry == this->end())
            continue;
        CartesianProduct multilinearity_indices(argument_sizes);
        for (auto arg_indices = multilinearity_indices.begin(); arg_indices != multilinearity_indices.end(); ++arg_indices) // Multi-linearity
        {
            size_t total_order = n;
            for (size_t i = 0; i != arguments.size(); ++i)
                total_order += (*arg_indices)[i];
            if (total_order > practical_precision)
                continue;
            std::vector< KontsevichGraphSum<T> > args(arguments.size());
            for (size_t i = 0; i != arguments.size(); ++i)
                args[i] = arguments[i][(*arg_indices)[i]];
            result[total_order] += entry->second(args);
        }
    }
    return result;
}

template <class T>
KontsevichGraphSeries<T> KontsevichGraphSeries<T>::skew_symmetrization() const
{
    KontsevichGraphSeries<T> total;
    total.precision(this->precision());
    for (auto current_term = this->begin(); current_term != this->end(); ++current_term)
        total[current_term->first] = current_term->second.skew_symmetrization();
    return total;
}

template <class T>
KontsevichGraphSeries<T>& KontsevichGraphSeries<T>::operator+=(const KontsevichGraphSeries<T>& rhs)
{
    size_t practical_precision = this->rbegin()->first;
    for (size_t n = 0; n <= practical_precision; ++n)
    {
        try {
            (*this)[n] += rhs.at(n);
        }
        catch (std::out_of_range) {}
    }
    return *this;
}

template <class T>
KontsevichGraphSeries<T> operator+(KontsevichGraphSeries<T> lhs, const KontsevichGraphSeries<T> &rhs)
{
    lhs += rhs;
    return lhs;
}

template <class T>
KontsevichGraphSeries<T>& KontsevichGraphSeries<T>::operator-=(const KontsevichGraphSeries<T>& rhs)
{
    size_t practical_precision = this->rbegin()->first;
    for (size_t n = 0; n <= practical_precision; ++n)
    {
        try {
            (*this)[n] -= rhs.at(n);
        }
        catch (std::out_of_range) { }
    }
    return *this;
}

template <class T>
KontsevichGraphSeries<T> operator-(KontsevichGraphSeries<T> lhs, const KontsevichGraphSeries<T>& rhs)
{
    lhs -= rhs;
    return lhs;
}

template <class T>
KontsevichGraphSeries<T> KontsevichGraphSeries<T>::from_istream(std::istream& is, std::function<T(std::string)> const& parser, std::function<bool(KontsevichGraph, size_t)> const& filter)
{
    KontsevichGraphSeries<T> graph_series;
    KontsevichGraphSum<T> term;
    size_t order = 0;
    for (std::string line; getline(is, line); )
    {
        if (line.length() == 0)
            continue;
        if (line[0] == 'h')
        {
            graph_series[order] = term;
            term = KontsevichGraphSum<T>({ });
            order = stoi(line.substr(2));
        }
        else
        {
            KontsevichGraph graph;
            std::stringstream ss(line);
            ss >> graph;
            if (filter && !filter(graph, order))
                continue;
            std::string coefficient_str;
            ss >> coefficient_str;
            T coefficient = parser(coefficient_str);
            term += KontsevichGraphSum<T>({ { coefficient, graph } });
        }
    }
    graph_series[order] = term; // the last one
    graph_series.precision(order);
    return graph_series;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const KontsevichGraphSeries<T>& series)
{
    if (series.size() == 0)
        return os << "0";
    auto final_term = series.end();
    --final_term;
    for (auto term = series.begin(); term != series.end(); term++)
    {
        if (term->first == 0)
            os << term->second;
        else
            os << "(" << term->second << ")*h^" << term->first;
        if (term != final_term)
            os << " + ";
    }
    return os;
}
