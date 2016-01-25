#include <kontsevich_graph_series.hpp>

template<class T>
KontsevichGraphSeries<T>::KontsevichGraphSeries(std::map< size_t, KontsevichGraphSum<T> > terms)
: d_terms(terms)
{}

template<class T>
KontsevichGraphSum<T>& KontsevichGraphSeries<T>::operator[](size_t order)
{
    return d_terms[order];
}

template <class T>
std::ostream& operator<<(std::ostream& os, const KontsevichGraphSeries<T>& series)
{
    if (series.d_terms.size() == 0)
        return os << "0";
    auto final_term = series.d_terms.end();
    --final_term;
    for (auto term = series.d_terms.begin(); term != series.d_terms.end(); term++)
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

template class KontsevichGraphSeries<int>;
template std::ostream& operator<<(std::ostream& os, const KontsevichGraphSeries<int>& gs);
