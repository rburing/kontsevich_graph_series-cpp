#include "kontsevich_graph_series.hpp"

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

template class KontsevichGraphSeries<int>;
template std::ostream& operator<<(std::ostream& os, const KontsevichGraphSeries<int>& gs);
