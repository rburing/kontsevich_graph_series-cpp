#include <kontsevich_graph_sum.hpp>

template<class T>
KontsevichGraphSum<T>::KontsevichGraphSum(std::vector< std::pair<T, KontsevichGraph> > terms)
: d_terms(terms)
{}

template<class T>
size_t KontsevichGraphSum<T>::size() const
{
    return d_terms.size();
}

template <class T>
void KontsevichGraphSum<T>::reduce()
{
    auto current_term = d_terms.begin();
    current_term->first *= current_term->second.sign();
    current_term->second.sign(1);
    while (current_term < d_terms.end())
    {
        auto subsequent_term = current_term + 1;
        while (subsequent_term < d_terms.end())
        {
            if (subsequent_term->second.abs() == current_term->second.abs())
            {
                current_term->first += subsequent_term->first * subsequent_term->second.sign();
                subsequent_term = d_terms.erase(subsequent_term);
            }
            else
                subsequent_term++;
        }
        current_term++;
    }
}

template <class T>
std::ostream& operator<<(std::ostream& os, const KontsevichGraphSum<T>& gs)
{
    if (gs.size() == 0)
        return os << "0";
    for (auto &term : gs.d_terms)
    {
        os << term.first << "*(" << term.second << ")";
        if (&term != &gs.d_terms.back())
            os << " + ";
    }
    return os;
}

template class KontsevichGraphSum<int>;
template std::ostream& operator<<(std::ostream& os, const KontsevichGraphSum<int>& gs);
