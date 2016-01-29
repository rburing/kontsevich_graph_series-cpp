#include "kontsevich_graph_sum.hpp"

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
        current_term++;
    }
}

template <class T>
std::ostream& operator<<(std::ostream& os, const KontsevichGraphSum<T>& gs)
{
    if (gs.size() == 0)
        return os << "0";
    for (auto &term : gs)
    {
        os << term.first << "*(" << term.second << ")";
        if (&term != &gs.back())
            os << " + ";
    }
    return os;
}

template class KontsevichGraphSum<int>;
template std::ostream& operator<<(std::ostream& os, const KontsevichGraphSum<int>& gs);
