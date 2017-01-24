#ifndef INCLUDED_KONTSEVICH_GRAPH_SUM_H_
#define INCLUDED_KONTSEVICH_GRAPH_SUM_H_

#include <vector>
#include <utility>
#include <iostream>
#include "kontsevich_graph.hpp"

template<class T> class KontsevichGraphSum;
template<class T> std::ostream& operator<<(std::ostream&, const std::pair<T, KontsevichGraph>&);
template<class T> std::ostream& operator<<(std::ostream&, const KontsevichGraphSum<T>&);
template<class T> std::istream& operator>>(std::istream&, KontsevichGraphSum<T>&);

template<class T>
class KontsevichGraphSum : public std::vector< std::pair<T, KontsevichGraph> >
{
    using std::vector< std::pair<T, KontsevichGraph> >::vector; // inherit constructors
    using std::vector< std::pair<T, KontsevichGraph> >::operator[]; // inherit subscript operator

    public:
    typedef std::pair<T, KontsevichGraph> Term;
    KontsevichGraphSum<T> operator[](std::vector<size_t> indegrees) const;
    T operator[](KontsevichGraph) const;
    KontsevichGraphSum<T> operator()(std::vector< KontsevichGraphSum<T> > arguments) const;
    KontsevichGraphSum<T>& operator+=(const KontsevichGraphSum<T>& rhs);
    KontsevichGraphSum<T>& operator-=(const KontsevichGraphSum<T>& rhs);
    KontsevichGraphSum<T>& operator=(const KontsevichGraphSum<T>&) = default;
    std::vector< std::vector<size_t> > in_degrees(bool ascending = false) const;
    KontsevichGraphSum<T> skew_symmetrization() const;
    void reduce();
    void reduce_mod_permutations();
    bool operator==(const KontsevichGraphSum<T>& other) const;
    bool operator==(int other) const;
    bool operator!=(const KontsevichGraphSum<T>& other) const;

    friend std::ostream& operator<< <>(std::ostream& os, const KontsevichGraphSum<T>::Term& term);
    friend std::ostream& operator<< <>(std::ostream& os, const KontsevichGraphSum<T>& gs);
    friend std::istream& operator>> <>(std::istream& is, KontsevichGraphSum<T>& sum);
};

template <class T>
KontsevichGraphSum<T> operator+(KontsevichGraphSum<T> lhs, const KontsevichGraphSum<T>& rhs);
template <class T>
KontsevichGraphSum<T> operator-(KontsevichGraphSum<T> lhs, const KontsevichGraphSum<T>& rhs);
template <class T>
KontsevichGraphSum<T> operator*(T lhs, KontsevichGraphSum<T> rhs);

#include "kontsevich_graph_sum.tpp"

#endif
