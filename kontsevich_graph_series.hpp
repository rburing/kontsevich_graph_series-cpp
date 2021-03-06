#ifndef INCLUDED_KONTSEVICH_GRAPH_SERIES_H_
#define INCLUDED_KONTSEVICH_GRAPH_SERIES_H_

#include "kontsevich_graph_sum.hpp"
#include <ostream>
#include <map>
#include <limits>

template<class T> class KontsevichGraphSeries;
template<class T> std::ostream& operator<<(std::ostream&, const KontsevichGraphSeries<T>&);

template<class T>
class KontsevichGraphSeries : public std::map< size_t, KontsevichGraphSum<T> >
{
    size_t d_precision = std::numeric_limits<std::size_t>::max();

    using std::map< size_t, KontsevichGraphSum<T> >::map; // inherit constructors

    public:
    size_t precision() const;
    size_t precision(size_t new_precision);
    KontsevichGraphSeries<T> operator()(std::vector< KontsevichGraphSeries<T> > arguments) const;
    KontsevichGraphSeries<T>& operator+=(const KontsevichGraphSeries<T>& rhs);
    KontsevichGraphSeries<T>& operator-=(const KontsevichGraphSeries<T>& rhs);
    KontsevichGraphSeries<T> symmetrization() const;
    KontsevichGraphSeries<T> skew_symmetrization() const;
    KontsevichGraphSeries<T> inverse() const;
    KontsevichGraphSeries<T> gauge_transform(const KontsevichGraphSeries<T>& gauge);
    bool operator==(int other) const;
    bool operator!=(int other) const;
    void reduce_mod_skew();

    static KontsevichGraphSeries<T> from_istream(std::istream& is, std::function<T(std::string)> const& parser, std::function<bool(KontsevichGraph, size_t)> const& filter = nullptr);

    friend std::ostream& operator<< <>(std::ostream& os, const KontsevichGraphSeries<T>& series);
};

template <class T>
KontsevichGraphSeries<T> operator+(KontsevichGraphSeries<T> lhs, const KontsevichGraphSeries<T>& rhs);
template <class T>
KontsevichGraphSeries<T> operator-(KontsevichGraphSeries<T> lhs, const KontsevichGraphSeries<T>& rhs);
template <class T>
KontsevichGraphSeries<T> gerstenhaber_bracket(const KontsevichGraphSeries<T>& left, const KontsevichGraphSeries<T>& right);
template <class T>
KontsevichGraphSeries<T> schouten_bracket(const KontsevichGraphSeries<T>& left, const KontsevichGraphSeries<T>& right);

#include "kontsevich_graph_series.tpp"

#endif
