#ifndef INCLUDED_KONTSEVICH_GRAPH_SERIES_H_
#define INCLUDED_KONTSEVICH_GRAPH_SERIES_H_

#include "kontsevich_graph_sum.hpp"
#include <ostream>
#include <map>

template<class T> class KontsevichGraphSeries;
template<class T> std::ostream& operator<<(std::ostream&, const KontsevichGraphSeries<T>&);

template<class T>
class KontsevichGraphSeries : public std::map< size_t, KontsevichGraphSum<T> >
{
    using std::map< size_t, KontsevichGraphSum<T> >::map; // inherit constructors

    public:
    void reduce();

    friend std::ostream& operator<< <>(std::ostream& os, const KontsevichGraphSeries<T>& series);
};

#endif
