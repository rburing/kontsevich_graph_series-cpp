#ifndef INCLUDED_KONTSEVICH_GRAPH_SERIES_H_
#define INCLUDED_KONTSEVICH_GRAPH_SERIES_H_

#include <kontsevich_graph_sum.hpp>
#include <ostream>
#include <map>

template<class T> class KontsevichGraphSeries;
template<class T> std::ostream& operator<<(std::ostream&, const KontsevichGraphSeries<T>&);

template<class T>
class KontsevichGraphSeries
{
    std::map< size_t, KontsevichGraphSum<T> > d_terms;
    
    public:
    KontsevichGraphSeries(std::map< size_t, KontsevichGraphSum<T> > terms);
    KontsevichGraphSum<T>& operator[](size_t order);

    friend std::ostream& operator<< <>(std::ostream& os, const KontsevichGraphSeries<T>& series);
};

#endif
