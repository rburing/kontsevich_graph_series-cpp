#ifndef INCLUDED_LEIBNIZ_GRAPH_H_
#define INCLUDED_LEIBNIZ_GRAPH_H_

#include "kontsevich_graph.hpp"
#include "kontsevich_graph_sum.hpp"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <istream>

template<class T> class LeibnizGraph;
template<class T> std::istream& operator>>(std::istream& is, LeibnizGraph<T>& g);

template<class T>
class LeibnizGraph : KontsevichGraph
{
    std::vector<KontsevichGraph::VertexPair> d_jacobiators;
    bool d_skew;
    size_t d_max_jac_indegree;
    std::vector< std::set<KontsevichGraph::Vertex*> > d_jacobiator_targets;
    std::map<KontsevichGraph::Vertex*, size_t> d_leibniz_targets;

public:
    LeibnizGraph() {};
    LeibnizGraph(KontsevichGraph graph, std::vector<KontsevichGraph::VertexPair> jacobiators, bool skew = false);
    LeibnizGraph(const LeibnizGraph<T>& other);
    LeibnizGraph& operator=(const LeibnizGraph<T>& other);
    std::string encoding() const;
    static std::map< LeibnizGraph<T>, T> map_from_istream(std::istream& is, std::function<T(std::string)> const& parser = nullptr);
    static std::set< LeibnizGraph<T> > those_yielding_kontsevich_graph(KontsevichGraph& graph, bool skew_leibniz = false);
    size_t max_jac_indegree() const;
    bool operator<(const LeibnizGraph<T>& rhs) const;
    KontsevichGraphSum<T> expansion(T prefactor = 1);
    void normalize();

private:
    void set_jacobiator_and_leibniz_targets();
    friend std::istream& operator>> <>(std::istream& is, LeibnizGraph<T>& g);
};

#include "leibniz_graph.tpp"

#endif
