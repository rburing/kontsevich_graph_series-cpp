#ifndef INCLUDED_LEIBNIZ_GRAPH_H_
#define INCLUDED_LEIBNIZ_GRAPH_H_

#include "kontsevich_graph.hpp"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <istream>

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
    std::string encoding() const;
    template<class T> static std::map<LeibnizGraph, T> map_from_istream(std::istream& is, std::function<T(std::string)> const& parser = nullptr);
    static std::set<LeibnizGraph> those_yielding_kontsevich_graph(KontsevichGraph& graph, bool skew_leibniz = false);
    size_t max_jac_indegree() const;
    bool operator<(const LeibnizGraph& rhs) const;
    void normalize();
private:
    friend std::istream& operator>>(std::istream& is, LeibnizGraph& g);
};

#include "leibniz_graph.tpp"

#endif
