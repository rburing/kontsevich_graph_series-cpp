#ifndef INCLUDED_LEIBNIZ_GRAPH_H_
#define INCLUDED_LEIBNIZ_GRAPH_H_

#include "kontsevich_graph.hpp"
#include <vector>
#include <string>
#include <map>
#include <istream>

class LeibnizGraph : public std::pair<KontsevichGraph, std::vector<KontsevichGraph::VertexPair> >
{
    bool skew;

public:
    LeibnizGraph() {};
    LeibnizGraph(KontsevichGraph graph, std::vector<KontsevichGraph::VertexPair> jacobiators, bool skew_leibniz = false);
    std::string encoding() const;
    template<class T> static std::map<LeibnizGraph, T> map_from_istream(std::istream& is, std::function<T(std::string)> const& parser = nullptr);
private:
    friend std::istream& operator>>(std::istream& is, LeibnizGraph& g);
};

#include "leibniz_graph.tpp"

#endif
