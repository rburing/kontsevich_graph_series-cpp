#ifndef INCLUDED_LEIBNIZ_GRAPH_H_
#define INCLUDED_LEIBNIZ_GRAPH_H_

#include "kontsevich_graph.hpp"
#include <vector>
#include <string>
#include <sstream>

class LeibnizGraph : public std::pair<KontsevichGraph, std::vector<KontsevichGraph::VertexPair> >
{
    using std::pair<KontsevichGraph, std::vector<KontsevichGraph::VertexPair> >::pair;

public:
    std::string encoding() const;
};

#include "leibniz_graph.tpp"

#endif
