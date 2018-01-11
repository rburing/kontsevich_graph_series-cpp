#include "leibniz_graph.hpp"

std::string LeibnizGraph::encoding() const
{
    std::stringstream ss;
    ss << this->second.size();
    ss << "   ";
    ss << this->first.encoding();
    ss << "   ";
    for (KontsevichGraph::VertexPair jacobiator : this->second)
        ss << jacobiator.first << " " << jacobiator.second;
    return ss.str();
}
