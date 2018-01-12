#include "leibniz_graph.hpp"
#include <sstream>

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

std::istream& operator>>(std::istream& is, LeibnizGraph& g)
{
    size_t jacobiators;
    is >> jacobiators;
    is >> g.first;
    g.second.clear();
    std::pair<size_t, size_t> jacobiator;
    size_t jacobiator_count = 0;
    while (jacobiator_count++ < jacobiators && is >> jacobiator.first >> jacobiator.second)
        g.second.push_back(jacobiator);
    return is;
}

template<class T>
std::map<LeibnizGraph, T> LeibnizGraph::map_from_istream(std::istream& is, std::function<T(std::string)> const& parser)
{
    std::map<LeibnizGraph, T> result;
    if (parser == nullptr)
        return result;
    for (std::string line; getline(is, line);)
    {
        if (line.length() == 0 || line[0] == '#') // also skip comments
            continue;
        std::stringstream ss(line);
        LeibnizGraph g;
        ss >> g;
        std::string coefficient_str;
        ss >> coefficient_str;
        T coefficient = parser(coefficient_str);
        result[g] = coefficient;
    }
    return result;
}
