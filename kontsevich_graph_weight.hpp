#ifndef INCLUDED_KONTSEVICH_GRAPH_WEIGHT_
#define INCLUDED_KONTSEVICH_GRAPH_WEIGHT_

#include <ginac/ginac.h>
#include "kontsevich_graph.hpp"
#include <vector>
#include <string>

GiNaC::matrix jacobian(GiNaC::matrix f, GiNaC::lst vars)
{
    int rows = f.rows();
    int cols = vars.nops();
    GiNaC::matrix J(rows, cols);
    for (int i = 0; i != rows; ++i)
    {
        for (int j = 0; j != cols; ++j)
        {
            GiNaC::symbol var = GiNaC::ex_to<GiNaC::symbol>(vars[j]);
            J(i,j) = f(i, 0).diff(var);
        }
    }
    return J;
}

GiNaC::ex phi(std::pair<GiNaC::ex,GiNaC::ex> p1, std::pair<GiNaC::ex,GiNaC::ex> p2)
{
    GiNaC::ex a = p1.first, b = p1.second, x = p2.first, y = p2.second;
    return GiNaC::atan(2*b*(a-x)/(pow(a-x,2) + (y+b)*(y-b)));
}

GiNaC::matrix gauss_map_jacobian(KontsevichGraph graph)
{
    // assume two ground vertices
    std::vector< std::pair<GiNaC::ex, GiNaC::ex> > external_coordinates = { { 0, 0 }, { 1, 0 } };
    std::vector< std::pair<GiNaC::ex, GiNaC::ex> > internal_coordinates;
    GiNaC::lst internal_coordinates_lst;
    for (size_t idx = 0; idx != graph.internal(); ++idx)
    {
        std::string xy = { char('a' + 2*idx), char('a' + 2*idx + 1) };
        GiNaC::symbol x_coordinate(xy.substr(0,1)), y_coordinate(xy.substr(1,1));
        internal_coordinates.push_back({ x_coordinate, y_coordinate });
        internal_coordinates_lst.append(x_coordinate);
        internal_coordinates_lst.append(y_coordinate);
    }
    auto coordinates = [&](KontsevichGraph::Vertex v) -> std::pair<GiNaC::ex,GiNaC::ex>
    {
        return (v < (char)graph.external()) ? external_coordinates[v] : internal_coordinates[v - graph.external()];
    };
    GiNaC::matrix Phi(2*graph.internal(), 1);
    for (KontsevichGraph::Vertex v : graph.internal_vertices())
    {
        KontsevichGraph::VertexPair targets = graph.targets(v);
        Phi(2*(v - (char)graph.external()),     0) = phi(coordinates(v), coordinates(targets.first));
        Phi(2*(v - (char)graph.external()) + 1, 0) = phi(coordinates(v), coordinates(targets.second));
    }
    return jacobian(Phi, internal_coordinates_lst);
}

GiNaC::ex weight_integrand(KontsevichGraph graph)
{
    return gauss_map_jacobian(graph).determinant(); // TODO: too slow
}

#endif
