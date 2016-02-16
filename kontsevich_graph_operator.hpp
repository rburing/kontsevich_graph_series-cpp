#ifndef INCLUDED_KONTSEVICH_GRAPH_OPERATOR_
#define INCLUDED_KONTSEVICH_GRAPH_OPERATOR_

#include <ginac/ginac.h>
#include "kontsevich_graph_series.hpp"
#include "util/cartesian_product.hpp"

struct PoissonStructure
{
    std::vector<GiNaC::symbol> coordinates;
    std::vector< std::vector<GiNaC::ex> > bivector;
};

GiNaC::ex operator_from_graph(KontsevichGraph graph, PoissonStructure poisson, std::vector<GiNaC::ex> arguments)
{
    GiNaC::ex result = 0;
    size_t dimension = poisson.coordinates.size();
    // edge labels run from 1 to dimension
    std::vector<size_t> max_index(2*graph.internal(), dimension);
    CartesianProduct index_product(max_index);
    for (auto indices = index_product.begin(); indices != index_product.end(); ++indices)
    {
        // TODO: test if graph.external() == arguments.size()
        GiNaC::ex summand = graph.sign();
        for (size_t n = 0; n != graph.vertices(); ++n)
        {
            GiNaC::ex factor;
            if (n < graph.external())
                factor = arguments[n];
            else
                factor = poisson.bivector[(*indices)[2*(n-graph.external())]][(*indices)[2*(n-graph.external()) + 1]];
            for (size_t j : graph.neighbors_in(n))
            {
                size_t incoming_index = (graph.targets(j).first == n) ? (*indices)[2*(j-graph.external())] : (*indices)[2*(j-graph.external()) + 1];
                factor = diff(factor, poisson.coordinates[incoming_index]);
            }
            summand *= factor;
        }
        result += summand;
    }
    return result;
}

GiNaC::ex evaluate(KontsevichGraphSum<GiNaC::ex> terms, PoissonStructure poisson, std::vector<GiNaC::ex> arguments)
{
    GiNaC::ex total = 0;
    for (auto& term : terms)
    {
        total += term.first * operator_from_graph(term.second, poisson, arguments);
    }
    return total;
}

#endif
