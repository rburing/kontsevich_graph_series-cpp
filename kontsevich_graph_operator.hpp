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
                size_t incoming_index = ((size_t)graph.targets(j).first == n) ? (*indices)[2*(j-graph.external())] : (*indices)[2*(j-graph.external()) + 1];
                factor = diff(factor, poisson.coordinates[incoming_index]);
            }
            summand *= factor;
        }
        result += summand;
    }
    return result;
}

// TODO: refactor the code duplicated above and below
void operator_coefficient_from_graph(KontsevichGraph graph, PoissonStructure poisson, GiNaC::ex coefficient, std::map< std::vector< std::multiset<size_t> >, GiNaC::ex >& accumulator)
{
    std::vector< std::multiset<size_t> > external_indices_template(graph.external());
    for (size_t n = 0; n != graph.external(); ++n)
        for (size_t j : graph.neighbors_in(n))
            external_indices_template[n].insert( ((size_t)graph.targets(j).first == n) ? 2*(j-graph.external()) : 2*(j-graph.external()) + 1 );

    GiNaC::ex result = 0;
    size_t dimension = poisson.coordinates.size();
    std::vector<size_t> max_index(2*graph.internal(), dimension);
    CartesianProduct index_product(max_index);
    for (auto indices = index_product.begin(); indices != index_product.end(); ++indices)
    {
        std::vector< std::multiset<size_t> > external_indices(graph.external());
        for (size_t n = 0; n != graph.external(); ++n)
            for (size_t j : external_indices_template[n])
                external_indices[n].insert((*indices)[j]);

        GiNaC::ex summand = coefficient * graph.sign();
        for (size_t n : graph.internal_vertices())
        {
            GiNaC::ex factor = poisson.bivector[(*indices)[2*(n-graph.external())]][(*indices)[2*(n-graph.external()) + 1]];
            for (size_t j : graph.neighbors_in(n))
            {
                size_t incoming_index = ((size_t)graph.targets(j).first == n) ? (*indices)[2*(j-graph.external())] : (*indices)[2*(j-graph.external()) + 1];
                factor = diff(factor, poisson.coordinates[incoming_index]);
            }
            summand *= factor;
        }
        accumulator[external_indices] += summand;
    }
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

std::map< std::vector< std::multiset<size_t> >, GiNaC::ex > evaluate_coefficients(KontsevichGraphSum<GiNaC::ex> terms, PoissonStructure poisson)
{
    std::map< std::vector< std::multiset<size_t> >, GiNaC::ex > accumulator;
    for (auto& term : terms)
    {
        operator_coefficient_from_graph(term.second, poisson, term.first, accumulator);
    }
    return accumulator;
}

#endif
