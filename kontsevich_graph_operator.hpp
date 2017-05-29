#ifndef INCLUDED_KONTSEVICH_GRAPH_OPERATOR_
#define INCLUDED_KONTSEVICH_GRAPH_OPERATOR_

#include <ginac/ginac.h>
#include "kontsevich_graph_series.hpp"
#include "util/cartesian_product.hpp"
#include "util/poisson_structure.hpp"

typedef std::multiset<size_t> multi_index;
typedef std::vector<multi_index> multi_indexes;

void map_operator_coefficients_from_graph(KontsevichGraph graph, PoissonStructure& poisson, std::function<void(multi_indexes, GiNaC::ex)> fun)
{
    multi_indexes external_indices_template(graph.external());
    for (size_t n = 0; n != graph.external(); ++n)
        for (size_t j : graph.neighbors_in(n))
            external_indices_template[n].insert( ((size_t)graph.targets(j).first == n) ? 2*(j-graph.external()) : 2*(j-graph.external()) + 1 );

    GiNaC::ex result = 0;
    size_t dimension = poisson.coordinates.size();
    std::vector<size_t> max_index(2*graph.internal(), dimension);
    CartesianProduct index_product(max_index);
    for (auto indices = index_product.begin(); indices != index_product.end(); ++indices)
    {
        multi_indexes external_indices(graph.external());
        for (size_t n = 0; n != graph.external(); ++n)
            for (size_t j : external_indices_template[n])
                external_indices[n].insert((*indices)[j]);

        GiNaC::ex summand = graph.sign();
        for (size_t n : graph.internal_vertices())
        {
            std::pair<size_t, size_t> bivector_indices { (*indices)[2*(n-graph.external())], (*indices)[2*(n-graph.external()) + 1] };
            multi_index partial_derivatives;
            for (size_t j : graph.neighbors_in(n))
            {
                size_t incoming_index = ((size_t)graph.targets(j).first == n) ? (*indices)[2*(j-graph.external())] : (*indices)[2*(j-graph.external()) + 1];
                partial_derivatives.insert(incoming_index);
            }
            try {
                summand *= poisson.bivector_derivatives_cache.at({ bivector_indices, partial_derivatives });
            }
            catch (std::out_of_range)
            {
                GiNaC::ex factor = poisson.bivector[bivector_indices.first][bivector_indices.second];
                for (size_t j : partial_derivatives)
                {
                    factor = diff(factor, poisson.coordinates[j]);
                }
                summand *= factor;
                poisson.bivector_derivatives_cache[{ bivector_indices, partial_derivatives }] = factor;
            }
        }
        fun(external_indices, summand);
    }
}

GiNaC::ex operator_from_graph(KontsevichGraph graph, PoissonStructure& poisson, std::vector<GiNaC::ex> arguments)
{
    // TODO: check if graph.external() == arguments.size()
    GiNaC::ex result = 0;
    map_operator_coefficients_from_graph(graph, poisson, [&result, &poisson, &arguments](multi_indexes derivatives, GiNaC::ex summand) {
            GiNaC::ex tail = 1;
            for (size_t m = 0; m != arguments.size(); ++m)
            {
                GiNaC::ex factor = arguments[m];
                for (size_t idx : derivatives[m])
                    factor = factor.diff(poisson.coordinates[idx]);
                tail *= factor;
            }
            result += summand * tail;
    });
    return result;
}

void operator_coefficient_from_graph(KontsevichGraph graph, PoissonStructure& poisson, GiNaC::ex coefficient, std::map<multi_indexes, GiNaC::ex>& accumulator)
{
    map_operator_coefficients_from_graph(graph, poisson, [&coefficient, &accumulator](multi_indexes derivatives, GiNaC::ex summand) {
        accumulator[derivatives] += coefficient  * summand;
    });
}

GiNaC::ex evaluate(KontsevichGraphSum<GiNaC::ex> terms, PoissonStructure& poisson, std::vector<GiNaC::ex> arguments)
{
    GiNaC::ex total = 0;
    for (auto& term : terms)
    {
        total += term.first * operator_from_graph(term.second, poisson, arguments);
    }
    return total;
}

std::map<multi_indexes, GiNaC::ex> evaluate_coefficients(KontsevichGraphSum<GiNaC::ex> terms, PoissonStructure& poisson)
{
    std::map<multi_indexes, GiNaC::ex> accumulator;
    for (auto& term : terms)
    {
        operator_coefficient_from_graph(term.second, poisson, term.first, accumulator);
    }
    return accumulator;
}

#endif
