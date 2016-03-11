#include "../kontsevich_graph_series.hpp"
#include "../kontsevich_graph_operator.hpp"
#include "../util/partitions.hpp"
#include "../util/factorial.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits>
using namespace std;
using namespace GiNaC;

size_t order = 4;

int main()
{
    map< size_t, set<KontsevichGraph> > relevants;
    map<KontsevichGraph, ex> weights;
    cout << "Reading in known graphs and their weights...";
    ifstream weights_file("data/known_weights.txt");
    weights_file.ignore(numeric_limits<streamsize>::max(), '\n'); // ignore first line with column headers
    KontsevichGraph graph;
    int numerator, denominator;
    char slash;
    while (weights_file >> graph >> numerator >> slash >> denominator)
    {
        ex weight = numerator/(ex)denominator;
        weights[graph] = weight;
        relevants[graph.internal()].insert(graph);
    }
    cout << " done." << endl;
    cout << "Listing relevant graphs (prime, positive differential order, nonzero, modulo mirror images) at order " << order << "...";
    cout.flush();
    relevants[order] = KontsevichGraph::graphs(order, 2, true, true,
                    [](KontsevichGraph g) -> bool
                    {
                        return g.positive_differential_order() && g.is_prime() && !g.is_zero();
                    });
    cout << " done." << endl;
    cout << "Making a table of primes and weights...";
    map< size_t, vector<KontsevichGraph> > primes;
    size_t weight_count = 0;
    lst weight_vars;
    for (size_t n = 1; n <= order; ++n)
    {
        primes[n].reserve(2*relevants[n].size());
        for (KontsevichGraph g : relevants[n])
        {
            if (n == order)
            {
                symbol weight("w_" + to_string(weight_count++));
                weights[g] = weight;
                weight_vars.append(weight);
            }
            primes[n].push_back(g);
            KontsevichGraph mirror = g.mirror_image();
            if (g.abs() != mirror.abs())
            {
                weights[mirror] = weights[g];
                if (n % 2 == 1)
                    weights[mirror] *= -1;
                primes[n].push_back(mirror);
            }
        }
    }
    cout << " done." << endl;
    cout << "Constructing star product at order";
    KontsevichGraphSeries<ex> star_product;
    star_product[0] = KontsevichGraphSum<ex> { { 1, KontsevichGraph(0, 2, {}) } };
    for (size_t n = 1; n <= order; ++n)
    {
        cout << " " << n << ",";
        cout.flush();
        ex major_coeff = 1/(ex)factorial(n);
        Partitions partitions(n);
        for (auto partition = partitions.begin(); partition != partitions.end(); ++partition)
        {
            // count multiplicities of parts
            std::map<size_t, std::vector<size_t> > multiplicity;
            for (size_t i = 0; i != (*partition).size(); ++i)
                multiplicity[(*partition)[i]].push_back(i);
            // ignore multiplicity 1
            for (auto part : *partition)
                if (multiplicity[part].size() == 1)
                    multiplicity.erase(part);

            std::vector<size_t> prime_sizes((*partition).size());
            for (size_t i = 0; i != (*partition).size(); ++i)
                prime_sizes[i] = primes[(*partition)[i]].size();
            CartesianProduct decompositions(prime_sizes);
            for (auto decomposition = decompositions.begin(); decomposition != decompositions.end(); ++decomposition)
            {
                bool accept = true;
                for (auto& part : multiplicity)
                {
                    size_t prev = 0;
                    for (auto& idx : part.second) // when drawing from multiple equal sets
                    {
                        size_t current = (*decomposition)[idx];
                        if (current < prev) // accept only non-decreasing indices, to avoid duplicates
                            accept = false;
                        prev = current;
                    }
                }
                if (!accept)
                    continue;

                KontsevichGraph composite(0, 2, {});
                ex coeff = major_coeff;
                for (size_t i = 0; i != (*partition).size(); ++i)
                {
                    composite *= primes[(*partition)[i]][(*decomposition)[i]];
                    coeff *= weights[primes[(*partition)[i]][(*decomposition)[i]]];
                }
                coeff *= composite.multiplicity();
                star_product[n] += KontsevichGraphSum<ex>({ { coeff, composite} });
            }
        }
        // TODO: n is not yet considered a partition of n, so treat it separately for now:
        for (KontsevichGraph prime : primes[n])
        {
            star_product[n] += KontsevichGraphSum<ex>({ { major_coeff * weights[prime] * prime.multiplicity(), prime } });
        }
    }
    cout << " done:" << endl;
    cout << star_product << "\n";
    cout << "Number of terms in star product:\n";
    for (size_t n = 0; n <= order; ++n)
    {
        cout << "h^" << n << ": " << star_product[n].size() << " total\n";
        for (std::vector<size_t> indegrees : star_product[n].in_degrees())
        {
            for (size_t j = 0; j != indegrees.size(); ++j)
                cout << indegrees[j] << " ";
            cout << ": " << star_product[n][indegrees].size() << "\n";
        }
    }
    cout << "Computing associator...";
    cout.flush();
    KontsevichGraphSeries<ex> arg = { { 0, { { 1, KontsevichGraph(0, 1, {}) } }} };
    KontsevichGraphSeries<ex> assoc = star_product({ star_product, arg }) - star_product({ arg, star_product });
    cout << endl;
    cout << "Reducing associator...";
    cout.flush();
    assoc.reduce();
    cout << endl;
    cout << "Number of terms in associator:\n";
    lst weight_system;
    for (size_t n = 0; n <= order; ++n)
    {
        cout << "h^" << n << ":\n";
        cout << assoc[n].size() << " total\n";
        for (std::vector<size_t> indegrees : assoc[n].in_degrees())
        {
            for (size_t j = 0; j != indegrees.size(); ++j)
                cout << indegrees[j] << " ";
            cout << ": " << assoc[n][indegrees].size() << "\n";

            symbol x("x"), y("y"), z("z");
            ex h = 1/4 * pow(x,5)*pow(y,3)*pow(z,4) + pow(y,5)*z*x + pow(z,5)*pow(x,2);

            std::vector<PoissonStructure> poisson_structures {
                { { x, y, z }, { {0, h.diff(z), -h.diff(y)},
                                 {-h.diff(z), 0, h.diff(x) },
                                 { h.diff(y), -h.diff(x), 0 } } },
            };
            map< vector<symbol>, vector<ex> > arguments {
                { { x, y, z }, { pow(x,3)*z*pow(y,2), pow(y,3) + y*pow(x,10), pow(z,4) * y * pow(x,2) } }
            };

            for (PoissonStructure& poisson : poisson_structures)
            {
                ex result = evaluate(assoc[n][indegrees], poisson, arguments[poisson.coordinates]).expand();
                vector<size_t> degrees;
                for (symbol& variable : poisson.coordinates)
                    degrees.push_back(result.degree(variable));
                CartesianProduct monomialdegrees_list(degrees);
                for (auto monomialdegrees = monomialdegrees_list.begin(); monomialdegrees != monomialdegrees_list.end(); ++monomialdegrees)
                {
                    ex result2 = result;
                    for (size_t i = 0; i != poisson.coordinates.size(); ++i)
                    {
                        result2 = result2.coeff(poisson.coordinates[i], (*monomialdegrees)[i]).expand();
                    }
                    if (result2 != 0)
                    {
                        weight_system.append(result2 == 0);
                        cout << result2 << " == 0\n";
                    }
                }
            }
        }
    }
    cout << "Got system of " << weight_system.nops() << " linear equations in " << weight_vars.nops() << " unknowns:\n";
    for (ex eq : weight_system)
        if (eq.lhs() != eq.rhs()) // not a tautology
            cout << eq << endl;
    cout << "Solving it...\n";
    for (ex eq : lsolve(weight_system, weight_vars))
        if (eq.lhs() != eq.rhs()) // not a tautology
            cout << eq << endl;
}
