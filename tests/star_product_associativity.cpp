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

DECLARE_FUNCTION_3P(phi)
REGISTER_FUNCTION(phi, dummy())
DECLARE_FUNCTION_3P(u)
REGISTER_FUNCTION(u, dummy())

int main()
{
    map< size_t, set<KontsevichGraph> > relevants;
    map<KontsevichGraph, ex> weights;
    cout << "Reading in known graphs and their weights...";
    ifstream weights_file("data/known_weights.txt");
    weights_file.ignore(numeric_limits<streamsize>::max(), '\n'); // ignore first line with column headers
    KontsevichGraph graph;
    string weight_str;
    parser reader;
    while (weights_file >> graph >> weight_str)
    {
        ex weight = reader(weight_str);
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
    symtab weights_table;
    for (size_t n = 1; n <= order; ++n)
    {
        primes[n].reserve(2*relevants[n].size());
        for (KontsevichGraph g : relevants[n])
        {
            if (n == order)
            {
                symbol weight("w_" + to_string(weight_count++));
                weights_table[weight.get_name()] = weight;
                if (weights.find(g) == weights.end()) // if not already known
                {
                    weights[g] = weight;
                    weight_vars.append(weight);
                }
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
    cout << "Reading in known weight relations...";
    parser weights_reader(weights_table);
    ifstream weight_relations_file("data/weight_relations.txt");
    lst weight_relations;
    for (string lhs, rhs; getline(weight_relations_file, lhs, '=') && weight_relations_file.ignore(1) && getline(weight_relations_file, rhs); )
    {
        weight_relations.append(weights_reader(lhs) == weights_reader(rhs));
    }
    for (auto& some_pair: weights)
    {
        some_pair.second = some_pair.second.subs(weight_relations);
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
    star_product.reduce();
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
    set<ex, ex_is_less> weight_system;
    cout << "Computing cyclic linear weight relations...\n";
    for (auto& term : star_product[order])
    {
        KontsevichGraph graph = term.second;
        // coefficient = multiplicity * weight / n! => weight = n! * coefficient / multiplicity
        ex linear_combination = (graph.internal() % 2 ? 1 : -1) * star_product[order][graph] / graph.multiplicity();
        // iterate over subsets of set of edges: 2^(2n) elements
        vector<size_t> selector(2*graph.internal(), 2);
        CartesianProduct edge_selector(selector);
        for (auto edge_selection = edge_selector.begin(); edge_selection != edge_selector.end(); ++edge_selection)
        {
            std::vector<KontsevichGraph::VertexPair> targets = graph.targets();

            // skip those where an edge is already connected to the first ground vertex
            bool skip = false;
            for (size_t choice = 0; choice != 2*graph.internal(); ++choice)
                if ((*edge_selection)[choice] && ((choice % 2 == 0 && targets[choice/2].first == 0) || (choice % 2 == 1 && targets[choice/2].second == 0)))
                    skip = true;
            if (skip)
                continue;

            // replace edges by edges to the first ground vertex:
            for (size_t choice = 0; choice != 2*graph.internal(); ++choice)
            {
                if ((*edge_selection)[choice]) // replace edge
                {
                    if (choice % 2 == 0)
                        targets[choice/2].first = 0;
                    else
                        targets[choice/2].second = 0;
                }
            }
            KontsevichGraph modified_graph(graph.internal(), graph.external(), targets, graph.sign());
            // build linear relation
            linear_combination -= (modified_graph.in_degree(0) % 2 ? 1 : -1) * star_product[order][modified_graph] / modified_graph.multiplicity();
        }
        weight_system.insert(linear_combination == 0);
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
    for (size_t n = 0; n <= order; ++n)
    {
        cout << "h^" << n << ":\n";
        cout << assoc[n].size() << " total\n";
        for (std::vector<size_t> indegrees : assoc[n].in_degrees())
        {
            for (size_t j = 0; j != indegrees.size(); ++j)
                cout << indegrees[j] << " ";
            cout << ": " << assoc[n][indegrees].size() << "\n";
            cout.flush();

            symbol x("x"), y("y"), z("z");

            std::vector<PoissonStructure> poisson_structures {
                { { x, y, z }, { {0, u(x,y,z)*phi(x,y,z).diff(z), -u(x,y,z)*phi(x,y,z).diff(y)},
                                 {-u(x,y,z)*phi(x,y,z).diff(z), 0, u(x,y,z)*phi(x,y,z).diff(x) },
                                 { u(x,y,z)*phi(x,y,z).diff(y), -u(x,y,z)*phi(x,y,z).diff(x), 0 } } },
            };

            for (PoissonStructure& poisson : poisson_structures)
            {
                typedef std::vector< std::multiset<size_t> > multi_index;
                map< multi_index, map<ex, ex, ex_is_less> > coefficients;
                for (auto& term : assoc[n][indegrees])
                {
                    map_operator_coefficients_from_graph(term.second, poisson, [&coefficients, &term](multi_index arg_derivatives, GiNaC::ex summand) {
                            ex result = (term.first * summand).expand();
                            if (result == 0)
                                return;
                            if (!is_a<add>(result))
                                result = lst(result);
                            for (auto term : result)
                            {
                                ex coefficient = 1;
                                ex derivatives = 1;
                                if (!is_a<mul>(term))
                                    term = lst(term);
                                for (auto factor : term)
                                {
                                    if  (is_a<GiNaC::function>(factor) || is_a<fderivative>(factor))
                                        derivatives *= factor;
                                    else if (is_a<numeric>(factor) || is_a<symbol>(factor))
                                        coefficient *= factor;
                                    else if (is_a<power>(factor))
                                    {
                                        if (is_a<GiNaC::function>(factor.op(0)) || is_a<fderivative>(factor.op(0)))
                                            derivatives *= factor;
                                        else if (is_a<numeric>(factor.op(0)) || is_a<symbol>(factor.op(0)))
                                            coefficient *= factor;
                                    }
                                    else
                                    {
                                        cout << "What the hell is " << factor << "?\n";
                                        // TODO: return 1;
                                    }
                                }
                                coefficients[arg_derivatives][derivatives] += coefficient;
                            }
                    });
                    cout << ".";
                    cout.flush();
                }
                cout << "\n";

                for (auto pair : coefficients)
                {
                    for (auto pair2 : pair.second)
                    {
                        ex result2 = pair2.second;
                        if (result2 == 0)
                            continue;
                        for (ex var : weight_vars)
                        {
                            if (result2.coeff(var) != 0) // first weight variable occurring in expression
                            {
                                result2 /= result2.coeff(var); // divide by its coefficient
                                break;
                            }
                        }
                        weight_system.insert(result2 == 0);
                    }
                }
            }
        }
    }
    cout << "Got system of " << weight_system.size() << " linear equations in " << weight_vars.nops() << " unknowns:\n";
    lst weight_system_lst;
    for (ex eq : weight_system)
    {
        if (eq.lhs() != eq.rhs()) // not a tautology
        {
            cout << eq << endl;
            weight_system_lst.append(eq);
        }
    }
    cout << "Solving it...\n";
    for (ex eq : lsolve(weight_system_lst, weight_vars))
        if (eq.lhs() != eq.rhs()) // not a tautology
            cout << eq << endl;
}
