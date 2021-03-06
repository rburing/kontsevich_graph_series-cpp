#include "../kontsevich_graph_series.hpp"
#include "../kontsevich_graph_operator.hpp"
#include "../util/partitions.hpp"
#include "../util/factorial.hpp"
#include "../util/poisson_structure.hpp"
#include "../util/poisson_structure_examples.hpp" // for poisson_structures
#include <ginac/ginac.h>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;
using namespace GiNaC;

void equations_from_particular_poisson(KontsevichGraphSum<ex> graph_sum, PoissonStructure& poisson, set<ex, ex_is_less>& linear_system, lst& unknowns, vector<size_t> point)
{
    lst point_substitution;
    for (size_t i = 0; i != poisson.coordinates.size(); ++i)
        point_substitution.append(poisson.coordinates[i] == point[i]);
    typedef std::vector< std::multiset<size_t> > multi_index;
    map< multi_index, ex > coefficients;
    size_t count = 0;
    for (auto& term : graph_sum)
    {
        map_operator_coefficients_from_graph(term.second, poisson, [&coefficients, &term, &point_substitution](multi_index arg_derivatives, GiNaC::ex summand) {
            ex result = (term.first * summand).subs(point_substitution).expand();
            coefficients[arg_derivatives] += result;
        });
        cerr << "\r" << ++count << " / " << graph_sum.size();
    }
    for (auto& entry : coefficients)
    {
        if (entry.second == 0)
            continue;
        ex result = entry.second;
        // TODO: normalize equation
        pair< set<ex,ex_is_less>::const_iterator, bool > insertion = linear_system.insert(result == 0);
        if (insertion.second) // new equation
            cout << result << "==0\n";
    }
    cout.flush();
    cerr << "\n";
}

void equations_from_polynomial_poisson(KontsevichGraphSum<ex> graph_sum, PoissonStructure& poisson, set<ex, ex_is_less>& linear_system, lst& unknowns)
{
    typedef std::vector< std::multiset<size_t> > multi_index;
    map< multi_index, ex > coefficients;
    size_t count = 0;
    for (auto& term : graph_sum)
    {
        map_operator_coefficients_from_graph(term.second, poisson, [&coefficients, &term](multi_index arg_derivatives, GiNaC::ex summand) {
            ex result = (term.first * summand).expand();
            coefficients[arg_derivatives] += result;
        });
        cerr << "\r" << ++count << " / " << graph_sum.size();
    }
    for (auto& entry : coefficients)
    {
        if (entry.second == 0)
            continue;
        ex result = entry.second;
        vector<size_t> degrees;
        for (symbol const& variable : poisson.coordinates)
            degrees.push_back(result.degree(variable));
        CartesianProduct monomialdegrees_list(degrees);
        for (auto monomialdegrees = monomialdegrees_list.begin(); monomialdegrees != monomialdegrees_list.end(); ++monomialdegrees)
        {
            ex result2 = result;
            for (size_t i = 0; i != poisson.coordinates.size(); ++i)
                result2 = result2.coeff(poisson.coordinates[i], (*monomialdegrees)[i]).expand();
            if (result2 == 0)
                continue;
            // TODO: normalize equation
            pair< set<ex,ex_is_less>::const_iterator, bool > insertion = linear_system.insert(result2 == 0);
            if (insertion.second) // new equation
                cout << result2 << "==0\n";
        }
    }
    cout.flush();
    cerr << "\n";
}

void equations_from_generic_poisson(KontsevichGraphSum<ex> graph_sum, PoissonStructure& poisson, set<ex, ex_is_less>& linear_system, lst& unknowns)
{
    typedef std::vector< std::multiset<size_t> > multi_index;
    map< multi_index, map<ex, ex, ex_is_less> > coefficients;
    size_t count = 0;
    for (auto& term : graph_sum)
    {
        map_operator_coefficients_from_graph(term.second, poisson, [&coefficients, &term](multi_index arg_derivatives, GiNaC::ex summand) {
                ex result = (term.first * summand).expand();
                if (result == 0)
                    return;
                if (!is_a<add>(result))
                    result = lst({ result });
                for (auto term : result)
                {
                    ex coefficient = 1;
                    ex derivatives = 1;
                    if (!is_a<mul>(term))
                        term = lst({ term });
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
                            cerr << "What the hell is " << factor << "?\n";
                        }
                    }
                    coefficients[arg_derivatives][derivatives] += coefficient;
                }
        });
        cerr << "\r" << ++count << " / " << graph_sum.size();
    }
    for (auto pair : coefficients)
    {
        for (auto pair2 : pair.second)
        {
            ex result2 = pair2.second;
            if (result2 == 0)
                continue;
            // TODO: normalize equation
            auto insertion = linear_system.insert(result2 == 0);
            if (insertion.second) // new equation
                cout << result2 << "==0\n";
        }
    }
    cout.flush();
    cerr << "\n";
}

int main(int argc, char* argv[])
{
    if ((argc != 3 && argc != 4) || poisson_structures.find(argv[2]) == poisson_structures.end())
    {
        cerr << "Usage: " << argv[0] << " <graph-series-filename> <poisson-structure> [--linear-solve]\n\n"
             << "Poisson structures can be chosen from the following list:\n";
        for (auto const& entry : poisson_structures)
        {
            cerr << "- " << entry.first << "\n";
        }
        return 1;
    }
    bool solve = argc == 4 && string(argv[3]) == "--linear-solve";

    PoissonStructure& poisson = poisson_structures[argv[2]];

    // Reading in graph series:
    string graph_series_filename(argv[1]);
    ifstream graph_series_file(graph_series_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> graph_series = KontsevichGraphSeries<ex>::from_istream(graph_series_file, [&coefficient_reader](std::string s) -> ex { return coefficient_reader(s); });
    size_t order = graph_series.precision();

    // Make a list of unknowns
    lst unknowns;
    for (std::pair<string, ex> pair : coefficient_reader.get_syms())
        unknowns.append(pair.second);

    // Point for particular Poisson structure:
    vector<size_t> point(poisson.coordinates.size());
    for (size_t i = 0; i != poisson.coordinates.size(); ++i)
        point[i] = i+1;
    // Right now, only the point (1, 2, ..., dim). TODO: more points, options

    cerr << "Number of terms:\n";
    set<ex, ex_is_less> linear_system;
    for (size_t n = 0; n <= order; ++n)
    {
        cerr << "h^" << n << ":\n";
        cerr << graph_series[n].size() << " total\n";
        for (std::vector<size_t> indegrees : graph_series[n].in_degrees(true))
        {
            for (size_t j = 0; j != indegrees.size(); ++j)
                cerr << indegrees[j] << " ";
            cerr << ": " << graph_series[n][indegrees].size() << "\n";
            
            switch (poisson.type)
            {
                case PoissonStructure::Type::Polynomial:
                    equations_from_polynomial_poisson(graph_series[n][indegrees], poisson, linear_system, unknowns);
                    break;
                case PoissonStructure::Type::Generic:
                    equations_from_generic_poisson(graph_series[n][indegrees], poisson, linear_system, unknowns);
                    break;
                case PoissonStructure::Type::Particular:
                    equations_from_particular_poisson(graph_series[n][indegrees], poisson, linear_system, unknowns, point);
                    break;
            }
        }
    }
    if (solve)
    {
        cout << "Got system of " << linear_system.size() << " linear equations in " << unknowns.nops() << " unknowns.\n";
        lst linear_system_lst;
        for (ex eq : linear_system)
            if (eq.lhs() != eq.rhs()) // not a tautology
                linear_system_lst.append(eq);
        cout << "Solving it...\n";
        for (ex eq : lsolve(linear_system_lst, unknowns))
            if (eq.lhs() != eq.rhs()) // not a tautology
                cout << eq << endl;
    }
}
