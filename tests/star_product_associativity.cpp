#include "../kontsevich_graph_series.hpp"
#include "../kontsevich_graph_operator.hpp"
#include "../util/partitions.hpp"
#include "../util/factorial.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits>
#include <sstream>
using namespace std;
using namespace GiNaC;

DECLARE_FUNCTION_3P(phi)
REGISTER_FUNCTION(phi, dummy())
DECLARE_FUNCTION_3P(u)
REGISTER_FUNCTION(u, dummy())

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " <star-product-filename>\n";
        return 1;
    }

    // Reading in star product:
    string star_product_filename(argv[1]);
    ifstream star_product_file(star_product_filename);
    parser coefficient_reader;
    KontsevichGraphSeries<ex> star_product;
    KontsevichGraphSum<ex> term;
    size_t order = 0;
    for (string line; getline(star_product_file, line); )
    {
        if (line.length() == 0)
            continue;
        if (line[0] == 'h')
        {
            star_product[order] = term;
            term = KontsevichGraphSum<ex>({ });
            order = stoi(line.substr(2));
        }
        else
        {
            KontsevichGraph graph;
            stringstream ss(line);
            ss >> graph;
            string coefficient_str;
            ss >> coefficient_str;
            ex coefficient = coefficient_reader(coefficient_str);
            term += KontsevichGraphSum<ex>({ { coefficient, graph } });
        }
    }
    star_product[order] = term; // the last one

    // Make a list of weight variables
    lst weight_vars;
    for (std::pair<string, ex> pair : coefficient_reader.get_syms())
        weight_vars.append(pair.second);

    star_product.reduce();

    cout << "Number of terms in star product:\n";
    for (size_t n = 0; n <= order; ++n)
    {
        cout << "h^" << n << ": " << star_product[n].size() << " total\n";
        for (std::vector<size_t> indegrees : star_product[n].in_degrees(true))
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
    set<ex, ex_is_less> weight_system;
    for (size_t n = 0; n <= order; ++n)
    {
        cout << "h^" << n << ":\n";
        cout << assoc[n].size() << " total\n";
        for (std::vector<size_t> indegrees : assoc[n].in_degrees(true))
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
