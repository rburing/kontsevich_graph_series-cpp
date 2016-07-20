#ifndef INCLUDED_POISSON_STRUCTURE_
#define INCLUDED_POISSON_STRUCTURE_

#include <ginac/ginac.h>
#include <map>
#include <string>

struct PoissonStructure
{
    std::vector<GiNaC::symbol> coordinates;
    std::vector< std::vector<GiNaC::ex> > bivector;
};

/* Examples: */

DECLARE_FUNCTION_3P(phi)
REGISTER_FUNCTION(phi, dummy())
DECLARE_FUNCTION_3P(u)
REGISTER_FUNCTION(u, dummy())
GiNaC::symbol x("x"), y("y"), z("z");

std::map<std::string, PoissonStructure> poisson_structures {
    {"3d-generic", { { x, y, z }, { {0, u(x,y,z)*phi(x,y,z).diff(z), -u(x,y,z)*phi(x,y,z).diff(y)},
                                    {-u(x,y,z)*phi(x,y,z).diff(z), 0, u(x,y,z)*phi(x,y,z).diff(x) },
                                    { u(x,y,z)*phi(x,y,z).diff(y), -u(x,y,z)*phi(x,y,z).diff(x), 0 } } } },
};

#endif
