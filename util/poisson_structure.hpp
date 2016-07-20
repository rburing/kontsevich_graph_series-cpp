#ifndef INCLUDED_POISSON_STRUCTURE_
#define INCLUDED_POISSON_STRUCTURE_

#include <ginac/ginac.h>

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

std::vector<PoissonStructure> poisson_structures {
    { { x, y, z }, { {0, u(x,y,z)*phi(x,y,z).diff(z), -u(x,y,z)*phi(x,y,z).diff(y)},
                     {-u(x,y,z)*phi(x,y,z).diff(z), 0, u(x,y,z)*phi(x,y,z).diff(x) },
                     { u(x,y,z)*phi(x,y,z).diff(y), -u(x,y,z)*phi(x,y,z).diff(x), 0 } } },
};

#endif
