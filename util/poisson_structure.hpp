#ifndef INCLUDED_POISSON_STRUCTURE_
#define INCLUDED_POISSON_STRUCTURE_

#include <ginac/ginac.h>
#include <map>
#include <string>

struct PoissonStructure
{
    std::vector<GiNaC::symbol> coordinates;
    std::vector< std::vector<GiNaC::ex> > bivector;
    enum class Type { Generic, Polynomial, Particular };
    Type type;
};

/* Examples: */

DECLARE_FUNCTION_3P(phi)
REGISTER_FUNCTION(phi, dummy())
DECLARE_FUNCTION_3P(u)
REGISTER_FUNCTION(u, dummy())
GiNaC::symbol r("r"), t("t");
GiNaC::symbol x("x"), y("y"), z("z");
GiNaC::symbol a("a"), b("b"), l("l");

std::map<std::string, PoissonStructure> poisson_structures {
    {"2d-polar",   { { r, t }, { { 0, 1/r },
                                 { -1/r, 0 } },
                     PoissonStructure::Type::Particular } },
    {"3d-generic", { { x, y, z }, { {0, u(x,y,z)*phi(x,y,z).diff(z), -u(x,y,z)*phi(x,y,z).diff(y)},
                                    {-u(x,y,z)*phi(x,y,z).diff(z), 0, u(x,y,z)*phi(x,y,z).diff(x) },
                                    { u(x,y,z)*phi(x,y,z).diff(y), -u(x,y,z)*phi(x,y,z).diff(x), 0 } },
                     PoissonStructure::Type::Generic } },
    {"3d-polynomial", { { x, y, z }, { { 0, x*y*z, x*y*z},
                                       {-x*y*z, 0, x*y*z},
                                       {-x*y*z, -x*y*z, 0} },
                        PoissonStructure::Type::Polynomial } },
    {"3d-trig", { { a, b, l }, { { 0, cos(l), -sin(l)/b },
                                 {-cos(l), 0, sin(l)/a},
                                 {sin(l)/b, -sin(l)/a, 0} },
                  PoissonStructure::Type::Particular } },
};

#endif
