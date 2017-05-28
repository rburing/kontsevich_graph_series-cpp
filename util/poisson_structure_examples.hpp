#ifndef INCLUDED_POISSON_STRUCTURE_EXAMPLES_
#define INCLUDED_POISSON_STRUCTURE_EXAMPLES_

#include <ginac/ginac.h>
#include <map>
#include <string>
#include "poisson_structure.hpp"

/* Examples: */

DECLARE_FUNCTION_3P(phi)
REGISTER_FUNCTION(phi, dummy())
DECLARE_FUNCTION_3P(u)
REGISTER_FUNCTION(u, dummy())
GiNaC::symbol r("r"), t("t");
GiNaC::symbol x("x"), y("y"), z("z"), w("w");
GiNaC::symbol a("a"), b("b"), l("l");
GiNaC::symbol b01("b01"), b02("b02"), b03("b03"), b12("b12"), b13("b13"), b23("b23");
GiNaC::symbol x11("x11"), x12("x12"), x13("x13"), x21("x21"), x22("x22"), x23("x23"), x31("x31"), x32("x32"), x33("x33");
GiNaC::symbol x1("x1"), x2("x2"), x3("x3"), x4("x4");
DECLARE_FUNCTION_4P(f1)
REGISTER_FUNCTION(f1, dummy())
DECLARE_FUNCTION_4P(f2)
REGISTER_FUNCTION(f2, dummy())
DECLARE_FUNCTION_4P(f3)
REGISTER_FUNCTION(f3, dummy())

std::map<std::string, PoissonStructure> poisson_structures {
    {"2d-polar",   { { r, t }, { { 0, 1/r },
                                 { -1/r, 0 } },
                     PoissonStructure::Type::Particular } },
    {"3d-generic", { { x, y, z }, { {0, u(x,y,z)*phi(x,y,z).diff(z), -u(x,y,z)*phi(x,y,z).diff(y)},
                                    {-u(x,y,z)*phi(x,y,z).diff(z), 0, u(x,y,z)*phi(x,y,z).diff(x) },
                                    { u(x,y,z)*phi(x,y,z).diff(y), -u(x,y,z)*phi(x,y,z).diff(x), 0 } },
                     PoissonStructure::Type::Generic } },
    {"3d-determinant", { { x, y, z }, { {0, phi(x,y,z).diff(z), -phi(x,y,z).diff(y)},
                                    {-phi(x,y,z).diff(z), 0, phi(x,y,z).diff(x) },
                                    { phi(x,y,z).diff(y), -phi(x,y,z).diff(x), 0 } },
                     PoissonStructure::Type::Generic } },
    {"4d-determinant", { { x1, x2, x3, x4 }, { 
{0, -(f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x3)), (f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x2)), -(f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x3)  - f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x2))},
{(f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x3)), 0, -(f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x1)), (f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x3) - f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x1))},
{-(f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x2)), (f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x1)), 0, -(f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x2) - f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x1))},
{(f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x3)  - f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x2)), -(f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x3) - f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x1)), (f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x2) - f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x1)), 0 }
                                                },
                     PoissonStructure::Type::Generic } },
    {"4d-rank2", { { x1, x2, x3, x4 }, { 
{0, -f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x3)), f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x2)), -f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x3)  - f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x2))},
{f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x3)), 0, -f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x1)), f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x3) - f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x1))},
{-f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x2)), f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x4) - f2(x1,x2,x3,x4).diff(x4)*f3(x1,x2,x3,x4).diff(x1)), 0, -f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x2) - f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x1))},
{f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x3)  - f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x2)), -f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x3) - f2(x1,x2,x3,x4).diff(x3)*f3(x1,x2,x3,x4).diff(x1)), f1(x1,x2,x3,x4)*(f2(x1,x2,x3,x4).diff(x1)*f3(x1,x2,x3,x4).diff(x2) - f2(x1,x2,x3,x4).diff(x2)*f3(x1,x2,x3,x4).diff(x1)), 0 }
                                                },
                     PoissonStructure::Type::Generic } },
    {"3d-polynomial", { { x, y, z }, { { 0, x*y*z, x*y*z},
                                       {-x*y*z, 0, x*y*z},
                                       {-x*y*z, -x*y*z, 0} },
                        PoissonStructure::Type::Polynomial } },
    {"3d-trig", { { a, b, l }, { { 0, cos(l), -sin(l)/b },
                                 {-cos(l), 0, sin(l)/a},
                                 {sin(l)/b, -sin(l)/a, 0} },
                  PoissonStructure::Type::Particular } },
    {"4d-quadratic", { { x, y, z, w }, { { 0, b01*x*y,  b02*x*z,  b03*w*x },
                                         { -b01*x*y, 0, b12*y*z, b13*w*y },
                                         { -b02*x*z, -b12*y*z, 0, b23*w*z },
                                         { -b03*w*x, -b13*w*y, -b23*w*z, 0} },
                       PoissonStructure::Type::Polynomial } },
    {"9d-rank6", { { x11, x12, x13, x21, x22, x23, x31, x32, x33 }, { 
{0, x11*x11*x12 + x12*x12*x21, x11*x11*x13 + 2*x12*x13*x21 + x13*x13*x31, x11*x11*x21 + x12*x21*x21, 2*x11*x12*x21 + 2*x12*x21*x22, 2*x11*x13*x21 + 2*x13*x21*x22 + x12*x21*x23 + x13*x23*x31, x13*x31*x31 + (x11*x11 + 2*x12*x21)*x31, 2*(x11*x12 + x12*x22)*x31 + (x12*x21 + x13*x31)*x32, 2*x13*x21*x32 + 2*x13*x31*x33 + 2*(x11*x13 + x12*x23)*x31},
{-x11*x11*x12 - x12*x12*x21, 0, 2*x12*x13*x22 - x12*x12*x23 + x13*x13*x32, 0, x12*x12*x21 + x12*x22*x22, x12*x13*x21 + 2*x13*x22*x22 + x13*x23*x32, (x12*x21 + x13*x31)*x32, x12*x12*x31 + 2*x12*x22*x32 + x13*x32*x32, x12*x13*x31 + 2*x13*x32*x33 + (2*x13*x22 + x12*x23)*x32},
{-x11*x11*x13 - 2*x12*x13*x21 - x13*x13*x31, -2*x12*x13*x22 + x12*x12*x23 - x13*x13*x32, 0, -x12*x21*x23 - x13*x23*x31, x12*x13*x21 - x13*x23*x32, x13*x13*x21 + 2*x13*x22*x23 - x12*x23*x23, 0, x12*x13*x31 + x12*x23*x32, x13*x13*x31 + 2*x13*x23*x32 + x13*x33*x33},
{-x11*x11*x21 - x12*x21*x21, 0, x12*x21*x23 + x13*x23*x31, 0, x12*x21*x21 + x21*x22*x22, x13*x21*x21 + 2*x21*x22*x23 + x23*x23*x31, 2*x21*x22*x31 + x23*x31*x31 - x21*x21*x32, x23*x31*x32 + (x12*x21 + 2*x22*x22)*x31, x21*x23*x32 + 2*x23*x31*x33 + (x13*x21 + 2*x22*x23)*x31},
{-2*x11*x12*x21 - 2*x12*x21*x22, -x12*x12*x21 - x12*x22*x22, -x12*x13*x21 + x13*x23*x32, -x12*x21*x21 - x21*x22*x22, 0, x22*x22*x23 + x23*x23*x32, -x12*x21*x31 + x23*x31*x32, x22*x22*x32 + x23*x32*x32, 2*x22*x23*x32 + 2*x23*x32*x33},
{-2*x11*x13*x21 - 2*x13*x21*x22 - x12*x21*x23 - x13*x23*x31, -x12*x13*x21 - 2*x13*x22*x22 - x13*x23*x32, -x13*x13*x21 - 2*x13*x22*x23 + x12*x23*x23, -x13*x21*x21 - 2*x21*x22*x23 - x23*x23*x31, -x22*x22*x23 - x23*x23*x32, 0, -x13*x21*x31 - x21*x23*x32, 0, x23*x23*x32 + x23*x33*x33},
{-x13*x31*x31 - (x11*x11 + 2*x12*x21)*x31, -(x12*x21 + x13*x31)*x32, 0, -2*x21*x22*x31 - x23*x31*x31 + x21*x21*x32, x12*x21*x31 - x23*x31*x32, x13*x21*x31 + x21*x23*x32, 0, x12*x31*x31 + 2*x22*x31*x32 - x21*x32*x32, x13*x31*x31 + 2*x23*x31*x32 + x31*x33*x33},
{-2*(x11*x12 + x12*x22)*x31 - (x12*x21 + x13*x31)*x32, -x12*x12*x31 - 2*x12*x22*x32 - x13*x32*x32, -x12*x13*x31 - x12*x23*x32, -x23*x31*x32 - (x12*x21 + 2*x22*x22)*x31, -x22*x22*x32 - x23*x32*x32, 0, -x12*x31*x31 - 2*x22*x31*x32 + x21*x32*x32, 0, x23*x32*x32 + x32*x33*x33},
{-2*x13*x21*x32 - 2*x13*x31*x33 - 2*(x11*x13 + x12*x23)*x31, -x12*x13*x31 - 2*x13*x32*x33 - (2*x13*x22 + x12*x23)*x32, -x13*x13*x31 - 2*x13*x23*x32 - x13*x33*x33, -x21*x23*x32 - 2*x23*x31*x33 - (x13*x21 + 2*x22*x23)*x31, -2*x22*x23*x32 - 2*x23*x32*x33, -x23*x23*x32 - x23*x33*x33, -x13*x31*x31 - 2*x23*x31*x32 - x31*x33*x33, -x23*x32*x32 - x32*x33*x33, 0}
                                                                }, PoissonStructure::Type::Particular } },
};

#endif
