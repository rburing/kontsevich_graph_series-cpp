#ifndef INCLUDED_POISSON_STRUCTURE_
#define INCLUDED_POISSON_STRUCTURE_

#include <vector>
#include <ginac/ginac.h>

struct PoissonStructure
{
    std::vector<GiNaC::symbol> coordinates;
    std::vector< std::vector<GiNaC::ex> > bivector;
    // used in poisson_make_vanish:
    enum class Type { Generic, Polynomial, Particular };
    Type type;
};

#endif
