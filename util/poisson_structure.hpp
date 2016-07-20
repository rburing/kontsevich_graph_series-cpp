#ifndef INCLUDED_POISSON_STRUCTURE_
#define INCLUDED_POISSON_STRUCTURE_

#include <ginac/ginac.h>

struct PoissonStructure
{
    std::vector<GiNaC::symbol> coordinates;
    std::vector< std::vector<GiNaC::ex> > bivector;
};

#endif
