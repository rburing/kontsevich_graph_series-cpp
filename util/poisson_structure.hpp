#ifndef INCLUDED_POISSON_STRUCTURE_
#define INCLUDED_POISSON_STRUCTURE_

#include <vector>
#include <map>
#include <ginac/ginac.h>

typedef std::multiset<size_t> multi_index;

struct PoissonStructure
{
    std::vector<GiNaC::symbol> coordinates;
    std::vector< std::vector<GiNaC::ex> > bivector;
    // used in poisson_make_vanish:
    enum class Type { Generic, Polynomial, Particular };
    Type type;
    std::map< std::pair< std::pair<size_t, size_t>, multi_index >, GiNaC::ex> bivector_derivatives_cache;
};

#endif
