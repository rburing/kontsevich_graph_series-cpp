#include <vector>
#include <utility>
#include <cstddef>
#include <ostream>

class KontsevichGraph
{
    size_t d_internal;
    size_t d_external;
    std::vector< std::pair<size_t, size_t> > d_targets;

    public:
    KontsevichGraph(size_t internal, size_t external, std::vector< std::pair<size_t, size_t> > targets);

    private:
    friend std::ostream& operator<<(std::ostream &os, const KontsevichGraph& g);
};
