#ifndef INCLUDED_PERMUTATIONS_H_
#define INCLUDED_PERMUTATIONS_H_

#include <vector>
#include <algorithm>

template <class T>
inline int parity(std::vector<T> permutation)
{
    size_t swaps = 0;
    for (size_t i = 0; i < permutation.size(); ++i) {
        while (i != (size_t)permutation[i]) {
            ++swaps;
            std::swap(permutation[i], permutation[permutation[i]]);
        }
    }
    return swaps % 2 == 0 ? 1 : -1;
}

#endif
