#ifndef INCLUDED_SORT_PAIRS_H_
#define INCLUDED_SORT_PAIRS_H_

#include <vector>

inline std::pair<char, char> exchange_pair(std::pair<char, char> p)
{
    return { p.second, p.first };
}

inline size_t sort_pairs(std::vector< std::pair<char, char> >::iterator pairs_begin, std::vector< std::pair<char, char> >::iterator pairs_end)
{
    size_t exchanges = 0;
    for (auto pair = pairs_begin; pair != pairs_end; pair++)
    {
        std::pair<char, char> exchanged = exchange_pair(*pair);
        if (exchanged < *pair)
        {
            *pair = exchanged;
            exchanges++;
        }
    }
    return exchanges;
}

#endif
