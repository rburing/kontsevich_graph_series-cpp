#ifndef INCLUDED_SORT_PAIRS_H_
#define INCLUDED_SORT_PAIRS_H_

inline std::pair<size_t, size_t> exchange_pair(std::pair<size_t, size_t> p)
{
    return { p.second, p.first };
}

inline size_t sort_pairs(std::vector< std::pair<size_t, size_t> >& pairs)
{
    size_t exchanges = 0;
    for (auto pair = pairs.begin(); pair != pairs.end(); pair++)
    {
        std::pair<size_t, size_t> exchanged = exchange_pair(*pair);
        if (exchanged < *pair)
        {
            *pair = exchanged;
            exchanges++;
        }
    }
    return exchanges;
}

#endif
