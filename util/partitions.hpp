#ifndef INCLUDED_PARTITIONS_H_
#define INCLUDED_PARTITIONS_H_

#include <vector>
#include <algorithm>

class Partitions : public std::iterator<std::input_iterator_tag, std::vector<size_t> >
{
    size_t d_number;
    std::vector<size_t> d_partition;
    size_t d_parts;

    public:
    Partitions(size_t number)
    : d_number(number), d_partition(number+2), d_parts(1)
    {
        d_partition[0] = 0;
        d_partition[1] = number;
        ++(*this);
    }

    Partitions begin() const
    {
        return Partitions(d_number);
    }

    Partitions end() const
    {
        Partitions end = *this;
        end.d_parts = 0;
        return end;
    }

    bool operator==(const Partitions& rhs) const
    {
        return d_number == rhs.d_number && d_parts == rhs.d_parts;
    }

    bool operator!=(const Partitions& rhs) const
    {
        return !(*this == rhs);
    }

    Partitions& operator++()
    {
        size_t x = d_partition[d_parts - 1] + 1;
        size_t y = d_partition[d_parts] - 1;
        --d_parts;
        while (x <= y)
        {
            d_partition[d_parts] = x;
            y -= x;
            ++d_parts;
        }
        d_partition[d_parts] = x + y;
        return *this;
    }

    std::vector<size_t> operator*()
    {
        std::vector<size_t> partition(d_partition.begin(), d_partition.begin() + d_parts + 1);
        return partition;
    }
};

#endif
