#ifndef INCLUDED_CARTESIAN_PRODUCT_H_
#define INCLUDED_CARTESIAN_PRODUCT_H_

#include <vector>
#include <algorithm>

class CartesianProduct : public std::iterator<std::input_iterator_tag, std::vector<size_t> >
{
    std::vector<size_t>& d_ends;
    long long d_length;
    long long d_position;

    public:
    CartesianProduct(std::vector<size_t>& ends) : d_ends(ends)
    {
      auto product = [](long long a, size_t b) { return a*b; };
      d_length = accumulate(ends.begin(), ends.end(), 1LL, product);
      d_position = 0;
    }

    CartesianProduct(std::vector<size_t>& ends, long long length, long long position)
    : d_ends(ends), d_length(length), d_position(position)
    {}

    CartesianProduct begin()
    {
        return CartesianProduct(d_ends, d_length, 0);
    }

    CartesianProduct end()
    {
        return CartesianProduct(d_ends, d_length, d_length);
    }

    bool operator==(const CartesianProduct& rhs)
    {
        return d_position == rhs.d_position;
    }

    bool operator!=(const CartesianProduct& rhs)
    {
        return d_position != rhs.d_position;
    }

    CartesianProduct& operator++()
    {
        d_position++;
        return *this;
    }

    CartesianProduct& operator++(int)
    {
        d_position++;
        return *this;
    }

    std::vector<size_t> operator*()
    {
        std::vector<size_t> u(d_ends.size());
        lldiv_t q { d_position, 0 };
        for( long long i=d_ends.size()-1 ; 0<=i ; --i ) {
          q = std::div( q.quot, d_ends[i] );
          u[i] = q.rem;
        }
        return u;
    }
};

#endif
