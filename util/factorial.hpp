#ifndef INCLUDED_FACTORIAL_H_
#define INCLUDED_FACTORIAL_H_

inline size_t factorial(size_t n)
{
    size_t factorial = 1;
    for (size_t i = n; i != 1; --i)
        factorial *= i;
    return factorial;
}

#endif
