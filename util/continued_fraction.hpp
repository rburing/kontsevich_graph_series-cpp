#ifndef CONTINUED_FRACTION_H_
#define CONTINUED_FRACTION_H_

#include <vector>
#include <cmath>
#include <cln/cln.h>
#include <ginac/ginac.h>

std::vector<GiNaC::ex> convergents(double number, double eps)
{
    int p1 = 1, p2 = 0, q1 = 0, q2 = 1, tmp;
    std::vector<GiNaC::ex> result;
    while (true)
    {
        int a = std::floor(number);
        tmp = p1;
        p1 = a*p1 + p2; p2 = tmp;
        tmp = q1;
        q1 = a*q1 + q2; q2 = tmp;
        result.push_back((GiNaC::ex)p1 / q1);
        if (fabs(number - std::floor(number)) < eps)
            break;
        number = 1/(number - std::floor(number));
    }
    return result;
}

GiNaC::ex best_rational_approximation(double number, double eps)
{
    std::vector<GiNaC::ex> convergents_list = convergents(number, eps);
    return convergents_list[convergents_list.size()-1];
}

#endif
