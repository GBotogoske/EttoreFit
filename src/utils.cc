#include "utils.hh"
#include <algorithm>
#include <iostream>

void roll_vector(std::vector<double>& v, int shift)
{
    int n = v.size();
    if (n == 0) return;

    shift = ((shift % n) + n) % n;

    std::rotate(v.begin(), v.end() - shift, v.end());
}