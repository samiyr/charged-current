#ifndef UTILITY_H
#define UTILITY_H

#include <vector>

#define POW2(x) (x) * (x)
#define POW4(x) (x) * (x) * (x) * (x)


template <typename T>
constexpr std::vector<T> vector_intersection(std::vector<T> v1, std::vector<T> v2) {
    std::vector<T> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(),v1.end(), v2.begin(),v2.end(), back_inserter(v3));
    return v3;
}

#endif