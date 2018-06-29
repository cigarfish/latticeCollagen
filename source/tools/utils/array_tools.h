#ifndef ARRAY_TOOLS_H
#define ARRAY_TOOLS_H

#include <array>
#include <utility>
#include <iostream>
#include <iterator>

template <class T, size_t... Is, size_t N>
constexpr std::array<T, N> multiply(std::array<T, N> const &src,
                                    std::index_sequence<Is...>,
                                    T const& mul)
{
    return std::array<T, N>{{(src[Is] * mul)...}};
}

template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr)
{
    std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(os, " "));
    return os;
}

#endif
