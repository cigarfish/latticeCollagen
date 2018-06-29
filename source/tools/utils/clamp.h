#ifndef CLAMP_H
#define CLAMP_H

#if __cplusplus <= 201402L

#include <algorithm>
#include <functional>

// from the C++17 standard

namespace std {

template<class T>
    constexpr const T& clamp( const T& v, const T& lo, const T& hi )
    {
        return clamp( v, lo, hi, std::less<>() );
    }

template<class T, class Compare>
    constexpr const T& clamp( const T& v, const T& lo, const T& hi, Compare comp )
    {
        return assert( !comp(hi, lo) ),
               comp(v, lo) ? lo : comp(hi, v) ? hi : v;
    }

} // end namespace

#endif

#endif
