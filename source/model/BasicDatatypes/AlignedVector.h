////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  AlignedVector.h                                               //
//                                                                            //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                //
//    Created:  2016-04-20                                                    //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIGNED_VECTOR_H
#define ALIGNED_VECTOR_H

#include <vector>
#include <ostream>

#if defined(_MSC_VER)
    #include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
    #include <x86intrin.h>
#endif

#if defined(_MSC_VER)
    #define ALIGNED_(x) __declspec(align(x))
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
    #define ALIGNED_(x) __attribute__ ((aligned(x)))
#endif

#define ALIGNED_TYPE_(t, x) typedef ALIGNED_(x) t

#include "tools/allocator/AlignmentAllocator.h"

using Memory::AlignedAllocator;
using Memory::Alignment;

#ifdef __AVX__
    #ifdef DEBUG
        #pragma message "Selecting AVX 32 byte alignment"
    #endif

    static constexpr std::size_t alignment = to_integral(Alignment::AVX);

    template<typename T>
    using AligndVector = typename std::vector<T, AlignedAllocator<T, Alignment::AVX>>;

#elif __SSE__
    #ifdef DEBUG
        #pragma message "Selecting SSE 16 byte alignment"
    #endif

    template<typename T>
    using AligndVector = typename std::vector<T, AlignedAllocator<T, Alignment::SSE>>;

    static constexpr std::size_t alignment = to_integral(Alignment::SSE);

#else
    #ifdef DEBUG
        #pragma message "Selecting no forced alignment"
    #endif

    static constexpr std::size_t alignment = to_integral(Alignment::Normal);

    template<typename T>
    using AligndVector = typename std::vector<T>;

#endif // end alignment switch

// Create aligned types
ALIGNED_TYPE_(double,           alignment) align_double;
ALIGNED_TYPE_(float,            alignment) align_float;
ALIGNED_TYPE_(int,              alignment) align_int;
ALIGNED_TYPE_(unsigned int,     alignment) align_uint;
ALIGNED_TYPE_(long,             alignment) align_long;
ALIGNED_TYPE_(unsigned long,    alignment) align_ulong;

ALIGNED_TYPE_(double,           alignment) aligned_double;
ALIGNED_TYPE_(float,            alignment) aligned_float;
ALIGNED_TYPE_(int,              alignment) aligned_int;
ALIGNED_TYPE_(unsigned int,     alignment) aligned_uint;
ALIGNED_TYPE_(long,             alignment) aligned_long;
ALIGNED_TYPE_(unsigned long,    alignment) aligned_ulong;

template<typename T>
std::ostream& operator<<(std::ostream& os, const AligndVector<T>& vec)
{
    const std::size_t n {vec.size()};
    os << "[";
    if (n == 0)
    {
        os << "]" << std::endl;
        return os;
    }
    for (std::size_t i = 0; i < n-1; i++)
    {
       os << vec[i] << ", ";
    }
    os << vec[n-1] << "]" << std::endl;
    return os;
}

#endif
