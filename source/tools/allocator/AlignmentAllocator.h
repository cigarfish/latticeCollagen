////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  AlignmentAllocator.h                                          //
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

#ifndef ALIGNED_ALLOCATOR_H
#define ALIGNED_ALLOCATOR_H

//
// TODO: This header file currently does two different things! That needs to be
// sorted out.
//

#include <stdlib.h>
#include <malloc.h>
#include <memory>
#include <type_traits>

#include "tools/concepts/concepts.h"

#define MALLOC __attribute__((malloc))
/*
#if __has_attribute(alloc_size)
    #define ALLOC_SIZE(args...) __attribute__((alloc_size(args)))
#else
    #define ALLOC_SIZE(args...)
#endif
*/
#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
    #define RESTRICT __restrict__
#endif
/*
#if __has_attribute(alloc_align) || (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 9)))
    #define ALLOC_ALIGN(arg) __attribute__((alloc_align(arg)))
#else
    #define ALLOC_ALIGN(arg)
#endif

#if (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))) || defined(__clang__)
    #define ASSUME_ALIGNED(ptr, alignment) __builtin_assume_aligned(ptr, alignment)
// TODO check for icc
//#elif __has_attribute(__assume_aligned)
//    #define ASSUME_ALIGNED(ptr, alignment) __assume_aligned(ptr, alignment)
#else
    #define ASSUME_ALIGNED(ptr, alignment)
#endif
*/
using std::size_t;

namespace Memory {

static inline
bool is_aligned(const void * RESTRICT ptr, size_t alignment)
{
    return (uintptr_t)ptr % alignment == 0;
}

enum class Alignment : short
{
    Normal  = sizeof(void*),
    SSE     = 16,
    AVX     = 32,
};

void * aligned_alloc(size_t alignment, size_t size);// ALLOC_ALIGN(1) ALLOC_SIZE(2) MALLOC;
void * aligned_realloc(void * ptr, size_t alignment, size_t size);// ALLOC_ALIGN(2) ALLOC_SIZE(3);
void aligned_free(void * ptr) noexcept;

void* allocate_aligned_memory(size_t align, size_t size);// ALLOC_ALIGN(1) ALLOC_SIZE(2) MALLOC;
void deallocate_aligned_memory(void* ptr) noexcept;
size_t align_on(size_t value, size_t alignment) noexcept;

struct aligned_memory_header
{
    size_t offset;
    size_t allocated_size;
};

template<typename T, Alignment Align = Alignment::AVX>
class AlignedAllocator;

template<Alignment Align>
class AlignedAllocator<void, Align>
{
    public:
        using pointer       = void*;
        using const_pointer = const void*;
        using value_type    = void;

        template <class U> struct rebind { typedef AlignedAllocator<U, Align> other; };
};

template<typename T, Alignment Align>
class AlignedAllocator
{
    public:
        using value_type        = T;
        using pointer           = T*;
        using const_pointer     = const T*;
        using reference         = T&;
        using const_reference   = const T&;
        using size_type         = size_t;
        using difference_type   = ptrdiff_t;

        using propagate_on_container_move_assignment = std::true_type;

        template<class U>
        struct rebind { using other = AlignedAllocator<U, Align>; };

    public:
        AlignedAllocator() noexcept
        {}

        template<class U>
        AlignedAllocator(const AlignedAllocator<U, Align>&) noexcept
        {}

        size_type max_size() const noexcept
        { return (size_type(~0) - size_type(Align)) / sizeof(T); }

        pointer address(reference x) const noexcept
        { return std::addressof(x); }

        const_pointer address(const_reference x) const noexcept
        { return std::addressof(x); }

        pointer allocate(size_type n, typename AlignedAllocator<void, Align>::const_pointer = 0)
        {
            const size_type alignment = static_cast<size_type>( Align );
            void * ptr = allocate_aligned_memory(alignment, n * sizeof(T));
            if (ptr == nullptr)
                throw std::bad_alloc();

            return reinterpret_cast<pointer>(ptr);
        }

        void deallocate(pointer p, size_type) noexcept
        { return deallocate_aligned_memory(p);}

        template<class U, class ...Args>
        void construct(U* p, Args&&... args)
        { ::new(reinterpret_cast<void*>(p)) U(std::forward<Args>(args)...); }

        void destroy(pointer p)
        { p->~T(); }
};

template<typename T, Alignment Align>
class AlignedAllocator<const T, Align>
{
    public:
        using value_type        = T;
        using pointer           = const T*;
        using const_pointer     = const T*;
        using reference         = const T&;
        using const_reference   = const T&;
        using size_type         = size_t;
        using difference_type   = ptrdiff_t;

        using propagate_on_container_move_assignment = std::true_type;

        template<class U>
        struct rebind { using other = AlignedAllocator<U, Align>; };

    public:
        AlignedAllocator() noexcept
        {}

        template<class U>
        AlignedAllocator(const AlignedAllocator<U, Align>&) noexcept
        {}

        size_type max_size() const noexcept
        { return (size_type(~0) - size_type(Align)) / sizeof(T); }

        const_pointer address(const_reference x) const noexcept
        { return std::addressof(x); }

        pointer allocate(size_type n, typename AlignedAllocator<void, Align>::const_pointer = 0)
        {
            const size_type alignment = static_cast<size_type>( Align );
            void * ptr = allocate_aligned_memory(alignment, n * sizeof(T));
            if (ptr == nullptr)
                throw std::bad_alloc();

            return reinterpret_cast<pointer>(ptr);
        }

        void deallocate(pointer p, size_type) noexcept
        { return deallocate_aligned_memory(p);}

        template<class U, class ...Args>
        void construct(U* p, Args&&... args)
        { ::new(reinterpret_cast<void*>(p)) U(std::forward<Args>(args)...); }

        void destroy(pointer p)
        { p->~T(); }
};


template<typename T, Alignment TAlign, typename U, Alignment UAlign>
inline bool operator==(const AlignedAllocator<T, TAlign>&, const AlignedAllocator<U, UAlign>&) noexcept
{ return TAlign == UAlign; }

template<typename T, Alignment TAlign, typename U, Alignment UAlign>
inline bool operator!=(const AlignedAllocator<T, TAlign>&, const AlignedAllocator<U, UAlign>&) noexcept
{ return TAlign != UAlign; }

} // end namespace Memory

#endif
