
#include "AlignmentAllocator.h"

#include <cassert>
#include <cstring>
#include <iostream>

using namespace std;

using std::size_t;

namespace Memory {

    template<typename T>
    constexpr bool is_power_of_two(T x)
    {
        return x && ((x & (x-1)) == 0);
    }

    size_t align_on(size_t value, size_t alignment) noexcept
    {
        return (value+alignment-1) & ~(alignment-1);
    }

    void* allocate_aligned_memory(size_t align, size_t size)
    {
        assert(align >= sizeof(void*));
        assert(is_power_of_two(align));

        if (size == 0)
            return nullptr;

        void* ptr = nullptr;
        int rc = posix_memalign(&ptr, align, size);

        if (rc != 0)
            return nullptr;

        return ptr;
    }

    void deallocate_aligned_memory(void * ptr) noexcept
    {
        return free(ptr);
    }

    // these functions implement an aligned allocation
    // with support for reallocation
    //
    // this always wastes the size of the alignment + size of the header
    // but it is fast
    //
    // These functions are based on ideas in eigen and netlib
    //
    void * aligned_alloc(size_t alignment, size_t size)
    {
        assert(alignment >= sizeof(void*));
        assert(is_power_of_two(alignment));

        if (size == 0)
            return nullptr;

        void * ptr = malloc(size + sizeof(aligned_memory_header) + alignment);
        if (ptr == nullptr)
            return nullptr;

        size_t offset = align_on(reinterpret_cast<size_t>(ptr) +
                             sizeof(aligned_memory_header), alignment)
                            - reinterpret_cast<size_t>(ptr);

        aligned_memory_header * header = reinterpret_cast<aligned_memory_header*>(
            static_cast<char*>(ptr) + offset) - 1;

        header->offset = offset;
        header->allocated_size = size + sizeof(aligned_memory_header) + alignment;

        return static_cast<char*>(ptr) + offset;
    }

    void * aligned_realloc(void * ptr, size_t alignment, size_t size)
    {
        assert(alignment >= sizeof(void*));
        assert(is_power_of_two(alignment));

        aligned_memory_header* header = nullptr;

        if (ptr == nullptr)
            return aligned_alloc(alignment, size);

        header = reinterpret_cast<aligned_memory_header*>(ptr) - 1;

        size_t new_size = size ? size + sizeof(aligned_memory_header)
            + alignment : 0u;

        size_t old_size = header ? header->allocated_size : 0u;
        size_t old_offset = header ? header->offset : 0u;

        if (new_size == old_size)
            return ptr;

        if (new_size == 0)
        {
            aligned_free(ptr);
            return nullptr;
        }

        void * new_ptr = realloc(static_cast<char*>(ptr) - old_offset, new_size);
        if (new_ptr == nullptr)
            return nullptr;

        size_t offset = align_on(reinterpret_cast<size_t>(new_ptr) +
                                 sizeof(aligned_memory_header), alignment)
                      - reinterpret_cast<size_t>(new_ptr);

        if (offset != old_offset)
            memmove(static_cast<char*>(new_ptr) + offset,
                    static_cast<char*>(new_ptr) + old_offset,
                    std::min(size, old_size));

        header = reinterpret_cast<aligned_memory_header*>(static_cast<char*>(new_ptr)
                                                          + offset) - 1;
        header->offset = offset;
        header->allocated_size = new_size + sizeof(aligned_memory_header)
            + alignment;

        return static_cast<char*>(new_ptr) + offset;
    }

    void aligned_free(void * ptr) noexcept
    {
        aligned_memory_header * header;

        if (ptr == nullptr)
            return;

        header = static_cast<aligned_memory_header*>(ptr) - 1;
        size_t offset = header->offset;

        free(static_cast<char*>(ptr) - offset);
    }


} // end namespace Memory
