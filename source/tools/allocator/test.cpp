#include <iostream>
#include <memory>
#include <cstdlib>

using namespace std;

static inline bool is_aligned(const void *__restrict__ pointer, size_t byte_count)
{ return (uintptr_t)pointer % byte_count == 0; }

struct aligned_memory_header
{
	size_t offset;
	size_t allocated_size;
};

int main()
{

    size_t alignment = 32;
	size_t size = sizeof(double) * 1021 * 1023;
	cout <<"Allocating " << size+ sizeof(aligned_memory_header) + alignment
	    << " bytes.\n";
    void * ptr = malloc(size + sizeof(aligned_memory_header) + alignment);
    //void * ptr;
    //int rc = posix_memalign(&ptr, alignment, 100);
    //if (rc != 0) ptr = nullptr;

    cout << "ptr:"<<ptr <<endl;
    size_t value = reinterpret_cast<std::size_t>(ptr)
        + sizeof(aligned_memory_header);
    size_t offset = (value+alignment-1) & ~(alignment-1);
    size_t offset2 = offset - reinterpret_cast<size_t>(ptr);

    aligned_memory_header* header = reinterpret_cast<aligned_memory_header*>
        (static_cast<char*>(ptr) + offset2) - 1;
    cout << "header:"<< header << endl;
    header->offset = offset2;
    header->allocated_size = size +sizeof(aligned_memory_header) + alignment;

    cout << "value:" << value << " offset:" << offset << endl;
    cout << "offset2:" << offset2 << endl;
    cout << "header_offset:" << header->offset << endl;
    cout << "alloc_sz:" << header->allocated_size << endl;

    if (ptr) free(ptr);

    if (is_aligned(ptr, alignment))
        cout << "memory aligned" << endl;
    else
        cout << "memory not aligned" << endl;

    void * newptr = static_cast<char*>(ptr) + offset2;
    cout << "and the shifted memory?" << endl;
    if (is_aligned(newptr, alignment))
        cout << "memory aligned" << endl;
    else
        cout << "memory not aligned" << endl;


    return 0;
}
