
#include "type_info.h"

#if defined(__GNUC__)

#include <cxxabi.h>

std::string demangle(const char * name)
{
    int status = -4;
    std::unique_ptr<char, void(*)(void*)> res {
        abi::__cxa_demangle(name, nullptr, nullptr, &status),
        std::free
    };
    return (status==0) ? res.get() : name;
}

#else

std::string demangle(const char * name)
{
    return name;
}

#endif
