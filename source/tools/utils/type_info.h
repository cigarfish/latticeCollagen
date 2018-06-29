#ifndef TYPE_INFO_H
#define TYPE_INFO_H

#include <string>
#include <typeinfo>
#include <memory>
#include <cstdlib>

std::string demangle(const char * name);

template <typename... Args>
struct Impl;

template <typename First, typename... Args>
struct Impl<First, Args...>
{
    static std::string name()
    {
        return std::string(demangle(typeid(First).name())) + ", "
            + Impl<Args...>::name();
    }
};

template<>
struct Impl<>
{
    static std::string name()
    {
        return "";
    }
};

template<typename... Args>
std::string type_name()
{
    return Impl<Args...>::name().c_str();
}

#endif
