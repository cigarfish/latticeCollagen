
#pragma once

#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include "tools/model/CSModelTools.h"

using std::string;

template <class T>
struct TypeInfo
{
    using Return    = const T&;
    using SReturn   = T;
    using Parameter = const T&;
    using Reference = T&;
    static SReturn fromString(const string& val)
    {
        std::stringstream s(val);
        T ret;
        s >> ret;
        // this seems to not give reliable data when loading vectors?
        //if (s.fail())
        //{
        //    throw string("Unable to read value from string \'") + val + "'";
        //}
        return ret;
    }
    static const string name();
};


template<>
struct TypeInfo<double>
{
    using Return    = double;
    using SReturn   = double;
    using Parameter = double;
    using Reference = double&;
    static SReturn fromString(const string& val)
    {
        std::stringstream s(val);
        double ret;
        s >> ret;
        if (s.fail())
        {
            throw string("Unable to read value from string \'") + val + "'";
        }
        return ret;
    }
    static const string name()
    {
        return "Double";
    }
};


template<>
struct TypeInfo<float>
{
    using Return    = float;
    using SReturn   = float;
    using Parameter = float;
    using Reference = float&;
    static SReturn fromString(const string& val)
    {
        std::stringstream s(val);
        float ret;
        s >> ret;
        if (s.fail())
        {
            throw string("Unable to read value from string\' ") + val + "'";
        }
        return ret;
    }
    static const string name()
    {
        return "Float";
    }
};


template<>
struct TypeInfo<bool>
{
    using Return    = bool;
    using SReturn   = bool;
    using Parameter = bool;
    using Reference = bool&;
    static SReturn fromString(string val)
    {
        lower_case(val);
        if ( val == "true" )
            return true;
        else if ( val == "false")
            return false;
        else
            throw string("Unable to read value from string\' ") + val + "'";
    }
    static const string name()
    {
        return "Bool";
    }
};


template<>
struct TypeInfo<string>
{
    using Return    = const string&;
    using SReturn   = string;
    using Parameter = const string&;
    using Reference = string&;
    static SReturn fromString(const string& val) { return val; }
    static const string name()
    {
        return "String";
    }
};

