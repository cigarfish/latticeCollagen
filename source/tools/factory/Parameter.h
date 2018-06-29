
#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <cstring>

#include "Types.h"
#include "BiologyLink.h"
#include "tools/utils/exceptions.h"
//#include "xmlParser/xmlParser.h"

using std::string;
using namespace std;

class RequiredPolicy
{
public:
    RequiredPolicy() {}
    bool isRequired() const { return true; }
    bool isMissing() const { return false; }
    const string& stringVal() const { return string_val; }
    const string& unitVal() const { return unit_val; }
protected:
    void assertDefined() const {}
    void setMissing() const { throw string("Required parameter is not set"); }
    void setStringVal (const string& val) { string_val = val; }
    void setUnitVal (const string& unit) { unit_val = unit; }
    ~RequiredPolicy() {}
private:
    string string_val, unit_val;
};

class OptionalPolicy
{
public:
    OptionalPolicy() {}
    bool isRequired() const { return false; }
    bool isMissing() const { return is_missing; }
    const string& stringVal() const { return string_val; }
    const string& unitVal() const { return unit_val; }
protected:
    void assertDefined() const
    {
        if (is_missing)
            throw string("Optional parameter queried although it is not defined!");
    }
    void setMissing() { is_missing=true; }
    void setStringVal (const string& val) { string_val = val; is_missing = false; }
    void setUnitVal (const string& unit) { unit_val = unit; }
    ~OptionalPolicy() {}
private:
    bool is_missing;
    string string_val, unit_val;
};

class DefaultValPolicy
{
public:
    DefaultValPolicy() : default_defined(false) {}
    bool isRequired() const { return false; }
    bool isMissing() const { return false; }
    const string& stringVal() const { return string_val; }
    const string& unitVal() const { return unit_val; }
    void setDefault(string val) { default_value = val; default_defined = true; }
    // this needs to be public so we can set a default unit!
    void setUnitVal (const string& unit) { unit_val = unit; }
protected:
    void assertDefined() const {}
    void setMissing()
    {
        if (!default_defined)
            throw string("No default value provided");
        string_val = default_value;
    }
    void setStringVal (const string& val) { string_val = val; }
    ~DefaultValPolicy() {}
private:
    bool default_defined;
    string default_value, string_val, unit_val;
};

class ParameterBase
{
public:
    virtual void loadFromXML(XMLNode node) = 0;
    virtual void init() = 0;
    virtual void read(string value) = 0;
    virtual string XMLPath() const = 0;
};

template <class ValueType, class RequirementPolicy>
class Reader : public RequirementPolicy
{
public:
    typename TypeInfo<ValueType>::SReturn operator()() const
    {
        RequirementPolicy::assertDefined();
        return const_val;
    }

    typename TypeInfo<ValueType>::SReturn get() const
    {
        RequirementPolicy::assertDefined();
        return const_val;
    }

    // HACK
    void set(typename TypeInfo<ValueType>::SReturn value)
    {
        const_val = value;
    }

protected:
    Reader() {}
    ~Reader() {}

    void read(const string& string_val)
    {
        const_val = TypeInfo<ValueType>::fromString(string_val);
    }

    void init() {}

private:
    ValueType const_val;
};

// TODO this can be merged with the previous as long as biolink understands
// that empty units are non changed!
template <class ValueType, class RequirementPolicy>
class QuantityReader : public RequirementPolicy
{
public:
    typename TypeInfo<ValueType>::SReturn operator()() const
    {
        RequirementPolicy::assertDefined();
        return const_val;
    }

    typename TypeInfo<ValueType>::SReturn get() const
    {
        RequirementPolicy::assertDefined();
        return const_val;
    }

protected:
    QuantityReader() {}
    ~QuantityReader() {}

    void read(const string& string_val)
    {
        const_val = TypeInfo<ValueType>::fromString(string_val);
        // apply biolink to transform to internal units
        auto instance = BiologyLink::Instance();
        //std::cout << "Transforming to internal " << string_val << std::endl;
        const_val = instance->scaleBiologyToInternal(const_val,
                                                     RequirementPolicy::unitVal());
    }

    void init() {}

private:
    ValueType const_val;
};

template <class ValueType, class RequirementPolicy>
class PointerQuantityReader : public RequirementPolicy
{
public:
    void setPointer(ValueType * ptr)
    {
        const_val = ptr;
    }

    typename TypeInfo<ValueType>::SReturn operator()() const
    {
        RequirementPolicy::assertDefined();
        assert(const_val);
        return * const_val;
    }

    typename TypeInfo<ValueType>::SReturn get() const
    {
        RequirementPolicy::assertDefined();
        assert(const_val);
        return * const_val;
    }

protected:
    PointerQuantityReader() {}
    ~PointerQuantityReader() {}

    void read(const string& string_val)
    {
        assert(const_val);
        (*const_val) = TypeInfo<ValueType>::fromString(string_val);
        auto instance = BiologyLink::Instance();
        (*const_val) = instance->scaleBiologyToInternal(*const_val,
                                                        RequirementPolicy::unitVal());
    }

    void init() {}

private:
    ValueType * const_val = nullptr;
};

template <class ValueType, class RequirementPolicy>
class PointerReader : public RequirementPolicy
{
public:
    void setPointer(ValueType * ptr)
    {
        const_val = ptr;
    }

    typename TypeInfo<ValueType>::SReturn operator()() const
    {
        RequirementPolicy::assertDefined();
        assert(const_val);
        return * const_val;
    }

    typename TypeInfo<ValueType>::SReturn get() const
    {
        RequirementPolicy::assertDefined();
        assert(const_val);
        return * const_val;
    }

protected:
    PointerReader() {}
    ~PointerReader() {}

    void read(const string& string_val)
    {
        assert(const_val);
        (*const_val) = TypeInfo<ValueType>::fromString(string_val);
    }

    void init() {}

private:
    ValueType * const_val = nullptr;
};

vector<string> tokenize(const string& str, const string& delimeters = " ");
XMLNode getXMLNode(const XMLNode XML_base, const string& path);
string strip_last_token(string& s, const string& del);
//string remove_last_token(string& s, const string& del);
XMLNode getXMLNode(const XMLNode XML_base, const string& path,
                   const string& value);
XMLNode getParameterNode(const XMLNode& XML_base, const string& name);
XMLNode getParameterNode(const XMLNode& XML_base, XMLNode& node,
                         const string& name);

template <class T>
bool getXMLAttribute(const XMLNode XML_base, string path, T& value,
                     bool verbose=false)
{
    //std::cout << "getXMLAttribute(" << XML_base.getName() << ", " << path << ", "
    //    << value << ", " << verbose << ")" << std::endl;
    //std::cout << "YOO" << std::endl;
    string attribute;
    attribute = strip_last_token(path,"/");
    //std::cout << "Attribute:" << attribute << std::endl;
    if (attribute.empty()) return false;

    if (verbose)
    {
        cout << "getXMLAttribute: seeking for: ";
        cout << "path:" << path << std::endl;
        if (path.length()==0) cout << XML_base.getName();
        else cout << path;
        cout << "->" << attribute;
    }

    XMLNode xNode;
    if (!path.empty())
    {
        xNode = getXMLNode(XML_base, path);
        if (verbose)
            std::cout << "1 Trying to access node with " << path << std::endl;
        //std::cout << "Node:" << xNode.getName() << std::endl;
    }
    else
    {
        xNode = XML_base;
    }

    //if (xNode.isAttributeSet(attribute.c_str()))
    //    std::cout << "Attribute is defined" << std::endl;
    //else
    //    std::cout << "Attribute is not defined" << std::endl;

    if (xNode.isEmpty())
    { if (verbose) cout << " .. not found (xNode empty)" << endl; return false; }

    string str_val;
    if (lower_case(attribute) == "text")
    {
        if (xNode.nText())
            str_val = xNode.getText();
        else
        {
            if (verbose) cout << " .. not found (text)" << endl;
            return false;
        }
    }
    else if (xNode.isAttributeSet(attribute.c_str()))
        str_val = xNode.getAttribute(attribute.c_str());
    else
    {
        if (verbose)
            cout << " .. not found (attribute)" << endl;
        return false;
    }

    //std::cout <<"str_val:" << str_val << std::endl;

    std::stringstream k(str_val);
    T tmp;
    k >> tmp;
    if (k.fail())
    {
        if (verbose)
            cout << "failure during value conversion" << endl;
        return false;
    }
    if (path.length()==0)
        path = XML_base.getName();
    if (verbose)
        cout << ": " << tmp << endl;
    value=tmp;
    return true;
}

template <>
bool getXMLAttribute<string>(const XMLNode XML_base, string path, string& value,
                             bool verbose);

template<>
bool getXMLAttribute<bool>(const XMLNode XML_base, string path, bool& value,
                           bool verbose);

template <class T, template <class S, class R> class XMLValueInterpreter = Reader, class RequirementPolicy = RequiredPolicy>
struct Parameter : public ParameterBase, public XMLValueInterpreter<T, RequirementPolicy>
{
public:
    using ValueType = T;

    Parameter()
        : ParameterBase(), xml_path("")
    {}

    void setXMLPath(const string& path) { xml_path = path; }
    string XMLPath() const override { return xml_path; }

    void loadFromXML(XMLNode node) override
    {
        //std::cout << "Loading Parameter from " << node.getName()
        //    << std::endl;

        //std::cout << "XML Path: " << xml_path << std::endl;

        stored_node = node;
        try {
            string raw_string, raw_unit;

            if (!getXMLAttribute(node, xml_path, raw_string, false))
                XMLValueInterpreter<T, RequirementPolicy>::setMissing();
            else
                XMLValueInterpreter<T, RequirementPolicy>::setStringVal(raw_string);

            if (! XMLValueInterpreter<T, RequirementPolicy>::isMissing())
                XMLValueInterpreter<T, RequirementPolicy>::read(XMLValueInterpreter<T, RequirementPolicy>::stringVal());
        }
        catch (string e)
        {
            std::cerr << e << " in " << xml_path << std::endl;
        }
    }

    void init() override
    {
        try {
            XMLValueInterpreter<T, RequirementPolicy>::init();
        }
        catch (string e)
        {
            std::cerr << e << std::endl;
        }
    }

    void read(string value) override
    {
        XMLValueInterpreter<T, RequirementPolicy>::read(value);
    }

    bool isDefined() const
    {
        return ! XMLValueInterpreter<T, RequirementPolicy>::isMissing();
    }

private:
    string xml_path;
    XMLNode stored_node;
};

//
// Improve this class
//
// make sure that it looks for parameters deeper in the XML tree
// so that I don't have to write -> phi0/value
//
// so essentially a class that can take several of these values
//
// so a class that has two XMLValueInterpreters one for the unit and the other
// for the value
//
//
template <class T, template <class S, class R> class XMLValueInterpreter = Reader, class RequirementPolicy = RequiredPolicy>
struct Parameter2 : public ParameterBase, public XMLValueInterpreter<T, RequirementPolicy>
{
public:
    using ValueType = T;

    Parameter2()
        : ParameterBase(), xml_path("")
    {}

    void setXMLPath(const string& path) { xml_path = path; }
    string XMLPath() const override { return xml_path; }

    void loadFromXML(XMLNode node) override
    {
        //std::cout << "Loading Parameter from " << node.getName()
        //    << std::endl;

        //std::cout << "XML Path: " << xml_path << std::endl;

        stored_node = node;
        try {
            string raw_string, raw_unit;

            // can we do this better?
            XMLNode parameter_node = getParameterNode(node, xml_path.c_str());
            //if (parameter_node.isEmpty())
            //    throw ParameterNotFound(xml_path);

            if (!getXMLAttribute(parameter_node, "value", raw_string, false))
                XMLValueInterpreter<T, RequirementPolicy>::setMissing();
            else
                XMLValueInterpreter<T, RequirementPolicy>::setStringVal(raw_string);

            if (getXMLAttribute(parameter_node, "unit", raw_unit, false))
                XMLValueInterpreter<T, RequirementPolicy>::setUnitVal(raw_unit);

            if (! XMLValueInterpreter<T, RequirementPolicy>::isMissing())
                XMLValueInterpreter<T, RequirementPolicy>::read(XMLValueInterpreter<T, RequirementPolicy>::stringVal());
        }
        catch (string e)
        {
            std::cerr << e << " in " << xml_path << std::endl;
        }
    }

    void init() override
    {
        try {
            XMLValueInterpreter<T, RequirementPolicy>::init();
        }
        catch (string e)
        {
            std::cerr << e << std::endl;
        }
    }

    void read(string value) override
    {
        XMLValueInterpreter<T, RequirementPolicy>::read(value);
    }

    bool isDefined() const
    {
        return ! XMLValueInterpreter<T, RequirementPolicy>::isMissing();
    }

private:
    string xml_path;
    XMLNode stored_node;
};

