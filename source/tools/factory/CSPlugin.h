////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSPlugin.h                                                    //
//                                                                            //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                //
//    Created:  2016-10-20                                                    //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <sstream>
#include <string>
#include <vector>

#include "Parameter.h"
#include "CSRegistrar.h"
#include "xmlParser/xmlParser.h"

// Logger support
#include "spdlog/spdlog.h"

using std::string;
using std::vector;

class CSModel;
class ModelElement;
class BoundingBoxList;

class CSPlugin
{
public:
    CSPlugin()
        : mpModel(nullptr)
    {}

    virtual ~CSPlugin() {}

    CSPlugin(XMLNode& node,
             std::stringstream& errors,
             std::stringstream& warnings)
    {
        createFromXML(node, errors, warnings);
    }

    CSPlugin(const CSPlugin& other)
    {}

    CSPlugin(CSPlugin&& other)
    {}

    virtual void createFromXML(XMLNode& node,
                               std::stringstream & warnings,
                               std::stringstream & errors)
    {
        stored_node = node;
        plugin_name = node.getName();

        for (std::size_t i = 0; i < parameters.size(); ++i)
            parameters[i]->loadFromXML(node);
    }

    virtual std::string toString() const
    {
        return "CSPlugin";
    }

    virtual void Reset() = 0;
    virtual void initialize() {}

    void registerParameter(ParameterBase& parameter)
    {
        parameters.push_back(&parameter);
    }

    string getFullName() const
    {
        return plugin_name;
    }

    bool setParameter(string xml_path, string value)
    {
        return true;
    }

    virtual void Update()
    {}

    virtual unsigned int Priority() const
    {
        return 0;
    }

    void setModelPointer(CSModel * model)
    {
        mpModel = model;
    }

    virtual void output(unsigned int outidx) const
    {}

    virtual void report() const
    {}

    virtual void print(std::ostream& stream) const
    {
        stream << "CSPlugin";
    }

    void setLogger(std::shared_ptr<spdlog::logger> new_console)
    {
        console = new_console;
    }

protected:
    bool verbose = false;
    string plugin_name;
    XMLNode stored_node;
    CSModel * mpModel;

    // logger console
    std::shared_ptr<spdlog::logger> console;

private:
    vector<ParameterBase*> parameters;
};


template <class CellPointerType>
class CSModelPlugin : public CSPlugin
{
public:

    CSModelPlugin()
        : CSPlugin()
    {}

    virtual ~CSModelPlugin() {}

    CSModelPlugin(XMLNode& node,
             std::stringstream& errors,
             std::stringstream& warnings)
        : CSPlugin(node, errors, warnings)
    {}

    CSModelPlugin(const CSModelPlugin& other)
        : CSPlugin(other)
    {}

    CSModelPlugin(CSModelPlugin&& other)
        : CSPlugin(other)
    {}

    virtual std::string toString() const
    {
        return "CSModelPlugin";
    }

    virtual void print(std::ostream& stream) const
    {
        stream << "CSModelPlugin";
    }

    // Main interface change!
    virtual void update(CellPointerType cell) = 0;
};

// allow sorting of CSPlugin
struct less_than_key
{
    inline bool operator()(const std::unique_ptr<CSPlugin>& plugin1,
                           const std::unique_ptr<CSPlugin>& plugin2) const
    {
        return (plugin1->Priority() < plugin2->Priority());
    }
};


// Factory method
template <typename... Args>
auto createPlugin(const string& name, Args&&... args)
{
    return create<CSPlugin, Args...>(name, std::forward<Args>(args)...);
}




