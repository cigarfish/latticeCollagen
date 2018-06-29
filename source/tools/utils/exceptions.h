////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  exceptions.h                                                  //
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

#ifndef TISIM_EXCEPTIONS_H
#define TISIM_EXCEPTIONS_H

#include <exception>
#include <system_error>
#include <stdexcept>
#include <string>

class PluginConstructionFailedException : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    PluginConstructionFailedException(const std::string& msg)
        : runtime_error("Plugin " + msg + " construction failed!")
    {}
};

class PluginNotFound : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    PluginNotFound(const std::string& msg)
        : runtime_error("Plugin " + msg + " not found!")
    {}
};

class InternalPluginError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    InternalPluginError(const std::string& pluginName, const std::string& msg)
        : runtime_error(pluginName + ":" + msg + "!")
    {}
};

class ModelCellPopulationNotFound : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    ModelCellPopulationNotFound(const std::string& msg)
        : runtime_error("Cell Population " + msg + " not found!")
    {}
};

class ModelConstructionFailedException : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    ModelConstructionFailedException(const std::string & msg)
        : runtime_error("Model " + msg + " construction failed!")
    {}
};

class NotImplementedException : public std::logic_error
{
public:
    using std::logic_error::logic_error;

    NotImplementedException(const std::string & msg)
        : logic_error("Function " + msg + " not implemented")
    {}
};

class NullPtrDereference : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    NullPtrDereference(const std::string & msg)
        : runtime_error("Pointer " + msg + " is NULL!")
    {}
};

class FailureToConverge : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    FailureToConverge(const std::string & method_name,
                      const std::size_t steps)
        : runtime_error("Failure to converge in " + std::to_string(steps) +
                        " steps in method " + method_name + "!")
    {}
};

class SolverFailure : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    SolverFailure(const std::string & error_msg)
        : runtime_error(error_msg)
    {}

    SolverFailure(const std::string & method_name,
                  const double displacement)
        : runtime_error("Displacement in current step is " +
                        std::to_string(displacement) +
                        " steps in method " + method_name + "!")
    {}
};

class NotFiniteNumber : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    NotFiniteNumber(const std::string & msg,
                    const double number)
        : runtime_error("Number in " + msg + " has invalid value "
                        + std::to_string(number) + ".")
    {}
};

class FileNotFound : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    FileNotFound(const std::string & filename)
        : runtime_error("File " + filename + " does not exist.")
    {}
};

class ParameterNotFound : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    ParameterNotFound(const std::string & parameter)
        : runtime_error("Parameter " + parameter + " not found.")
    {}
};

class ExceededTolerance : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;

    ExceededTolerance(const std::string & where)
        : runtime_error("Exceeded tolerance in " + where + " .")
    {}
};

class CellSysInputOutputFailure : public std::system_error
{
    public:
        using std::system_error::system_error;

    CellSysInputOutputFailure(std::error_code ec, const std::string & what)
        : system_error(ec, what)
    {}
};

#endif
