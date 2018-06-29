////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  Function.h                                                    //
//                                                                            //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                //
//    Created:  2016-10-06                                                    //
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

#include "SimulationObject.h"
//#include "CSRegistrar.h"

class CSParameterContext;

class StrainFunction : public SimulationObject
{
public:
    StrainFunction();
    //StrainFunction(XMLNode& functionNode,
    //               std::stringstream& errors,
    //               std::stringstream& warnings);

    virtual ~StrainFunction();

    virtual double operator()(const double strain) const = 0;
    virtual StrainFunction* clone() const = 0;
    virtual double stiffness() = 0;
    //virtual void print(std::ostream& stream) const = 0;
};

// factory
//template <typename... Args>
//auto createFunction(const std::string& name, Args&&... args)
//{
//    return create<StrainFunction, Args...>(name, std::forward<Args>(args)...);
//}

