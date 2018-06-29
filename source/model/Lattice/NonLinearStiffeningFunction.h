////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  NonLinearFunction.h                                           //
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

#include "SimulationObject.h"
#include "Function.h"
//#include "Parameter.h"

//
// This force calculation is done according to Steinwachs 2015 et. al.
//
//

class NonLinearStrainFunction final : public StrainFunction
{
public:
    NonLinearStrainFunction();

    //NonLinearStrainFunction(XMLNode& functionNode,
    //               std::stringstream& errors,
    //               std::stringstream& warnings);

    NonLinearStrainFunction(const NonLinearStrainFunction& other);
    NonLinearStrainFunction(NonLinearStrainFunction&& other);

    NonLinearStrainFunction& operator=(const NonLinearStrainFunction& other) = delete;
    NonLinearStrainFunction& operator=(NonLinearStrainFunction&& other) = delete;

    virtual ~NonLinearStrainFunction();

    double operator()(const double strain) const override;

    NonLinearStrainFunction* clone() const override
    {
        return new NonLinearStrainFunction(*this);
    }

    double stiffness() override;

    void print(std::ostream& stream) const override;
    bool ready() const override;

private:
    //Parameter2<double, QuantityReader, RequiredPolicy> k0;
    //Parameter2<double, QuantityReader, RequiredPolicy> d0;
    //Parameter2<double, QuantityReader, RequiredPolicy> ls;
    //Parameter2<double, QuantityReader, RequiredPolicy> ds;
	double k0;
	double d0;
	double ls;
	double ds;

    double ExpMinusOne(const double value) const;
};

//static CSRegistrar<StrainFunction, NonLinearStrainFunction>
//    NonLinearFunctionFactoryReg("NonLinearSpringFunction");

//static CSRegistrar<StrainFunction, NonLinearStrainFunction, XMLNode&,
//    std::stringstream&, std::stringstream&>
//    NonLinearFunctionFactoryReg2("NonLinearSpringFunction");

