////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  LinearFunction.h                                              //
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

class LinearStrainFunction final : public StrainFunction
{
public:
    LinearStrainFunction();

    //LinearStrainFunction(XMLNode& functionNode,
    //               std::stringstream& errors,
    //               std::stringstream& warnings);

    LinearStrainFunction(const LinearStrainFunction& other);
    LinearStrainFunction(LinearStrainFunction&& other);

    LinearStrainFunction& operator=(const LinearStrainFunction& other) = delete;
    LinearStrainFunction& operator=(LinearStrainFunction&& other) = delete;

    virtual ~LinearStrainFunction();

    double operator()(const double strain) const override;

    LinearStrainFunction* clone() const override
    {
        return new LinearStrainFunction(*this);
    }

    double stiffness() override;

    void print(std::ostream& stream) const override;
    bool ready() const override;

	// added by Jieling
	void setK(double K) { k = K; }

private:
    // FIXME and combine again
    //Parameter2<double, QuantityReader, RequiredPolicy> k;
	double k;
};

//static CSRegistrar<StrainFunction, LinearStrainFunction>
//    LinearFunctionFactoryReg("LinearSpringFunction");

//static CSRegistrar<StrainFunction, LinearStrainFunction, XMLNode&,
//    std::stringstream&, std::stringstream&>
//    LinearFunctionFactoryReg2("LinearSpringFunction");
