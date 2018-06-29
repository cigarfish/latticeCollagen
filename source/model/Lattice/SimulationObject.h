////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  SimulationObject.h                                            //
//                                                                            //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                //
//    Created:  2016-08-02                                                    //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

//#include "tools/factory/Parameter.h"
//#include "3rdparty/xmlParser/xmlParser.h"
#include "../BasicDatatypes/Color.h" // added by Jieling

class SimulationObject
{
public:
    SimulationObject();
    SimulationObject(/*XMLNode& node,*/
                     std::stringstream& errors,
                     std::stringstream& warnings);

    SimulationObject(const SimulationObject& other);
    SimulationObject(SimulationObject&& other);

    virtual ~SimulationObject();

    virtual void print(std::ostream& stream) const;
    virtual bool ready() const;
    virtual void initialize(void);
    //virtual void update(void) = 0;

    virtual void load(/*XMLNode& node,*/
                      std::stringstream& errors,
                      std::stringstream& warnings);

    std::string getName() const;

    // TODO make this Required again
    //Parameter<unsigned int, Reader, OptionalPolicy> id;
	unsigned int id;

	// added by Jieling
	ARGBColor color;
	void SetColor(double r, double g, double b, double alpha) { color.red = r; color.green = g; color.blue = b; color.alpha = alpha; }

protected:
    // maybe make this a parameter?
    bool verbose = false;
    std::string object_name;
    //XMLNode stored_node;

    //void registerParameter(ParameterBase& parameter);

private:
    //vector<ParameterBase*> parameters;
};

std::ostream& operator<<(std::ostream& stream, const SimulationObject& obj);
