////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  SimulationObject.cpp                                          //
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

#include "SimulationObject.h"
#include <iostream>

std::ostream& operator<<(std::ostream& stream, const SimulationObject& obj)
{
    obj.print(stream);
    return stream;
}

SimulationObject::SimulationObject()
{
    //id.setXMLPath("Id");
    //registerParameter(id);
}

SimulationObject::SimulationObject(/*XMLNode& node,*/
                                  std::stringstream& errors,
                                  std::stringstream& warnings)
    : SimulationObject()
{
    load(/*node,*/ errors, warnings);
}

SimulationObject::SimulationObject(const SimulationObject& other)
    : id(other.id)
{}

SimulationObject::SimulationObject(SimulationObject&& other)
    : id(other.id)
{}

SimulationObject::~SimulationObject()
{}

void SimulationObject::load(/*XMLNode& node,*/
                            std::stringstream& errors,
                            std::stringstream& warnings)
{
    //stored_node = node;
    //object_name = node.getName();

    //for (std::size_t i = 0; i < parameters.size(); ++i)
    //    parameters[i]->loadFromXML(node);
}

void SimulationObject::print(std::ostream& stream) const
{
    stream << "SimulationObject().";
}

bool SimulationObject::ready() const
{
    return true;
}

void SimulationObject::initialize()
{}

//void SimulationObject::registerParameter(ParameterBase& parameter)
//{
    //parameters.push_back(&parameter);
//}

std::string SimulationObject::getName() const
{
    return object_name;
}

