////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  Spring.cpp                                                    //
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

#include "LatticeSpring.h"
#include "LinearFunction.h"

LatticeSpring::LatticeSpring()
    : SimulationObject(),
    lattice(nullptr)
{
    //mEquilibrium.setXMLPath("mEquilibrium");
    mEquilibrium = true;

	fibre = NULL;
    //registerParameter(mEquilibrium);
}

void LatticeSpring::load(/*XMLNode& springNode,*/
                  std::stringstream& errors,
                  std::stringstream& warnings)
{
    SimulationObject::load(/*springNode,*/ errors, warnings);
}

LatticeSpring::LatticeSpring(const LatticeSpring& spring)
    : SimulationObject(spring)
{
    //msStiffness.reset(spring.msStiffness->clone());
}

LatticeSpring::LatticeSpring(LatticeSpring&& spring)
    : SimulationObject(spring),
    msStiffness(std::move(spring.msStiffness))
{
    //spring->msStiffness.reset(nullptr);
}

LatticeSpring::~LatticeSpring()
{
    lattice = nullptr;
}

double LatticeSpring::stiffness()
{
    return msStiffness->stiffness();
}

void LatticeSpring::setLattice(Lattice* l) { lattice = l; }

