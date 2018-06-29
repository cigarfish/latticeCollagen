////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  Spring.h                                                      //
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

#include <vector>
#include <string>
#include <memory>

#include "SimulationObject.h"
//#include "CSRegistrar.h"
//#include "CSPlugin.h"
#include "Function.h"

class ModelElementLatticeNode;
class ModelElementFibre;
class Lattice;

class LatticeSpring : public SimulationObject // changed by Jieling, from Spring to LatticeSpring
{
public:

	enum SpringType {
		TypeLinearSpring=0,
		TypeRotationalSpring
	};

    LatticeSpring();
    //Spring(XMLNode& node,
    //       std::stringstream& errors,
    //       std::stringstream& warnings);

    LatticeSpring(const LatticeSpring& spring);
    LatticeSpring(LatticeSpring&& spring);
    virtual ~LatticeSpring();

    virtual void update(double timeStep) = 0;
    virtual bool isRelaxed() const = 0;

    virtual double eq_position() const = 0;
    virtual double stiffness();

	double mYoung;
	double mPoisson; 
	double mRadius; // added by Jieling

	ModelElementFibre * fibre;

    // TODO can we remove this?
    void setLattice(Lattice * l);

    virtual std::vector<unsigned int> nodeIds() const = 0;
    virtual std::vector<ModelElementLatticeNode*> nodes() const = 0;
    virtual unsigned int numberOfNodes() const = 0;
	// added by Jieling
	//virtual void changeNode(ModelElementLatticeNode* N, int index) const = 0;

    virtual void load(/*XMLNode& springNode,*/
                      std::stringstream& errors,
                      std::stringstream& warnings) override;

	SpringType mSpringType;

protected:
    // indicates whether the spring is in equilibrium state when read
    //Parameter<bool, Reader, DefaultValPolicy> mEquilibrium;
	bool mEquilibrium;

    //std::unique_ptr<StrainFunction> msStiffness;
	StrainFunction *msStiffness;

    // TEMP
    Lattice * lattice;
};

struct find_spring_by_node_id : std::unary_function<const std::unique_ptr<LatticeSpring>&, bool>
{
    unsigned long id;
    find_spring_by_node_id(unsigned long _id): id(_id) {}
    //bool operator()(const std::unique_ptr<LatticeSpring>& spring)
    //{
    //    auto nodeIds = spring->nodeIds();
    //    return std::find(nodeIds.begin(), nodeIds.end(), id) != nodeIds.end();
    //}
};

// factory
//template <typename... Args>
//auto createSpring(const std::string& name, Args&&... args)
//{
//    return create<LatticeSpring, Args...>(name, std::forward<Args>(args)...);
//}

