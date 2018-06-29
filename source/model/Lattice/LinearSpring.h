////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  LinarSpring.h                                                 //
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
#include <utility>
#include <string>
#include <vector>
#include <sstream>
#include <memory>
//#include <QtXml>

#include "Vector.h"
#include "SimulationObject.h"
#include "LatticeSpring.h"
#include "ModelElementLatticeNode.h"
#include "model/Elements/ModelElementSphere.h"
//#include "CSParameterContext.h"
#include "Function.h"
//#include "CSRegistrar.h"
//#include "xmlParser/xmlParser.h"
//#include "tools/factory/Parameter.h"
#include "gui/GLTools/CSGLBar.h"

using std::vector;

/*
 * Structure for linear spring
 *
 * Need node connections
 *
 * Suppose I am dealing with ecm cell i
 * then this cell has a map connections (ptr to neighbour, spring struct).
 *
 *
 *
 * Structure for angular springs
 *
 * The involved nodes are (b, i, b')
 * where the angular spring is located at node i
 * the angle is formed by bi and ib'
 *
 *
 * For the force calculation we require
 *
 * 1) beta
 * 2) deviation from the original angle of the angular spring
 * 3)
 *
 * n_b = c_b / | c_b |
 *
 *
 */

class Lattice;

class LinearSpring final : public LatticeSpring
{
public:
    LinearSpring();
	// added by Jieling
	LinearSpring(ModelElementLatticeNode *N1, ModelElementLatticeNode *N2);
    //LinearSpring(XMLNode& node, std::stringstream& errors,
    //             std::stringstream& warnings);

    LinearSpring(const LinearSpring& other);
    LinearSpring(LinearSpring&& other);

    LinearSpring& operator=(const LinearSpring& other) = delete;
    LinearSpring& operator=(LinearSpring&& other) = delete;

    virtual ~LinearSpring();

private:
    // member variables
    ModelElementLatticeNode *n1, *n2;

    // Parameters that are loaded from XML
    //Parameter2<double> l0;
	double l0;

    // Parameter for the nodes of the spring
    //Parameter<unsigned int> Node1_Id, Node2_Id;
	unsigned int Node1_Id;
	unsigned int Node2_Id;

    // Equilibrium length of the spring
    double length0;

public:
    virtual void print(std::ostream& stream) const override;
    virtual void initialize() override;
    virtual void update(double timeStep) override; // update the force
    virtual bool ready() const override;
    virtual bool isRelaxed() const override;

    virtual double eq_position() const override;

    virtual std::vector<unsigned int> nodeIds() const override;
    virtual std::vector<ModelElementLatticeNode*> nodes() const override;
    virtual unsigned int numberOfNodes() const override;
	// added by Jieling
	BoundingBox * boundingBox();
	CSGLObject * GLObject() { return mpGLObject; }
	double getLength() { return (n1->position - n2->position).Norm(); }
	void strainTestForce();
	void unboundCheck();
	// for plasticity
	bool harden; // if it is harden: permanently elongated
	double hardenRatio; // post-yield modulus ratio
	double strainHarden; // the strain for hardening
	void elongation();
	double getStrainY();
	// for unbinding
	double n1_k_off;
	double n2_k_off; // initial: 0; 
	bool n1_unbound;
	bool n2_unbound; // false: bound; true: unbound
	void changeN1(ModelElementLatticeNode* N);
	void changeN2(ModelElementLatticeNode* N);
	double getStrain();

private:
    // parameter
    //Parameter2<bool, Reader, DefaultValPolicy> mDensityDependency;
	bool mDensityDependency;

    Vector3f position(void) const;
    double displacement(void) const;
    double strain(void) const;

// added by Jieling
public:
	BoundingBox mBoundingBox;
	CSGLObject * mpGLObject;
	double initialL0;
	double strainY;
	double initialLY;
	double curStrain;
	double curStrainStaining;
	unsigned int index;
};

//static CSRegistrar<Spring, LinearSpring> LinearSpringFactoryReg("LinearSpring");
//static CSRegistrar<Spring, LinearSpring, XMLNode&, std::stringstream&,
//    std::stringstream&> LinearSpringFactoryReg2("LinearSpring");

