////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  RotationalSpring.h                                            //
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

#include "Vector.h"
#include "SimulationObject.h"
#include "LatticeSpring.h"
#include "ModelElementLatticeNode.h"
#include "model/Elements/ModelElementSphere.h"
#include "Function.h"
#include "tools/utils/exceptions.h"
//#include "CSRegistrar.h"
//#include "xmlParser/xmlParser.h"
//#include "tools/factory/Parameter.h"

using std::vector;

class Lattice;

//
// Figure out how to deal with rotations of the rotational axis of this spring
//
// This implements a rotational spring element. The three nodes of the
// rotational spring are given denotes by N1, N2, NC. The torsional spring is
// placed at NC between the vectors formed by NCN1 and NCN2.
//
//

class RotationalSpring final : public LatticeSpring
{
public:
    RotationalSpring();
    //RotationalSpring(XMLNode& node, std::stringstream& errors,
    //                 std::stringstream& warnings);
	RotationalSpring(ModelElementLatticeNode *N1, ModelElementLatticeNode *N2, ModelElementLatticeNode *Center);

    RotationalSpring(const RotationalSpring& other);
    RotationalSpring(RotationalSpring&& other);

    RotationalSpring& operator=(const RotationalSpring& other) = delete;
    RotationalSpring& operator=(RotationalSpring&& other) = delete;

    virtual ~RotationalSpring();

private:
    // lattice nodes
    ModelElementLatticeNode *n1, *n2, *center;

    // Parameters - this is provided if the current is different
    //Parameter2<double> Phi0;
	double Phi0;

    // The nodes of the spring
    //Parameter<unsigned int> Node1_Id, Node2_Id, NodeC_Id;
	unsigned Node1_Id;
	unsigned Node2_Id;
	unsigned NodeC_Id;

    // Holds the relaxed angle
    double phi0;

	// added by Jieling
	double curv0; // the default curvature
	double fibreLength;
	Vector3f n10;
	Vector3f n20;
	Vector3f center0; // the original position of the nodes, to calculate the change of the angle

public:
    virtual void print(std::ostream& stream) const override;
    virtual void initialize() override;
    virtual void update(double timeStep) override;
    virtual bool ready() const override;
    virtual bool isRelaxed() const override;

    virtual double eq_position() const override;

    virtual std::vector<unsigned int> nodeIds() const override;
    virtual std::vector<ModelElementLatticeNode*> nodes() const override;
    virtual unsigned int numberOfNodes() const override;
	// added by Jieling
	//virtual const Vector3f& GetEndPoints(int index) const override;
	//virtual void setEndPoints() override;
	ModelElementLatticeNode *getN1() { return n1; }
	ModelElementLatticeNode *getN2() { return n2; }
	ModelElementLatticeNode *getCenter() { return center; }
	void changeN1(ModelElementLatticeNode* N);
	void changeN2(ModelElementLatticeNode* N);
	void changeCenter(ModelElementLatticeNode* N);
	double getLength1C() { return (n1->position - center->position).Norm(); }
	double getLength2C() { return (n2->position - center->position).Norm(); }
	double getPhi0() { return phi0; }
	double getCurvature();
	double getBendingForce();
	double getEI();
	void setPhi0(double p) { phi0 = p; }
	void reset();
	void setFibreLength(double l) { fibreLength = l; reset(); }
	void setFibreLengthOnly();
	void strainTestForce();
};

// factory
//static CSRegistrar<Spring, RotationalSpring>
//    RotationalSpringFactoryReg("RotationalSpring");
//static CSRegistrar<Spring, RotationalSpring, XMLNode&, std::stringstream&,
//    std::stringstream&> RotationalSpringFactoryReg2("RotationalSpring");

