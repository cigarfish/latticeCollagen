////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  ModelElementLatticeNote.h                                     //
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

#include <string>
#include <functional>

#include "SimulationObject.h"
#include "model/Elements/ModelElementSphere.h"
#include "LatticeSpring.h"
#include "ModelElementFibre.h"
#include "model/Elements/ModelElementVesselSphere.h" // added by Jieling
//#include "CSRegistrar.h"
//#include "xmlParser/xmlParser.h"

// added by Jieling
namespace H5 { class CompType; };

struct ModelElementLatticeNode : public SimulationObject, public ModelElementSphere
{
    ModelElementLatticeNode();
    ModelElementLatticeNode(double x, double y, double z, double r);
    ModelElementLatticeNode(/*XMLNode& node,*/
                            std::stringstream& errors,
                            std::stringstream& warnings);

    ModelElementLatticeNode(const ModelElementLatticeNode& element);
    ModelElementLatticeNode(ModelElementLatticeNode&& element);

    ModelElementLatticeNode& operator=(const ModelElementLatticeNode&) = delete;
    ModelElementLatticeNode& operator=(ModelElementLatticeNode&&) = delete;

	// added by Jieling
	int latticeIndex;
	std::vector<ModelElementLatticeNode*> nodes; // neighbors
	std::vector<ModelElementFibre*> fibres; // fibres
	std::vector<LatticeSpring*> springs; // only for LinearSpring
	std::vector<LatticeSpring*> Rsprings; // only for RotationalSpring
	std::vector<LatticeSpring*> RspringsToDelete; // for RotationalSpring to delete
	ModelElementVesselSphere *vesselNeighbor;
	ModelElementVesselSphere *alongVesselNeighbor;
	static void HDF5DataFormat(H5::CompType &);
	static H5::CompType ParseHDF5DataFormat(H5::CompType & inputTypeDefinition,
		std::stringstream &,
		std::stringstream & warnings);
	void addNeighbor(ModelElementLatticeNode* E);
	void removeNeighbor(ModelElementLatticeNode* E);
	void addSpring(LatticeSpring *E);
	void removeSpring(LatticeSpring *E);
	void removeRspring(LatticeSpring *E);
	bool addFibre(ModelElementFibre *F);
	bool removeFibre(ModelElementFibre *F);

	unsigned int mIndex;

    virtual ~ModelElementLatticeNode();

    void Reset() override;

    void load(/*XMLNode& springNode,*/ std::stringstream& errors,
              std::stringstream& warnings) override;

    void print(std::ostream& stream) const override;
    bool ready() const override;

    void decay(const double rate);

    // density of decay
    double density;

    // Parameters for Pauli repulsion force!
    double hard_magnitude;
    double hard_force_crit;

    // some vars used to debug
    bool pulling;
    double pulling_force;

    bool isAnchor() const;

    // these maybe should be defined in ModelElement?
    //Parameter2<double, PointerQuantityReader, DefaultValPolicy> mYoung;
    //Parameter2<double, PointerReader, DefaultValPolicy> mPoisson;
	double mYoung;
	double mPoisson;

	// added by Jieling
	bool touched; // force applied on
	bool fibreEnd; // if it is the endnode of the fibre
	bool top; // for test in the cubic lattice
	bool bottom; // for test in the cubic lattice
	bool free; // if it is a free node (without binding to another fibre)
	Vector3f mLinearForce;
	Vector3f mRotationalForce;
	Vector3f mStrainTestForce;
	Vector3f mStrainTestForceLoaded;

    //Parameter2<double, PointerReader, DefaultValPolicy> mForceHardMag;
    //Parameter2<double, PointerReader, DefaultValPolicy> mForceHardCrit;
	double mForceHardMag;
	double mForceHardCrit;

    //Parameter<bool, Reader, DefaultValPolicy> mDecayAble;
	bool mDecayAble;

    void setFriction(double coefficient);
};

struct find_id : std::unary_function<ModelElementLatticeNode, bool>
{
    unsigned long id;
    find_id(unsigned long _id) : id(_id) {}
    bool operator()(const ModelElementLatticeNode* node) const
    {
        return node->id == id;
    }
};

// TEMP FIXME
//static CSRegistrar<ModelElementLatticeNode, ModelElementLatticeNode,
//    XMLNode&, std::stringstream&, std::stringstream&>
//    ModelElementLatticeNodeFactoryReg("Node");

//template <typename... Args>
//auto createLatticeNode(const std::string& name, Args&&... args)
//{
//    return create<ModelElementLatticeNode, Args...>
//        (name, std::forward<Args>(args)...);
//}
