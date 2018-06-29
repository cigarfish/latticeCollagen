////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  Lattice.h                                                     //
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

#include <sstream>
#include <vector>
#include <array>
#include <memory>
#include <map>

//#include "3rdparty/xmlParser/xmlParser.h"
//#include "tools/factory/CSPlugin.h"
//#include "CSRegistrar.h"
#include "SimulationObject.h"
#include "LatticeSpring.h"
//#include "tools/factory/Parameter.h"
#include "model/Model/CSModel.h"
#include "model/Lattice/ModelElementLatticeNode.h"
// added by Jieling
#include "tools/math/delaunay.h"

class ModelElementLatticeNode;
class ModelElementFibre;
class BoundingBoxList;

using std::vector;

// This class really implements a network as defined by Stein et. al 2008
//
// A network is a triple N = (V, X, F)
//
// where V = list of vertices (1, 2, ..., n)
//       X = coordinates of the verticies (u1, u2, ..., un)
//       F = vectors of the fibres in the network
//
// F = (f1, f2, ..., fn) where fk = (v_1^k, ...., v_pk^k) a vector of the
// indeces of the coordinates of the vertices which belong to a given fibre
//
class Lattice
{
public:

    Lattice();
    ~Lattice();

    std::size_t size() const;
    std::size_t NumberOfNodes() const;
    std::size_t NumberOfSprings() const;

    vector<ModelElementLatticeNode*>&  getNodes();
    const vector<ModelElementLatticeNode*>&  getNodes() const;
    //vector<std::unique_ptr<Spring>>& getSprings();
	// added by Jieling
	vector<LatticeSpring*> &getSprings() { return springs; };
	vector<ModelElementFibre*> &getFibres() { return fibres; }

    // Function to remove nodes from the lattice
	bool addNode(ModelElementLatticeNode * to_add);
	bool removeNode(ModelElementLatticeNode * to_remove);
    bool removeNodes(const std::vector<ModelElementLatticeNode * >& to_remove); // DO NOT USE
	// added by Jieling
	// the two basic functions to remove a fibre/spring from the lattice system
	bool removeSpring(LatticeSpring * to_remove);
	void removeNodeRSprings(ModelElementLatticeNode * N); // remove all rotational springs from the node
	void removeNodeRSpringsL(ModelElementLatticeNode * N1, ModelElementLatticeNode * N2); // remove the rotational spring regarding N1 and N2
	bool removeFibre(ModelElementFibre * to_remove);
	void freeFibre(ModelElementFibre *F, int eIndex); // eIndex: 1/2 for node 1/2
	void createNetwork(double cutoff, int subDivide);
	// cutoff: only spring with distance < cutoff
	// subDivide: divide the spring into subDivide segements
	void updateFibreBound(); // fibre unbound if strain > th

	void initialize();// override;
	void Update();// override;
	void Reset();// override;
	void report() const;// override;
    unsigned int Priority() const// override
    {
        return 25;
    }

    std::string toString() const// override
    {
        return "Lattice";
    }

	void setArena(CSGLArena *Arena) { mpArena = Arena; }
	void setBoundingBoxList(BoundingBoxList *B);
    // Interface to access vertices and fibres in the network
    ModelElementLatticeNode* getNode(const std::size_t idx);
    const ModelElementLatticeNode* getNode(const std::size_t idx) const;

    // Get Nodes by node id
    ModelElementLatticeNode* getNodeById(const std::size_t id);
    const ModelElementLatticeNode* getNodeById(const std::size_t id) const;

    // get coordinate of node
    Vector3f& getNodeCoordinate(const std::size_t idx);
    //const Vector3f& getNodeCoordinate(const std::size_t idx) const;

    // get fibres
    // what do i want to get?

    // methods to check dimension of the lattice
    std::array<double, 3> LatticeDimension() const;
    //void output(unsigned int outidx) const override;

private:
    // the nodes in the lattice
    const std::string nodeName {"Node"};
    vector<ModelElementLatticeNode*> nodes;
	// added by Jieling
	vector<ModelElementLatticeNode*> mNodesToDelete;
	vector<LatticeSpring*> mSpringsToDelete;

    // node density threshold for removal
    static constexpr double mDensityThreshold = 1.e-8;

    // Count number of non-anchor nodes
    unsigned int getNumberOfNodes() const;

    // the springs in the lattice here so we can update them
    const std::string springName {"Spring"};
    //vector<std::unique_ptr<Spring>> springs;
	vector<LatticeSpring*> springs;

    // the fibres in the network
    const std::string fibreName {"Fibre"};
    vector<ModelElementFibre*> fibres; // revised by Jieling

    // Lattice parameters
    //Parameter2<std::string, Reader, RequiredPolicy> mPathName;
    //Parameter2<std::string, Reader, DefaultValPolicy> mCellPopulationName;
	std::string mCellPopulationName;
    //Parameter2<bool, Reader, DefaultValPolicy> mIncludeAnchorsInBoundingBox;
	bool mIncludeAnchorsInBoundingBox;
    // TODO make this available for all plugins
    //Parameter2<bool, Reader, DefaultValPolicy> mVerbose;
    //Parameter2<bool, Reader, DefaultValPolicy> mCheckLatticeDimension;
	bool mCheckLatticeDimension;
    //std::string mConnectionsFile, mH5PartFileNameECM, mH5PartFileNameCon;

    // friction coefficient between nodes
    //Parameter2<double, QuantityReader, RequiredPolicy> mGamma;
	double mGamma;

    // nodes for output
    std::size_t numberOfNodesOutput;

    // enlargement tolerance
    static constexpr double mEnlargementFactor = 5;

    // lattice dimensions
    std::array<double, 3> originalDimensions;
    std::array<double, 3> mMaxAllowedLattice;

public:
    // pointer to the appropriate boundingboxlist
    BoundingBoxList * mpBoundingBoxList;
	// pointer to the arena
	CSGLArena *mpArena;

    // output details
    //void writeH5PartOutput(unsigned int outputTime) const;
    void SaveECMConnectionsRaw(unsigned int idx) const;
    //void writeH5PartOutputECMConnections(unsigned int outputTime) const;

protected:
	CSModel *mpModel;
};

//static CSRegistrar<CSPlugin, Lattice> LatticeFactory("Lattice");
//static CSRegistrar<CSPlugin, Lattice, XMLNode&,
//    std::stringstream&, std::stringstream&> LatticeFactory2("Lattice");


