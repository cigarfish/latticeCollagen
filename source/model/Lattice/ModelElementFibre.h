////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  ModelElementFibre.h                                           //
//                                                                            //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                //
//    Created:  2017-01-02                                                    //
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
#include "LatticeSpring.h"
//#include "CSRegistrar.h"
//#include "xmlParser/xmlParser.h"

using std::vector;

class ModelElementLatticeNode;
//class LinearSpring;

struct ModelElementFibre : public SimulationObject
{
    using NodeIds = std::vector<std::size_t>;

    ModelElementFibre();
    //ModelElementFibre(XMLNode& node, std::stringstream& errors,
    //                 std::stringstream& warnings);

    ModelElementFibre& operator=(const ModelElementFibre&) = delete;
    ModelElementFibre& operator=(ModelElementFibre&&) = delete;

    virtual ~ModelElementFibre();

    //void addNode(std::size_t NodeId);
    //void removeNode(std::size_t NodeId);
	// revised by Jieling
	void addNode(ModelElementLatticeNode* N);
	void removeNode(ModelElementLatticeNode* N);
	bool checkNode(ModelElementLatticeNode *N);
	void addSpring(LatticeSpring *S);
	void removeSpring(LatticeSpring* S);
	bool checkSpring(LatticeSpring *S);
	void unboundE(int end, ModelElementLatticeNode* N); // detach the fibre from end node 1/2 and bind to N
	void unboundM(ModelElementLatticeNode* O, ModelElementLatticeNode *N); // detach the fibre from middle node O and bind to N
	LatticeSpring* locateSpring(ModelElementLatticeNode* e1, ModelElementLatticeNode* e2);
	vector<ModelElementLatticeNode*>& getNodes() { return nodes; }
    const vector<ModelElementLatticeNode*>& getNodes() const;
	// added by Jieling
	vector<LatticeSpring*> &getSprings() { return springs; }
	bool locateNodeSpring(ModelElementLatticeNode *node1, ModelElementLatticeNode *node2, std::vector<ModelElementLatticeNode*> &ns, std::vector<LatticeSpring*> &ls);
	int relativePosition(ModelElementLatticeNode *node1, ModelElementLatticeNode *node2, 
		std::vector<ModelElementLatticeNode*> &n1s, std::vector<LatticeSpring*> &l1s, 
		std::vector<ModelElementLatticeNode*> &n2s, std::vector<LatticeSpring*> &l2s);
	bool relativeSinglePosition(ModelElementLatticeNode *node, 
		std::vector<ModelElementLatticeNode*> &ns, std::vector<LatticeSpring*> &ls);
	ModelElementLatticeNode* splitNode(ModelElementLatticeNode *node1, ModelElementLatticeNode *node2);
	bool selfMerge();
	void addVirtualNodes();
	double getEndBendingForce(int endType);
	
	// added by Jieling
	ModelElementLatticeNode* n1;
	ModelElementLatticeNode* n2; // the two end-nodes of the fiber
	vector<ModelElementLatticeNode*> nodes; // nodes within the fibre
	vector<LatticeSpring*> springs; // LinearSprings within the fibre
	// bending type: 
	// -1: not able to bend
	// 0: fixed end node 
	bool bending; // if this fibre has bending energy
	Vector3f mVirtualN1; // a virtual node next to n1
	Vector3f mVirtualN2; // a virtual node next to n2
	double length; // fibre length
	bool n1Bound;
	bool n2Bound; // if these two end-nodes are bound to the other fibres
	LatticeSpring* getE1();
	LatticeSpring* getE2(); // the two end-node springs
	
    //virtual void load(XMLNode& node,
    //                  std::stringstream& errors,
    //                  std::stringstream& warnings) override;

    virtual void print(std::ostream& stream) const override;

    virtual bool ready() const override;
    virtual void initialize(void) override;
    // make abstract TODO
    void update(void) ;

//private:
public:
    // list of nodes that belong to the fibre
    // list of node ids? or should I put pointers here
    //NodeIds nodes;
    //Parameter<unsigned int> fibre_id;
	unsigned int fibre_id;
};

// TEMP FIXME
//static CSRegistrar<ModelElementFibre, ModelElementFibre, XMLNode&,
//    std::stringstream&, std::stringstream&> ModelElementFibreFactoryReg("Fibre");

//template <typename... Args>
//auto createFibre(const std::string& name, Args&&... args)
//{
//    return create<ModelElementFibre, Args...>(name, std::forward<Args>(args)...);
//}
