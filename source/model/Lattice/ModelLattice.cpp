////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  Lattice.cpp                                                   //
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

#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <tuple>

//#include "H5PartHelpers.h"
//#include <H5Cpp.h>
//#include "H5Part.h"
//#include "H5Block.h"
//#include "H5BlockReadWrite.h"

#include "LatticeSpring.h"
#include "ModelLattice.h"
#include "ModelElementLatticeNode.h"
//#include "CSParameterContextTemporary.h"
#include "ModelElementFibre.h"

#include "tools/factory/Types.h"
#include "tools/utils/array_tools.h"
#include "tools/utils/exceptions.h"
#include "LinearSpring.h"
#include "RotationalSpring.h"
#include "CSModel.h"
#include "tools/model/BoundingBoxList.h"
#include "gui/CSGLArena.h"
#include "gui/GLTools/CSGLBar.h"

Lattice::Lattice()
    : numberOfNodesOutput(0)
{
    // set plugin name
    //plugin_name = "Lattice";

    //mPathName.setXMLPath("path");
    //registerParameter(mPathName);

    mGamma = 0.5;
    //registerParameter(mGamma);

    mCellPopulationName = "Lattice";
    //registerParameter(mCellPopulationName);

    mIncludeAnchorsInBoundingBox = false;
    //registerParameter(mIncludeAnchorsInBoundingBox);

    //mVerbose.setXMLPath("Verbose");
    //mVerbose.setDefault("false");
    //registerParameter(mVerbose);

    mCheckLatticeDimension = false;
    //registerParameter(mCheckLatticeDimension);
}

Lattice::~Lattice()
{
    // de-register cell population from CSModel
    //mpModel->removeCellPopulation(mCellPopulationName);

    //const std::string anchorPopulationName = mCellPopulationName + "_anchors";
    //mpModel->removeCellPopulation(anchorPopulationName);
}

void Lattice::setBoundingBoxList(BoundingBoxList *B)
{
	mpBoundingBoxList = B;
}

ModelElementLatticeNode* Lattice::getNode(const std::size_t idx)
{
    // this does range checking
    return nodes.at(idx);
}

const ModelElementLatticeNode* Lattice::getNode(const std::size_t idx) const
{
    return nodes.at(idx);
}

ModelElementLatticeNode* Lattice::getNodeById(const std::size_t id)
{
    auto NodeIt = std::find_if(nodes.begin(), nodes.end(), find_id(id));
    if (NodeIt == nodes.end())
    {
        //if (verbose)
        //    console->error("Node number {0:d} not found.", id);
		std::cout << "	-> Node " << id << " is not found." << std::endl;

        return nullptr;
    }
    return *NodeIt;
}

const ModelElementLatticeNode* Lattice::getNodeById(const std::size_t id) const
{
    auto NodeIt = std::find_if(nodes.begin(), nodes.end(), find_id(id));
    if (NodeIt == nodes.end())
    {
        //if (verbose)
        //    console->error("Node number {0:d} not found.", id);
		std::cout << "	-> Node " << id << " is not found." << std::endl;

        return nullptr;
    }
    return *NodeIt;
}

Vector3f& Lattice::getNodeCoordinate(const std::size_t idx)
{
	return getNode(idx)->position;//getNode(idx)->coordinates();
}

void Lattice::initialize()
{
    //console->info("Initializing lattice.");

    //if (mIncludeAnchorsInBoundingBox())
    //    console->warn("Including lattice anchors in the bounding box for "
    //                  "collision detection! This may have some undesired "
    //                  "consequences.  Are you sure you want to continue?");

    //if (mCheckLatticeDimension())
    //    console->info("Enabling lattice dimension check!");
    //else
    //    console->warn("Lattice dimension check is not enabled!");

    //mH5PartFileNameECM = mpModel->mOutputPrefix() + "_ECM_h5.h5part";
    //mH5PartFileNameCon = mpModel->mOutputPrefix() + "_ECMConnections_h5.h5part";

    // move somewhere else!
    //console->info("Creating h5part lattice output files!");

    //H5PartFile * H5PartOutECM =
    //    TiSimOpenH5PartFile(mH5PartFileNameECM.c_str(), H5PART_WRITE);

    //H5PartFile * H5PartOutCon =
    //    TiSimOpenH5PartFile(mH5PartFileNameCon.c_str(), H5PART_WRITE);

    //TiSimCloseH5PartFile(H5PartOutECM);
    //TiSimCloseH5PartFile(H5PartOutCon);

    //console->info("Nodes... ( {0:d} ).", nodes.size());

    // Move to a function?
    //console->info("Creating lattice bounding box.");
    //if (mpBoundingBoxList != nullptr)
    //    delete mpBoundingBoxList;
    //mpBoundingBoxList = new BoundingBoxList(mpModel->dimension);

    // register the cell population
    //mpModel->registerCellPopulation(mCellPopulationName());

    // needs to be here apparently
    const std::string anchorPopulationName = mCellPopulationName + "_anchors";

    for (auto& node : nodes)
    {
        node->initialize();
        node->setFriction(mGamma);
        assert(node->ready());

		if (node->isAnchor() && mIncludeAnchorsInBoundingBox) {
			//mpModel->registerCell(static_cast<ModelElement*>(node), anchorPopulationName);
		}
		else {
			//mpModel->registerCell(static_cast<ModelElement*>(node), mCellPopulationName);
		}
    }

    // Compute the number of nodes we should output
    numberOfNodesOutput = getNumberOfNodes();

    //if (mVerbose())
    //{
    //    std::cout << "Nodes..." << std::endl;
    //    for (const auto& node : nodes)
    //        std::cout << *node << std::endl;
    //}

    //console->info("Springs... ({0:d}).", springs.size());
    for (const auto& spring : springs)
        spring->initialize();

    //console->info("fECM: {0:f}.", mGamma());
    for (const auto& spring : springs)
    {
		if (!spring->isRelaxed())
			//console->warn("Spring {0:d} is not relaxed.", spring->id());
			std::cout << "	-> Spring " << spring->id << " is not relaxed." << std::endl;
        //std::cout << *spring << std::endl;
    }

    //console->info("Fibres... ({0:d}).", fibres.size());
    for (const auto& fibre : fibres)
        fibre->initialize();

    //console->info("Lattice initialization completed.");

    // Verify construction!
    //console->info("Lattice verification of initialization.");
    bool ready = true;
    for (const auto& node : nodes)
        if (!node->ready())
            ready = false;

    for (const auto& spring : springs)
        if (!spring->ready())
            ready = false;

    for (const auto& fibre : fibres)
        if (!fibre->ready())
            ready = false;

    // store the original dimensions
    originalDimensions = LatticeDimension();
    mMaxAllowedLattice = multiply(originalDimensions,
                                  std::make_index_sequence<3>{},
                                  mEnlargementFactor);

    // just do this here for the moment
    SaveECMConnectionsRaw(0);

    //if (mVerbose())
    //    print(std::cout);

	if (ready)
		//console->info("Lattice initialization successful!");
		std::cout << "	-> Lattice initialization successful!" << std::endl;
	else
		//console->error("Something went wrong during lattice initialization.");
		std::cout << "	-> Something went wrong during lattice initialization." << std::endl;
}

void Lattice::Update()
{
    // first check any nodes that may be removed
    std::vector<ModelElementLatticeNode *> disappeared;
    for (const auto& node : nodes)
    {
        if (node->density < mDensityThreshold)
            disappeared.push_back(node);
    }

    // Remove nodes!
    removeNodes(disappeared);

    //std::cout << "Lattice update" << std::endl;
    for (auto& spring : springs)
        spring->update(1.);

    // Check lattice dimensions
    if (mCheckLatticeDimension)
    {
        auto currentDimensions = LatticeDimension();
        if (currentDimensions > mMaxAllowedLattice)
        {
            std::cerr << "Error Lattice has enlarged above its tolerance.\n"
                << "\tMaximum allowed lattice size is " << mMaxAllowedLattice << ".\n"
                << "\tCurrent lattice size is " << currentDimensions << ".\n"
                << std::endl;
            throw ExceededTolerance("Lattice");
        }

        for (auto& node : nodes)
        {
            //if (node->attributes["displacement"] > 0.002)
            //    std::cout << *node << std::endl;

            assert(node->density <= 1.);
        }

    } // End Lattice-Verify
}

void Lattice::Reset()
{
    for (auto& node : nodes)
    {
        node->Reset();
        // TODO should this really be resetted each time?
        node->setFriction(mGamma);
    }
}

std::size_t Lattice::size() const { return (nodes.size() + springs.size()); }
std::size_t Lattice::NumberOfNodes() const { return nodes.size(); }
std::size_t Lattice::NumberOfSprings() const { return springs.size(); }

vector<ModelElementLatticeNode*>& Lattice::getNodes()
{
    return nodes;
}

const vector<ModelElementLatticeNode*>& Lattice::getNodes() const
{
    return nodes;
}

std::array<double, 3> Lattice::LatticeDimension() const
{
    double min_x = 0., max_x = 0.;
    double min_y = 0., max_y = 0.;
    double min_z = 0., max_z = 0.;

    for (const auto& node : nodes)
    {
        auto pos = node->position;

        if (min_x > pos.x)
            min_x = pos.x;
        else if (max_x < pos.x)
            max_x = pos.x;

        if (min_y > pos.y)
            min_y = pos.y;
        else if (max_y < pos.y)
            max_y = pos.y;

        if (min_z > pos.z)
            min_z = pos.z;
        else if (max_z < pos.z)
            max_z = pos.z;
    }

    return {max_x - min_x, max_y - min_y, max_z - min_z};
}

void Lattice::report() const
{}

void Lattice::SaveECMConnectionsRaw(unsigned int idx) const
{
    // TODO broken
    //if (mConnectionsFile == "")
    //    return;

    if (!springs.size())
        return;

    using saveData_type = std::tuple<ModelElementLatticeNode*,
          ModelElementLatticeNode*, double, double>;

    std::vector<saveData_type> verts;
    unsigned int found = 0;

    for ( const auto& spring : springs)
    {
        if (spring->numberOfNodes() != 2)
            continue;

        found++;
        auto nodes = spring->nodes();
        //verts.emplace_back(make_tuple(nodes[0], nodes[1], spring->eq_position(),
        //                              spring->stiffness()));
    }

    //console->info("Saving {0:d} ecm cells, and {1:d} connections.", found, verts.size());
    //ofstream myfile;
    //std::string mConnectionsFile = mpModel->mOutputPrefix()
    //    +"_ECMconnections_"+to_string(idx)+".dat";
    //myfile.open (mConnectionsFile.c_str());
    //myfile << "#\tx0\ty0\tz0\t#\tx1\ty1\tz1\tl0\tk";
    //cout << "writing" << std::endl;
    //for (const auto& v : verts)
    //{
    //    myfile << "\n" << get<0>(v)->mSimObjId << "\t" << get<0>(v)->position.x
    //        << "\t" << get<0>(v)->position.y << "\t" << get<0>(v)->position.z;
    //    myfile << "\t" << get<1>(v)->mSimObjId << "\t" << get<1>(v)->position.x
    //        << "\t" << get<1>(v)->position.y << "\t" << get<1>(v)->position.z;
    //    myfile << "\t" << get<2>(v) << "\t" << get<3>(v);
    //}
    //myfile.close();
}

bool Lattice::addNode(ModelElementLatticeNode * to_add)
{
	//if (to_add != NULL)
	//	std::cout << "	-> Adding a node" << std::endl;

	nodes.push_back(to_add);

	mpBoundingBoxList->add(to_add);

	return true;
}

bool Lattice::removeNode(ModelElementLatticeNode * to_remove)
{
	if (to_remove != NULL)
		std::cout << "	-> Removing a node " << to_remove->mGlobalIndex << std::endl;
	else
		return false;

	std::vector<ModelElementLatticeNode *>::iterator found
		= std::find(nodes.begin(), nodes.end(), to_remove);

	// remove from the nodes list
	if (found != nodes.end())
		nodes.erase(found);
	else
		std::cerr << "	-> Error in ModelLattice: Attempt to remove unregistered LatticeNode!\n" << std::endl;

	// remove from the bounding box
	mpBoundingBoxList->remove(to_remove);

	// remove from the arena
	mpArena->removeObject(to_remove->GLObject());

	// delete node
	to_remove->~ModelElementLatticeNode();

	/*for (int i = 0; i < (int)to_remove->springs.size(); i++)
	{
		// edge to remove
		LatticeSpring *sI = to_remove->springs.at(i);
		std::vector<LatticeSpring*>::iterator foundI
			= std::find(to_remove->springs.begin(), to_remove->springs.end(), sI);
		if (foundI != to_remove->springs.end())
			to_remove->springs.erase(foundI);
		else
			std::cerr << "	-> Error in ModelLattice: Attempt to remove unregistered spring!" << std::endl;
		
		// also to remove the spring from other LatticeNode
		for (int j = 0; j < sI->nodes().size(); j++)
		{
			ModelElementLatticeNode *nJ = sI->nodes().at(j);
			if (to_remove != nJ) // not to_remove node
			{
				std::vector<LatticeSpring*>::iterator foundJ
					= std::find(nJ->springs.begin(), nJ->springs.end(), sI);
				if (foundJ != nJ->springs.end())
					nJ->springs.erase(foundJ);
				else
					std::cerr << "	-> Error in ModelLattice: Attempt to remove unregistered spring!" << std::endl;
			}
		}

		if (sI->mSpringType == LatticeSpring::TypeLinearSpring)
		{
		}
		else
		{
		}
	}*/

	return true;
}

bool Lattice::removeNodes(const std::vector<ModelElementLatticeNode * >& to_remove)
{
	// DO NOT USE, UNDER TEST
	/*
	if (to_remove.size())
		//console->info("Removing springs and nodes!");
		std::cout << "	-> Removing springs and nodes!" << std::endl;

    for (auto node : to_remove)
    {
        springs.erase(std::remove_if(springs.begin(), springs.end(),
                                     find_spring_by_node_id(node->mSimObjId)),
                      springs.end());
    }

    for (auto node : to_remove)
    {
        // make sure we dont remove any anchors
        assert(!node->isAnchor());
		//Jieling, TODO
        //mpModel->removeCell(node, mCellPopulationName);
        nodes.erase(std::remove(nodes.begin(), nodes.end(), node), nodes.end());
    }

    // make sure this is reduced!
    numberOfNodesOutput = getNumberOfNodes();
	*/
    return true;
}

bool Lattice::addLinearSpring(LatticeSpring * to_add)
{
	springs.push_back(to_add);
	return true;
}

bool Lattice::addRotationalSpring(LatticeSpring * to_add)
{
	springs.push_back(to_add);
	return true;
}

void Lattice::mergeTwoNodes(ModelElementLatticeNode * n1, ModelElementLatticeNode * n2)
{
	// remove all rotational springs containing n1 and n2
	while (!n1->Rsprings.empty())
	{
		LatticeSpring *rsi = n1->Rsprings.back();
		removeLatticeSpring(rsi);
		n1->Rsprings.pop_back();
	}
	while (!n2->Rsprings.empty())
	{
		LatticeSpring *rsi = n2->Rsprings.back();
		removeLatticeSpring(rsi);
		n2->Rsprings.pop_back();
	}
	// if there is a spring between n1 and n2, get the fibre
	ModelElementFibre *fn1n2 = NULL;
	for (int i = 0; i < n1->springs.size(); i++)
	{
		if (n1->springs.at(i)->nodes().at(0) == n2 ||
			n1->springs.at(i)->nodes().at(1) == n2)
		{
			fn1n2 = n1->springs.at(i)->fibre;
			break;
		}
	}
	if (fn1n2 == NULL)
	{
		// check if there is fibre forms a circle between n1 and n2
		// if so, cut it into two halves
		vector<ModelElementFibre*> cFn1n2;
		circleFibre(n1, n2, cFn1n2);
		if (cFn1n2.size() > 0)
		{
			for (int i = 0; i < cFn1n2.size(); i++)
			{
				ModelElementLatticeNode *cFsn = cFn1n2.at(i)->splitNode(n1, n2);
				if (cFsn != NULL)
				{
					cutFibre(cFn1n2.at(i), cFsn);
				}
			}
		}
	}
	// if n1 or n2 is a fibre end node
	bool fEnd = n1->fibreEnd || n2->fibreEnd;
	// the new position in the middle of n1 and n2
	double cn1n2x = (n1->position.x + n2->position.x) / 2;
	double cn1n2y = (n1->position.y + n2->position.y) / 2;
	double cn1n2z = (n1->position.z + n2->position.z) / 2;
	// check the fibre list of both n1 and n2 to find all common nodes
	// and the fibres which share these common nodes
	vector<ModelElementLatticeNode*> commonNodes;
	vector<ModelElementFibre*> connectedFibres;
	//vector<int> commonType; // -1: single common node for the pair of fibres; >=0: two common nodes for the pair of fibres
	//vector<bool> commonUse; // this is only for the case with two common nodes
	for (int i = 0; i < n1->fibres.size(); i++)
	{
		if (n1->fibres.at(i) == fn1n2)
			continue;
		ModelElementLatticeNode *n1n1 = n1->fibres.at(i)->n1;
		for (int j = 0; j < n2->fibres.size(); j++)
		{
			if (n2->fibres.at(j) == fn1n2)
				continue;
			if (n2->fibres.at(j)->checkNode(n1n1))
			{
				bool addc = true;
				for (int k = 0; k < commonNodes.size(); k++)
				{
					if (commonNodes.at(k) == n1n1)
					{
						addc = false;
						break;
					}
				}
				if (addc)
				{
					commonNodes.push_back(n1n1);
					connectedFibres.push_back(n1->fibres.at(i));
					connectedFibres.push_back(n2->fibres.at(j));
				}
			}
		}
		ModelElementLatticeNode *n1n2 = n1->fibres.at(i)->n2;
		for (int j = 0; j < n2->fibres.size(); j++)
		{
			if (n2->fibres.at(j) == fn1n2)
				continue;
			if (n2->fibres.at(j)->checkNode(n1n2))
			{
				bool addc = true;
				for (int k = 0; k < commonNodes.size(); k++)
				{
					if (commonNodes.at(k) == n1n2)
					{
						addc = false;
						break;
					}
				}
				if (addc)
				{
					commonNodes.push_back(n1n2);
					connectedFibres.push_back(n1->fibres.at(i));
					connectedFibres.push_back(n2->fibres.at(j));
				}
			}
		}
		for (int i1 = 0; i1 < n1->fibres.at(i)->nodes.size(); i1++)
		{
			if (n1->fibres.at(i) == fn1n2)
				continue;
			for (int j1 = 0; j1 < n2->fibres.size(); j1++)
			{
				if (n2->fibres.at(j1) == fn1n2)
					continue;
				if (n2->fibres.at(j1)->checkNode(n1->fibres.at(i)->nodes.at(i1)))
				{
					bool addc = true;
					for (int k = 0; k < commonNodes.size(); k++)
					{
						if (commonNodes.at(k) == n1->fibres.at(i)->nodes.at(i1))
						{
							addc = false;
							break;
						}
					}
					if (addc)
					{
						commonNodes.push_back(n1->fibres.at(i)->nodes.at(i1));
						connectedFibres.push_back(n1->fibres.at(i));
						connectedFibres.push_back(n2->fibres.at(j1));
					}
				}
			}
		}
	}
	// if the distance between two nodes are too close, we just merge them into one
	if (n1->checkNeighbor(n2))
	{
		// first remove/update the fibre
		if ((fn1n2->n1 == n1 && fn1n2->n2 == n2) ||
			(fn1n2->n1 == n2 && fn1n2->n2 == n1))
		{
			// this fibre will be removed from the system
			removeFibre(fn1n2);
		}
		else if ((fn1n2->n1 == n1 && fn1n2->n2 != n2) ||
				 (fn1n2->n2 == n1 && fn1n2->n1 != n2))
		{
			// this fibre has additional springs only starting from n2
			//   o----o-----o
			//   n1   n2
			// remove n2 from the fibre's node list
			// remove spring n1-n2 from the fibre's spring list 
			fn1n2->removeNode(n2);
			LatticeSpring *lsn1n2 = n1->locateSpring(n2);
			if (lsn1n2 == NULL)
				std::cout << "	-> Error: the spring n1-n2 is not found!" << std::endl;
			else
				fn1n2->removeSpring(lsn1n2);
			// change the end node of the fibre
			if (fn1n2->n1 == n1)
				fn1n2->n1 = n2;
			else if (fn1n2->n2 == n1)
				fn1n2->n2 = n2;
			n1->removeFibre(fn1n2);
			// update neighboring relationship between n1 and n2
			n1->removeNeighbor(n2);
			n2->removeNeighbor(n1);
			n1->removeSpring(lsn1n2);
			n2->removeSpring(lsn1n2);
			removeLatticeSpring(lsn1n2); // delete the spring n1-n2
		}
		else if ((fn1n2->n1 != n1 && fn1n2->n2 == n2) ||
				 (fn1n2->n2 != n1 && fn1n2->n1 == n2))
		{
			// this fibre has additional springs only starting from n1
			//   o-----o-----o
			//         n1    n2
			// remove n1 from the fibre's node list and change n1 to the end node
			// remove spring n1-n2 from the fibre's spring list
			fn1n2->removeNode(n1);
			LatticeSpring *lsn1n2 = n1->locateSpring(n2);
			if (lsn1n2 == NULL)
				std::cout << "	-> Error: the spring n1-n2 is not found!" << std::endl;
			else
				fn1n2->removeSpring(lsn1n2);
			// change the end node of the fibre
			if (fn1n2->n1 == n2)
				fn1n2->n1 = n1;
			else if (fn1n2->n2 == n2)
				fn1n2->n2 = n1;
			n2->removeFibre(fn1n2);
			// update neighboring relationship between n1 and n2
			n1->removeNeighbor(n2);
			n2->removeNeighbor(n1);
			n1->removeSpring(lsn1n2);
			n2->removeSpring(lsn1n2);
			removeLatticeSpring(lsn1n2); // delete the spring n1-n2
		}
		else if ((fn1n2->n1 != n1 && fn1n2->n2 != n2) ||
				 (fn1n2->n2 != n1 && fn1n2->n1 != n2))
		{
			// this fibre has additional springs starting from n1 and n2
			//  o-----o------o-------o
			//        n1     n2
			// remove n2 from the fibre's node list
			// remove spring n1-n2 from the fibre's spring list
			fn1n2->removeNode(n2);
			LatticeSpring *lsn1n2 = n1->locateSpring(n2);
			if (lsn1n2 == NULL)
				std::cout << "	-> Error: the spring n1-n2 is not found!" << std::endl;
			else
				fn1n2->removeSpring(lsn1n2);
			// update neighboring relationship between n1 and n2
			n1->removeNeighbor(n2);
			n2->removeNeighbor(n1);
			n1->removeSpring(lsn1n2);
			n2->removeSpring(lsn1n2);
			removeLatticeSpring(lsn1n2); // delete the spring n1-n2
		}
	}
	else
	{
		// there is no spring between n1 and n2
	}
	// update the coordination for n1 and n2
	n1->position.x = cn1n2x;
	n1->position.y = cn1n2y;
	n1->position.z = cn1n2z;
	n2->position.x = cn1n2x;
	n2->position.y = cn1n2y;
	n2->position.z = cn1n2z;

	// start to update/delete the fibre containing n2 
	for (int i = 0; i < commonNodes.size(); i++)
	{
		ModelElementLatticeNode *cNi = commonNodes.at(i);
		ModelElementFibre *cFn1i = connectedFibres.at(i * 2);
		ModelElementFibre *cFn2i = connectedFibres.at(i * 2 + 1);
		std::vector<LatticeSpring*> lsFn1i;
		std::vector<ModelElementLatticeNode*> lnFn1i;
		std::vector<LatticeSpring*> lsFn2i;
		std::vector<ModelElementLatticeNode*> lnFn2i;
		// locate the nodes and springs in each of the two fibres 
		// lsFn1i contains all springs between n1 and cNi
		// lnFn1i contains all nodes between n1 and cNi (including cNi but not n1)
		// the same for lsFn2i/lnFn2i
		cFn1i->locateNodeSpring(n1, cNi, lnFn1i, lsFn1i);
		cFn2i->locateNodeSpring(n2, cNi, lnFn2i, lsFn2i);
		// then remove all nodes and springs in fibre containing n2
		if (lnFn1i.size() == 1 && lsFn1i.size() == 1 && 
			lnFn2i.size() == 1 && lsFn2i.size() == 1)
		{
			// only single spring in the fibre to be merged
			// remove neighboring node relationship
			cNi->removeNeighbor(n2);
			n2->removeNeighbor(cNi);
			// remove neighboring spring relationship
			LatticeSpring *nls = lsFn2i.at(0); // only single spring
			double nlsY = nls->mYoung;
			double nlsP = nls->mPoisson;
			double nlsR = nls->mRadius;
			cNi->removeSpring(nls);
			n2->removeSpring(nls);
			// now cut fibre cFn2i into pieces (because the commone springs are merged into cFn1i)
			if ((cFn2i->n1 == cNi && cFn2i->n2 == n2) ||
				(cFn2i->n2 == cNi && cFn2i->n1 == n2))
			{
				// if fibre cFn2i has only one spring between cNi and n2
				// remove the fibre from the system
				removeFibre(cFn2i); // delete the fibre/springs from the system
			}
			else if ((cFn2i->n1 == cNi && cFn2i->n2 != n2) ||
				 	 (cFn2i->n2 == cNi && cFn2i->n1 != n2))
			{
				// fibre has additional springs only starting from n2
				// need to change the fibre's end node to n2
				if (cFn2i->n1 == cNi)
				{
					cFn2i->n1 = n2;
				}
				else
				{
					cFn2i->n2 = n2;
				}
				cNi->removeFibre(cFn2i);
				// then remove the merged spring from the fibre
				cFn2i->removeSpring(nls);
				cFn2i->removeNode(n2);
				// delete the merged spring
				removeLatticeSpring(nls);
				// delete the corresponding rotational springs
				removeNodeRSpringsL(cNi, n2);
				removeNodeRSpringsL(n2, cNi);
				// add this fibre into n1
				//n1->addFibre(cFn2i);
			}
			else if ((cFn2i->n1 != cNi && cFn2i->n2 == n2) ||
					 (cFn2i->n2 != cNi && cFn2i->n1 == n2))
			{
				// fibre has additional springs only starting from cNi
				// need to change the fibre's end node to cNi
				if (cFn2i->n2 == n2)
				{
					cFn2i->n2 = cNi;
				}
				else
				{
					cFn2i->n1 = cNi;
				}
				n2->removeFibre(cFn2i);
				// then remove the merged spring from the fibre
				cFn2i->removeSpring(nls);
				cFn2i->removeNode(cNi);
				// delete the merged springs
				removeLatticeSpring(nls);
				// delete the corresponding rotational springs
				removeNodeRSpringsL(cNi, n2);
				removeNodeRSpringsL(n2, cNi);
			}
			else if ((cFn2i->n1 != cNi && cFn2i->n2 != n2) ||
					 (cFn2i->n2 != cNi && cFn2i->n1 != n2))
			{
				// fibre has additional springs starting from cNi and n2
				// need to add one new fibre
				vector<ModelElementLatticeNode*> n1s;
				vector<LatticeSpring*> l1s;
				vector<ModelElementLatticeNode*> n2s;
				vector<LatticeSpring*> l2s;
				// l1s takes all springs between fibre_n1 and the inquired node
				// n1s takes all nodes between fibre_n1 and the inquired node (excluding fibre_n1, including the inquired node)
				// the same for l2s and n2s (between fibre_n2 and the inquired node)
				int orderType = cFn2i->relativePosition(n2, cNi, n1s, l1s, n2s, l2s); // determine the end node of the fibre
				if (orderType == 0)
				{
					// fibre_n1 - n2 - cNi - fibre_n2
					cFn2i->removeSpring(nls);
					cFn2i->removeNode(n2);
					cFn2i->removeNode(cNi);
					// delete the merged spring
					removeLatticeSpring(nls);
					// delete the corresponding rotational springs
					removeNodeRSpringsL(cNi, n2);
					removeNodeRSpringsL(n2, cNi);
					// create a new fibre
					ModelElementFibre *newF = new ModelElementFibre();
					newF->n1 = cNi;
					newF->n2 = cFn2i->n2;
					newF->n1Bound = false;
					newF->n2Bound = false;
					// add nodes/springs into the new fibre
					for (int i0 = 0; i0 < n2s.size(); i0++)
					{
						if (n2s.at(i0) != cNi)
						{
							newF->addNode(n2s.at(i0));
							cFn2i->removeNode(n2s.at(i0));
							n2s.at(i0)->removeFibre(cFn2i); // remove the old fibre
							n2s.at(i0)->addFibre(newF); // add the new fibre
						}
					}
					cNi->removeFibre(cFn2i);
					cNi->addFibre(newF);
					cFn2i->n2->removeFibre(cFn2i);
					cFn2i->n2->addFibre(newF);
					for (int i0 = 0; i0 < l2s.size(); i0++)
					{
						cFn2i->removeSpring(l2s.at(i0));
						newF->addSpring(l2s.at(i0));
						l2s.at(i0)->setFibre(newF);
					}
					fibres.push_back(newF); // add the new fibre into the system
					newF->fibre_id = fibres.size();
					// change the end node for the original fibre
					cFn2i->n2 = n2;
				}
				else if (orderType == 1)
				{
					// fibre_n1 - cNi - n2 - fibre_n2
					cFn2i->removeSpring(nls);
					cFn2i->removeNode(cNi);
					cFn2i->removeNode(n2);
					// delete the merged spring
					removeLatticeSpring(nls);
					// delete the corresponding rotational springs
					removeNodeRSpringsL(cNi, n2);
					removeNodeRSpringsL(n2, cNi);
					// create a new fibre
					ModelElementFibre *newF = new ModelElementFibre();
					newF->n1 = cFn2i->n1;
					newF->n2 = cNi;
					newF->n1Bound = false;
					newF->n2Bound = false;
					// add nodes/springs into the new fibre
					for (int i0 = 0; i0 < n1s.size(); i0++)
					{
						if (n1s.at(i0) != cNi)
						{
							newF->addNode(n1s.at(i0));
							cFn2i->removeNode(n1s.at(i0));
							n1s.at(i0)->removeFibre(cFn2i); // remove the old fibre
							n1s.at(i0)->addFibre(newF); // add the new fibre
						}
					}
					cFn2i->n1->removeFibre(cFn2i);
					cFn2i->n1->addFibre(newF);
					cNi->removeFibre(cFn2i);
					cNi->addFibre(newF);
					for (int i0 = 0; i0 < l1s.size(); i0++)
					{
						cFn2i->removeSpring(l1s.at(i0));
						newF->addSpring(l1s.at(i0));
						l1s.at(i0)->setFibre(newF);
					}
					fibres.push_back(newF); // add the new fibre into the system
					newF->fibre_id = fibres.size();
					// change the end node for the original fibre
					cFn2i->n1 = n2;
				}
				else
				{
					std::cout << "	-> Error: the two inquired nodes are not in the fibre!" << std::endl;
				}
			}
			// the merged spring will enhance the stiffness of the springs between n1 and cNi
			int mnn = (int)lnFn1i.size();
			for (int i0 = 0; i0 < lnFn1i.size(); i0++)
			{
				if (lnFn1i.at(i0) != cNi)
				{
					double ncxi = (cNi->position.x * (i0 + 1) + n1->position.x * (mnn - i0 - 1)) / mnn;
					double ncyi = (cNi->position.y * (i0 + 1) + n1->position.y * (mnn - i0 - 1)) / mnn;
					double nczi = (cNi->position.z * (i0 + 1) + n1->position.z * (mnn - i0 - 1)) / mnn;
					lnFn1i.at(i0)->position.x = ncxi;
					lnFn1i.at(i0)->position.y = ncyi;
					lnFn1i.at(i0)->position.z = nczi;
				}
			}
			for (int i0 = 0; i0 < lsFn1i.size(); i0++)
			{
				double lsiY = lsFn1i.at(i0)->mYoung;
				double lsiP = lsFn1i.at(i0)->mPoisson;
				double lsiR = lsFn1i.at(i0)->mRadius;
				double lsiPnew = (nlsP + lsiP) / 2;
				double lsiRnew = (nlsR + lsiR) / 2;
				// Young's modulus is recalculated for the merged springs 
				// k = pi*R^2*Y/L, then Y_new s.t. pi*R_new^2*Y_new = pi*R_1^2*Y_1 + pi*R_2^2*Y_2
				double lsiYnew = (lsiR * lsiR * lsiY + nlsR * nlsR * nlsY) / (lsiRnew * lsiRnew);
				lsFn1i.at(i0)->updateYoung(lsiYnew);
				lsFn1i.at(i0)->nodes().at(0)->mYoung = lsiYnew;
				lsFn1i.at(i0)->nodes().at(1)->mYoung = lsiYnew; // update Y for the spring and its nodes
				lsFn1i.at(i0)->updatePoisson(lsiPnew);
				lsFn1i.at(i0)->updateRadius(lsiRnew);
				(static_cast<LinearSpring*>(lsFn1i.at(i0)))->reset(); // reset the stiffness of the spring
			}
		}
		else
		{
			// more than one spring in the fibre, this would be disregarded
		}
	}
	// update the neighboring relationship for n1 and n2
	// move n2's neighboring nodes/springs into n1
	while (!n2->nodes.empty())
	{
		// n2's neighboring nodes, move them into n1
		ModelElementLatticeNode *n2nj = n2->nodes.back();
		// update all the rotational springs at n2nj
		for (int i = 0; i < n2nj->Rsprings.size(); i++)
		{
			RotationalSpring *ri = static_cast<RotationalSpring*>(n2nj->Rsprings.at(i));
			if (ri->nodes().at(0) == n2)
				ri->changeN1(n1);
			else if (ri->nodes().at(2) == n2)
				ri->changeN2(n1);
			ri->reset();
		}
		n2->removeNeighbor(n2nj);
		n2nj->removeNeighbor(n2); // remove neighboring nodes
		n1->addNeighbor(n2nj);
		n2nj->addNeighbor(n1); // add neighboring nodes
	}
	while (!n2->springs.empty())
	{
		// n2's neighboring springs, move them into n1
		LatticeSpring *n2sj = n2->springs.back();
		n2->removeSpring(n2sj); // remove neighboring springs
		// spring's end node also changes to n1 from n2
		if (n2sj->nodes().at(0) == n2)
		{
			(static_cast<LinearSpring*>(n2sj))->changeN1(n1);
			ModelElementEdge *n2sje = (static_cast<LinearSpring*>(n2sj))->getBoundEdge();
			n2sje->setPoint(&n1->position, 0);
			n2sje->e1Index = n1->mGlobalIndex;
		}
		else if (n2sj->nodes().at(1) == n2)
		{
			(static_cast<LinearSpring*>(n2sj))->changeN2(n1);
			ModelElementEdge *n2sje = (static_cast<LinearSpring*>(n2sj))->getBoundEdge();
			n2sje->setPoint(&n1->position, 1);
			n2sje->e2Index = n1->mGlobalIndex;
		}
		n1->addSpring(n2sj); // add neighboring springs
	}
	// move all n2's fibre into n1
	while (!n2->fibres.empty())
	{
		// n2's fibres, move them into n1
		ModelElementFibre *n2fj = n2->fibres.back();
		if (n2fj->n1 == n2)
		{
			n2fj->n1 = n1;
			n1->addFibre(n2fj);
		}
		else if (n2fj->n2 == n2)
		{
			n2fj->n2 = n1;
			n1->addFibre(n2fj);
		}
		else
		{
			if (n2fj->checkNode(n1) && n2fj->checkNode(n2))
			{
				std::cout << "	-> Error: there is another fibre with spring n1-n2!" << std::endl;
			}
			if (n2fj->checkNode(n2))
			{
				n2fj->removeNode(n2);
				n2fj->addNode(n1);
				n1->addFibre(n2fj);
			}
		}
		n2->removeFibre(n2fj);
	}
	// delete n2
	removeNode(n2);
	// add new rotational springs at n1
	for (int i = 0; i < n1->fibres.size(); i++)
	{
		addRotationalSpringNode(n1, n1->fibres.at(i));
	}
	/*if (fEnd)
	{
		std::vector<LatticeSpring*> tempLS;
		// rotational spring is constructed between fibres with n1 as end node
		for (int i = 0; i < n1->springs.size(); i++)
		{
			if (n1->springs.at(i)->fibre->n1 == n1 ||
				n1->springs.at(i)->fibre->n2 == n1)
			{
				// the end node of the fibre containing this spring is n1
				tempLS.push_back(n1->springs.at(i));
			}
		}
		if (tempLS.size() >= 2)
		{
			for (int i = 0; i < tempLS.size() - 1; i++)
			{
				ModelElementLatticeNode * tempn1 = tempLS.at(i)->nodes().at(0);
				if (tempn1 == n1)
					tempn1 = tempLS.at(i)->nodes().at(1);
				for (int j = i + 1; j < tempLS.size(); j++)
				{
					ModelElementLatticeNode * tempn2 = tempLS.at(j)->nodes().at(0);
					if (tempn2 == n1)
						tempn2 = tempLS.at(j)->nodes().at(1);
					RotationalSpring *rsij = new RotationalSpring(tempn1, tempn2, n1);
					n1->Rsprings.push_back(rsij);
					springs.push_back(rsij);
				}
			}
		}
	}*/
}

bool Lattice::addEdge(ModelElementEdge * to_add)
{
	edges.push_back(to_add);
	return true;
}

bool Lattice::removeEdge(ModelElementEdge * to_remove)
{
	if (to_remove != NULL)
		std::cout << "	-> Removing an edge" << std::endl;
	else
		return false;

	std::vector<ModelElementEdge *>::iterator found
		= std::find(edges.begin(), edges.end(), to_remove);

	// remove from the edges list
	if (found != edges.end())
		edges.erase(found);
	else
		std::cerr << "	-> Error in ModelLattice: Attempt to remove unregistered ModelElementEdge!\n" << std::endl;

	// remove from the bounding box
	mpBoundingBoxList->remove(to_remove);

	// delete edge
	to_remove->~ModelElementEdge();

	return true;
}

bool Lattice::removeLatticeSpring(LatticeSpring * to_remove)
{
	if (to_remove != NULL)
		std::cout << "	-> Removing an edge" << std::endl;
	else
		return false;

	std::vector<LatticeSpring *>::iterator found
		= std::find(springs.begin(), springs.end(), to_remove);

	if (found != springs.end())
		springs.erase(found);
	else
		std::cerr << "	-> Error in ModelLattice: Attempt to remove unregistered Spring!\n" << std::endl;

	// remove the boundingbox edge
	if (to_remove->mSpringType == LatticeSpring::TypeLinearSpring)
	{
		LinearSpring *lstr = static_cast<LinearSpring*>(to_remove);
		mpArena->removeObject(lstr->GLObject()); // remove from the arena

		ModelElementEdge *remove_edge = to_remove->getBoundEdge();
		std::vector<ModelElementEdge*>::iterator foundE
			= std::find(edges.begin(), edges.end(), remove_edge);
		if (foundE != edges.end())
			edges.erase(foundE);
		else
			std::cerr << "	-> Error in ModelLattice: Attempt to remove unregistered Edge!\n" << std::endl;
		// remove from the bounding box
		mpBoundingBoxList->remove(remove_edge);
		// delete edge
		remove_edge->~ModelElementEdge();
	}

	// delete spring
	to_remove->~LatticeSpring();

	return true;
}

bool Lattice::removeSpring(LatticeSpring * to_remove)
{
	// remove a spring in all structures: node, spring, fibre
	if (to_remove != NULL)
		std::cout << "	-> Removing an edge" << std::endl;
	else
		return false;

	// remove from the springs list
	std::vector<LatticeSpring *>::iterator found
		= std::find(springs.begin(), springs.end(), to_remove);

	if (found != springs.end())
		springs.erase(found);
	else
		std::cerr << "	-> Error in ModelLattice: Attempt to remove unregistered Spring!\n" << std::endl;

	if (to_remove->mSpringType == LatticeSpring::TypeLinearSpring)
	{
		// linear spring
		for (int i = 0; i < to_remove->nodes().size(); i++)
		{
			// remove spring from the nodes
			ModelElementLatticeNode * rNi = to_remove->nodes().at(i);
			rNi->removeSpring(to_remove);
		}
		// remove spring nodes from the nodes
		to_remove->nodes().at(1)->removeNeighbor(to_remove->nodes().at(0));
		to_remove->nodes().at(0)->removeNeighbor(to_remove->nodes().at(1));
	}
	else
	{
		// rotational spring
		ModelElementLatticeNode * rNc = to_remove->nodes().at(1); // only the center node
		rNc->removeRspring(to_remove);
	}

	// remove spring from fibre
	if (to_remove->mSpringType == LatticeSpring::TypeLinearSpring)
	{
		LinearSpring *to_removeL = (LinearSpring*)to_remove;
		ModelElementFibre *toF = to_removeL->fibre;
		// check the spring is one the end or inside the fibre
		if (to_removeL->nodes().at(0) == toF->n1 ||
			to_removeL->nodes().at(0) == toF->n2 ||
			to_removeL->nodes().at(1) == toF->n1 ||
			to_removeL->nodes().at(1) == toF->n2)
		{
			// on the two ends of the fibre
			toF->removeSpring(to_remove); // remove the spring from the fibre
			to_remove->fibre = NULL;
			if (to_removeL->nodes().at(0) == toF->n1)
			{
				ModelElementLatticeNode * newN1 = to_removeL->nodes().at(1);
				toF->removeNode(newN1); // remove the new end node from the fibre
				toF->n1 = newN1;
				toF->n1Bound = false;
			}
			else if (to_removeL->nodes().at(1) == toF->n1)
			{
				ModelElementLatticeNode * newN1 = to_removeL->nodes().at(0);
				toF->removeNode(newN1); // remove the new end node from the fibre
				toF->n1 = newN1;
				toF->n1Bound = false;
			}
			else if (to_removeL->nodes().at(0) == toF->n2)
			{
				ModelElementLatticeNode * newN2 = to_removeL->nodes().at(1);
				toF->removeNode(newN2); // remove the new end node from the fibre
				toF->n2 = newN2;
				toF->n2Bound = false;
			}
			else if (to_removeL->nodes().at(1) == toF->n2)
			{
				ModelElementLatticeNode * newN2 = to_removeL->nodes().at(0);
				toF->removeNode(newN2); // remove the new end node from the fibre
				toF->n2 = newN2;
				toF->n2Bound = false;
			}
		}
		else
		{
			// inside the fibre
			toF->removeSpring(to_remove);
			to_remove->fibre = NULL;
			// first check which fibre end-node can be reached from node0 of the spring
			ModelElementLatticeNode *removeN0next = to_remove->nodes().at(0);
			ModelElementLatticeNode *removeN0nextback = to_remove->nodes().at(0);
			bool hit = false;
			while (!hit)
			{
				for (int i0 = 0; i0 < removeN0next->nodes.size(); i0++)
				{
					if (removeN0next->nodes.at(i0) != removeN0nextback)
					{
						bool onFibre = false;
						for (int j0 = 0; j0 < removeN0next->nodes.at(i0)->fibres.size(); j0++)
						{
							if (removeN0next->nodes.at(i0)->fibres.at(j0) == toF)
								onFibre = true;
						}
						if (onFibre)
						{
							removeN0nextback = removeN0next;
							removeN0next = removeN0next->nodes.at(i0);
							break;
						}
					}
				}
				if (removeN0next == toF->n1 ||
					removeN0next == toF->n2)
					hit = true;
			}
			if (removeN0next == toF->n1)
			{
				/**********************
				    removed spring
				 o----o  -  o--------o
				 n1   N0    N1       n2
				       fibre

				 1. update fibre n1-n2
					to n1-N0
				 2. add one new fibre
					N1-n2
				 3. update lists of n2
				**********************/
				ModelElementFibre *newFibre = new ModelElementFibre(); // add a new fibre
				newFibre->n1 = to_remove->nodes().at(1);
				newFibre->n2 = toF->n2;
				newFibre->n1Bound = false;
				newFibre->n2Bound = true;
				// update the end-nodes of the old fibre
				toF->removeNode(to_remove->nodes().at(0)); // remove N0 from the old fibre
				toF->n2 = to_remove->nodes().at(0);
				// start from N1 to n2
				ModelElementLatticeNode *moveN1next = to_remove->nodes().at(1);
				ModelElementLatticeNode *moveN1nextback = to_remove->nodes().at(1);
				hit = false;
				while (!hit)
				{
					// move nodes & springs between N1 and n2 from the old fiber to the new fiber
					for (int i0 = 0; i0 < moveN1next->nodes.size(); i0++)
					{
						if (moveN1next->nodes.at(i0) != moveN1nextback)
						{
							bool onFibre = false;
							for (int j0 = 0; j0 < moveN1next->nodes.at(i0)->fibres.size(); j0++)
							{
								if (moveN1next->nodes.at(i0)->fibres.at(j0) == toF)
									onFibre = true;
							}
							if (onFibre)
							{
								moveN1nextback = moveN1next;
								moveN1next = moveN1next->nodes.at(i0);
								LatticeSpring* to_move = toF->locateSpring(moveN1nextback, moveN1next);
								toF->removeNode(moveN1nextback);
								toF->removeSpring(to_move); // remove node/spring from the old fibre
								if (moveN1nextback != to_remove->nodes().at(1)) // exclude the end-node
								{
									newFibre->addNode(moveN1nextback);
									moveN1nextback->removeFibre(toF);
									moveN1nextback->addFibre(newFibre); // renew the fibre flag
								}
								newFibre->addSpring(to_move); // add node/spring into the new fibre
								to_move->fibre = newFibre; // renew the fibre pointer
								break;
							}
						}
					}
					if (moveN1next == toF->n2)
						hit = true;
				}
				// renew the fibre flag for the end-nodes of the new fibre
				newFibre->n1->removeFibre(toF);
				newFibre->n1->addFibre(newFibre);
				newFibre->n2->removeFibre(toF);
				newFibre->n2->addFibre(newFibre);
			}
			else
			{
				/**********************
				removed spring
				o----o  -  o--------o
				n1   N1    N0       n2
				      fibre
				 1. update fibre n1-n2
					to n1-N1
				 2. add one new fibre 
					N0-n2
				 3. update lists of n2
				**********************/
				ModelElementFibre *newFibre = new ModelElementFibre(); // add a new fibre
				newFibre->n1 = to_remove->nodes().at(0);
				newFibre->n2 = toF->n2;
				newFibre->n1Bound = false;
				newFibre->n2Bound = true;
				// update the end-nodes of the old fibre
				toF->removeNode(to_remove->nodes().at(1)); // remove N1 from the old fibre
				toF->n2 = to_remove->nodes().at(1);
				// start from N0 to n2
				ModelElementLatticeNode *moveN0next = to_remove->nodes().at(1);
				ModelElementLatticeNode *moveN0nextback = to_remove->nodes().at(1);
				hit = false;
				while (!hit)
				{
					// move nodes & springs between N0 and n2 from the old fiber to the new fiber
					for (int i0 = 0; i0 < moveN0next->nodes.size(); i0++)
					{
						if (moveN0next->nodes.at(i0) != moveN0nextback)
						{
							bool onFibre = false;
							for (int j0 = 0; j0 < moveN0next->nodes.at(i0)->fibres.size(); j0++)
							{
								if (moveN0next->nodes.at(i0)->fibres.at(j0) == toF)
									onFibre = true;
							}
							if (onFibre)
							{
								moveN0nextback = moveN0next;
								moveN0next = moveN0next->nodes.at(i0);
								LatticeSpring* to_move = toF->locateSpring(moveN0nextback, moveN0next);
								toF->removeNode(moveN0nextback); // remove node/spring from the old fibre
								if (moveN0nextback != to_remove->nodes().at(0)) // exclude the end-node
								{
									newFibre->addNode(moveN0nextback);
									moveN0nextback->removeFibre(toF);
									moveN0nextback->addFibre(newFibre); // renew the fibre flag
								}
								newFibre->addSpring(to_move); // add node/spring into the new fibre
								to_move->fibre = newFibre; // renew the fibre pointer
								break;
							}
						}
					}
					if (moveN0next == toF->n2)
						hit = true;
				}
				// renew the fibre flag for the end-nodes of the new fibre
				newFibre->n1->removeFibre(toF);
				newFibre->n1->addFibre(newFibre);
				newFibre->n2->removeFibre(toF);
				newFibre->n2->addFibre(newFibre);
			}
		}
	}

	return true;
}

void Lattice::removeNodeRSprings(ModelElementLatticeNode * N)
{
	while ( !N->Rsprings.empty())
	{
		LatticeSpring *Si = N->Rsprings.back();
		// remove the spring from the lattice system
		std::vector<LatticeSpring*>::iterator found
			= std::find(springs.begin(), springs.end(), Si);
		if (found != springs.end())
			springs.erase(found);
		// remove the spring from the system
		removeLatticeSpring(Si);
		// remove the spring from the node list
		N->Rsprings.pop_back();
	}
}

void Lattice::removeNodeRSpringsL(ModelElementLatticeNode * N1, ModelElementLatticeNode * N2)
{
	// the center lattice node is N1
	for (int i = 0; i < N1->Rsprings.size(); i++)
	{
		RotationalSpring * Ri = static_cast<RotationalSpring *>(N1->Rsprings.at(i));
		if (Ri->nodes().at(1) != N1)
		{
			std::cout << "	-> Error: the center lattice node of the RotationalSpring is wrong!" << std::endl;
			continue;
		}
		if (Ri->nodes().at(0) == N2 || Ri->nodes().at(2) == N2)
		{
			// any RotationalSpring involving N2 which is to be detached
			N1->RspringsToDelete.push_back(Ri);
		}
	}
	while (!N1->RspringsToDelete.empty())
	{
		vector<LatticeSpring*>::iterator found
			= std::find(springs.begin(), springs.end(), N1->RspringsToDelete.back());
		if (found != springs.end())
			springs.erase(found);
		vector<LatticeSpring*>::iterator foundI
			= std::find(N1->Rsprings.begin(), N1->Rsprings.end(), N1->RspringsToDelete.back());
		if (foundI != N1->Rsprings.end())
			N1->Rsprings.erase(foundI);
		removeLatticeSpring(N1->RspringsToDelete.back());
		N1->RspringsToDelete.pop_back();
	}
}

bool Lattice::addFibre(ModelElementFibre * to_add)
{
	fibres.push_back(to_add);
	to_add->fibre_id = fibres.size();
	return true;
}

void Lattice::circleFibre(LatticeSpring * ls, vector<ModelElementFibre*> & fs)
{
	if (ls->mSpringType == LatticeSpring::TypeLinearSpring)
	{
		ModelElementFibre *lsF = ls->fibre;
		ModelElementLatticeNode *lsN1 = ls->nodes().at(0);
		ModelElementLatticeNode *lsN2 = ls->nodes().at(1);
		for (int i = 0; i < lsN1->fibres.size(); i++)
		{
			if (lsN1->fibres.at(i) == lsF)
				continue;

			for (int j = 0; j < lsN2->fibres.size(); j++)
			{
				if (lsN2->fibres.at(j) == lsF)
					continue;

				if (lsN1->fibres.at(i) == lsN2->fibres.at(j))
				{
					bool addf = true;
					for (int k = 0; k < fs.size(); k++)
					{
						if (fs.at(k) == lsN1->fibres.at(i))
						{
							addf = false;
							break;
						}
					}
					if (addf)
						fs.push_back(lsN1->fibres.at(i));
				}
			}
		}
	}
}

void Lattice::circleFibre(ModelElementLatticeNode *node1, ModelElementLatticeNode *node2, vector<ModelElementFibre*> & fs)
{
	for (int i = 0; i < node1->fibres.size(); i++)
	{
		for (int j = 0; j < node2->fibres.size(); j++)
		{
			if (node1->fibres.at(i) == node2->fibres.at(j))
			{
				bool addf = true;
				for (int k = 0; k < fs.size(); k++)
				{
					if (fs.at(k) == node1->fibres.at(i))
					{
						addf = false;
						break;
					}
				}
				if (addf)
					fs.push_back(node1->fibres.at(i));
			}
		}
	}
}

void Lattice::cutFibre(ModelElementFibre * to_cut, ModelElementLatticeNode * splitN)
{
	if (to_cut->nodes.size() == 0)
		return;

	if (!to_cut->checkNode(splitN))
		std::cout << "	-> Error: the cutting point node is not found in the fibre!" << std::endl;

	// remove all rotational springs in splitN
	while (!splitN->Rsprings.empty())
	{
		LatticeSpring *rsn = splitN->Rsprings.back();
		removeLatticeSpring(rsn);
		splitN->Rsprings.pop_back();
	}

	vector<ModelElementLatticeNode*> ns;
	vector<LatticeSpring*> ls;
	// ls takes all springs between fibre_n1 and the inquired node splitN
	// ns takes all nodes between fibre_n1 and splitN (excluding fibre_n1, including splitN)
	to_cut->relativeSinglePosition(splitN, ns, ls);
	ModelElementLatticeNode * fN1 = to_cut->n1;
	ModelElementLatticeNode * fN2 = to_cut->n2;
	// fN1 - splitN - fN2
	to_cut->removeNode(splitN);
	// create a new fibre
	ModelElementFibre * newF = new ModelElementFibre();
	newF->n1 = fN1;
	newF->n2 = splitN;
	newF->n1Bound = false;
	newF->n2Bound = false;
	// add nodes/springs into the new fibre
	for (int i0 = 0; i0 < ns.size(); i0++)
	{
		if (ns.at(i0) != splitN)
		{
			newF->addNode(ns.at(i0));
			to_cut->removeNode(ns.at(i0));
			ns.at(i0)->removeFibre(to_cut); // remove the old fibre
			ns.at(i0)->addFibre(newF); // add the new fibre
		}
	}
	fN1->removeFibre(to_cut);
	fN1->addFibre(newF);
	splitN->addFibre(newF);
	for (int i0 = 0; i0 < ls.size(); i0++)
	{
		to_cut->removeSpring(ls.at(i0));
		newF->addSpring(ls.at(i0));
		ls.at(i0)->setFibre(newF);
	}
	fibres.push_back(newF);
	newF->fibre_id = fibres.size();
	//  change the end node for the original fibre
	to_cut->n1 = splitN;
}

bool Lattice::addRotationalSpringNode(ModelElementLatticeNode * nQ, ModelElementFibre * fQ)
{
	if (!nQ->checkFibre(fQ)) // nQ is not in the fibre
		return false;

	if (fQ->n1 == nQ || fQ->n2 == nQ) // nQ is the end node of the fibre
		return false;

	ModelElementLatticeNode *node1 = NULL;
	ModelElementLatticeNode *node2 = NULL;
	bool n1take = false;
	bool n2take = false;

	for (int i = 0; i < nQ->nodes.size(); i++)
	{
		if (nQ->nodes.at(i)->checkFibre(fQ))
		{
			if (!n1take && !n2take)
			{
				node1 = nQ->nodes.at(i);
				n1take = true;
			}
			else if (n1take && !n2take)
			{
				node2 = nQ->nodes.at(i);
				n2take = true;
			}
			else if (n1take && n2take)
			{
				break;
			}
		}
	}

	if (node1 == NULL || node2 == NULL)
	{
		std::cout << "	-> Error: neighboring nodes of node " << nQ->mGlobalIndex << " are missing!" << std::endl;
		return false;
	}

	bool exist = true;
	for (int i = 0; i < nQ->Rsprings.size(); i++)
	{
		RotationalSpring *rI = static_cast<RotationalSpring*>(nQ->Rsprings.at(i));
		if ((rI->nodes().at(0) == node1 && rI->nodes().at(2) == node2 && rI->nodes().at(1) == nQ) ||
			(rI->nodes().at(0) == node2 && rI->nodes().at(2) == node1 && rI->nodes().at(1) == nQ))
		{
			rI->reset();
			exist = false;
		}
	}
	if (exist)
	{
		RotationalSpring *rQ = new RotationalSpring(node1, node2, nQ);
		nQ->Rsprings.push_back(rQ);
		addRotationalSpring(rQ);
	}

	return true;
}

bool Lattice::removeFibre(ModelElementFibre * to_remove)
{
	// remove a fibre in all structure: node, spring
	if (to_remove != NULL)
		std::cout << "	-> Removing a fibre " << to_remove->n1->mGlobalIndex << "-" << to_remove->n2->mGlobalIndex << " from the lattice system" << std::endl;
	else
		return false;

	// remove the springs from the fibre
	for (int i = 0; i < to_remove->springs.size(); i++)
	{
		LatticeSpring * Si = to_remove->springs.at(i);
		Si->nodes().at(0)->removeSpring(Si);
		Si->nodes().at(1)->removeSpring(Si); // remove spring from the two end-nodes
		Si->nodes().at(0)->removeNeighbor(Si->nodes().at(1));
		Si->nodes().at(1)->removeNeighbor(Si->nodes().at(0)); // remove neighbor-nodes from the two end-nodes
		//to_remove->removeSpring(Si); // remove spring from the fibre
		mSpringsToDelete.push_back(Si);
		Si->fibre = NULL; // remove the fibre flag
	}
	while (!mSpringsToDelete.empty())
	{
		// remove spring from the fibre
		to_remove->removeSpring(mSpringsToDelete.back());
		// remove the corresponding rotational springs
		removeNodeRSpringsL(mSpringsToDelete.back()->nodes().at(0), mSpringsToDelete.back()->nodes().at(1));
		removeNodeRSpringsL(mSpringsToDelete.back()->nodes().at(1), mSpringsToDelete.back()->nodes().at(0));
		// remove the linear spring from the system
		removeLatticeSpring(mSpringsToDelete.back());
		mSpringsToDelete.pop_back();
	}
	for (int i = 0; i < to_remove->nodes.size(); i++)
	{
		ModelElementLatticeNode * Ni = to_remove->nodes.at(i);
		if (Ni->fibres.size() == 1)
			mNodesToDelete.push_back(Ni); // delete the fibre nodes having no crosslink with other fibres 
		else if (Ni->fibres.size() > 1)
		{
			//ModelElementFibre * Fj = Ni->fibres.at(j);
			Ni->removeFibre(to_remove); // remove the fibre from the nodes
		}
		to_remove->removeNode(Ni); // remove the nodes from the fibre
	}
	to_remove->n1->removeFibre(to_remove);
	to_remove->n2->removeFibre(to_remove); // remove the fibre from the two end-nodes of the fibre

	// remove the single nodes from the system
	while (!mNodesToDelete.empty())
	{
		removeNode(mNodesToDelete.back());
		mNodesToDelete.pop_back();
	}

	std::vector<ModelElementFibre*>::iterator foundF
		= std::find(fibres.begin(), fibres.end(), to_remove);
	fibres.erase(foundF);

	// delete fibre
	to_remove->~ModelElementFibre();

	return true;
}

void Lattice::freeFibre(ModelElementFibre *F, int eIndex)
{
	if (eIndex == 1)
	{
		// n1
	}
	else
	{
		// n2
	}
}

void Lattice::createNetwork(double cutoff, int subDivide)
{
	if (subDivide < 1) // the fibre segment is at least 1 
		subDivide = 1;
	// Preliminary test: only linear springs first
	// delaunay triangulation of the Lattice Nodes
	CD3DW *Del = new CD3DW();
	int nsize = (int)nodes.size();
	for (int i = 0; i < nsize; i++)
	{
		// add Lattice Nodes into the system one-by-one
		Del->add(nodes.at(i)->position.x,
				nodes.at(i)->position.y,
				nodes.at(i)->position.z,
				nodes.at(i)->mRadius,
				i,
				&nodes.at(i));
	}
	CComplex *Complex = new CComplex(Del);
	// add LinearSpring
	int Pn = (int)Complex->P.size();
	for (int i = 0; i < Pn; i++)
	{
		CSimplex *p = Complex->P[i];
		int size = p->size();
		if (size == 2) // edges
		{
			C1Simplex *p1 = (C1Simplex*)p;
			int n1 = p1->A->n();
			int n2 = p1->B->n();
			Vector3f v12 = nodes.at(n1)->position - nodes.at(n2)->position;
			double dist12 = v12.Norm();
			if (dist12 > cutoff) continue;
			// add LinearSpring
			//double toss = ((double)rand() / (RAND_MAX));
			//if (toss > 0.8) continue;
			
			// check first if this edge has overlap with existing edges
			// these are stored into two categories: springs and nodes
			// for the crossed spring, if the crossed point is too close 
			// to the end node ( < fibreSize * mergeFract), the end node
			// would be taken as a overlapSpringEndNodes
			bool addFibre = true;
			// crossed springs
			vector<LinearSpring*> overlapSprings; // if any overlapped spring found
			vector<Vector3f> overlapSpringsCrossedPoints; // crossed points of the overlapped spring list above
			vector<double> overlapSpringsFraction; // the fraction of the crossed points
			vector<bool> overlapSpringsUse; // if we do add a cross-link node with this spring
			// crossed nodes
			vector<ModelElementLatticeNode*> overlapSpringEndNodes; // the crossed points are the end nodes of other springs
			vector<Vector3f> overlapSpringEndNodesPoints; // the crossed points
			// end node crosses with other spring
			vector<LinearSpring*> overlapEndNodeSpring; // the end node crosses with other springs
			vector<Vector3f> overlapEndNodesPoints; // the crossed points
			vector<int> overlapEndNodes; // 0: n1; 1: n2
			vector<bool> overlapEndNodeSpringUse; // if we do add a cross-link node with this spring
			// bounding box limits
			double noder = (nodes.at(n1)->mRadius + nodes.at(n2)->mRadius) / 2;
			double nodex1 = nodes.at(n1)->position.x;
			double nodey1 = nodes.at(n1)->position.y;
			double nodez1 = nodes.at(n1)->position.z;
			double nodex2 = nodes.at(n2)->position.x;
			double nodey2 = nodes.at(n2)->position.y;
			double nodez2 = nodes.at(n2)->position.z;
			// bounding box of the edge under check
			BoundingBox bbox;
			bbox.xmin = (nodex1 < nodex2 ? nodex1 : nodex2) - noder;
			bbox.xmax = (nodex1 > nodex2 ? nodex1 : nodex2) + noder;
			bbox.ymin = (nodey1 < nodey2 ? nodey1 : nodey2) - noder;
			bbox.ymax = (nodey1 > nodey2 ? nodey1 : nodey2) + noder;
			bbox.zmin = (nodez1 < nodez2 ? nodez1 : nodez2) - noder;
			bbox.zmax = (nodez1 > nodez2 ? nodez1 : nodez2) + noder;
			ModelElement::Type edgeType = ModelElement::TypeCellLatticeSpring;
			ModelElement::Type vesselType = ModelElement::TypeVesselSphere;
			ModelElement::Type cellType = ModelElement::TypeCellSpherical;
			CSListContainer<ModelElement *> possiblyCrossedCells;
			unsigned long index = mpBoundingBoxList->mStart;
			/************************
			 if the inquired bbox
			 is completely burried 
			 by a box
			 o-------------o
			     o-------o bbox
			************************/
			bool largeCover = false;
			ModelElement *largeCE = NULL;
			while (*mpBoundingBoxList->mppLimits[index] < bbox.xmin)
			{
				// added by Jieling
				BoundingBox *testBB = mpBoundingBoxList->mppElements[index]->getBoundingBox();
				if (testBB->xmax > bbox.xmax)
				{
					largeCover = true;
					largeCE = mpBoundingBoxList->mppElements[index];
					break;
				}

				if ((index = mpBoundingBoxList->mpForward[index]) == mpBoundingBoxList->mEnd)
					break;
			}
			if (index != mpBoundingBoxList->mEnd)
			{
				// find all limits up to bbox.xmax
				do
				{
					if (mpBoundingBoxList->mppElements[index]->mType != edgeType &&
						mpBoundingBoxList->mppElements[index]->mType != vesselType &&
						mpBoundingBoxList->mppElements[index]->mType != cellType)
						continue; // if it is not hepatocyte or linear spring or vessel sphere, revised by Jieling

					if (mpBoundingBoxList->mppElements[index]->mType == edgeType)
					{
						ModelElementEdge *edgeTest = static_cast<ModelElementEdge*>(mpBoundingBoxList->mppElements[index]);
						LinearSpring *boundSpring = static_cast<LinearSpring*>(edgeTest->getBoundSpring());
						if (boundSpring->nodes().at(0) == nodes.at(n1) ||
							boundSpring->nodes().at(1) == nodes.at(n1) ||
							boundSpring->nodes().at(0) == nodes.at(n2) ||
							boundSpring->nodes().at(1) == nodes.at(n2))
							continue; // if it is the spring directly connected with the spring		
					}

					if (*mpBoundingBoxList->mppLimits[index] > bbox.xmax)
					{
						if (largeCover)
						{
							if (mpBoundingBoxList->mppElements[index] == largeCE)
								break;
						}
						else
						{
							break;
						}
					}

					// test if found element's other bb limitis are in range
					BoundingBox *testBB = mpBoundingBoxList->mppElements[index]->getBoundingBox();

					if (testBB->xmin > bbox.xmax || bbox.xmin > testBB->xmax)
						continue;

					if (testBB->ymin > bbox.ymax || bbox.ymin > testBB->ymax)
						continue;

					if (mpBoundingBoxList->mDimensions == 3)
						if (testBB->zmin > bbox.zmax || bbox.zmin > testBB->zmax)
							continue;

					if (mpBoundingBoxList->mppElements[index]->mType == cellType ||
						mpBoundingBoxList->mppElements[index]->mType == edgeType ||
						mpBoundingBoxList->mppElements[index]->mType == vesselType)
					{
						bool addCC = true;
						for (ModelElement **chIter = possiblyCrossedCells.begin();
							chIter != possiblyCrossedCells.end();
							++chIter)
						{
							if ((*chIter) == mpBoundingBoxList->mppElements[index])
							{
								addCC = false;
								break;
							}
						}
						if (addCC)
							possiblyCrossedCells.push_back(mpBoundingBoxList->mppElements[index]);
					}

				} while ((index = mpBoundingBoxList->mpForward[index]) != mpBoundingBoxList->mEnd);
			}
			// end of searching crossed cell
			for (ModelElement **iter = possiblyCrossedCells.begin();
				iter != possiblyCrossedCells.end();
				++iter)
			{
				if ((*iter)->mType == cellType) // cell-spring interaction test
				{
					CellSpherical *cellI = static_cast<CellSpherical*>(*iter);
					if (CSModelTools::edgeIntersectSpere(nodes.at(n1)->position, nodes.at(n2)->position, 0.0848,
						cellI->position, cellI->mRadius)) // temporarily using 0.0848 as spring radius
					{
						addFibre = false;
						break; // overlapped with existing cells, do not add this fibre
					}
				}
				else if ((*iter)->mType == vesselType) // vessel-spring interaction test
				{
					ModelElementVesselSphere * vesselI = static_cast<ModelElementVesselSphere*>(*iter);
					if (CSModelTools::edgeIntersectSpere(nodes.at(n1)->position, nodes.at(n2)->position, 0.0848,
						vesselI->position, vesselI->mRadius)) // temporarily using 0.0848 as spring radius
					{
						addFibre = false;
						break; // overlapped with existing sinusoid, do not add this fibre
					}
				}
				else if ((*iter)->mType == edgeType) // spring-spring interaction test
				{
					ModelElementEdge * edgeI = static_cast<ModelElementEdge*>(*iter);
					LatticeSpring *springI = edgeI->getBoundSpring();
					LinearSpring *lspringI = static_cast<LinearSpring*>(springI);
					double lspringL = lspringI->getLength();
					Vector3f c1, c2, dir2; // unit direction from potential crossed spring to the spring under test
					double sc1, tc2;
					double distI = CSModelTools::GetDistanceTwoLines(nodes.at(n1)->position, nodes.at(n2)->position,
						springI->nodes().at(0)->position, springI->nodes().at(1)->position, dir2, sc1, tc2, c1, c2); // distance between two springs
					if (distI < lspringI->mRadius * 2) // if these two springs are overlapped
					{
						addFibre = true;
						// check if nodes need to be merged using both fraction and absolute distance
						// if the distance between two nodes are too close, even if the fraction is okay
						// we still merge them
						double sc1D = sc1 * dist12;
						double tc2D = tc2 * lspringL;
						double Db = fibreSize * mergeFract;
						int sStatus = 0; // 0: close to n1; 1: between n1 and n2; 2: close to n2
						int tStatus = 0; // 0: close to n1; 1: between n1 and n2; 2: close to n2
						if (sc1 < mergeFract)
						{	
							sStatus = 0;
						}
						else if (sc1 > mergeFract && sc1 < 1. - mergeFract)
						{
							sStatus = 1;
							if (sc1D < Db)
							{
								sStatus = 0;
							}
							else if (sc1D > dist12 - Db)
							{
								sStatus = 2;
							}
						}
						else if (sc1 > 1. - mergeFract)
						{
							sStatus = 2;
						}
						if (tc2 < mergeFract)
						{
							tStatus = 0;
						}
						else if (tc2 > mergeFract && tc2 < 1. - mergeFract)
						{
							tStatus = 1;
							if (tc2D < Db)
							{
								tStatus = 0;
							}
							else if (tc2D > lspringL - Db)
							{
								tStatus = 2;
							}
						}
						else if (tc2 > 1. - mergeFract)
						{
							tStatus = 2;
						}
						
						if (sStatus == 1 &&
							tStatus == 1)
						{
							bool addN = true;
							// regular one: crossed springs
							for (int i0 = 0; i0 < overlapSpringEndNodes.size(); i0++)
							{
								if (lspringI->nodes().at(0) == overlapSpringEndNodes.at(i0) ||
									lspringI->nodes().at(1) == overlapSpringEndNodes.at(i0))
								{
									addN = false;
									break;
								}
							}
							if (addN)
							{
								overlapSprings.push_back(lspringI); // add one overlapped spring
								overlapSpringsCrossedPoints.push_back(c1);
								overlapSpringsCrossedPoints.push_back(c2); // add overlapped points
								overlapSpringsFraction.push_back(sc1);
								overlapSpringsFraction.push_back(tc2); // add fraction of overlapped points
								overlapSpringsUse.push_back(true); // default: true
							}
						}
						else if (sStatus == 1 &&
								 tStatus == 0)
						{
							// particular one: crossed with end node n1
							std::vector<ModelElementLatticeNode*>::iterator foundLN =
								std::find(overlapSpringEndNodes.begin(), overlapSpringEndNodes.end(), lspringI->nodes().at(0));
							if (foundLN == overlapSpringEndNodes.end())
							{
								// check neighbor node relationship
								if (lspringI->nodes().at(0)->checkNeighbor(nodes.at(n1)) ||
									lspringI->nodes().at(0)->checkNeighbor(nodes.at(n2)))
								{
									addFibre = false;
								}
								for (int i0 = 0; i0 < overlapSprings.size(); i0++)
								{
									if (overlapSprings.at(i0)->nodes().at(0) == lspringI->nodes().at(0) ||
										overlapSprings.at(i0)->nodes().at(1) == lspringI->nodes().at(0))
									{
										// the crossed node is also the end node of one crossed spring
										overlapSpringsUse.at(i0) = false; // then do not consider this crossed spring
										break;
									}
								}
								for (int i0 = 0; i0 < overlapEndNodeSpring.size(); i0++)
								{
									if (overlapEndNodeSpring.at(i0)->nodes().at(0) == lspringI->nodes().at(0) ||
										overlapEndNodeSpring.at(i0)->nodes().at(1) == lspringI->nodes().at(0))
									{
										// the crossed node is also the end node of one crossed spring
										overlapEndNodeSpringUse.at(i0) = false; // then do not consider this crossed spring
									}
								}
								overlapSpringEndNodes.push_back(lspringI->nodes().at(0));
								overlapSpringEndNodesPoints.push_back(c1);
							}
							else
							{
								// the node is already a crossed one 
								// check neighbor node relationship
								if (lspringI->nodes().at(0)->checkNeighbor(nodes.at(n1)) ||
									lspringI->nodes().at(0)->checkNeighbor(nodes.at(n2)))
								{
									addFibre = false;
								}
							}
						}
						else if (sStatus == 1 &&
								 tStatus == 2)
						{
							// particular one: crossed with end node n2
							std::vector<ModelElementLatticeNode*>::iterator foundLN =
								std::find(overlapSpringEndNodes.begin(), overlapSpringEndNodes.end(), lspringI->nodes().at(1));
							if (foundLN == overlapSpringEndNodes.end())
							{
								// check neighbor node relationship
								if (lspringI->nodes().at(1)->checkNeighbor(nodes.at(n1)) ||
									lspringI->nodes().at(1)->checkNeighbor(nodes.at(n2)))
								{
									addFibre = false;
								}
								for (int i0 = 0; i0 < overlapSprings.size(); i0++)
								{
									if (overlapSprings.at(i0)->nodes().at(0) == lspringI->nodes().at(1) ||
										overlapSprings.at(i0)->nodes().at(1) == lspringI->nodes().at(1))
									{
										// the crossed node is also the end node of one crossed spring
										overlapSpringsUse.at(i0) = false; // then do not consider this crossed spring
										break;
									}
								}
								for (int i0 = 0; i0 < overlapEndNodeSpring.size(); i0++)
								{
									if (overlapEndNodeSpring.at(i0)->nodes().at(0) == lspringI->nodes().at(1) ||
										overlapEndNodeSpring.at(i0)->nodes().at(1) == lspringI->nodes().at(1))
									{
										// the crossed node is also the end node of one crossed spring
										overlapEndNodeSpringUse.at(i0) = false; // then do not consider this crossed spring
									}
								}
								overlapSpringEndNodes.push_back(lspringI->nodes().at(1));
								overlapSpringEndNodesPoints.push_back(c1);
							}
							else
							{
								// the node is already a crossed one 
								// check neighbor node relationship
								if (lspringI->nodes().at(1)->checkNeighbor(nodes.at(n1)) ||
									lspringI->nodes().at(1)->checkNeighbor(nodes.at(n2)))
								{
									addFibre = false;
								}
							}
						}
						else if (sStatus == 0 &&
								 tStatus == 1)
						{
							// check the neighbor node relationship
							if (lspringI->nodes().at(0)->checkNeighbor(nodes.at(n1)) ||
								lspringI->nodes().at(1)->checkNeighbor(nodes.at(n1)))
							{
								addFibre = false;
							}
							// particular one: end node n1 crosses with the spring
							bool addN = true;
							for (int i0 = 0; i0 < overlapSpringEndNodes.size(); i0++)
							{
								if (lspringI->nodes().at(0) == overlapSpringEndNodes.at(i0) ||
									lspringI->nodes().at(1) == overlapSpringEndNodes.at(i0))
								{
									addN = false;
									break;
								}
							}
							if (addN)
							{
								overlapEndNodeSpring.push_back(lspringI); // add the overlapped spring
								overlapEndNodesPoints.push_back(c2); // add overlapped point
								overlapEndNodes.push_back(0); // n1 is the cross-link node
								overlapEndNodeSpringUse.push_back(true);
							}
						}
						else if (sStatus == 2 &&
								 tStatus == 1)
						{
							// check the neighbor node relationship
							if (lspringI->nodes().at(0)->checkNeighbor(nodes.at(n2)) ||
								lspringI->nodes().at(1)->checkNeighbor(nodes.at(n2)))
							{
								addFibre = false;
							}
							// particular one: end node n2 crosses with the spring
							bool addN = true;
							for (int i0 = 0; i0 < overlapSpringEndNodes.size(); i0++)
							{
								if (lspringI->nodes().at(0) == overlapSpringEndNodes.at(i0) ||
									lspringI->nodes().at(1) == overlapSpringEndNodes.at(i0))
								{
									addN = false;
									break;
								}
							}
							if (addN)
							{
								overlapEndNodeSpring.push_back(lspringI); // add the overlapped spring
								overlapEndNodesPoints.push_back(c2); // add overlapped point
								overlapEndNodes.push_back(1); // n2 is the cross-link node
								overlapEndNodeSpringUse.push_back(true);
							}
						}
						else if ((sStatus == 0 && tStatus == 0) ||
							(sStatus == 0 && tStatus == 2) ||
							(sStatus == 2 && tStatus == 0) ||
							(sStatus == 2 && tStatus == 2))
						{
							std::cout << "	-> Error: the spring end nodes are too close to each other! Try resampling the lattice nodes!" << std::endl;
						}
					}
				}
			}
			if (addFibre)
			{
				if (overlapSprings.size() > 0 || 
					overlapSpringEndNodes.size() > 0 ||
					overlapEndNodeSpring.size() > 0)
				{
					// first create a fibre and the basic spring for the new fibre
					ModelElementFibre *fNew = new ModelElementFibre();
					fNew->n1 = nodes.at(n1);
					fNew->n2 = nodes.at(n2); // add the two nodes into the fibre
					fNew->n1Bound = true;
					fNew->n2Bound = true;
					nodes.at(n1)->addFibre(fNew);
					nodes.at(n2)->addFibre(fNew);
					fibres.push_back(fNew);
					fNew->fibre_id = fibres.size();
					// create the first spring and add into the fibre
					LinearSpring *ls1 = new LinearSpring(nodes.at(n1), nodes.at(n2));
					ls1->index = springs.size();
					ls1->n1_unbound = false;
					ls1->n2_unbound = false;
					nodes.at(n1)->addSpring(ls1);
					nodes.at(n2)->addSpring(ls1); // add neighboring springs
					nodes.at(n1)->addNeighbor(nodes.at(n2));
					nodes.at(n2)->addNeighbor(nodes.at(n1)); // add neighbors
					springs.push_back(ls1);
					// add spring object into the arena
					mpArena->addObject(ls1->GLObject());
					fNew->addSpring(ls1); // add the spring to the fibre
					// add one bounding box edge to the bounding box
					ModelElementEdge *e11 = new ModelElementEdge(nodes.at(n1)->position.x,
						nodes.at(n1)->position.y,
						nodes.at(n1)->position.z);
					e11->setPoint(&nodes.at(n1)->position, 0);
					e11->setPoint(&nodes.at(n2)->position, 1);
					e11->e1Index = nodes.at(n1)->mGlobalIndex;
					e11->e2Index = nodes.at(n2)->mGlobalIndex;
					e11->setRadius((nodes.at(n1)->mRadius + nodes.at(n2)->mRadius) / 2);
					e11->setBoundingBox(e11->mRadius);
					e11->setBoundSpring(ls1);
					mpBoundingBoxList->add(e11); // add boundingbox edge to the boundingbox
					ls1->changeBoundEdge(e11);
					ls1->setFibre(fNew); // set the fibre for the edge
					edges.push_back(e11); // add boundingbox edge into the edges
					// add the fibre, and add cross-link to all overlapped springs
					for (int i0 = 0; i0 < overlapSprings.size(); i0++)
					{
						if (!overlapSpringsUse.at(i0)) // skip the crossed spring which is not considered
							continue;
						// divide one fibre into two
						// remove neighboring relationship
						// renew rotational spring orientation
						LinearSpring *lspringI = overlapSprings.at(i0);
						LatticeSpring *springI = static_cast<LatticeSpring*>(lspringI);
						ModelElementEdge *edgeI = springI->getBoundEdge();
						ModelElementFibre *fibreI = springI->fibre;
						ModelElementLatticeNode *n1I = lspringI->nodes().at(0);
						ModelElementLatticeNode *n2I = lspringI->nodes().at(1);
						// add a new lattice node, the coordinate for the crosslink
						// i0*2 for c1 and i0*2+1 for c2
						double c2x = (overlapSpringsCrossedPoints.at(i0 * 2).x + overlapSpringsCrossedPoints.at(i0 * 2 + 1).x) / 2;
						double c2y = (overlapSpringsCrossedPoints.at(i0 * 2).y + overlapSpringsCrossedPoints.at(i0 * 2 + 1).y) / 2;
						double c2z = (overlapSpringsCrossedPoints.at(i0 * 2).z + overlapSpringsCrossedPoints.at(i0 * 2 + 1).z) / 2;
						double c2r = (overlapSprings.at(i0)->nodes().at(0)->mRadius + overlapSprings.at(i0)->nodes().at(1)->mRadius) / 2;
						ModelElementLatticeNode *nodeC = new ModelElementLatticeNode(c2x, c2y, c2z, c2r);
						// now we divide the original spring n1I - n2I into two 
						// springs: n1I - nodeC, nodeC - n2I
						// the original fibre fibreI has only one spring n1I - n2I
						// now fibreI has two springs: n1I - nodeC and nodeC - n2I
						nodeC->mYoung = (n1I->mYoung + n2I->mYoung) / 2;
						nodeC->mPoisson = (n1I->mPoisson + n2I->mPoisson) / 2;
						nodeC->fibreEnd = false; // node inside the fibre
						addNode(nodeC);
						nodeC->latticeIndex = nodes.size();
						// remove the neighboring node relationship
						n1I->removeNeighbor(n2I);
						n2I->removeNeighbor(n1I);
						// n1 still keeps the spring
						n2I->removeSpring(springI);
						lspringI->changeN2(nodeC); // change the end node of the crossed spring
						edgeI->setPoint(&nodeC->position, 1); // change the limits of the edge boundingbox
						edgeI->e2Index = nodeC->mGlobalIndex;
						fibreI->addNode(nodeC); // add the crossed node into the fibre
						nodeC->addFibre(fibreI);
						// add one new spring for the fibre
						LinearSpring *lsNI = new LinearSpring(nodeC, n2I);
						lsNI->SetColor(0., 1., 0., 1.); // green
						lsNI->index = springs.size();
						lsNI->n1_unbound = false;
						lsNI->n2_unbound = false; // not free nodes
						springs.push_back(lsNI);
						// add neighboring spring relationship
						nodeC->addSpring(springI);
						nodeC->addSpring(lsNI);
						n2I->addSpring(lsNI); // add the new springs into n2I and nodeC
						// add neighboring node relationship
						n1I->addNeighbor(nodeC);
						nodeC->addNeighbor(n1I);
						n2I->addNeighbor(nodeC);
						nodeC->addNeighbor(n2I);
						// add spring object into the arena
						mpArena->addObject(lsNI->GLObject());
						fibreI->addSpring(lsNI); // add spring to the fibre
						lsNI->fibre = fibreI;
						// add one boundingbox edge to the bounding box
						ModelElementEdge *eNI = new ModelElementEdge(nodeC->position.x, nodeC->position.y, nodeC->position.z);
						eNI->setPoint(&nodeC->position, 0);
						eNI->setPoint(&n2I->position, 1);
						eNI->e1Index = nodeC->mGlobalIndex;
						eNI->e2Index = n2I->mGlobalIndex;
						eNI->setRadius((nodeC->mRadius + n2I->mRadius) / 2);
						eNI->setBoundingBox(eNI->mRadius);
						eNI->setBoundSpring(lsNI);
						lsNI->changeBoundEdge(eNI);
						mpBoundingBoxList->add(eNI); // add boundingbox edge to the boundingbox
						edges.push_back(eNI);
						// cut the new fibre into several spring segment
						LatticeSpring *cutLs = NULL; // the spring to be cut
						for (int i1 = 0; i1 < fNew->springs.size(); i1++)
						{
							// c1 is on fNew
							double fs1x = fNew->springs.at(i1)->nodes().at(0)->position.x;
							double fs1y = fNew->springs.at(i1)->nodes().at(0)->position.y;
							double fs1z = fNew->springs.at(i1)->nodes().at(0)->position.z;
							double fs2x = fNew->springs.at(i1)->nodes().at(1)->position.x;
							double fs2y = fNew->springs.at(i1)->nodes().at(1)->position.y;
							double fs2z = fNew->springs.at(i1)->nodes().at(1)->position.z;
							double c1x = overlapSpringsCrossedPoints.at(i0 * 2).x;
							double c1y = overlapSpringsCrossedPoints.at(i0 * 2).y;
							double c1z = overlapSpringsCrossedPoints.at(i0 * 2).z;
							double fractI = CSModelTools::fractionPointLine(c1x, c1y, c1z,
								fs1x, fs1y, fs1z,
								fs2x, fs2y, fs2z); // the fraction of crossed point on spring_i1 of fNew
							if (fractI >= 0 && fractI <= 1)
							{
								cutLs = fNew->springs.at(i1);
								break;
							}
						}
						if (cutLs == NULL)
						{
							std::cout << "	Error: the spring to be crossed is not found!" << std::endl;
						}
						else
						{
							// now cut spring cutLs
							LinearSpring *cutLsl = static_cast<LinearSpring*>(cutLs);
							ModelElementEdge *cutEI = cutLsl->getBoundEdge();
							ModelElementLatticeNode *cutN1I = cutLsl->nodes().at(0);
							ModelElementLatticeNode *cutN2I = cutLsl->nodes().at(1);
							// remove the neighboring node relationship
							cutN1I->removeNeighbor(cutN2I);
							cutN2I->removeNeighbor(cutN1I);
							// cutN1I still keeps the spring
							cutN2I->removeSpring(cutLs);
							cutLsl->changeN2(nodeC); // change the end node of the crossed spring
							cutEI->setPoint(&nodeC->position, 1); // change the limits of the edge boundingbox
							cutEI->e2Index = nodeC->mGlobalIndex;
							fNew->addNode(nodeC); // add the crossed node in the fibre fNew
							nodeC->addFibre(fNew);
							// add one new spring for the fibre fNew
							LinearSpring *lsNC = new LinearSpring(nodeC, cutN2I);
							lsNC->SetColor(0., 1., 0., 1.); // green
							lsNC->index = springs.size();
							lsNC->n1_unbound = false;
							lsNC->n2_unbound = false; // not free nodes
							springs.push_back(lsNC);
							// add neighboring spring relationship
							nodeC->addSpring(cutLs);
							nodeC->addSpring(lsNC);
							cutN2I->addSpring(lsNC); // add the new springs into cutN2I and nodeC
							// add neighboring node relationship
							cutN1I->addNeighbor(nodeC);
							nodeC->addNeighbor(cutN1I);
							cutN2I->addNeighbor(nodeC);
							nodeC->addNeighbor(cutN2I);
							// add spring object into the arena
							mpArena->addObject(lsNC->GLObject());
							fNew->addSpring(lsNC); // add spring to the fibre
							lsNC->fibre = fNew;
							// add one boundingbox edge to the bounding box
							ModelElementEdge *eNC = new ModelElementEdge(nodeC->position.x, nodeC->position.y, nodeC->position.z);
							eNC->setPoint(&nodeC->position, 0);
							eNC->setPoint(&cutN2I->position, 1);
							eNC->e1Index = nodeC->mGlobalIndex;
							eNC->e2Index = cutN2I->mGlobalIndex;
							eNC->setRadius((nodeC->mRadius + cutN2I->mRadius) / 2);
							eNC->setBoundingBox(eNC->mRadius);
							eNC->setBoundSpring(lsNC);
							lsNC->changeBoundEdge(eNC);
							mpBoundingBoxList->add(eNC); // add boundingbox edge to the boundingbox
							edges.push_back(eNC);
						}
					}
					for (int i0 = 0; i0 < overlapSpringEndNodes.size(); i0++)
					{
						double c2x = (overlapSpringEndNodesPoints.at(i0).x + overlapSpringEndNodes.at(i0)->position.x) / 2;
						double c2y = (overlapSpringEndNodesPoints.at(i0).y + overlapSpringEndNodes.at(i0)->position.y) / 2;
						double c2z = (overlapSpringEndNodesPoints.at(i0).z + overlapSpringEndNodes.at(i0)->position.z) / 2;
						// cut the new fibre into several spring segments
						LatticeSpring *cutLs = NULL; // the spring to be cut
						for (int i1 = 0; i1 < fNew->springs.size(); i1++)
						{
							// c1 is on fNew
							double fs1x = fNew->springs.at(i1)->nodes().at(0)->position.x;
							double fs1y = fNew->springs.at(i1)->nodes().at(0)->position.y;
							double fs1z = fNew->springs.at(i1)->nodes().at(0)->position.z;
							double fs2x = fNew->springs.at(i1)->nodes().at(1)->position.x;
							double fs2y = fNew->springs.at(i1)->nodes().at(1)->position.y;
							double fs2z = fNew->springs.at(i1)->nodes().at(1)->position.z;
							double c1x = overlapSpringEndNodesPoints.at(i0).x;
							double c1y = overlapSpringEndNodesPoints.at(i0).y;
							double c1z = overlapSpringEndNodesPoints.at(i0).z;
							double fractI = CSModelTools::fractionPointLine(c1x, c1y, c1z,
								fs1x, fs1y, fs1z,
								fs2x, fs2y, fs2z); // the fraction of the crossed point on spring_i1 of fNew
							if (fractI >= 0 && fractI <= 1)
							{
								cutLs = fNew->springs.at(i1);
								break;
							}
						}
						if (cutLs == NULL)
						{
							std::cout << "	Error: the spring to be crossed is not found!" << std::endl;
						}
						else
						{
							// renew the coordinate
							overlapSpringEndNodes.at(i0)->position.x = c2x;
							overlapSpringEndNodes.at(i0)->position.y = c2y;
							overlapSpringEndNodes.at(i0)->position.z = c2z;
							// new cut spring cutLs
							LinearSpring *cutLsl = static_cast<LinearSpring*>(cutLs);
							ModelElementEdge *cutEI = cutLsl->getBoundEdge();
							ModelElementLatticeNode *cutN1I = cutLsl->nodes().at(0);
							ModelElementLatticeNode *cutN2I = cutLsl->nodes().at(1);
							// remove the neighboring node relationship
							cutN1I->removeNeighbor(cutN2I);
							cutN2I->removeNeighbor(cutN1I);
							// cutN1I still keeps the spring
							cutN2I->removeSpring(cutLs);
							cutLsl->changeN2(overlapSpringEndNodes.at(i0)); // change the end node of the newly added spring
							cutEI->setPoint(&(overlapSpringEndNodes.at(i0)->position), 1); // change the limits of the edge boundingbox
							cutEI->e2Index = overlapSpringEndNodes.at(i0)->mGlobalIndex;
							fNew->addNode(overlapSpringEndNodes.at(i0)); // add the crossed node in the fibre fNew
							overlapSpringEndNodes.at(i0)->addFibre(fNew);
							// add one new spring for the fibre fNew
							LinearSpring *lsNC = new LinearSpring(overlapSpringEndNodes.at(i0), cutN2I);
							lsNC->SetColor(0., 1., 0., 1.); // green
							lsNC->index = springs.size();
							lsNC->n1_unbound = false;
							lsNC->n2_unbound = false; // not free nodes
							springs.push_back(lsNC);
							// add neighboring spring relationship
							overlapSpringEndNodes.at(i0)->addSpring(cutLs);
							overlapSpringEndNodes.at(i0)->addSpring(lsNC);
							cutN2I->addSpring(lsNC); // add the new springs into cutN2I and crossed node
							// add neighboring node relationship
							cutN1I->addNeighbor(overlapSpringEndNodes.at(i0));
							overlapSpringEndNodes.at(i0)->addNeighbor(cutN1I);
							cutN2I->addNeighbor(overlapSpringEndNodes.at(i0));
							overlapSpringEndNodes.at(i0)->addNeighbor(cutN2I);
							// add spring object into the arena
							mpArena->addObject(lsNC->GLObject());
							fNew->addSpring(lsNC); // add spring to the fibre
							lsNC->fibre = fNew;
							// add one boundingbox edge to the bounding box
							ModelElementEdge *eNC = new ModelElementEdge(overlapSpringEndNodes.at(i0)->position.x, 
								overlapSpringEndNodes.at(i0)->position.y, overlapSpringEndNodes.at(i0)->position.z);
							eNC->setPoint(&(overlapSpringEndNodes.at(i0)->position), 0);
							eNC->setPoint(&cutN2I->position, 1);
							eNC->e1Index = overlapSpringEndNodes.at(i0)->mGlobalIndex;
							eNC->e2Index = cutN2I->mGlobalIndex;
							eNC->setRadius(overlapSpringEndNodes.at(i0)->mRadius);
							eNC->setBoundingBox(eNC->mRadius);
							eNC->setBoundSpring(lsNC);
							lsNC->changeBoundEdge(eNC);
							mpBoundingBoxList->add(eNC); // add boundingbox edge to the boundingbox
							edges.push_back(eNC);
						}
					}
					for (int i0 = 0; i0 < overlapEndNodeSpring.size(); i0++)
					{
						if (!overlapEndNodeSpringUse.at(i0)) // not use the crossed spring which is not considered
							continue;
						LinearSpring *lspringI = overlapEndNodeSpring.at(i0);
						LatticeSpring *springI = static_cast<LatticeSpring*>(lspringI);
						ModelElementEdge *edgeI = springI->getBoundEdge();
						ModelElementFibre *fibreI = springI->fibre;
						ModelElementLatticeNode *n1I = lspringI->nodes().at(0);
						ModelElementLatticeNode *n2I = lspringI->nodes().at(1);
						// add a new lattice node, the coordinate for the crosslink
						int cType = overlapEndNodes.at(i0);
						ModelElementLatticeNode *nodeC = nodes.at(n1);
						if (cType == 1)
							nodeC = nodes.at(n2);
						// remove the neighboring node relationship
						n1I->removeNeighbor(n2I);
						n2I->removeNeighbor(n1I);
						// n1 still keeps the spring
						n2I->removeSpring(springI);
						lspringI->changeN2(nodeC); // change the end node of the crossed spring
						edgeI->setPoint(&nodeC->position, 1); // change the limits of the edge boundingbox
						edgeI->e2Index = nodeC->mGlobalIndex;
						fibreI->addNode(nodeC); // add the crossed node into the fibre
						nodeC->addFibre(fibreI);
						// add one new spring for the fibre
						LinearSpring *lsNI = new LinearSpring(nodeC, n2I);
						lsNI->SetColor(0., 1., 0., 1.); // green
						lsNI->index = springs.size();
						lsNI->n1_unbound = false;
						lsNI->n2_unbound = false; // not free nodes
						springs.push_back(lsNI);
						// add neighboring spring relationship
						nodeC->addSpring(springI);
						nodeC->addSpring(lsNI);
						n2I->addSpring(lsNI); // add the new springs into n2I and nodeC
						// add neighboring node relationship
						n1I->addNeighbor(nodeC);
						nodeC->addNeighbor(n1I);
						n2I->addNeighbor(nodeC);
						nodeC->addNeighbor(n2I);
						// add spring object into the arena
						mpArena->addObject(lsNI->GLObject());
						fibreI->addSpring(lsNI); // add spring to the fibre
						lsNI->fibre = fibreI;
						// add one boundingbox edge to the bounding box
						ModelElementEdge *eNI = new ModelElementEdge(nodeC->position.x, nodeC->position.y, nodeC->position.z);
						eNI->setPoint(&nodeC->position, 0);
						eNI->setPoint(&n2I->position, 1);
						eNI->e1Index = nodeC->mGlobalIndex;
						eNI->e2Index = n2I->mGlobalIndex;
						eNI->setRadius((nodeC->mRadius + n2I->mRadius) / 2);
						eNI->setBoundingBox(eNI->mRadius);
						eNI->setBoundSpring(lsNI);
						lsNI->changeBoundEdge(eNI);
						mpBoundingBoxList->add(eNI); // add boundingbox edge to the boundingbox
						edges.push_back(eNI);
					}
					if (subDivide > 1) // one fibre with more than 1 spring segment
					{
						int diffSub = subDivide - (int)overlapSprings.size();
						if (diffSub > 0)
						{
							double basicLength = (nodes.at(n1)->position - nodes.at(n2)->position).Norm() / subDivide;
							for (int i2 = 0; i2 < fNew->springs.size(); i2++)
							{
								double LengthI2 = (fNew->springs.at(i2)->nodes().at(0)->position -
									fNew->springs.at(i2)->nodes().at(1)->position).Norm();
								if (LengthI2 >= basicLength * 2)
								{
									LatticeSpring *flsI = fNew->springs.at(i2);
									LinearSpring *flslI = static_cast<LinearSpring*>(flsI);
									ModelElementEdge *feI = flslI->getBoundEdge();
									ModelElementLatticeNode *flN1I = flslI->nodes().at(0);
									ModelElementLatticeNode *flN2I = flslI->nodes().at(1);
									double fcx1 = (flN1I->position.x + flN2I->position.x) / 2;
									double fcy1 = (flN1I->position.y + flN2I->position.y) / 2;
									double fcz1 = (flN1I->position.z + flN2I->position.z) / 2;
									double fcr1 = (flN1I->mRadius + flN2I->mRadius) / 2;
									// add a new node
									ModelElementLatticeNode *nodeFI = new ModelElementLatticeNode(fcx1, fcy1, fcz1, fcr1);
									nodeFI->mYoung = (flN1I->mYoung + flN2I->mYoung) / 2;
									nodeFI->mPoisson = (flN1I->mPoisson + flN2I->mPoisson) / 2;
									nodeFI->fibreEnd = false; // node inside the fibre
									addNode(nodeFI);
									nodeFI->latticeIndex = nodes.size();
									// remove the neighboring node relationship
									flN1I->removeNeighbor(flN2I);
									flN2I->removeNeighbor(flN1I);
									// flN1I still keeps the spring
									flN2I->removeSpring(flsI);
									flslI->changeN2(nodeFI); // change the end node of the spring
									feI->setPoint(&nodeFI->position, 1); // change the limits of the edge boundingbox
									feI->e2Index = nodeFI->mGlobalIndex;
									fNew->addNode(nodeFI); // add the new node into the fibre
									nodeFI->addFibre(fNew);
									// add one new spring for the fibre
									LinearSpring *lsNF = new LinearSpring(nodeFI, flN2I);
									lsNF->SetColor(0., 1., 0., 1.); // green
									lsNF->index = springs.size();
									lsNF->n1_unbound = false;
									lsNF->n2_unbound = false; // not free nodes
									springs.push_back(lsNF);
									// add neighboring spring relationship
									nodeFI->addSpring(flsI);
									nodeFI->addSpring(lsNF);
									flN2I->addSpring(lsNF); // add the new springs into flN2I and nodeFI
									// add neighboring node relationship
									flN1I->addNeighbor(nodeFI);
									nodeFI->addNeighbor(flN1I);
									flN2I->addNeighbor(nodeFI);
									nodeFI->addNeighbor(flN2I);
									// add spring object into the arena
									mpArena->addObject(lsNF->GLObject());
									fNew->addSpring(lsNF); // add spring to the fibre
									lsNF->fibre = fNew;
									// add one boundingbox edge to the bounding box
									ModelElementEdge *eNF = new ModelElementEdge(nodeFI->position.x, nodeFI->position.y, nodeFI->position.z);
									eNF->setPoint(&nodeFI->position, 0);
									eNF->setPoint(&flN2I->position, 1);
									eNF->e1Index = nodeFI->mGlobalIndex;
									eNF->e2Index = flN2I->mGlobalIndex;
									eNF->setRadius((nodeFI->mRadius + flN2I->mRadius) / 2);
									eNF->setBoundingBox(eNF->mRadius);
									eNF->setBoundSpring(lsNF);
									lsNF->changeBoundEdge(eNF);
									mpBoundingBoxList->add(eNF); // add boundingbox edge to the boundingbox
									edges.push_back(eNF);
								}
							}
						}
					}
					// check if spring merge is needed
					double fNewL = (fNew->n1->position - fNew->n2->position).Norm();
					for (int i1 = 0; i1 < fNew->springs.size(); i1++)
					{
						LatticeSpring *fNewsi = fNew->springs.at(i1);
						ModelElementLatticeNode *fNewsin1 = fNewsi->nodes().at(0);
						ModelElementLatticeNode *fNewsin2 = fNewsi->nodes().at(1);
						double fNewsiL = (fNewsin1->position - fNewsin2->position).Norm();
						if (fNewsiL < fNewL / 5)
						{
							// temporary test: if the spring < 1/5 of the fibre's length
							mergeTwoNodes(fNewsin1, fNewsin2); // merge the two nodes in the spring
						}
					}
					// update the bounding box after the spring division
					mpBoundingBoxList->update();
				}
				else if (overlapSprings.size() == 0 &&
						 overlapSpringEndNodes.size() == 0 &&
						 overlapEndNodeSpring.size() == 0) // no overlap spring found
				{
					if (subDivide == 1) // only one spring segment for the fibre
					{
						// no overlapped spring is found
						// now add one new fibre
						LinearSpring *lsN = new LinearSpring(nodes.at(n1), nodes.at(n2));
						lsN->SetColor(0., 1., 0., 1.); // green
						lsN->index = springs.size();
						lsN->n1_unbound = false;
						lsN->n2_unbound = false;
						nodes.at(n1)->springs.push_back(lsN);
						nodes.at(n2)->springs.push_back(lsN);
						nodes.at(n1)->nodes.push_back(nodes.at(n2));
						nodes.at(n2)->nodes.push_back(nodes.at(n1));
						springs.push_back(lsN);
						// add spring object into the arena
						mpArena->addObject(lsN->GLObject());
						// add spring into the fibre
						ModelElementFibre *f1 = new ModelElementFibre();
						f1->n1 = nodes.at(n1);
						f1->n2 = nodes.at(n2); // add n1 and n2 into the fibre
						f1->n1Bound = true;
						f1->n2Bound = true;
						f1->springs.push_back(lsN); // add spring into the fibre
						nodes.at(n1)->fibres.push_back(f1);
						nodes.at(n2)->fibres.push_back(f1); // add fibre into n1 and n2
						fibres.push_back(f1); // add fibre into the fibres
						f1->fibre_id = fibres.size();
						// add one boundingbox edge to the bounding box
						ModelElementEdge *e1 = new ModelElementEdge(nodes.at(n1)->position.x,
							nodes.at(n1)->position.y,
							nodes.at(n1)->position.z);
						e1->setPoint(&nodes.at(n1)->position, 0);
						e1->setPoint(&nodes.at(n2)->position, 1);
						e1->e1Index = nodes.at(n1)->mGlobalIndex;
						e1->e2Index = nodes.at(n2)->mGlobalIndex;
						e1->setRadius((nodes.at(n1)->mRadius + nodes.at(n2)->mRadius) / 2);
						e1->setBoundingBox(e1->mRadius);
						e1->setBoundSpring(lsN);
						mpBoundingBoxList->add(e1); // add boundingbox edge to the boundingbox
						lsN->changeBoundEdge(e1);
						lsN->setFibre(f1); // set the fibre for the edge
						edges.push_back(e1); // add boundingbox edge into the edges
					}
					else // one fibre with more than 1 spring segment
					{
						// add spring into the fibre
						ModelElementFibre *f1 = new ModelElementFibre();
						f1->n1 = nodes.at(n1);
						f1->n2 = nodes.at(n2); // add n1 and n2 into the fibre
						f1->n1Bound = true;
						f1->n2Bound = true;
						nodes.at(n1)->fibres.push_back(f1);
						nodes.at(n2)->fibres.push_back(f1); // add fibre into n1 and n2
						// add each node and spring into the fibre
						for (int j = 1; j < subDivide; j++)
						{
							// subdivide the linear spring into segments
							// first add subDivide - 1 new nodes here
							double xj = (nodes.at(n2)->position.x * j + nodes.at(n1)->position.x * (subDivide - j)) / subDivide;
							double yj = (nodes.at(n2)->position.y * j + nodes.at(n1)->position.y * (subDivide - j)) / subDivide;
							double zj = (nodes.at(n2)->position.z * j + nodes.at(n1)->position.z * (subDivide - j)) / subDivide;
							double rj = (nodes.at(n2)->mRadius * j + nodes.at(n1)->mRadius * (subDivide - j)) / subDivide;
							ModelElementLatticeNode *nodeJ = new ModelElementLatticeNode(xj, yj, zj, rj);
							nodeJ->mYoung = nodes.at(n1)->mYoung;
							nodeJ->mPoisson = nodes.at(n2)->mPoisson;
							nodeJ->fibreEnd = false; // node inside the fibre
							if (nodes.at(n1)->top && nodes.at(n2)->top)
								nodeJ->top = true;
							if (nodes.at(n1)->bottom && nodes.at(n2)->bottom)
								nodeJ->bottom = true;
							addNode(nodeJ);
							nodeJ->latticeIndex = nodes.size(); // add node to the fibre
							f1->nodes.push_back(nodeJ);
							nodeJ->fibres.push_back(f1); // add fibre into the new node
							//mpArena->addObject(nodeJ->GLObject()); // just for test
							if (j == 1)
							{
								LinearSpring *lsNJ = new LinearSpring(nodes.at(n1), nodeJ);
								lsNJ->SetColor(0., 1., 0., 1.); // green
								lsNJ->index = springs.size();
								lsNJ->n1_unbound = false;
								lsNJ->n2_unbound = false;
								// add spring to nodes
								nodes.at(n1)->springs.push_back(lsNJ);
								nodeJ->springs.push_back(lsNJ);
								// add nodes to each other
								nodes.at(n1)->nodes.push_back(nodeJ);
								nodeJ->nodes.push_back(nodes.at(n1));
								// add spring to network
								springs.push_back(lsNJ);
								// add spring object into the arena
								mpArena->addObject(lsNJ->GLObject());
								f1->springs.push_back(lsNJ); // add spring to the fibre
								lsNJ->fibre = f1;
								// add one boundingbox edge to the bounding box
								ModelElementEdge *e1j = new ModelElementEdge(nodes.at(n1)->position.x,
									nodes.at(n1)->position.y,
									nodes.at(n1)->position.z);
								e1j->setPoint(&nodes.at(n1)->position, 0);
								e1j->setPoint(&nodeJ->position, 1);
								e1j->e1Index = nodes.at(n1)->mGlobalIndex;
								e1j->e2Index = nodeJ->mGlobalIndex;
								e1j->setRadius((nodes.at(n1)->mRadius + nodeJ->mRadius) / 2);
								e1j->setBoundingBox(e1j->mRadius);
								e1j->setBoundSpring(lsNJ);
								mpBoundingBoxList->add(e1j); // add boundingbox edge to the boundingbox
								lsNJ->changeBoundEdge(e1j);
								edges.push_back(e1j); // add boundingbox edge into the edges
							}
							else
							{
								ModelElementLatticeNode *nodeJ1 = nodes.at(nodes.size() - 2);
								LinearSpring *lsNJ = new LinearSpring(nodeJ1, nodeJ);
								lsNJ->SetColor(0., 1., 0., 1.); // green
								lsNJ->index = springs.size();
								lsNJ->n1_unbound = false;
								lsNJ->n2_unbound = false;
								nodeJ1->springs.push_back(lsNJ);
								nodeJ->springs.push_back(lsNJ);
								nodeJ1->nodes.push_back(nodeJ);
								nodeJ->nodes.push_back(nodeJ1);
								springs.push_back(lsNJ);
								// add spring object into the arena
								mpArena->addObject(lsNJ->GLObject());
								f1->springs.push_back(lsNJ); // add spring to the fibre
								lsNJ->fibre = f1;
								// add one boundingbox edge to the bounding box
								ModelElementEdge *e1j = new ModelElementEdge(nodeJ1->position.x,
									nodeJ1->position.y,
									nodeJ1->position.z);
								e1j->setPoint(&nodeJ1->position, 0);
								e1j->setPoint(&nodeJ->position, 1);
								e1j->e1Index = nodeJ1->mGlobalIndex;
								e1j->e2Index = nodeJ->mGlobalIndex;
								e1j->setRadius((nodeJ1->mRadius + nodeJ->mRadius) / 2);
								e1j->setBoundingBox(e1j->mRadius);
								e1j->setBoundSpring(lsNJ);
								mpBoundingBoxList->add(e1j); // add boundingbox edge to the boundingbox
								lsNJ->changeBoundEdge(e1j);
								edges.push_back(e1j); // add boundingbox edge into the edges
							}
						}
						// add the last segment
						ModelElementLatticeNode *nodeJ1 = nodes.at(nodes.size() - 1);
						LinearSpring *lsNJ = new LinearSpring(nodeJ1, nodes.at(n2));
						lsNJ->SetColor(0., 1., 0., 1.); // green
						lsNJ->index = springs.size();
						lsNJ->n1_unbound = false;
						lsNJ->n2_unbound = false;
						nodeJ1->springs.push_back(lsNJ);
						nodes.at(n2)->springs.push_back(lsNJ);
						nodeJ1->nodes.push_back(nodes.at(n2));
						nodes.at(n2)->nodes.push_back(nodeJ1);
						springs.push_back(lsNJ);
						// add spring object into the arena
						mpArena->addObject(lsNJ->GLObject());
						f1->springs.push_back(lsNJ); // add spring to the fibre
						lsNJ->fibre = f1;
						fibres.push_back(f1); // add fibre into the fibres
						f1->fibre_id = fibres.size();
						// add one boundingbox edge to the bounding box
						ModelElementEdge *e1 = new ModelElementEdge(nodeJ1->position.x,
							nodeJ1->position.y,
							nodeJ1->position.z);
						e1->setPoint(&nodeJ1->position, 0);
						e1->setPoint(&nodes.at(n2)->position, 1);
						e1->e1Index = nodeJ1->mGlobalIndex;
						e1->e2Index = nodes.at(n2)->mGlobalIndex;
						e1->setRadius((nodeJ1->mRadius + nodes.at(n2)->mRadius) / 2);
						e1->setBoundingBox(e1->mRadius);
						e1->setBoundSpring(lsNJ);
						mpBoundingBoxList->add(e1); // add boundingbox edge to the boundingbox
						lsNJ->changeBoundEdge(e1);
						edges.push_back(e1); // add boundingbox edge into the edges
					}
					mpBoundingBoxList->update();
				}
				nodes.at(n1)->fibreEnd = true;
				nodes.at(n2)->fibreEnd = true; // the fibre end-node
			}
		}
	}
	// add RotationalSpring
	// it is along the fibre
	for (int i = 0; i < fibres.size(); i++)
	{
		ModelElementFibre * fI = fibres.at(i);
		for (int j = 0; j < fI->nodes.size(); j++)
		{
			ModelElementLatticeNode * nJ = fI->nodes.at(j);
			ModelElementLatticeNode * n1 = NULL;
			ModelElementLatticeNode * n2 = NULL;
			bool n1take = false;
			bool n2take = false;
			for (int k = 0; k < nJ->nodes.size(); k++)
			{
				ModelElementLatticeNode * nK = nJ->nodes.at(k);
				if (nK->checkFibre(fI))
				{
					if (!n1take && !n2take)
					{
						n1 = nK;
						n1take = true;
					}
					else if (n1take && !n2take)
					{
						n2 = nK;
						n2take = true;
					}
					else if (n1take && n2take)
					{
						break;
					}
				}
			}
			if (n1 == NULL || n2 == NULL)
			{
				std::cout << "	-> The neighboring nodes for node " << nJ->mGlobalIndex << " are missed!" << std::endl;
				break;
			}
			// add one rotational spring: n1 - nJ - n2
			RotationalSpring *rJ = new RotationalSpring(n1, n2, nJ);
			nJ->Rsprings.push_back(rJ);
			springs.push_back(rJ);
		}
	}
	// delete delaunay triangulation
	Del->~CD3DW();
	// delete delaunay simplices
	Complex->~CComplex();
}

void Lattice::createNetworkSimple(double cutoff, int subDivide)
{
	// has to be called from the Model*.cpp
	// simply add fibres regardless if there is spring-spring overlap
	// then call the function growTest to do the spring-spring collision detect
	if (subDivide < 1) // the fibre segment is at least 1 
		subDivide = 1;
	CD3DW *Del = new CD3DW();
	int nsize = (int)nodes.size();
	for (int i = 0; i < nsize; i++)
	{
		// add Lattice Nodes into the system one-by-one
		Del->add(nodes.at(i)->position.x,
			nodes.at(i)->position.y,
			nodes.at(i)->position.z,
			nodes.at(i)->mRadius,
			i,
			&nodes.at(i));
	}

	CComplex *Complex = new CComplex(Del);
	// add LinearSpring
	int Pn = (int)Complex->P.size();
	for (int i = 0; i < Pn; i++)
	{
		CSimplex *p = Complex->P[i];
		int size = p->size();
		if (size == 2) // edges
		{
			C1Simplex *p1 = (C1Simplex*)p;
			int n1 = p1->A->n();
			int n2 = p1->B->n();
			Vector3f v12 = nodes.at(n1)->position - nodes.at(n2)->position;
			double dist12 = v12.Norm();
			if (dist12 > cutoff) continue;

			if (subDivide == 1) // only one spring segment for the fibre
			{
				// no overlapped spring is found
				// now add one new fibre
				LinearSpring *lsN = new LinearSpring(nodes.at(n1), nodes.at(n2));
				lsN->SetColor(0., 1., 0., 1.); // green
				lsN->index = springs.size();
				lsN->n1_unbound = false;
				lsN->n2_unbound = false;
				nodes.at(n1)->springs.push_back(lsN);
				nodes.at(n2)->springs.push_back(lsN);
				nodes.at(n1)->nodes.push_back(nodes.at(n2));
				nodes.at(n2)->nodes.push_back(nodes.at(n1));
				springs.push_back(lsN);
				// add spring object into the arena
				mpArena->addObject(lsN->GLObject());
				// add spring into the fibre
				ModelElementFibre *f1 = new ModelElementFibre();
				f1->n1 = nodes.at(n1);
				f1->n2 = nodes.at(n2); // add n1 and n2 into the fibre
				f1->n1Bound = true;
				f1->n2Bound = true;
				f1->springs.push_back(lsN); // add spring into the fibre
				nodes.at(n1)->fibres.push_back(f1);
				nodes.at(n2)->fibres.push_back(f1); // add fibre into n1 and n2
				fibres.push_back(f1); // add fibre into the fibres
				f1->fibre_id = fibres.size();
				// add one boundingbox edge to the bounding box
				ModelElementEdge *e1 = new ModelElementEdge(nodes.at(n1)->position.x,
					nodes.at(n1)->position.y,
					nodes.at(n1)->position.z);
				e1->setPoint(&nodes.at(n1)->position, 0);
				e1->setPoint(&nodes.at(n2)->position, 1);
				e1->e1Index = nodes.at(n1)->mGlobalIndex;
				e1->e2Index = nodes.at(n2)->mGlobalIndex;
				e1->setRadius((nodes.at(n1)->mRadius + nodes.at(n2)->mRadius) / 2);
				e1->setBoundingBox(e1->mRadius);
				e1->setBoundSpring(lsN);
				mpBoundingBoxList->add(e1); // add boundingbox edge to the boundingbox
				lsN->changeBoundEdge(e1);
				lsN->setFibre(f1); // set the fibre for the edge
				edges.push_back(e1); // add boundingbox edge into the edges
			}
			else // one fibre with more than 1 spring segment
			{
				// add spring into the fibre
				ModelElementFibre *f1 = new ModelElementFibre();
				f1->n1 = nodes.at(n1);
				f1->n2 = nodes.at(n2); // add n1 and n2 into the fibre
				f1->n1Bound = true;
				f1->n2Bound = true;
				nodes.at(n1)->fibres.push_back(f1);
				nodes.at(n2)->fibres.push_back(f1); // add fibre into n1 and n2
													// add each node and spring into the fibre
				for (int j = 1; j < subDivide; j++)
				{
					// subdivide the linear spring into segments
					// first add subDivide - 1 new nodes here
					double xj = (nodes.at(n2)->position.x * j + nodes.at(n1)->position.x * (subDivide - j)) / subDivide;
					double yj = (nodes.at(n2)->position.y * j + nodes.at(n1)->position.y * (subDivide - j)) / subDivide;
					double zj = (nodes.at(n2)->position.z * j + nodes.at(n1)->position.z * (subDivide - j)) / subDivide;
					double rj = (nodes.at(n2)->mRadius * j + nodes.at(n1)->mRadius * (subDivide - j)) / subDivide;
					ModelElementLatticeNode *nodeJ = new ModelElementLatticeNode(xj, yj, zj, rj);
					nodeJ->mYoung = nodes.at(n1)->mYoung;
					nodeJ->mPoisson = nodes.at(n2)->mPoisson;
					nodeJ->fibreEnd = false; // node inside the fibre
					if (nodes.at(n1)->top && nodes.at(n2)->top)
						nodeJ->top = true;
					if (nodes.at(n1)->bottom && nodes.at(n2)->bottom)
						nodeJ->bottom = true;
					addNode(nodeJ);
					nodeJ->latticeIndex = nodes.size(); // add node to the fibre
					f1->nodes.push_back(nodeJ);
					nodeJ->fibres.push_back(f1); // add fibre into the new node
												 //mpArena->addObject(nodeJ->GLObject()); // just for test
					if (j == 1)
					{
						LinearSpring *lsNJ = new LinearSpring(nodes.at(n1), nodeJ);
						lsNJ->SetColor(0., 1., 0., 1.); // green
						lsNJ->index = springs.size();
						lsNJ->n1_unbound = false;
						lsNJ->n2_unbound = false;
						// add spring to nodes
						nodes.at(n1)->springs.push_back(lsNJ);
						nodeJ->springs.push_back(lsNJ);
						// add nodes to each other
						nodes.at(n1)->nodes.push_back(nodeJ);
						nodeJ->nodes.push_back(nodes.at(n1));
						// add spring to network
						springs.push_back(lsNJ);
						// add spring object into the arena
						mpArena->addObject(lsNJ->GLObject());
						f1->springs.push_back(lsNJ); // add spring to the fibre
						lsNJ->fibre = f1;
						// add one boundingbox edge to the bounding box
						ModelElementEdge *e1j = new ModelElementEdge(nodes.at(n1)->position.x,
							nodes.at(n1)->position.y,
							nodes.at(n1)->position.z);
						e1j->setPoint(&nodes.at(n1)->position, 0);
						e1j->setPoint(&nodeJ->position, 1);
						e1j->e1Index = nodes.at(n1)->mGlobalIndex;
						e1j->e2Index = nodeJ->mGlobalIndex;
						e1j->setRadius((nodes.at(n1)->mRadius + nodeJ->mRadius) / 2);
						e1j->setBoundingBox(e1j->mRadius);
						e1j->setBoundSpring(lsNJ);
						mpBoundingBoxList->add(e1j); // add boundingbox edge to the boundingbox
						lsNJ->changeBoundEdge(e1j);
						edges.push_back(e1j); // add boundingbox edge into the edges
					}
					else
					{
						ModelElementLatticeNode *nodeJ1 = nodes.at(nodes.size() - 2);
						LinearSpring *lsNJ = new LinearSpring(nodeJ1, nodeJ);
						lsNJ->SetColor(0., 1., 0., 1.); // green
						lsNJ->index = springs.size();
						lsNJ->n1_unbound = false;
						lsNJ->n2_unbound = false;
						nodeJ1->springs.push_back(lsNJ);
						nodeJ->springs.push_back(lsNJ);
						nodeJ1->nodes.push_back(nodeJ);
						nodeJ->nodes.push_back(nodeJ1);
						springs.push_back(lsNJ);
						// add spring object into the arena
						mpArena->addObject(lsNJ->GLObject());
						f1->springs.push_back(lsNJ); // add spring to the fibre
						lsNJ->fibre = f1;
						// add one boundingbox edge to the bounding box
						ModelElementEdge *e1j = new ModelElementEdge(nodeJ1->position.x,
							nodeJ1->position.y,
							nodeJ1->position.z);
						e1j->setPoint(&nodeJ1->position, 0);
						e1j->setPoint(&nodeJ->position, 1);
						e1j->e1Index = nodeJ1->mGlobalIndex;
						e1j->e2Index = nodeJ->mGlobalIndex;
						e1j->setRadius((nodeJ1->mRadius + nodeJ->mRadius) / 2);
						e1j->setBoundingBox(e1j->mRadius);
						e1j->setBoundSpring(lsNJ);
						mpBoundingBoxList->add(e1j); // add boundingbox edge to the boundingbox
						lsNJ->changeBoundEdge(e1j);
						edges.push_back(e1j); // add boundingbox edge into the edges
					}
				}
				// add the last segment
				ModelElementLatticeNode *nodeJ1 = nodes.at(nodes.size() - 1);
				LinearSpring *lsNJ = new LinearSpring(nodeJ1, nodes.at(n2));
				lsNJ->SetColor(0., 1., 0., 1.); // green
				lsNJ->index = springs.size();
				lsNJ->n1_unbound = false;
				lsNJ->n2_unbound = false;
				nodeJ1->springs.push_back(lsNJ);
				nodes.at(n2)->springs.push_back(lsNJ);
				nodeJ1->nodes.push_back(nodes.at(n2));
				nodes.at(n2)->nodes.push_back(nodeJ1);
				springs.push_back(lsNJ);
				// add spring object into the arena
				mpArena->addObject(lsNJ->GLObject());
				f1->springs.push_back(lsNJ); // add spring to the fibre
				lsNJ->fibre = f1;
				fibres.push_back(f1); // add fibre into the fibres
				f1->fibre_id = fibres.size();
				// add one boundingbox edge to the bounding box
				ModelElementEdge *e1 = new ModelElementEdge(nodeJ1->position.x,
					nodeJ1->position.y,
					nodeJ1->position.z);
				e1->setPoint(&nodeJ1->position, 0);
				e1->setPoint(&nodes.at(n2)->position, 1);
				e1->e1Index = nodeJ1->mGlobalIndex;
				e1->e2Index = nodes.at(n2)->mGlobalIndex;
				e1->setRadius((nodeJ1->mRadius + nodes.at(n2)->mRadius) / 2);
				e1->setBoundingBox(e1->mRadius);
				e1->setBoundSpring(lsNJ);
				mpBoundingBoxList->add(e1); // add boundingbox edge to the boundingbox
				lsNJ->changeBoundEdge(e1);
				edges.push_back(e1); // add boundingbox edge into the edges
			}
			nodes.at(n1)->fibreEnd = true;
			nodes.at(n2)->fibreEnd = true; // the fibre end-node
		}
	}
	// add RotationalSpring
	// it is along the fibre
	for (int i = 0; i < fibres.size(); i++)
	{
		ModelElementFibre * fI = fibres.at(i);
		for (int j = 0; j < fI->nodes.size(); j++)
		{
			ModelElementLatticeNode * nJ = fI->nodes.at(j);
			ModelElementLatticeNode * n1 = NULL;
			ModelElementLatticeNode * n2 = NULL;
			bool n1take = false;
			bool n2take = false;
			for (int k = 0; k < nJ->nodes.size(); k++)
			{
				ModelElementLatticeNode * nK = nJ->nodes.at(k);
				if (nK->checkFibre(fI))
				{
					if (!n1take && !n2take)
					{
						n1 = nK;
						n1take = true;
					}
					else if (n1take && !n2take)
					{
						n2 = nK;
						n2take = true;
					}
					else if (n1take && n2take)
					{
						break;
					}
				}
			}
			if (n1 == NULL || n2 == NULL)
			{
				std::cout << "	-> The neighboring nodes for node " << nJ->mGlobalIndex << " are missed!" << std::endl;
				break;
			}
			// add one rotational spring: n1 - nJ - n2
			RotationalSpring *rJ = new RotationalSpring(n1,n2,nJ);
			nJ->Rsprings.push_back(rJ);
			springs.push_back(rJ);
		}
	}
	// delete delaunay triangulation
	Del->~CD3DW();
	// delete delaunay simplices
	Complex->~CComplex();
}

void Lattice::updateFibreBound()
{
	std::string outputFileName1 = "debugFile.txt";
	std::ofstream FO;
	FO.open(outputFileName1.c_str(), std::ios::out | std::ios::app);
	FO << "ModelLattice.cpp / updateFibreBound()	starts" << endl;
	// first check the strain for all linear springs
	for (int i = 0; i < springs.size(); i++)
	{
		if (springs.at(i)->mSpringType == LatticeSpring::TypeLinearSpring)
		{
			LinearSpring *Si = static_cast<LinearSpring*>(springs.at(i));
			Si->unboundCheck();
		}
	}
	// for all unbound fibres, check the next binding node
	for (int i = 0; i < fibres.size(); i++)
	{
		ModelElementFibre * Fi = fibres.at(i);
		FO << "	-> fibre " << i << " out of " << fibres.size() << std::endl;
		FO << "	-> Fi->n1 " << Fi->n1->mGlobalIndex << std::endl;
		FO << "	-> Fi->n2 " << Fi->n2->mGlobalIndex << std::endl;
		LinearSpring *e1 = static_cast<LinearSpring*>(Fi->getE1());
		LinearSpring *e2 = static_cast<LinearSpring*>(Fi->getE2());
		FO << "		-> e1->0 " << e1->nodes().at(0)->mGlobalIndex << std::endl;
		FO << "		-> e1->1 " << e1->nodes().at(1)->mGlobalIndex << std::endl;
		FO << "		-> e2->0 " << e2->nodes().at(0)->mGlobalIndex << std::endl;
		FO << "		-> e2->1 " << e2->nodes().at(1)->mGlobalIndex << std::endl;
		// first check e1
		if (e1->nodes().at(0) == Fi->n1)
		{
			if (e1->n1_unbound && !e1->nodes().at(0)->free) // for nodes(0)
			{
				FO << "		-> e1-Fi->n1 " << Fi->n1->mGlobalIndex << " is unbound" << std::endl;
				// fibre n1 is unbound
				int nbIndex = -1; // node id for the next binding
				double minLength = 1000000.;
				double fLength = (Fi->n2->position - Fi->n1->position).Norm();
				for (int j = 0; j < Fi->n1->nodes.size(); j++)
				{
					ModelElementLatticeNode *n1j = Fi->n1->nodes.at(j);
					bool outFibre = true;
					for (int k = 0; k < n1j->fibres.size(); k++)
					{
						if (n1j->fibres.at(k) == Fi)
						{
							outFibre = false;
							break;
						}
					}
					if (outFibre) // the node must be not on the same fibre
					{
						double dLength = (Fi->n2->position - n1j->position).Norm();
						if (dLength < fLength && dLength < minLength)
						{
							minLength = dLength;
							nbIndex = j;
						}
					}
				}
				FO << "		-> nbIndex: " << nbIndex << std::endl;
				//nbIndex = -1; // test to disable sliding
				if (nbIndex >= 0)
				{
					// fibre slides to the next node
					Fi->unboundE(1, Fi->n1->nodes.at(nbIndex)); // slide from n1 to node-nbindex
				}
				else
				{
					// fibre with free end-node n1
					// add one free node for the fibre
					ModelElementLatticeNode *e10 = e1->nodes().at(0);
					ModelElementLatticeNode *e11 = e1->nodes().at(1);
					e10->removeNeighbor(e11);
					e11->removeNeighbor(e10); // remove the neighbors
					e10->removeSpring(e1); // remove the spring from the end-node
					e10->removeFibre(Fi); // remove the fibre from the end-node
					// add one new lattice node here
					double e1L = (e10->position - e11->position).Norm();
					double e1IL = e1->initialL0;
					double nnx = e11->position.x + (e10->position.x - e11->position.x) * e1IL / e1L;
					double nny = e11->position.y + (e10->position.y - e11->position.y) * e1IL / e1L;
					double nnz = e11->position.z + (e10->position.z - e11->position.z) * e1IL / e1L;
					double nnr = e11->mRadius;
					ModelElementLatticeNode *nn = new ModelElementLatticeNode(nnx, nny, nnz, nnr);
					nn->mYoung = e11->mYoung;
					nn->mPoisson = e11->mPoisson;
					nn->addNeighbor(e11);
					e11->addNeighbor(nn); // add neighbors to the new node
					nn->addSpring(e1); // add spring to the new node
					nn->addFibre(Fi); // add fibre to the new node
					this->addNode(nn);
					mpArena->addObject(nn->GLObject());
					//e1->nodes().at(0) = nn;
					e1->changeN1(nn);
					e1->n1_unbound = true;
					e1->nodes().at(0)->free = true;
					Fi->n1 = nn;
					Fi->n1Bound = false;
					
				}
			}
		}
		else // e1->nodes(1) == Fi->n1
		{
			if (e1->n2_unbound && !e1->nodes().at(1)->free) // for nodes(1)
			{
				FO << "		-> e1-Fi->n1 " << Fi->n1->mGlobalIndex << " is unbound" << std::endl;
				// fibre n1 is unbound
				int nbIndex = -1;
				double minLength = 1000000.;
				double fLength = (Fi->n2->position - Fi->n1->position).Norm();
				for (int j = 0; j < Fi->n1->nodes.size(); j++)
				{
					ModelElementLatticeNode *n1j = Fi->n1->nodes.at(j);
					bool outFibre = true;
					for (int k = 0; k < n1j->fibres.size(); k++)
					{
						if (n1j->fibres.at(k) == Fi)
						{
							outFibre = false;
							break;
						}
					}
					if (outFibre) // the node must be not on the same fibre
					{
						double dLength = (Fi->n2->position - n1j->position).Norm();
						if (dLength < fLength && dLength < minLength)
						{
							minLength = dLength;
							nbIndex = j;
						}
					}
				}
				FO << "		-> nbIndex: " << nbIndex << std::endl;
				//nbIndex = -1; // test to disable sliding
				if (nbIndex >= 0)
				{
					// fibre slides to the next node
					Fi->unboundE(2, Fi->n1->nodes.at(nbIndex)); // slide from n2 to node-nbindex
				}
				else
				{
					// fibre with free end-node n1
					// add one free node for the fibre
					ModelElementLatticeNode *e10 = e1->nodes().at(0);
					ModelElementLatticeNode *e11 = e1->nodes().at(1);
					e10->removeNeighbor(e11);
					e11->removeNeighbor(e10); // remove the neighbors
					e11->removeSpring(e1); // remove the spring from the end-node
					e11->removeFibre(Fi); // remove the fibre from the end-node
					// add one new lattice node here
					double e1L = (e10->position - e11->position).Norm();
					double e1IL = e1->initialL0;
					double nnx = e10->position.x + (e11->position.x - e10->position.x) * e1IL / e1L;
					double nny = e10->position.y + (e11->position.y - e10->position.y) * e1IL / e1L;
					double nnz = e10->position.z + (e11->position.z - e10->position.z) * e1IL / e1L;
					double nnr = e10->mRadius;
					ModelElementLatticeNode *nn = new ModelElementLatticeNode(nnx, nny, nnz, nnr);
					nn->mYoung = e10->mYoung;
					nn->mPoisson = e10->mPoisson;
					nn->addNeighbor(e10);
					e10->addNeighbor(nn); // add neighbors to the new node
					nn->addSpring(e1); // add spring to the new node
					nn->addFibre(Fi); // add fibre to the new node
					this->addNode(nn); // add the new node to the lattice
					mpArena->addObject(nn->GLObject());
					//e1->nodes().at(1) = nn;
					e1->changeN2(nn);
					e1->nodes().at(1)->free = true;
					Fi->n1 = nn;
					Fi->n1Bound = false;
				}
			}
		}
		// then check e2
		if (e2->nodes().at(0) == Fi->n2)
		{
			if (e2->n1_unbound && !e2->nodes().at(0)->free) // for nodes(0)
			{
				FO << "		-> e2-Fi->n2 " << Fi->n2->mGlobalIndex << " is unbound" << std::endl;
				// fibre n2 is unbound
				int nbIndex = -1; // node id for the next binding
				double minLength = 1000000.;
				double fLength = (Fi->n2->position - Fi->n1->position).Norm();
				for (int j = 0; j < Fi->n2->nodes.size(); j++)
				{
					ModelElementLatticeNode *n2j = Fi->n2->nodes.at(j);
					bool outFibre = true;
					for (int k = 0; k < n2j->fibres.size(); k++)
					{
						if (n2j->fibres.at(k) == Fi)
						{
							outFibre = false;
							break;
						}
					}
					if (outFibre) // the node must be not on the same fibre
					{
						double dLength = (Fi->n1->position - n2j->position).Norm();
						if (dLength < fLength && dLength < minLength)
						{
							minLength = dLength;
							nbIndex = j;
						}
					}
				}
				FO << "		-> nbIndex: " << nbIndex << std::endl;
				//nbIndex = -1; // test to disable sliding
				if (nbIndex >= 0)
				{
					// fibre slides to the next node
					Fi->unboundE(1, Fi->n2->nodes.at(nbIndex)); // slide from n1 to node-nbindex
				}
				else
				{
					// fibre with free end-node n1
					// add one free node for the fibre
					ModelElementLatticeNode *e20 = e2->nodes().at(0);
					ModelElementLatticeNode *e21 = e2->nodes().at(1);
					e20->removeNeighbor(e21);
					e21->removeNeighbor(e20); // remove the neighbors
					e20->removeSpring(e2); // remove the spring from the end-node
					e20->removeFibre(Fi); // remove the fibre from the end-node
					// add one new lattice node here
					double e2L = (e20->position - e21->position).Norm();
					double e2IL = e2->initialL0;
					double nnx = e21->position.x + (e20->position.x - e21->position.x) * e2IL / e2L;
					double nny = e21->position.y + (e20->position.y - e21->position.y) * e2IL / e2L;
					double nnz = e21->position.z + (e20->position.z - e21->position.z) * e2IL / e2L;
					double nnr = e21->mRadius;
					ModelElementLatticeNode *nn = new ModelElementLatticeNode(nnx, nny, nnz, nnr);
					nn->mYoung = e21->mYoung;
					nn->mPoisson = e21->mPoisson;
					nn->addNeighbor(e21);
					e21->addNeighbor(nn); // add neighbors to the new node
					nn->addSpring(e2); // add spring to the new node
					nn->addFibre(Fi); // add fibre to the new node
					this->addNode(nn);
					//e2->nodes().at(0) == nn;
					e2->changeN1(nn);
					e2->nodes().at(0)->free = true;
					Fi->n2 = nn;
					Fi->n2Bound = false;
				}
			}
		}
		else // e2->nodes(1) == Fi->n2
		{
			if (e2->n2_unbound && !e2->nodes().at(1)->free) // for nodes(1)
			{
				FO << "		-> e2-Fi->n2 " << Fi->n2->mGlobalIndex << " is unbound" << std::endl;
				// fibre n2 is unbound
				int nbIndex = -1; // node id for the next binding
				double minLength = 1000000.;
				double fLength = (Fi->n2->position - Fi->n1->position).Norm();
				for (int j = 0; j < Fi->n2->nodes.size(); j++)
				{
					ModelElementLatticeNode *n2j = Fi->n2->nodes.at(j);
					bool outFibre = true;
					for (int k = 0; k < n2j->fibres.size(); k++)
					{
						if (n2j->fibres.at(k) == Fi)
						{
							outFibre = false;
							break;
						}
					}
					if (outFibre) // the node must be not on the same fibre
					{
						double dLength = (Fi->n1->position - n2j->position).Norm();
						//if (Fi->n1->mGlobalIndex == 3 && Fi->n2->mGlobalIndex == 5)
						//{
						//	FO << "		-> step_" << j << ", dLength: " << dLength << ", fLength: " << fLength << std::endl;
						//}
						if (dLength < fLength && dLength < minLength)
						{
							minLength = dLength;
							nbIndex = j;
							//if (Fi->n1->mGlobalIndex == 3 && Fi->n2->mGlobalIndex == 5)
							//{
							//	FO << "		-> step_" << j << ", minLength: " << minLength << ", nbIndex: " << nbIndex << std::endl;
							//}
						}
					}
				}
				FO << "		-> nbIndex: " << nbIndex << std::endl;
				//nbIndex = -1; // test to disable sliding
				if (nbIndex >= 0)
				{
					// fibre slides to the next node
					Fi->unboundE(2, Fi->n2->nodes.at(nbIndex)); // slide from n2 to node-nbindex
				}
				else
				{
					// fibre with free end-node n2
					// ad one free node for the fibre
					ModelElementLatticeNode *e20 = e2->nodes().at(0);
					ModelElementLatticeNode *e21 = e2->nodes().at(1);
					e20->removeNeighbor(e21);
					e21->removeNeighbor(e20); // remove the neighbors
					e21->removeSpring(e2); // remove the spring from the end-node
					e21->removeFibre(Fi); // remove the fibre from the end-node
					// add one new lattice node here
					double e2L = (e20->position - e21->position).Norm();
					double e2IL = e2->initialL0;
					double nnx = e20->position.x + (e21->position.x - e20->position.x) * e2IL / e2L;
					double nny = e20->position.y + (e21->position.y - e20->position.y) * e2IL / e2L;
					double nnz = e20->position.z + (e21->position.z - e20->position.z) * e2IL / e2L;
					double nnr = e20->mRadius;
					ModelElementLatticeNode *nn = new ModelElementLatticeNode(nnx, nny, nnz, nnr);
					nn->mYoung = e20->mYoung;
					nn->mPoisson = e20->mPoisson;
					nn->addNeighbor(e20);
					e20->addNeighbor(nn); // add neighbors to the new node
					nn->addSpring(e2); // add spring to the new node
					nn->addFibre(Fi); // add fibre to the new node
					this->addNode(nn);
					//e2->nodes().at(1) == nn;
					e2->changeN2(nn);
					e2->nodes().at(1)->free = true;
					Fi->n2 = nn;
					Fi->n2Bound = false;
				}
			}
		}
	}
	FO.close();
}

unsigned int Lattice::getNumberOfNodes() const
{
    unsigned int numberOfNodes = 0;
    for (const auto node : nodes)
        if (!node->isAnchor())
            numberOfNodes++;

    return numberOfNodes;
}
