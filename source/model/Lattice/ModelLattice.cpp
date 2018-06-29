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
		std::cout << "	-> Removing a node" << std::endl;
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
		N1->RspringsToDelete.pop_back();
	}
}

bool Lattice::removeFibre(ModelElementFibre * to_remove)
{
	// remove a fibre in all structure: node, spring
	if (to_remove != NULL)
		std::cout << "	-> Removing a fibre from the lattice system" << std::endl;
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
		mSpringsToDelete.pop_back();
	}
	for (int i = 0; i < to_remove->nodes.size(); i++)
	{
		ModelElementLatticeNode * Ni = to_remove->nodes.at(i);
		if (Ni->fibres.size() == 1)
			mNodesToDelete.push_back(Ni); // delete the fibre nodes having no crosslink with other fibres 
		for (int j = 0; j < Ni->fibres.size(); j++)
		{
			ModelElementFibre * Fj = Ni->fibres.at(j);
			Ni->removeFibre(Fj); // remove the fibre from the nodes
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
			if (subDivide == 1)
			{
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
			}
			else
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
					double xj = (nodes.at(n2)->position.x * j + nodes.at(n1)->position.x * (subDivide - j))/subDivide;
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
					mpArena->addObject(nodeJ->GLObject()); // just for test
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
					}
					else
					{
						ModelElementLatticeNode *nodeJ1 = nodes.at(nodes.size() - 2);
						LinearSpring *lsNJ = new LinearSpring(nodeJ1,nodeJ);
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
					}
				}
				// add the last segment
				ModelElementLatticeNode *nodeJ1 = nodes.at(nodes.size() - 1);
				LinearSpring *lsNJ = new LinearSpring(nodeJ1,nodes.at(n2));
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
			}
			nodes.at(n1)->fibreEnd = true;
			nodes.at(n2)->fibreEnd = true; // the fibre end-node
		}
	}
	int nn = nodes.size();
	// add RotationalSpring
	for (int i = 0; i < nn; i++)
	{
		if (!nodes.at(i)->fibreEnd)
			continue;
		// only look for the end-nodes of fibres
		int sni = nodes.at(i)->springs.size();
		if (sni >= 2)
		{
			/***************************
			  o--o--o
				 |\
				 | \
				 o  o 
			 RotationalSpring is only 
			 stored in the center node
			***************************/
			ModelElementLatticeNode *ni = nodes.at(i);
			for (int j = 0; j < sni - 1; j++)
			{
				ModelElementLatticeNode *nj = ni->nodes.at(j);
				for (int k = j + 1; k < sni; k++)
				{
					ModelElementLatticeNode *nk = ni->nodes.at(k);
					RotationalSpring *rsN = new RotationalSpring(nj,nk,ni);
					ni->Rsprings.push_back(rsN);
					springs.push_back(rsN);
				}
			}
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
