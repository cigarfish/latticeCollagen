////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  ModelElementFibre.cpp                                         //
//                                                                            //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                //
//    Created:  2017-02-02                                                    //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "ModelElementFibre.h"
#include "ModelElementLatticeNode.h"
#include "LinearSpring.h"
// TODO this just doesnt seem right
#include "../../Core.h"

#include <algorithm>

ModelElementFibre::ModelElementFibre()
    : SimulationObject()
{
    //fibre_id.setXMLPath("Id");
    //registerParameter(fibre_id);
}

//ModelElementFibre::ModelElementFibre(XMLNode& node, std::stringstream& errors,
//                                     std::stringstream& warnings)
//    : ModelElementFibre()
//{
    // call function load here!
//}

ModelElementFibre::~ModelElementFibre()
{}

//void ModelElementFibre::load(XMLNode& node,
//                             std::stringstream& errors,
//                             std::stringstream& warnings)
//{
    // TODO IMPLEMENT
//}

//void ModelElementFibre::addNode(std::size_t NodeId)
// revised by Jieling
void ModelElementFibre::addNode(ModelElementLatticeNode *N)
{
    // check if it is already present
    if (std::find(std::begin(nodes), std::end(nodes), N) == nodes.end())
    {
        nodes.push_back(N);
    }
    else
    {
        std::cout << "Node is already in the fibre " << fibre_id << "." << std::endl;
    }
}

//void ModelElementFibre::removeNode(std::size_t NodeId)
// revised by Jieling
void ModelElementFibre::removeNode(ModelElementLatticeNode *N)
{
	std::vector<ModelElementLatticeNode*>::iterator foundN
		= std::find(std::begin(nodes), std::end(nodes), N);

	if ( foundN == nodes.end())
		std::cout << "Attempting to delete node but it can't be found in fibre " << fibre_id << "." << std::endl;
	else
		nodes.erase(foundN);
}

void ModelElementFibre::addSpring(LatticeSpring *S)
{
	std::vector<LatticeSpring*>::iterator foundS
		= std::find(springs.begin(), springs.end(), S);

	if (foundS != springs.end())
		std::cout << "Spring is already in the fibre" << std::endl;
	else
		springs.push_back(S);
}

void ModelElementFibre::removeSpring(LatticeSpring *S)
{
	std::vector<LatticeSpring*>::iterator foundS
		= std::find(springs.begin(), springs.end(), S);

	if (foundS == springs.end())
		std::cout << "Attempting to delete spring but it can't be found in fibre" << std::endl;
	else
		springs.erase(foundS);
}

void ModelElementFibre::unboundE(int end, ModelElementLatticeNode * N)
{
	if (end == 1)
	{
		// detach from node1
		ModelElementLatticeNode *n1O = n1;
		ModelElementLatticeNode *n1N = NULL;
		bool inFibre = false;
		for (int i = 0; i < n1->nodes.size(); i++)
		{
			ModelElementLatticeNode * Ni = n1->nodes.at(i);
			for (int j = 0; j < Ni->fibres.size(); j++)
			{
				if (Ni->fibres.at(j) == this)
				{
					inFibre = true;
					n1N = Ni;
					break;
				}
			}
			if (inFibre)
				break;
		}
		if (inFibre)
		{
			LatticeSpring * Si = locateSpring(n1, n1N);
			LinearSpring * Li = static_cast<LinearSpring*>(Si);
			n1N->removeNeighbor(n1);
			n1->removeNeighbor(n1N); // remove the neighboring nodes from the end-node
			n1->removeSpring(Si); // remove the spring from the end-node
			if (Si->nodes().at(0) == n1)
			{
				Li->changeN1(N);
				Li->n1_unbound = false;
				//Si->nodes().at(0) = N;
			}
			else
			{
				Li->changeN2(N);
				Li->n2_unbound = false;
				//Si->nodes().at(1) = N; // renew the spring's end-node
			}
			n1O->removeFibre(this); // remove the fibre from the old end-node
			N->addFibre(this); // add the fibre to the new end-node
			N->addSpring(Si); // add the spring to the new end-node
			n1N->addNeighbor(N);
			N->addNeighbor(n1N); // add neighbor node to the new end-node
			this->n1 = N; // change the end-node of the fibre
		}
	}
	else
	{
		// detach from node2
		ModelElementLatticeNode *n2O = n2;
		ModelElementLatticeNode *n2N = NULL;
		bool inFibre = false;
		for (int i = 0; i < n2->nodes.size(); i++)
		{
			ModelElementLatticeNode * Ni = n2->nodes.at(i);
			for (int j = 0; j < Ni->fibres.size(); j++)
			{
				if (Ni->fibres.at(j) == this)
				{
					inFibre = true;
					n2N = Ni;
					break;
				}
			}
			if (inFibre)
				break;
		}
		if (inFibre)
		{
			LatticeSpring * Si = locateSpring(n2, n2N);
			LinearSpring *Li = static_cast<LinearSpring*>(Si);
			n2N->removeNeighbor(n2);
			n2->removeNeighbor(n2N); // remove the neighboring nodes from the end-node
			n2->removeSpring(Si); // remove the spring from the end-node
			if (Si->nodes().at(0) == n2)
			{
				Li->changeN1(N);
				Li->n1_unbound = false;
				//Si->nodes().at(0) = N;
			}
			else
			{
				Li->changeN2(N);
				Li->n2_unbound = false;
				//Si->nodes().at(1) = N; // renew the spring's end-node
			}
			n2O->removeFibre(this); // remove the fibre from the old end-node
			N->addFibre(this); // add the fibre to the new end-node
			N->addSpring(Si); // add the spring to the new end-node
			n2N->addNeighbor(N);
			N->addNeighbor(n2N); // add neighbor node to the new end-node
			this->n2 = N; // change the end-node of the fibre
		}
	}
}

void ModelElementFibre::unboundM(ModelElementLatticeNode *O, ModelElementLatticeNode *N)
{
	// detach the neighbor nodes
	for (int i = 0; i < O->springs.size(); i++)
	{
		LatticeSpring *Si = O->springs.at(i);
		ModelElementLatticeNode *Ni = Si->nodes().at(0);
		if (Si->fibre != this)
			continue;
		// only for strings in the same fibre
		if (Ni == O)
		{ 
			// 0: O, 1: Ni 
			Ni = Si->nodes().at(1);
			Ni->removeNeighbor(O);
			O->removeNeighbor(Ni); // remove the neighbor relationship from the old node
			Ni->addNeighbor(N);
			N->addNeighbor(Ni);
			Si->nodes().at(0) = N; // move from O to N
		}
		else
		{
			// 1: O, 0: Ni
			Ni = Si->nodes().at(0);
			Ni->removeNeighbor(O);
			O->removeNeighbor(Ni); // remove the neighbor relationship from the old node
			Ni->addNeighbor(N);
			N->addNeighbor(Ni);
			Si->nodes().at(1) = N;
		}
	}
	O->removeFibre(this);
	this->removeNode(O); // remove the old node from the fibre
	N->addFibre(this);
	this->addNode(N); // add the new node into the fibre
}

LatticeSpring* ModelElementFibre::locateSpring(ModelElementLatticeNode* e1, ModelElementLatticeNode* e2)
{
	LatticeSpring *foundS = NULL;
	for (int i = 0; i < springs.size(); i++)
	{
		if ((springs.at(i)->nodes().at(0) == e1 && springs.at(i)->nodes().at(1) == e2) ||
			(springs.at(i)->nodes().at(1) == e1 && springs.at(i)->nodes().at(0) == e2))
		{
			foundS = springs.at(i);
			break;
		}
	}
	return foundS;
}

LatticeSpring* ModelElementFibre::getE1()
{
	for (int i = 0; i < n1->springs.size(); i++)
	{
		if (n1->springs.at(i)->fibre == this)
			return n1->springs.at(i);
	}
	return NULL;
}

LatticeSpring* ModelElementFibre::getE2()
{
	for (int i = 0; i < n2->springs.size(); i++)
	{
		if (n2->springs.at(i)->fibre == this)
			return n2->springs.at(i);
	}
	return NULL;
}

const vector<ModelElementLatticeNode*>& ModelElementFibre::getNodes() const
{
    return nodes;
}

void ModelElementFibre::initialize()
{}

void ModelElementFibre::update()
{}

bool ModelElementFibre::ready() const
{
    return true;
}

void ModelElementFibre::print(std::ostream& stream) const
{}





