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
	bending = false;
	length = 0.;
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

bool ModelElementFibre::checkNode(ModelElementLatticeNode *N)
{
	if (n1 == N)
		return true;
	else if (n2 == N)
		return true;
	else
	{
		for (int i = 0; i < nodes.size(); i++)
		{
			if (nodes.at(i) == N)
				return true;
		}
	}
	return false;
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

bool ModelElementFibre::checkSpring(LatticeSpring *S)
{
	std::vector<LatticeSpring*>::iterator foundS
		= std::find(springs.begin(), springs.end(), S);

	if (foundS != springs.end())
		return true;
	else
		return false;

	return false;
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

bool ModelElementFibre::locateNodeSpring(ModelElementLatticeNode *node1, ModelElementLatticeNode *node2, std::vector<ModelElementLatticeNode*> &ns, std::vector<LatticeSpring*> &ls)
{
	// return all nodes and springs between node1 and node2
	bool node1In = false, node2In = false;
	if (node1 == n1 || node1 == n2)
		node1In = true;
	if (node2 == n1 || node2 == n2)
		node2In = true;
	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes.at(i) == node1)
		{
			node1In = true;
			break;
		}
	}
	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes.at(i) == node2)
		{
			node2In = true;
			break;
		}
	}

	if (node1In && node2In)
	{
		ModelElementLatticeNode *startN = node1;
		LatticeSpring *startL = NULL;
		for (int i = 0; i < startN->springs.size(); i++)
		{
			if (startN->springs.at(i)->fibre == this)
			{
				startL = startN->springs.at(i); // the first spring contains startN
				break;
			}
		}
		if (startL == NULL)
			std::cout << "	-> Error: the spring containing the node is not found!" << std::endl;
		if (startN == n1)
		{
			// start from n1
			ModelElementLatticeNode* nextN = startN;
			if (startL->nodes().at(0) == startN)
				nextN = startL->nodes().at(1);
			else
				nextN = startL->nodes().at(0);
			ls.push_back(startL);
			ns.push_back(nextN);
			if (nextN == node2) // startN always starts from node1
			{
				// if node2 is the end node of the starting spring
				return true;
			}
			else
			{
				LatticeSpring *backL = startL;
				bool hit = false;
				while (!hit)
				{
					for (int i = 0; i < nextN->springs.size(); i++)
					{
						if (nextN->springs.at(i)->fibre == this)
						{
							// spring inside of this fibre
							if (nextN->springs.at(i) != backL)
							{
								startL = nextN->springs.at(i);
								break;
							}
						}
					}
					if (startL->nodes().at(0) == nextN)
					{
						nextN = startL->nodes().at(1);
						startN = startL->nodes().at(0);
					}
					else
					{
						nextN = startL->nodes().at(0);
						startN = startL->nodes().at(1);
					}
					backL = startL;
					ls.push_back(startL);
					ns.push_back(nextN);
					if (nextN == node2 || nextN == n2)
						hit = true;
				}
				return true;
			}
		}
		else if (startN == n2)
		{
			ModelElementLatticeNode * nextN = startN;
			if (startL->nodes().at(0) == startN)
				nextN = startL->nodes().at(1);
			else
				nextN = startL->nodes().at(0);
			ls.push_back(startL);
			ns.push_back(nextN);
			if (nextN == node2) // startN always starts from node1
			{
				// if node1 is the end node of the starting spring
				return true;
			}
			else
			{
				LatticeSpring *backL = startL;
				bool hit = false;
				while (!hit)
				{
					for (int i = 0; i < nextN->springs.size(); i++)
					{
						if (nextN->springs.at(i)->fibre == this)
						{
							// spring inside of this fibre
							if (nextN->springs.at(i) != backL)
							{
								startL = nextN->springs.at(i);
								break;
							}
						}
					}
					if (startL->nodes().at(0) == nextN)
					{
						nextN = startL->nodes().at(1);
						startN = startL->nodes().at(0);
					}
					else
					{
						nextN = startL->nodes().at(0);
						startN = startL->nodes().at(1);
					}
					backL = startL;
					ls.push_back(startL);
					ns.push_back(nextN);
					if (nextN == node2 || nextN == n1)
						hit = true;
				}
				return true;
			}
		}
		else
		{
			// startN (node1) is not the end node of the fibre
			// there are two possible directions to reach node2
			ModelElementLatticeNode* nextN = startN;
			int dir1Index = -1, dir2Index = -1; // trace the direction of node1, one of it will reach node2
			int dirIndex = -1; // the direction that node2 can be reached
			ModelElementLatticeNode *startN1 = startN, *startN2 = startN;
			ModelElementLatticeNode *nextN1 = startN, *nextN2 = startN;
			LatticeSpring *startL1 = NULL, *startL2 = NULL;
			for (int i = 0; i < startN->springs.size(); i++)
			{
				if (startN->springs.at(i)->fibre == this)
				{
					if (dir1Index == -1)
					{
						dir1Index = i;
						startL1 = startN->springs.at(i);
						if (startN->springs.at(i)->nodes().at(0) == startN)
							nextN1 = startN->springs.at(i)->nodes().at(1);
						else
							nextN1 = startN->springs.at(i)->nodes().at(0);
					}
					else if (dir1Index >= 0 && dir2Index == -1)
					{
						if (dir1Index != i)
						{
							dir2Index = i;
							startL2 = startN->springs.at(i);
							if (startN->springs.at(i)->nodes().at(0) == startN)
								nextN2 = startN->springs.at(i)->nodes().at(1);
							else
								nextN2 = startN->springs.at(i)->nodes().at(0);
							break;
						}
					}
				}
			}
			if (dir1Index < 0 || dir2Index < 0 ||
				(dir1Index >= 0 && dir2Index >= 0 && dir1Index == dir2Index))
			{
				std::cout << "	Error: node inside the fibre has wrong neighboring spring relationship!" << std::endl;
				return false;
			}
			else
			{
				// trace from startN -> nextN1
				if (nextN1 == node2)
				{
					dirIndex = dir1Index;
				}
				else
				{
					LatticeSpring *backL = startL1;
					bool hit = false;
					while (!hit)
					{
						bool foundNext = false;
						for (int i = 0; i < nextN1->springs.size(); i++)
						{
							if (nextN1->springs.at(i)->fibre == this)
							{
								// spring inside of this fibre
								if (nextN1->springs.at(i) != backL)
								{
									startL1 = nextN1->springs.at(i);
									foundNext = true;
									break;
								}
							}
						}
						if (foundNext)
						{
							if (startL1->nodes().at(0) == nextN1)
							{
								nextN1 = startL1->nodes().at(1);
								startN1 = startL1->nodes().at(0);
							}
							else
							{
								nextN1 = startL1->nodes().at(0);
								startN1 = startL1->nodes().at(1);
							}
							backL = startL1;
						}
						if (nextN1 == node2 || nextN1 == n1 || nextN1 == n2)
						{
							hit = true;
							if (nextN1 == node2)
								dirIndex = dir1Index;
						}
					}
				}
				// trace from startN -> nextN2
				if (nextN2 == node2)
				{
					dirIndex = dir2Index;
				}
				else
				{
					LatticeSpring *backL = startL2;
					bool hit = false;
					while (!hit)
					{
						bool foundNext = false;
						for (int i = 0; i < nextN2->springs.size(); i++)
						{
							if (nextN2->springs.at(i)->fibre == this)
							{
								// spring inside of this fibre
								if (nextN2->springs.at(i) != backL)
								{
									startL2 = nextN2->springs.at(i);
									foundNext = true;
									break;
								}
							}
						}
						if (foundNext)
						{
							if (startL2->nodes().at(0) == nextN2)
							{
								nextN2 = startL2->nodes().at(1);
								startN2 = startL2->nodes().at(0);
							}
							else
							{
								nextN2 = startL2->nodes().at(0);
								startN2 = startL2->nodes().at(1);
							}
							backL = startL2;
						}
						if (nextN2 == node2 || nextN2 == n1 || nextN2 == n2)
						{
							hit = true;
							if (nextN2 == node2)
								dirIndex = dir2Index;
						}
					}
				}
				// trace from node1 to node2
				startN = node1;
				startL = startN->springs.at(dirIndex);
				if (startL->nodes().at(0) == startN)
					nextN = startL->nodes().at(1);
				else
					nextN = startL->nodes().at(0);
				LatticeSpring *backL = startL;
				ls.push_back(startL);
				ns.push_back(nextN);
				if (nextN == node2) // startN always starts from node1
				{
					return true;
				}
				else
				{
					bool hit = false;
					while (!hit)
					{
						bool foundNext = false;
						for (int i = 0; i < nextN->springs.size(); i++)
						{
							if (nextN->springs.at(i)->fibre == this)
							{
								if (nextN->springs.at(i) != backL)
								{
									startL = nextN->springs.at(i);
									foundNext = true;
									break;
								}
							}
						}
						if (foundNext)
						{
							if (startL->nodes().at(0) == nextN)
							{
								nextN = startL->nodes().at(1);
								startN = startL->nodes().at(0);
							}
							else
							{
								nextN = startL->nodes().at(0);
								startN = startL->nodes().at(1);
							}
							backL = startL;
							ls.push_back(startL);
							ns.push_back(nextN);
						}
						if (nextN == node2)
						{
							hit = true;
						}
					}
					return true;
				}
			}
		}
		return true;
	}
	else
		return true;
}

int ModelElementFibre::relativePosition(ModelElementLatticeNode *node1, ModelElementLatticeNode *node2, 
	std::vector<ModelElementLatticeNode*> &n1s, std::vector<LatticeSpring*> &l1s,
	std::vector<ModelElementLatticeNode*> &n2s, std::vector<LatticeSpring*> &l2s)
{
	/**********************
	 get the relative position
	 of the input nodes and
	 the end nodes of the fibre

	 0: o-----o-----o-----o
	    n1  node1 node2   n2
	 1: o-----o-----o-----o
	    n1  node2 node1   n2
	 -1: node1/2 not in the fibre

	 n1s and l1s are the nodes/springs
	 between n1 and node1/2
	 n2s and l2s are the nodes/springs
	 between n2 and node2/1
	**********************/
	if (!node1->checkFibre(this) || !node2->checkFibre(this))
		return -1; // inquired nodes are not in the fibre

	if (n1 == node1 || n2 == node2)
		return 0;
	else if (n1 == node2 || n2 == node1)
		return 1;

	bool touch_node1 = false;
	bool touch_node2 = false;
	int type = 0;
	ModelElementLatticeNode *startN = n1;
	ModelElementLatticeNode *nextN = startN;
	LatticeSpring *startL = NULL;
	for (int i = 0; i < startN->springs.size(); i++)
	{
		if (startN->springs.at(i)->fibre == this)
		{
			startL = startN->springs.at(i); // the first spring containing startN
			break;
		}
	}
	if (startL == NULL)
		std::cout << "	-> Error: the spring containing the node is not found!" << std::endl;
	bool hit = false;
	if (startL->nodes().at(0) == startN)
		nextN = startL->nodes().at(1);
	else
		nextN = startL->nodes().at(0);
	LatticeSpring *backL = startL;
	l1s.push_back(startL);
	n1s.push_back(nextN);
	if (nextN == node1)
	{
		touch_node1 = true;
		hit = true;
	}
	else if (nextN == node2)
	{
		touch_node2 = true;
		hit = true;
	}
	else if (nextN == n2)
	{
		hit = true;
	}
	while (!hit)
	{
		bool foundNext = false;
		for (int i = 0; i < nextN->springs.size(); i++)
		{
			if (nextN->springs.at(i)->fibre == this)
			{
				// spring inside of this fibre
				if (nextN->springs.at(i) != backL)
				{
					startL = nextN->springs.at(i);
					foundNext = true;
					break;
				}
			}
		}
		if (foundNext)
		{
			if (startL->nodes().at(0) == nextN)
			{
				nextN = startL->nodes().at(1);
				startN = startL->nodes().at(0);
			}
			else
			{
				nextN = startL->nodes().at(0);
				startN = startL->nodes().at(1);
			}
			backL = startL;
			l1s.push_back(startL);
			n1s.push_back(nextN);
		}
		if (nextN == node1)
		{
			touch_node1 = true;
			hit = true;
		}
		else if (nextN == node2)
		{
			touch_node2 = true;
			hit = true;
		}
		else if (nextN == n2)
		{
			hit = true;
		}
	}
	//
	if (touch_node1 && !touch_node2)
		type = 0;
	else if (!touch_node1 && touch_node2)
		type = 1;
	else
		type = -1;
	// then run from n2 to n1
	startN = n2;
	nextN = startN;
	for (int i = 0; i < startN->springs.size(); i++)
	{
		if (startN->springs.at(i)->fibre == this)
		{
			startL = startN->springs.at(i); // the first spring containing startN
			break;
		}
	}
	if (startL == NULL)
		std::cout << "	-> Error: the spring containing the node is not found!" << std::endl;
	hit = false;
	if (startL->nodes().at(0) == startN)
		nextN = startL->nodes().at(1);
	else
		nextN = startL->nodes().at(0);
	backL = startL;
	l2s.push_back(startL);
	n2s.push_back(nextN);
	if (nextN == node2)
	{
		hit = true;
	}
	else if (nextN == node1)
	{
		hit = true;
	}
	else if (nextN == n1)
	{
		hit = true;
	}
	while (!hit)
	{
		bool foundNext = false;
		for (int i = 0; i < nextN->springs.size(); i++)
		{
			if (nextN->springs.at(i)->fibre == this)
			{
				// spring inside of this fibre
				if (nextN->springs.at(i) != backL)
				{
					startL = nextN->springs.at(i);
					foundNext = true;
					break;
				}
			}
		}
		if (foundNext)
		{
			if (startL->nodes().at(0) == nextN)
			{
				nextN = startL->nodes().at(1);
				startN = startL->nodes().at(0);
			}
			else
			{
				nextN = startL->nodes().at(0);
				startN = startL->nodes().at(1);
			}
			backL = startL;
			l2s.push_back(startL);
			n2s.push_back(nextN);
		}
		if (nextN == node2)
		{
			hit = true;
		}
		else if (nextN == node1)
		{
			hit = true;
		}
		else if (nextN == n1)
		{
			hit = true;
		}
	}
	//
	return type;
}

bool ModelElementFibre::relativeSinglePosition(ModelElementLatticeNode *node,
	std::vector<ModelElementLatticeNode*> &ns, std::vector<LatticeSpring*> &ls)
{
	/**********************
	get the relative position
	of the input nodes and
	the end nodes of the fibre

	   o-----o-----o
	   n1  node   n2

	ns and ls are the nodes/springs
	between n1 and node
	**********************/
	if (!node->checkFibre(this))
		return false; // inquired node is not in the fibre

	if (n1 == node || n2 == node)
		return false;

	bool touch_node = false;
	ModelElementLatticeNode *startN = n1;
	ModelElementLatticeNode *nextN = startN;
	LatticeSpring *startL = NULL;
	for (int i = 0; i < startN->springs.size(); i++)
	{
		if (startN->springs.at(i)->fibre == this)
		{
			startL = startN->springs.at(i); // the first spring containing startN
			break;
		}
	}
	if (startL == NULL)
		std::cout << "	-> Error: the spring containing the node is not found!" << std::endl;
	bool hit = false;
	if (startL->nodes().at(0) == startN)
		nextN = startL->nodes().at(1);
	else
		nextN = startL->nodes().at(0);
	LatticeSpring *backL = startL;
	ls.push_back(startL);
	ns.push_back(nextN);
	if (nextN == node)
	{
		touch_node = true;
		hit = true;
	}
	while (!hit)
	{
		bool foundNext = false;
		for (int i = 0; i < nextN->springs.size(); i++)
		{
			if (nextN->springs.at(i)->fibre == this)
			{
				// spring inside of this fibre
				if (nextN->springs.at(i) != backL)
				{
					startL = nextN->springs.at(i);
					foundNext = true;
					break;
				}
			}
		}
		if (foundNext)
		{
			if (startL->nodes().at(0) == nextN)
			{
				nextN = startL->nodes().at(1);
				startN = startL->nodes().at(0);
			}
			else
			{
				nextN = startL->nodes().at(0);
				startN = startL->nodes().at(1);
			}
			backL = startL;
			ls.push_back(startL);
			ns.push_back(nextN);
		}
		if (nextN == node)
		{
			touch_node = true;
			hit = true;
		}
		else if (nextN == n2)
		{
			hit = true;
		}
	}

	return true;
}

ModelElementLatticeNode* ModelElementFibre::splitNode(ModelElementLatticeNode *node1, ModelElementLatticeNode *node2)
{
	ModelElementLatticeNode *sN = NULL;
	/************************
	          o sN
	         / \
	  node1 o   o node2
           /    |
	************************/
	std::vector<LatticeSpring*> ls;
	std::vector<ModelElementLatticeNode*> ns;
	// find all nodes between node1 and node2
	// ls includes all springs between node1 and node2
	// ns includes all nodes between node1 and node2 (excluding node1, including node2)
	locateNodeSpring(node1, node2, ns, ls);
	if (ns.size() > 1)
	{
		int cutIndex = -1;
		int cutNeighborSize = 0;
		for (int i = 0; i < ns.size(); i++)
		{
			if (ns.at(i) != node2 && ns.at(i) != node1)
			{
				if (ns.at(i)->nodes.size() >= cutNeighborSize)
				{
					cutIndex = i;
					cutNeighborSize = ns.at(i)->nodes.size();
				}
			}
		}
		if (cutIndex < 0)
			std::cout << "	> Error: the cutting point node is not found in the fibre " << n1->mGlobalIndex << " - " << n2->mGlobalIndex << std::endl;
		sN = ns.at(cutIndex);
	}
	return sN;
}

bool ModelElementFibre::selfMerge()
{
	// if the fibre forms a circle
	if (nodes.size() > 0)
	{

	}
	return false;
}

void ModelElementFibre::addVirtualNodes()
{
	// after adding nodes for the fibre, add two virtual nodes on the two ends of the fibre
	ModelElementLatticeNode *n1n = NULL;
	ModelElementLatticeNode *n2n = NULL;
	LatticeSpring *ln1 = getE1();
	LatticeSpring *ln2 = getE2();
	if (ln1->nodes().at(0) == n1)
	{
		n1n = ln1->nodes().at(1);
	}
	else
	{
		n1n = ln1->nodes().at(0);
	}
	if (ln2->nodes().at(0) == n2)
	{
		n2n = ln2->nodes().at(1);
	}
	else
	{
		n2n = ln2->nodes().at(0);
	}
	mVirtualN1.x = n1->position.x + (n1->position.x - n1n->position.x);
	mVirtualN1.y = n1->position.y + (n1->position.y - n1n->position.y);
	mVirtualN1.z = n1->position.z + (n1->position.z - n1n->position.z);
	mVirtualN2.x = n2->position.x + (n2->position.x - n2n->position.x);
	mVirtualN2.y = n2->position.y + (n2->position.y - n2n->position.y);
	mVirtualN2.z = n2->position.z + (n2->position.z - n2n->position.z);
}

double ModelElementFibre::getEndBendingForce(int endType)
{
	if (endType == 1)
	{
		// n1
		// virtual1 - n1 - n1n
		ModelElementLatticeNode *n1n = NULL;
		LatticeSpring *n1s = getE1();
		if (n1s->nodes().at(0) == n1)
			n1n = n1s->nodes().at(1);
		else
			n1n = n1s->nodes().at(0);
		double vx = mVirtualN1.x;
		double vy = mVirtualN1.y;
		double vz = mVirtualN1.z;
		double dist = (n1n->position - n1n->oCoord).Norm();
		if (dist > 1e-9)
		{
			// treat n1 - n1no as x-axis
			double xx = n1n->oCoord.x - n1->position.x;
			double xy = n1n->oCoord.y - n1->position.y;
			double xz = n1n->oCoord.z - n1->position.z;
			double xd = sqrt(xx*xx + xy*xy + xz*xz);
			xx /= xd;
			xy /= xd;
			xz /= xd;
			// get the normal of plane n1-n1n and n1-n1no, treat it as z-axis
			double zx = xy * (n1n->position.z - n1->position.z) -
				xz * (n1n->position.y - n1->position.y);
			double zy = xz * (n1n->position.x - n1->position.x) -
				xx * (n1n->position.z - n1->position.z);
			double zz = xx * (n1n->position.y - n1->position.y) -
				xy * (n1n->position.x - n1->position.x);
			double zd = sqrt(zx*zx + zy*zy + zz*zz);
			zx /= zd;
			zy /= zd;
			zz /= zd;
			// get the normal of x-axis and z-axis, treat it as y-axis
			double yx = zy*xz - zz*xy;
			double yy = zz*xx - zx*xz;
			double yz = zx*xy - zy*xx;
			// project n1n-n1 onto the new coordinate system
			double n1nonx = (n1n->position.x - n1->position.x) * xx + (n1n->position.y - n1->position.y) * xy + (n1n->position.z - n1->position.z) * xz;
			double n1nony = (n1n->position.x - n1->position.x) * yx + (n1n->position.y - n1->position.y) * yy + (n1n->position.z - n1->position.z) * yz;
			double n1nonz = (n1n->position.x - n1->position.x) * zx + (n1n->position.y - n1->position.y) * zy + (n1n->position.z - n1->position.z) * zz;
			// project the virtual node into the new coordinate system
			/*
				| a1 b1 c1 | vx   e1
				| a2 b2 c2 | vy = e2, Av = e
				| a3 b3 c3 | vz   e3

							 | A D G |
				A^-1 = det(A)| B E H |
							 | C F I |
			*/
			double a1 = xx;
			double b1 = xy;
			double c1 = xz;
			double e1 = -n1nonx + n1->position.x * xx + n1->position.y * xy + n1->position.z * xz;
			double a2 = yx;
			double b2 = yy;
			double c2 = yz;
			double e2 = n1nony + n1->position.x * yx + n1->position.y * yy + n1->position.z * yz;
			double a3 = zx;
			double b3 = zy;
			double c3 = zz;
			double e3 = n1nonz + n1->position.x * zx + n1->position.y * zy + n1->position.z * zz;
			// the inverse of A
			double detA = a1 * (b2 * c3 - c2 * b3) -
				a2 * (b1 * c3 - c1 * b3) +
				a3 * (b1 * c2 - c1 * b2);
			double A = (b2 * c3 - c2 * b3);
			double D = -(b1 * c3 - c1 * b3);
			double G = (b1 * c2 - c1 * b2);
			double B = -(a2 * c3 - c2 * a3);
			double E = (a1 * c3 - c1 * a3);
			double H = -(a1 * c2 - c1 * a2);
			double C = (a2 * b3 - b2 * a3);
			double F = -(a1 * b3 - b1 * a3);
			double I = (a1 * b2 - b1 * a2);
			vx = (A * e1 + D * e2 + G * e3) / detA;
			vy = (B * e1 + E * e2 + H * e3) / detA;
			vz = (C * e1 + F * e2 + I * e3) / detA;
			double vd = sqrt(vx*vx + vy*vy + vz*vz);
			//
			double virtualod = (mVirtualN1 - n1->position).Norm();
			vx = n1->position.x + (vx - n1->position.x) * virtualod / vd;
			vy = n1->position.y + (vy - n1->position.y) * virtualod / vd;
			vz = n1->position.z + (vz - n1->position.z) * virtualod / vd;
		}
		//
		double lnx = n1n->position.x - n1->position.x;
		double lny = n1n->position.y - n1->position.y;
		double lnz = n1n->position.z - n1->position.z;
		double lnd = sqrt(lnx * lnx + lny * lny + lnz * lnz);
		//
		double lnm1x = n1->position.x - vx;
		double lnm1y = n1->position.y - vy;
		double lnm1z = n1->position.z - vz;
		double lnm1d = sqrt(lnm1x * lnm1x + lnm1y * lnm1y + lnm1z * lnm1z);
		//
		double aR = (n1->mRadius + n1n->mRadius) / 2;
		double aY = (n1->mYoung + n1n->mYoung) / 2;
		double aS = aY * aR * aR * aR * aR * 3.14159265 / 4.; // EI
		//
		double lnlnm1 = lnx * lnm1x + lny * lnm1y + lnz * lnm1z;
		//
		double dx = (lnx - lnm1x) / lnd / lnm1d + lnlnm1 / lnd / lnm1d * (lnx / lnd / lnd - lnm1x / lnm1d / lnm1d);
		double dy = (lny - lnm1y) / lnd / lnm1d + lnlnm1 / lnd / lnm1d * (lny / lnd / lnd - lnm1y / lnm1d / lnm1d);
		double dz = (lnz - lnm1z) / lnd / lnm1d + lnlnm1 / lnd / lnm1d * (lnz / lnd / lnd - lnm1z / lnm1d / lnm1d);
		//
		double scaleF = aS / lnd;
		// calculate the moment:
		double curv = sqrt(2 * (1. - lnlnm1 / lnd / lnm1d));
		double moment = aS * curv / lnd;
		n1->mMoment = moment;
		//
		dx = -1. / lnd / lnm1d * (lnlnm1 * lnx / lnd / lnd - lnm1x);
		dy = -1. / lnd / lnm1d * (lnlnm1 * lny / lnd / lnd - lnm1y);
		dz = -1. / lnd / lnm1d * (lnlnm1 * lnz / lnd / lnd - lnm1z);
		// this dx, dy, dz is cos(theta(a)) for ra of d(a-1) and d(a), for the n2
		// assign it to n2
		n1n->mRotationalForce.x += dx * scaleF;
		n1n->mRotationalForce.y += dy * scaleF;
		n1n->mRotationalForce.z += dz * scaleF;
		n1n->directedForce.x += dx * scaleF;
		n1n->directedForce.y += dy * scaleF;
		n1n->directedForce.z += dz * scaleF;

		return moment;
	}
	else
	{
		// n2
		ModelElementLatticeNode *n2n = NULL;
		LatticeSpring *n2s = getE2();
		if (n2s->nodes().at(0) == n2)
			n2n = n2s->nodes().at(1);
		else
			n2n = n2s->nodes().at(0);
		double vx = mVirtualN2.x;
		double vy = mVirtualN2.y;
		double vz = mVirtualN2.z;
		double dist = (n2n->position - n2n->oCoord).Norm();
		if (dist > 1e-9)
		{
			// treat n2 - n2no as x-axis
			double xx = n2n->oCoord.x - n2->position.x;
			double xy = n2n->oCoord.y - n2->position.y;
			double xz = n2n->oCoord.z - n2->position.z;
			double xd = sqrt(xx*xx + xy*xy + xz*xz);
			xx /= xd;
			xy /= xd;
			xz /= xd;
			// get the normal of plane n1-n1n and n1-n1no, treat it as z-axis
			double zx = xy * (n2n->position.z - n2->position.z) -
				xz * (n2n->position.y - n2->position.y);
			double zy = xz * (n2n->position.x - n2->position.x) -
				xx * (n2n->position.z - n2->position.z);
			double zz = xx * (n2n->position.y - n2->position.y) -
				xy * (n2n->position.x - n2->position.x);
			double zd = sqrt(zx*zx + zy*zy + zz*zz);
			zx /= zd;
			zy /= zd;
			zz /= zd;
			// get the normal of x-axis and z-axis, treat it as y-axis
			double yx = zy*xz - zz*xy;
			double yy = zz*xx - zx*xz;
			double yz = zx*xy - zy*xx;
			// project n2n-n2 onto the new coordinate system
			double n2nonx = (n2n->position.x - n2->position.x) * xx + (n2n->position.y - n2->position.y) * xy + (n2n->position.z - n2->position.z) * xz;
			double n2nony = (n2n->position.x - n2->position.x) * yx + (n2n->position.y - n2->position.y) * yy + (n2n->position.z - n2->position.z) * yz;
			double n2nonz = (n2n->position.x - n2->position.x) * zx + (n2n->position.y - n2->position.y) * zy + (n2n->position.z - n2->position.z) * zz;
			// project the virtual node into the new coordinate system
			/*
				| a1 b1 c1 | vx   e1
				| a2 b2 c2 | vy = e2, Av = e
				| a3 b3 c3 | vz   e3

							 | A D G |
				A^-1 = det(A)| B E H |
							 | C F I |
			*/
			double a1 = xx;
			double b1 = xy;
			double c1 = xz;
			double e1 = -n2nonx + n2->position.x * xx + n2->position.y * xy + n2->position.z * xz;
			double a2 = yx;
			double b2 = yy;
			double c2 = yz;
			double e2 = n2nony + n2->position.x * yx + n2->position.y * yy + n2->position.z * yz;
			double a3 = zx;
			double b3 = zy;
			double c3 = zz;
			double e3 = n2nonz + n2->position.x * zx + n2->position.y * zy + n2->position.z * zz;
			// the inverse of A
			double detA = a1 * (b2 * c3 - c2 * b3) -
				a2 * (b1 * c3 - c1 * b3) +
				a3 * (b1 * c2 - c1 * b2);
			double A = (b2 * c3 - c2 * b3);
			double D = -(b1 * c3 - c1 * b3);
			double G = (b1 * c2 - c1 * b2);
			double B = -(a2 * c3 - c2 * a3);
			double E = (a1 * c3 - c1 * a3);
			double H = -(a1 * c2 - c1 * a2);
			double C = (a2 * b3 - b2 * a3);
			double F = -(a1 * b3 - b1 * a3);
			double I = (a1 * b2 - b1 * a2);
			vx = (A * e1 + D * e2 + G * e3) / detA;
			vy = (B * e1 + E * e2 + H * e3) / detA;
			vz = (C * e1 + F * e2 + I * e3) / detA;
			double vd = sqrt(vx*vx + vy*vy + vz*vz);
			//
			double virtualod = (mVirtualN2 - n2->position).Norm();
			vx = n2->position.x + (vx - n2->position.x) * virtualod / vd;
			vy = n2->position.y + (vy - n2->position.y) * virtualod / vd;
			vz = n2->position.z + (vz - n2->position.z) * virtualod / vd;
		}
		//
		double lnx = vx - n2->position.x;
		double lny = vy - n2->position.y;
		double lnz = vz - n2->position.z;
		double lnd = sqrt(lnx * lnx + lny * lny + lnz * lnz);
		//
		double lnm1x = n2->position.x - n2n->position.x;
		double lnm1y = n2->position.y - n2n->position.y;
		double lnm1z = n2->position.z - n2n->position.z;
		double lnm1d = sqrt(lnm1x * lnm1x + lnm1y * lnm1y + lnm1z * lnm1z);
		//
		double lnlnm1 = lnx * lnm1x + lny * lnm1y + lnz * lnm1z;
		//
		double dx = (lnx - lnm1x) / lnd / lnm1d + lnlnm1 / lnd / lnm1d * (lnx / lnd / lnd - lnm1x / lnm1d / lnm1d);
		double dy = (lny - lnm1y) / lnd / lnm1d + lnlnm1 / lnd / lnm1d * (lny / lnd / lnd - lnm1y / lnm1d / lnm1d);
		double dz = (lnz - lnm1z) / lnd / lnm1d + lnlnm1 / lnd / lnm1d * (lnz / lnd / lnd - lnm1z / lnm1d / lnm1d);
		//
		double aR = (n2->mRadius + n2n->mRadius) / 2;
		double aY = (n2->mYoung + n2n->mYoung) / 2;
		double aS = aY * aR * aR * aR * aR * 3.14159265 / 4.;
		double scaleF = aS / lnd;
		// calculate the moment:
		double curv = sqrt(2 * (1. - lnlnm1 / lnd / lnm1d));
		double moment = aS * curv / lnd;
		n2->mMoment = moment;
		dx = 1. / lnd / lnm1d * (lnlnm1 / lnm1d / lnm1d * lnm1x - lnx);
		dy = 1. / lnd / lnm1d * (lnlnm1 / lnm1d / lnm1d * lnm1y - lny);
		dz = 1. / lnd / lnm1d * (lnlnm1 / lnm1d / lnm1d * lnm1z - lnz);
		// this dx, dy, dz is cos(theta(a+2)) for ra of d(a+1) and d(a+2), for the n1
		// assign it to n1
		n2n->mRotationalForce.x += dx * scaleF;
		n2n->mRotationalForce.y += dy * scaleF;
		n2n->mRotationalForce.z += dz * scaleF;
		n2n->directedForce.x += dx * scaleF;
		n2n->directedForce.y += dy * scaleF;
		n2n->directedForce.z += dz * scaleF;

		return moment;
	}
	return 0.;
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





