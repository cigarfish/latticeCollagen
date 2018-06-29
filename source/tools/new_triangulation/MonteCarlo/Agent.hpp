/*
 * Agent.h
 *
 *  Created on: May 4, 2012
 *      Author: jagiella
 */

#ifndef AGENT_HPP_
#define AGENT_HPP_

template <class AgentType>
class AgentList;

class Agent {
	template <class AgentType>
	friend class AgentList;
private:
	int index;
public:
	Agent();
	virtual ~Agent();

	int getIndex();
};

template <class AgentType = Agent>
class AgentList {
private:
	// attributes
	AgentType ** agents;
	int countAgents;
	int countActiveAgents;

public:
	AgentList( int countNewAgents = 0);
	~AgentList();

	// methodes
	AgentType* activateAgent();
	void deactivateAgent( AgentType* agent);

	int size();
	int sizeActive();
	AgentType* get( int i);
	AgentType* operator[]( int i);
	//void print();
	//void printActiveAgents();
	//void printToPovray( const char *filename, VoronoiDiagram *vd);
};

#include "Agent.ipp"

#endif /* AGENT_H_ */
