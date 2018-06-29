/*
 * Agent.ipp
 *
 *  Created on: May 7, 2012
 *      Author: jagiella
 */

#ifndef AGENT_IPP_
#define AGENT_IPP_

#include <stdio.h>
#include <stdlib.h>

template <class Agent>
AgentList<Agent>::AgentList( int countNewAgents) {
	// TODO Auto-generated constructor stub
	int i;

	agents = (Agent **) calloc( countNewAgents, sizeof(Agent *));
	for( i=0; i<countNewAgents; i++){
		agents[i] = new Agent();
		agents[i]->index = i;
	}

	countActiveAgents = 0;
	countAgents = countNewAgents;
}

template <class Agent>
Agent* AgentList<Agent>::activateAgent()
{
	if( countActiveAgents==countAgents){
		this->agents = (Agent **) realloc( this->agents, (countAgents+1000) * sizeof(Agent *));
		for( int i=countAgents; i<countAgents+1000; i++){
			this->agents[i] = new Agent();
			this->agents[i]->index = i;
		}
		countAgents += 1000;
	}

	if( countActiveAgents<countAgents){
		/*this->agents[ countActiveAgents]->initAgent();
		this->agents[ countActiveAgents]->cellCount = 0;
		this->agents[ countActiveAgents]->maxCellCount = CountCellsPerVoronoiCell;
		this->agents[ countActiveAgents]->growingTumorCellCount = 0;
		this->agents[ countActiveAgents]->dividingTumorCellCount = 0;
		this->agents[ countActiveAgents]->divide = 1;*/

		return this->agents[ countActiveAgents++];
	}else{
		fprintf( stderr, "ERROR in Agent* AgentList::activateAgent()\nAgentList is already full!\n");
		exit( 0);
	}

}

template <class Agent>
void AgentList<Agent>::deactivateAgent( Agent* agent)
{
	int tempIndex    = agent->index;
	//Agent* tempAgent = agent;

	// decrease number of active agents
	--(this->countActiveAgents);

	// replace deactivated agent by another active one
	this->agents[ tempIndex]        = this->agents[ this->countActiveAgents];
	this->agents[ tempIndex]->index = tempIndex;

	// move deactivated agent
	this->agents[ this->countActiveAgents] = agent;
	agent->index = this->countActiveAgents;
}

template <class Agent>
int AgentList<Agent>::size() {return countAgents;}

template <class Agent>
int AgentList<Agent>::sizeActive() {return countActiveAgents;}


#endif /* AGENT_IPP_ */
