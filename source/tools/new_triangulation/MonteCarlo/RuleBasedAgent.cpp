/*
 * RuleBasedAgent.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jagiella
 */

#include "RuleBasedAgent.hpp"

RuleBasedAgent::RuleBasedAgent() {
	// TODO Auto-generated constructor stub
	ruleCount=0;
	rules=0;
}

RuleBasedAgent::~RuleBasedAgent() {
	// TODO Auto-generated destructor stub
}

void RuleBasedAgent::add( Process *rule)
{
	rules = (Process**) realloc( rules, sizeof(Process*) * (ruleCount+1));
	rules[ruleCount++] = rule;
}

void RuleBasedAgent::remove()
{
	for( int i=0; i<ruleCount; i++)
		delete rules[i];
	ruleCount = 0;
	free( rules);
	rules = 0;
}

Process **RuleBasedAgent::getRules()
{
	return rules;
}

int RuleBasedAgent::countRules()
{
	return ruleCount;
}
