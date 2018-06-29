/*
 * RuleBasedAgent.hpp
 *
 *  Created on: May 7, 2012
 *      Author: jagiella
 */

#ifndef RULEBASEDAGENT_HPP_
#define RULEBASEDAGENT_HPP_

#include "Agent.hpp"
#include "Process.hpp"

class RuleBasedAgent : public Agent {
private:
	Process **rules;
	int       ruleCount;
public:
	RuleBasedAgent();
	virtual ~RuleBasedAgent();

	void add( Process *rule);
	void remove();
	Process **getRules();
	int     countRules();
};

#endif /* RULEBASEDAGENT_HPP_ */
