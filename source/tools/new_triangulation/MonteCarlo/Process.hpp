/*
 * Process.hpp
 *
 *  Created on: 29.04.2012
 *      Author: jagiella
 */

#ifndef PROCESS_HPP_
#define PROCESS_HPP_

typedef float RATE_T;

//template <class StateType>
class Process {
	friend class ProcessTree;
private:
	// action tree attributes
	Process*	top;     // next action in probability list

	Process*	next;     // next action in probability list
	double	rateSumNext;
	int	sizeNext;

	Process*	prev;     // previous action in probability list
	double	rateSumPrev;
	int	sizePrev;
protected:
	RATE_T rate;
public:
	// constructor
	Process( RATE_T rate = 0.);

	// common methods
	bool	active();
	RATE_T	updateRate();
	void 	setRate( RATE_T value) { rate=value;}

	// specific methods
	virtual RATE_T getActualRate() = 0;
	virtual void perform() = 0;
	virtual bool condition() = 0;
	virtual void print();
};


class ProcessTree {
public:
	ProcessTree();
	~ProcessTree();

	void	add(    Process* action);
	void	add(    Process**action, int size);
	void	remove( Process* action);
	void	remove( Process** action, int size);

	RATE_T	rate();

	void   	update();
	void	update( Process* action);
	void	update( Process **p, int size);
	void	updateRate( Process* action, RATE_T newRate);

	Process*selectAction( RATE_T &time);
	void 	performGillespieStep( RATE_T &time);
	int		depth( Process* action);

	void	print();

		// attributes
	//double	rateSum;
	int		size;	// count of possible actions
	Process*root;	// first element in probability list
};

#endif /* PROCESS_HPP_ */
