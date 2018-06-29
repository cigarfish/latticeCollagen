/*
 * Process.cpp
 *
 *  Created on: May 2, 2012
 *      Author: jagiella
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "Process.hpp"


const char *actionTypeToString( int type){
	switch( type){
		/*case GROWTH:
		return "GROWTH";

		case DIVISION:
		return "DIVISION";

		case MIGRATION:
		return "MIGRATION";

		case NECROSIS:
		return "NECROSIS";

		case LYSIS:
		return "LYSIS";

		case VESSEL_GROWTH:
		return "VESSEL_GROWTH";*/

		default:
		return "UNKNOWN";
	}
	return "";
}

void ProcessTree::update( Process **p, int size)
{
	for( int i=0; i<size; i++)
		update( p[i]);
}


void ProcessTree::update( Process *p)
{
	// condition
	if( p->condition()){
		if( !p->active()){
			// update rate first
			p->updateRate();

			// add to tree
			this->add( p);
		}
		else
			// update rate only
			updateRate( p, p->getActualRate());
	}else
		if( p->active())
			// remove from tree
			this->remove( p);
}

void ProcessTree::update()
{
	Process* action = this->root;
	if( this->size > 0){

		//Process *actionStack[ this->size ], *actualAction;
		Process **actionStack = (Process**) malloc( sizeof(Process*)*this->size ), *actualAction;
		//double  rateSumStack[ this->size];
		//double  rateDifferenceStack[ this->size];
		int stackSize = 0, lastStackSize = 0;
		int i;

		for( i=0; i<this->size; i++){
			actionStack[ i] = NULL;
			//rateDifferenceStack[i] = 0.;
			//rateSumStack[i] = 0.;
		}

		actionStack[stackSize++] = action;
		do{
			actualAction = actionStack[stackSize-1];
			// TEST
			if( actualAction->next == actualAction || actualAction->next == actualAction){
				fprintf( stderr, "ERROR: in ActionTree::actualizeAllRates()\nactualAction->next == actualAction || actualAction->next == actualAction\n");
				exit( 0);
			}
			// root condition
			if( actualAction->top == actualAction && actualAction!=this->root){
				fprintf( stderr, "ERROR: in ActionTree::actualizeAllRates()\nactualAction->top == actualAction && actualAction!=this->root\n");
				exit( 0);
			}
			// TEST

			if( lastStackSize < stackSize){
				// climbing tree -> prev
				lastStackSize = stackSize;
				if( actualAction->sizePrev > 0){
					//rateDifferenceStack[stackSize] = 0.;
					actionStack[stackSize++] = actualAction->prev;
					//fprintf( stderr, "CLIMBING LEFT: stackSize = %i\n", stackSize);
				}//else
					//fprintf( stderr, "RESTING: stackSize = %i\n", stackSize);

			}else{
				lastStackSize = stackSize;
				// climb tree -> next
				if( actualAction->sizeNext > 0 &&  actualAction->next != actionStack[stackSize]){
					//rateDifferenceStack[stackSize] = 0.;
					actionStack[stackSize++] = actualAction->next;
					//fprintf( stderr, "CLIMBING RIGHT: stackSize = %i\n", stackSize);
				}
				// decend tree
				else{
					// adapt own rate
					//double oldRate = actualAction->rate;
					//rateDifferenceStack[stackSize-1] += (actualAction->actualizeRate() - oldRate);
					actualAction->updateRate();

					stackSize--;
					//fprintf( stderr, "DECENDING: stackSize = %i, from node %i\n", stackSize, actualAction->type);

					// adapt parent node
					if( stackSize>0){
						if( actionStack[stackSize-1]->prev == actualAction)
							actionStack[stackSize-1]->rateSumPrev = actionStack[stackSize]->rateSumPrev + actionStack[stackSize]->rateSumNext + actionStack[stackSize]->rate;
							//actionStack[stackSize-1]->rateSumPrev += rateDifferenceStack[stackSize];
						else
							actionStack[stackSize-1]->rateSumNext = actionStack[stackSize]->rateSumPrev + actionStack[stackSize]->rateSumNext + actionStack[stackSize]->rate;
							//actionStack[stackSize-1]->rateSumNext += rateDifferenceStack[stackSize];

						//rateDifferenceStack[stackSize-1] += rateDifferenceStack[stackSize];
					}else{
						//this->rateSum += rateDifferenceStack[stackSize];
						//this->rateSum = actionStack[stackSize]->rateSumPrev + actionStack[stackSize]->rateSumNext + actionStack[stackSize]->rate;
					}

				}
			}
		}while( stackSize>0);

		free( actionStack);
	}

}


ProcessTree::ProcessTree()
{
	size	= 0;
	root	= 0;
}



Process::Process( RATE_T rate)
{
	//this->type	= type;
	//this->rateSingleCell	= rate;
	this->rate = rate;
	//this->state = state;
//	this->internalState	= 0;
//	this->originalCell	= NULL;
//	this->destinationCell	= NULL;

	this->top	= 0;

	this->prev	= 0;
	this->rateSumPrev	= 0.;
	this->sizePrev	= 0;

	this->next	= 0;
	this->rateSumNext	= 0.;
	this->sizeNext	= 0;
}

//void Process::perform()
//{}

inline RATE_T Process::getActualRate()
{return rate;}

RATE_T Process::updateRate()
{
	this->rate = this->getActualRate();

	return this->rate;
}

/*RATE_T Process::getActualRate()
{
	RATE_T newRate = 0.;

	switch( this->type){

	}

	return 0;
}*/

void Process::print()
{
	fprintf( stderr, "[rate=%lf] -> prevSize=%i | nextSize=%i\n",/* actionTypeToString( this->type),*/ this->rate, this->sizePrev, this->sizeNext);
	if( this->sizePrev>0){
		fprintf( stderr, "prev [size=%i, rateSum=%lf]:\n", this->sizePrev, this->rateSumPrev);
		this->prev->print();
	}
	if( this->sizeNext>0){
		fprintf( stderr, "next [size=%i, rateSum=%lf]:\n", this->sizeNext, this->rateSumNext);
		this->next->print();
	}
}

bool Process::active()
{
	return top!=0;
}

void ProcessTree::print()
{
	fprintf( stderr, "root [size=%i, rateSum=%lf]:\n", this->size, rate());
	if( this->size>0){
		this->root->print();
	}
}

RATE_T ProcessTree::rate()
{
	//return this->rateSum;
	return ( root ? root->rate + root->rateSumPrev + root->rateSumNext : 0);
}


void ProcessTree::add( Process** action, int size)
{
//	fprintf( stderr, "%i\n", sizeof(Process*));
//	fprintf( stderr, "%i\n", sizeof(action));
	for( int p=0; p<size; p++)
		add( action[p]);
}

void ProcessTree::remove( Process** action, int size)
{
//	fprintf( stderr, "%i\n", sizeof(Process*));
//	fprintf( stderr, "%i\n", sizeof(action));
	for( int p=0; p<size; p++)
		remove( action[p]);
}

void ProcessTree::add( Process* action)
{
	//fprintf( stderr, "ActionTree::addAction( %s)\n", actionTypeToString( action->type));

	if( action->top != 0){
		//fprintf( stderr, "ERROR: Action is already element of action tree!\n");
		return;
	}

	/********        add action to prob_list    **********/

	// init action
	//action->actualizeRate();
	//fprintf( stderr, "AddAction(%s): rate=%lf\n", actionTypeToString( action->type), action->rate);

	//action->top = ;     // next action in probability list

	action->next = 0;     // next action in probability list
	action->rateSumNext = 0.;
	action->sizeNext = 0;

	action->prev = 0;     // previous action in probability list
	action->rateSumPrev = 0.;
	action->sizePrev = 0;


	if( this->size == 0){
		this->root = action;
		this->root->top = this->root;
		//this->rateSum = action->rate;
		this->size++;
		return;
	}else{
		Process *actualAction = this->root;
		do{
			if( actualAction->sizePrev <= actualAction->sizeNext){
				actualAction->sizePrev++;
				/* SLOW */
				/*if( actualAction->sizePrev>1)
					actualAction->rateSumPrev = action->rate + actualAction->prev->rateSumPrev + actualAction->prev->rateSumNext + actualAction->prev->rate;
				else
					actualAction->rateSumPrev = action->rate;*/
				/* SLOW END */
				/* FAST */ actualAction->rateSumPrev += action->rate;
				if( actualAction->sizePrev > 1)
					actualAction = actualAction->prev;
				else{
					actualAction->prev = action;
					action->top = actualAction;
				}

			}else{
				actualAction->sizeNext++;
				/* FAST */ actualAction->rateSumNext += action->rate;
				/* SLOW */
				/*if( actualAction->sizeNext>1)
					actualAction->rateSumNext = action->rate + actualAction->next->rateSumPrev + actualAction->next->rateSumNext + actualAction->next->rate;
				else
					actualAction->rateSumNext = action->rate;*/
				/* SLOW END*/
				if( actualAction->sizeNext > 1)
					actualAction = actualAction->next;
				else{
					actualAction->next = action;
					action->top = actualAction;
				}
			}
		}while( action->top == 0);
		/* FAST */ //this->rateSum += action->rate;
		/* SLOW */ //this->rateSum = this->root->rateSumPrev + this->root->rateSumNext + this->root->rate;
		this->size++;
	}

}

void ProcessTree::remove( Process* action)
{
	//fprintf( stderr, "ActionTree::deleteAction( %s)\n", actionTypeToString( action->type));

	/********        add action to prob_list    **********/
	if( action->top != NULL){
		// tree until root
		/* SLOW */ //action->rate=0.;
		Process *actualAction = action;
		while( actualAction != this->root){
			if( actualAction == actualAction->top->prev){
				// actualize prev attributes
				actualAction->top->sizePrev--;
				/* FAST */ actualAction->top->rateSumPrev -= action->rate;
				/* SLOW */ //actualAction->top->rateSumPrev = actualAction->rateSumPrev + actualAction->rateSumNext + actualAction->rate;
			}
			else{
				// actualize next attributes
				actualAction->top->sizeNext--;
				/* FAST */ actualAction->top->rateSumNext -= action->rate;
				/* SLOW */ //actualAction->top->rateSumNext = actualAction->rateSumPrev + actualAction->rateSumNext + actualAction->rate;
			}
			actualAction = actualAction->top;
		}

		// rearange subtrees
		Process* shortSubtree, *longSubtree;
		double longRateSum;
		int    longSize;
		if( action->sizePrev <= action->sizeNext){
			shortSubtree = action->prev;
			longSubtree  = action->next;
			longRateSum	= action->rateSumNext;
			longSize	= action->sizeNext;
		}else{
			shortSubtree = action->next;
			longSubtree  = action->prev;
			longRateSum	= action->rateSumPrev;
			longSize	= action->sizePrev;
		}
		//fprintf( stderr, "longSubtree:%i, shortSubtree:%i\n",
		  //       (longSubtree!=NULL ? longSubtree->type : -1),
		    //     (shortSubtree!=NULL ? shortSubtree->type : -1)         );

		// actualize tree
		/* FAST */ //this->rateSum -= action->rate;
		//fprintf( stderr, "size = %i\n", this->size);
		if( --(this->size) == 0){
			// tree is empty
			this->root = NULL;
			/* SLOW */ //this->rateSum = 0.;
			//longSubtree->top = longSubtree;
		}else{
			/* SLOW */ //this->rateSum = this->root->rateSumPrev + this->root->rateSumNext + this->root->rate;

			if( shortSubtree == NULL){
				// only one subtree has to be added
				if( this->root == action)
					this->root = longSubtree;
				else{
					if( action == action->top->prev){
						// actualize prev attributes
						action->top->prev = longSubtree;
					}else{
						// actualize next attributes
						action->top->next = longSubtree;
					}
				}
				if( longSubtree != NULL)
					longSubtree->top = action->top;


			}else{
				// attach first subtree
				if( this->root == action)
					this->root = shortSubtree;
				else{
					if( action == action->top->prev){
						// actualize prev attributes
						action->top->prev = shortSubtree;
					}else{
						// actualize next attributes
						action->top->next = shortSubtree;
					}
				}
				shortSubtree->top = action->top;

				// attach seconde subtree
				Process* lastAction = actualAction = shortSubtree;

				while( actualAction!=NULL){
					lastAction = actualAction;
					if( actualAction->sizePrev <= actualAction->sizeNext){
						actualAction->sizePrev    += longSubtree->sizePrev    + longSubtree->sizeNext    + 1;
						actualAction->rateSumPrev += longSubtree->rateSumPrev + longSubtree->rateSumNext + longSubtree->rate;
						actualAction = actualAction->prev;
						if( actualAction==NULL)
							lastAction->prev = longSubtree;
					}else{
						actualAction->sizeNext    += longSubtree->sizePrev    + longSubtree->sizeNext    + 1;
						actualAction->rateSumNext += longSubtree->rateSumPrev + longSubtree->rateSumNext + longSubtree->rate;
						actualAction = actualAction->next;
						if( actualAction==NULL)
							lastAction->next = longSubtree;
					}
				}
				/*if( lastAction->sizePrev == 0){
					lastAction->prev = longSubtree;
				}else{
					lastAction->next = longSubtree;
				}*/
				longSubtree->top = lastAction;

				//if( this->root == action)
				//	this->root = shortSubtree;
			}

		}

		// detach from tree
		action->top = NULL;
		action->prev = NULL;
		action->next = NULL;
		action->sizePrev = 0;
		action->sizeNext = 0;
		action->rateSumPrev = 0.;
		action->rateSumNext = 0.;
	}
	else{
		//fprintf( stderr, "ERROR: Can't delete action from action tree. Action isn't element of action tree!\n");
		//exit( 0);
	}

}

#define RAND01 ((double)rand()/(RAND_MAX+1.))
void ProcessTree::performGillespieStep( RATE_T &time)
{
	//int i;
	double random, sum;

	// CHOOSE ACTION
	Process* p_elem = this->root;

	// randomly choose part of rate sum
	double temp_rand = RAND01;
	random = temp_rand * rate();//this->rateSum;

	//fprintf( stderr, "size = %i, rateSum = %lf, random = %lf\n", this->size, this->rateSum, random);

	// sum probabilities until randomly choosen part is reached
	sum = 0.;
	do{
		//fprintf( stderr, "sum = %lf, random = %lf -> p_elem->rateSumPrev = %lf, p_elem->rateSumNext = %lf, p_elem->rate = %lf => %lf\n", sum, random, p_elem->rateSumPrev, p_elem->rateSumNext, p_elem->rate, p_elem->rateSumPrev+p_elem->rateSumNext+p_elem->rate);
		if( p_elem->sizePrev && random <= p_elem->rateSumPrev + sum){
			// take left branch
			//fprintf( stderr, "-> take left branch\n");
			p_elem = p_elem->prev;
		}else{
			sum += p_elem->rateSumPrev;
			if( p_elem->sizeNext && random <= p_elem->rateSumNext + sum){
				// take right branch
				//fprintf( stderr, "-> take right branch\n");
				p_elem = p_elem->next;
			}else{
				sum += p_elem->rateSumNext;
				if( random <= p_elem->rate + sum){
					// take this action
					//fprintf( stderr, "-> take this action\n");
					// calculate passed time
					random = RAND01;
					time += -log(1 - random) / rate();//this->rateSum;
					//*time += 1. / this->rateSum;


					/*int its = 10000;
					double min = 0.000001;
					double max = 1.;
					double dist = max - min;
					double stepsize = dist / (double)(its-1);
					double av = 0.;
					for(int i=0; i<its; i++){
						av += -log(min+i*stepsize);
					}
					av /= (double)its;
					printf( "%lf %lf %lf %lf\n", log(0.), log(1.), exp(1.), av); exit( 0);
					*/
					p_elem->perform();
					return;

				}else{
					fprintf( stderr, "ERROR: Something is going wrong in selectAction()\n");
					fprintf(stderr, "INFO: actionList->sum_of_prob=%e, random=%lf (%lf), sum=%lf\n", this->rate(), random, temp_rand, sum);

					exit( 0);
				}
			}
		}
	}while( p_elem!=NULL);
	printf("Something is going wrong in selectAction()\n");
	printf("INFO: actionList->sum_of_prob=%lf, random=%lf (%lf), sum=%lf\n", this->rate(), random, temp_rand, sum);
	//printActionList( actionList);
	//time_temp = difftime(time(NULL),t);
	//Sim_time += time_temp;
	exit(0);

}

Process* ProcessTree::selectAction( RATE_T &random)
{
	RATE_T sum;

	// CHOOSE ACTION
	Process* p_elem = this->root;

	// sum rates until randomly choosen part is reached
	sum = 0.;
	do{
		//fprintf( stderr, "sum = %lf, random = %lf -> p_elem->rateSumPrev = %lf, p_elem->rateSumNext = %lf, p_elem->rate = %lf => %lf\n", sum, random, p_elem->rateSumPrev, p_elem->rateSumNext, p_elem->rate, p_elem->rateSumPrev+p_elem->rateSumNext+p_elem->rate);
		if( p_elem->sizePrev && random <= p_elem->rateSumPrev + sum){
			// take left branch
			//fprintf( stderr, "-> take left branch\n");
			p_elem = p_elem->prev;
		}else{
			sum += p_elem->rateSumPrev;
			if( p_elem->sizeNext && random <= p_elem->rateSumNext + sum){
				// take right branch
				//fprintf( stderr, "-> take right branch\n");
				p_elem = p_elem->next;
			}else{
				sum += p_elem->rateSumNext;
				if( random <= p_elem->rate + sum){
					return p_elem;

				}else{
					fprintf( stderr, "ERROR: Something is going wrong in selectAction()\n");
					fprintf(stderr, "INFO: actionList->sum_of_prob=%lf, random=%lf, sum=%lf\n", this->rate(), random, sum);

					exit( 0);
				}
			}
		}
	}while( p_elem!=NULL);
	printf("Something is going wrong in selectAction()\n");
	printf("INFO: actionList->sum_of_prob=%lf, random=%lf, sum=%lf\n", this->rate(), random, sum);
	//printActionList( actionList);
	//time_temp = difftime(time(NULL),t);
	//Sim_time += time_temp;
	exit(0);

}

void ProcessTree::updateRate( Process* action, RATE_T newRate)
{
	//fprintf( stderr, "ActionTree::actualizeRate( %s, %e)\n", actionTypeToString( action->type), newRate);
	/* FAST */ double rateDifference = newRate - action->rate;
	if( rateDifference==0)
		return;

	if(action->top == NULL)
	return;
	//fprintf(stderr, "actualizeRate(%s: %lf -> %lf)", actionTypeToString(action->type), action->rate, newRate);

	// actualize action
	action->rate = newRate;

	// actualize tree
	Process *actualAction = action,
	        *lastAction;


	while( actualAction != this->root){
		lastAction = actualAction;
		actualAction = actualAction->top;
		if( actualAction->prev == lastAction){
			/* SLOW */ //actualAction->rateSumPrev = lastAction->rate + lastAction->rateSumPrev + lastAction->rateSumNext;
			/* FAST */ actualAction->rateSumPrev += rateDifference;
		}else{
			if(actualAction->next == lastAction)
				/* SLOW */ //actualAction->rateSumNext = lastAction->rate + lastAction->rateSumPrev + lastAction->rateSumNext;
				/* FAST */ actualAction->rateSumNext += rateDifference;
			else{
				fprintf(stderr, "ERROR in ActionTree::actualizeRate()\n");
				exit( 0);
			}
		}
	}

	/* SLOW */ //this->rateSum = this->root->rate + this->root->rateSumPrev + this->root->rateSumNext;
	/* FAST */ //this->rateSum += rateDifference;
	//fprintf(stderr, "rateDifference = %lf, rateSum = %lf\n", rateDifference, this->rateSum);

}
