#ifndef CELL_H
#define CELL_H

#define CELLCYCLE_STATE_QUIESCENT 0

#include "../BasicDatatypes/Vector.h"

class CSModel;
class ModelElement;
class Random;
class QXmlStreamWriter;


//! Defines the most basic cell class.
class Cell
{
public:

    enum State {
        StateQuiescent = 1,
        StateDividing  = 2,
        StateSPhase    = 4,
        // StateOther  = 8,
        // StateNext   = 16, etc in powers of two
    };

    bool getState( State state ) const { return (cellcycleState & state); };
    void setState( State state, bool onOff=true )
    { (onOff) ? cellcycleState |= state : cellcycleState &= (~state); };

    // proliferation control

    // since there are different models for the transition from the G1 to the G0
    // phase, let's introduce a flag and name the paradigms with an enum:
    enum ProliferationControl {
        ProliferationControlNone = 0, // quiescent for all times
        ProliferationControlPressure, // pressure based G1->G0->G1 transitions
        ProliferationControlCustom   // proliferation is triggered by direct calls to WakeUp()
    };

    ProliferationControl mProliferationControl;

    //! Implementation of the Cell's division process.
    /*!   \param parentModel The model to which the cell belongs, provided for
                             accessing e.g. the random number generator.
          \returns The new cell.
    */
    virtual ModelElement * Divide( CSModel * parentModel ) =0;

#pragma region Lifecycle properties

#pragma region Cycletime

	//! Cell cycle time in dimensionless units
    double cycleTime;

	virtual void SetCycleTimeGaussClamped(double mean, double stddev);

#pragma endregion

	//! Generation in relation to the population
	//! 0: Initial (eden) cell at t=0, n: Cell is daughter cell of a cell of generation n-1
	unsigned int mGeneration;

	//! Generic variable to store cell cycle state
	//! In the inherited class 'CellSpherical' cellCycleState is used to store whether a cell is quiescent cellCycleState[0] = true or not cellCycleState[0] = false 

	// Model1D:
	// [0] .. cell is quiescent (true) or not (false) // CELLCYCLE_STATE_QUIESCENT
    unsigned char cellcycleState;

	//! Cell volume
    //double volume;

    //! Pointer to the Model's Random instance
    //  This has to be set explicitly by the Model::AddCell() method
    Random * mpRandom;

    //! A vector for retaining the Langevin contribution to the total force.
    //  Needed for rescaling the Langevin contribution when choosing a different
    //  simulation time step size.
    Vector3f mLangevinForce;

#pragma endregion

    //! Default constructor
    Cell();

    //! Writes out an XML representation of this cell to the given QXmlStreamWriter
    virtual void writeXML( QXmlStreamWriter * ) const =0;
};

#endif
