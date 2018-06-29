#ifndef CELLSPHERICAL_H
#define CELLSPHERICAL_H

#include "Cell.h"
#include "../Elements/ModelElementSphere.h"
#include "../BasicDatatypes/Vector.h"

#include <string>
#include <vector>

class CSModel;
class QXmlStreamWriter;
namespace H5 { class CompType; }

//! Defines the most basic cell class.
class CellSpherical : public Cell, public ModelElementSphere
{
public:

    #pragma region Radius

    //! Initial radius of the cell at t=0 and right after cell division
    double initialRadius;

    //! Radius at which the cell divides into two
    double divisionRadius;

    //! Radius by which the current radius is increased in each timestep
    double deltaRadius; // Prelim: public

    #pragma endregion

    #pragma region Biophysics

    // moved to Cell.h

    #pragma endregion

    double mQuiescentMin;
    double mQuiescentMax;

    #pragma region Dumbbell

    bool mUseDumbbell;
    bool mDaughterCell;
    Vector3f mDivisionAxis;
    double defaultDivisionDistance;
    double mDumbbellPeriod;
    double mDumbbellInitialRadius;
    double deltaRadiusDumbbell;
    double mElongationFactor;
    CellSpherical * mpDivisionPartner;
    // For writing and reading hdf5
    long int mDivisionPartnerIndex;
    CSModel *mpModel;

    #pragma endregion

	// Use for example in embedding tissues to distinguish between external and growing tissue
	// 0 ... Normal, growing cell
	// 1 ... Eternally quiscent cell
	int mSubType;

    // Flag for cells directly on a necrotic region around the central vein
    // (scenario Lobule).
    bool mLesionEdge;

    // Holding the value of
    //   2*d_lesion/(d_lesion+d_lobuleEdge),
    // where d_lesion is the distance to the lesion edge cell in the line of
    // sight from the cell to the central vein, and d_lobuleEdge is the distance
    // to the closest edge of the lobule.
    double mLayerIndicator;

    // parameters for the S-phase - time of start and end of S phase after division:
    double mSPhaseStart;
    double mSPhaseEnd;
    double mDetachTime;

public:

    //! Default constructor
    CellSpherical(double x=0, double y=0, double z=0);

    //! Updates the growth step of the radius for one timestep
    void UpdateDeltaRadius(double timeStep);

    //! Grows the radius of the cell for one timestep
    void Grow();

    //! Re-implemented to set s-phase timings.
    virtual void SetCycleTimeGaussClamped(double mean, double stddev);

    //! Determines whether a cell reached a sufficient size (currentRadius) for division
    bool CanDivide();

    void WakeUp();

	//! Sets the celltype of the cell
	void SetCellSubType(int newCellSubType);

    virtual ModelElement * Divide( CSModel * parentModel );

    virtual CellSpherical * Duplicate( double x, double y, double z );

    // function to be called to determine if the cell goes into quiescent state.
    virtual void RestrictionCheckPoint();

    // overwritten to implement SPhase state dependent on mSPhaseStart/End.
    bool getState( Cell::State );

    #pragma region xml

    //! Writes out an XML representation of this cell to the given QXmlStreamWriter
    //  kept here for mere debugging purposes - WILL BE REMOVED in the future
    void writeXML( QXmlStreamWriter * ) const;

    //! Adds the member variables to save into the H5::CompType
    static void HDF5DataFormat( H5::CompType & );

    //! Builds up the H5::CompType from input data
    static H5::CompType ParseHDF5DataFormat( H5::CompType & inputType,
                                             std::stringstream & errors,
                                             std::stringstream & warnings );

    #pragma endregion
};

#endif
