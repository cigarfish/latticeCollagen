
#include "../../Core.h"
#include "CellSpherical.h"

#include "../Model/CSModel.h"
// preliminary:
#include "../Model/ModelCellsSpherical/ModelCellsSpherical.h"

#include "../../gui/GLTools/CSGLSphere.h"
#include "../../tools/Tools.h"

#include <sstream>

#include <H5Cpp.h>
#include "../../tools/dataIO/hdf5/CSHDF5Helpers.h"

#include <QtCore>


CellSpherical::CellSpherical(double x, double y, double z)
    : Cell(),
      ModelElementSphere(x, y, z),
      mUseDumbbell(false),
      mDaughterCell(false),
      mDumbbellPeriod(1800),
      mpDivisionPartner(NULL),
      mLesionEdge(false)
{
// Further init is done at AddCell in Model
  mType = ModelElement::TypeCellSpherical;

  // Default: Type 0 (Normal, growing cell)
  mSubType = 0;

  mProliferationControl = Cell::ProliferationControlPressure;
  this->lastPressure = 0;
}


void CellSpherical::Grow()
{
	// No growth of division for cells of type 1
	if (mSubType == 1) return;

    if ( getState( StateQuiescent ) )
    {
        RestrictionCheckPoint();
        if ( this->getState(Cell::StateQuiescent) )
            return;
    }

    if ( this->getState(Cell::StateDividing) )
    {
        // everything is done by the mother cell!
        if ( mDaughterCell )
            return;

        if ( mpDivisionPartner == NULL )
        {
            // we are just one loop after detach time
            setState(Cell::StateDividing, false);

            RestrictionCheckPoint();

            if (!getState(Cell::StateQuiescent))
            {
                UpdateDeltaRadius(mpModel->timeStep);
                mRadius += deltaRadius;
            }

            return;
        }

        // shrink the cell (using the current timestep(!)).
        mRadius += deltaRadiusDumbbell * mpModel->timeStep;
        mpDivisionPartner->mRadius += deltaRadiusDumbbell * mpModel->timeStep;
        if ( mRadius < initialRadius )
        {
            mRadius = initialRadius;
            mpDivisionPartner->mRadius = initialRadius;
        }

        // get the distance vector with the division partner
        Vector3f distanceVector = this->position - mpDivisionPartner->position;

        double distance = distanceVector.Normalize();
        if (distance==0.)
            distanceVector = mDivisionAxis;

        // advance the distance by the elongation factor times timeStep;
        double elongation = mElongationFactor * mpModel->timeStep;
        distanceVector *= elongation;

        position += distanceVector;
        mpDivisionPartner->position -= distanceVector;

        // now the distance has actually grown by 2*elongation, compare this to
        // 2*defaultDivisionDistance:
        if ( distance*.5 + elongation >= defaultDivisionDistance )
        {
            mRadius = initialRadius;
            mpDivisionPartner->mRadius = initialRadius;

            this->mDetachTime = mpModel->time;
            mpDivisionPartner->mDetachTime = this->mDetachTime;

            // decouple the two division partners
            mpDivisionPartner->mDaughterCell = false;
            mpDivisionPartner->mpDivisionPartner = NULL;
            mpDivisionPartner = NULL;
        }
    }
    else
    {
        #pragma region Debug output

        if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << "Cell " << this->mGlobalIndex << " Grow(). Radius: " << this->mRadius << " -> ";

        #pragma endregion

        UpdateDeltaRadius(mpModel->timeStep);
        mRadius += deltaRadius;

        #pragma region Debug output

        if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << this->mRadius << "\n";

        #pragma endregion
    }

}


void CellSpherical::SetCycleTimeGaussClamped( double mean, double stddev )
{
    Cell::SetCycleTimeGaussClamped(mean, stddev);

    // For a cell with a cycletime of exactly 24h, the S-phase would start exactly
    // 10 hours after entering the cell cycle, which is at 10/24 of the full cell
    // cycle.  It would end after 18h after entering the cell cycle, which is 3/4
    // of the full cell cycle.
    mSPhaseStart = cycleTime * 5/12;
    mSPhaseEnd = cycleTime * 3/4;

}


bool CellSpherical::CanDivide()
{
     if ( getState(Cell::StateDividing) )
        return false;

    if (mRadius >= divisionRadius)
        return true;

    return false;
}


void CellSpherical::WakeUp()
{
    if ( !getState(Cell::StateQuiescent) )
        return;

    setState( Cell::StateQuiescent, false );
    mDetachTime = mpModel->time;
}


void CellSpherical::SetCellSubType(int newCellSubType)
{
	mSubType = newCellSubType;

	if (mSubType == 1)
	{
        // Celltype 1 cells are always quiescent
        setState( Cell::StateQuiescent, true);
	}
}


ModelElement * CellSpherical::Divide( CSModel * parentModel )
{
    // preliminary:  accessing model parameters by casting to ModelCellsSpherical
    // to be replaced by having parameters local to cells.
    ModelCellsSpherical * model = static_cast<ModelCellsSpherical *>(parentModel);

    Vector3f newCellPosition = this->position;
    if ( model->is2D )
    {
        mpRandom->GetRandomUnitVector(&mDivisionAxis.x, &mDivisionAxis.y);
        mDivisionAxis.z = 0;
    }
    else
    {
      if( parentModel->dimension == 1 ){
        mDivisionAxis.x = 1;
        mDivisionAxis.y = 0;
        mDivisionAxis.z = 0;
      }else
        mpRandom->GetRandomUnitVector(&mDivisionAxis.x, &mDivisionAxis.y, &mDivisionAxis.z);
    }

    if (!mUseDumbbell)
    {
        mDivisionAxis.x *= model->defaultDivisionDistance;
        mDivisionAxis.y *= model->defaultDivisionDistance;
        if (!model->is2D)
            mDivisionAxis.z *= model->defaultDivisionDistance;
        newCellPosition += mDivisionAxis;
        mRadius = initialRadius;
    }

    CellSpherical * newCell = Duplicate( newCellPosition.x,
                                         newCellPosition.y,
                                         newCellPosition.z );

    #pragma region Create daughter cell 2 (by modifying the existing cell)

    SetCycleTimeGaussClamped( model->defaultCellCycleTime,
                              model->defaultCellCycleTimeStandardDeviation );
    newCell->SetCycleTimeGaussClamped( model->defaultCellCycleTime,
                                       model->defaultCellCycleTimeStandardDeviation );

    if ( mUseDumbbell )
    {
        newCell->mDaughterCell = true;

        // The cell grows up to divisionRadius before it goes into dumbbell
        // division (next cell cycle!):
        cycleTime -= mDumbbellPeriod;
        newCell->cycleTime -= mDumbbellPeriod;

        // calculate shrink rate of the dumbbell partners and growth rate of
        // their distance from each other.
        mElongationFactor   = defaultDivisionDistance  / mDumbbellPeriod;

        deltaRadiusDumbbell = ( initialRadius - divisionRadius ) / mDumbbellPeriod;

        newCell->mElongationFactor = mElongationFactor;
        newCell->deltaRadiusDumbbell = deltaRadiusDumbbell;

        // couple this cell with the newCell:
        // ToDo:  use indices so that saving dumbbell cells is possible
        this->mpDivisionPartner = newCell;
        newCell->mpDivisionPartner = this;

        // switch on the 'dividing' state
        this->setState(Cell::StateDividing);
        newCell->setState(Cell::StateDividing);

        newCell->mDivisionAxis = this->mDivisionAxis;

        newCell->mUseDumbbell = true;
    }
    else
    {
        // Update position
        this->position.x -= mDivisionAxis.x;
        this->position.y -= mDivisionAxis.y;
        if ( ! model->is2D )
            this->position.z -= mDivisionAxis.z;

        // Update corresponding growth step
        UpdateDeltaRadius( mpModel->timeStep );
        newCell->UpdateDeltaRadius( mpModel->timeStep );

        #pragma region Prelim: State transition (Growing -> Quiescent)

        RestrictionCheckPoint();
        newCell->RestrictionCheckPoint();

        #pragma endregion

        if ( !this->getState(StateQuiescent) )
        {
            this->mDetachTime = mpModel->time;
            newCell->mDetachTime = mpModel->time;
        }

    }

    #pragma endregion

    return (ModelElement *) newCell;
}


CellSpherical *
CellSpherical::Duplicate(double x, double y, double z)
{
    CellSpherical * newCell = new CellSpherical(x,y,z);

    newCell->mpModel  = this->mpModel;
    newCell->mpRandom = this->mpRandom;

    newCell->mUseDumbbell = mUseDumbbell;

    newCell->initialRadius = this->initialRadius;
    newCell->divisionRadius = this->divisionRadius;
    newCell->mRadius = this->mRadius;

    newCell->youngModulus = this->youngModulus;
    newCell->poissonRatio = this->poissonRatio;

    newCell->defaultDivisionDistance = this->defaultDivisionDistance;

    newCell->mQuiescentMin = this->mQuiescentMin;
    newCell->mQuiescentMax = this->mQuiescentMax;
    
    newCell->lastPressure = this->lastPressure;

    newCell->mProliferationControl = this->mProliferationControl;

    if (mUseDumbbell)
    {
        newCell->mDumbbellPeriod = mDumbbellPeriod;
        newCell->mDetachTime = this->mDetachTime;
    }

    return newCell;
}


void CellSpherical::UpdateDeltaRadius(double timeStep)
{
	deltaRadius = (timeStep / cycleTime) * (divisionRadius - initialRadius); // Recalculate radius growth step

	#pragma region Debug output

    if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << "CellSpherical::UpdateDeltaRadius(). Set new radius-delta to " << deltaRadius << " \n";

    #pragma endregion
}


void CellSpherical::writeXML(QXmlStreamWriter * xml) const
{
    xml->writeStartElement("Cell");
    xml->writeAttribute( "type", "CellSpherical" );

    xml->writeAttribute( "x", QString("%1").arg( position.x ) );
    xml->writeAttribute( "y", QString("%1").arg( position.y ) );
    xml->writeAttribute( "z", QString("%1").arg( position.z ) );

    xml->writeAttribute( "radius", QString("%1").arg( mRadius ) );

    char state[9] = "00000000";
    for ( int i = 0; i < 8; ++i )
        if ( cellcycleState & 1<<i )
          state[i] = '1';
    xml->writeAttribute( "cellcycleState", QString(state) );

    xml->writeEndElement();
}


void
CellSpherical::HDF5DataFormat( H5::CompType & typeDefinition )
{
    H5_DOUBLE_ARRAY( positionType, 3 );
//    H5_DOUBLE_ARRAY( velocitiesType, 3 );

    typeDefinition.insertMember( "Position", HOFFSET(CellSpherical, position), positionType );
//    typeDefinition.insertMember( "Velocities", HOFFSET(CellSpherical, mVelocities), positionType );
    typeDefinition.insertMember( "Radius", HOFFSET(CellSpherical, mRadius), H5::PredType::NATIVE_DOUBLE );
    typeDefinition.insertMember( "Cell Cycle State", HOFFSET(CellSpherical, cellcycleState), H5::PredType::NATIVE_B8 );
    typeDefinition.insertMember( "Cell Cycle Time", HOFFSET(CellSpherical, cycleTime), H5::PredType::NATIVE_DOUBLE );
    typeDefinition.insertMember( "Dumbbell Division Partner Index", HOFFSET(CellSpherical, mDivisionPartnerIndex), H5::PredType::NATIVE_LONG );
    typeDefinition.insertMember( "Cell Subtype", HOFFSET(CellSpherical, mSubType), H5::PredType::NATIVE_INT );
}


H5::CompType
CellSpherical::ParseHDF5DataFormat( H5::CompType & inputTypeDefinition,
                                    std::stringstream & /*errors*/,
                                    std::stringstream & warnings )
{
    int numMembers = inputTypeDefinition.getNmembers();

    H5::CompType typeDefinition( sizeof(CellSpherical) );

    // build the data type:
    for ( int i=0; i<numMembers; ++i )
    {
        std::string fieldName = inputTypeDefinition.getMemberName(i);

        if ( fieldName == "Position" )
        {
            H5_DOUBLE_ARRAY(positionType, 3);
            typeDefinition.insertMember( "Position",
                                   HOFFSET(CellSpherical, position),
                                   positionType );
        }
        else if ( fieldName == "Radius" )
        {
            typeDefinition.insertMember("Radius",
                                  HOFFSET(CellSpherical, mRadius),
                                  H5::PredType::NATIVE_DOUBLE);
        }
        else if ( fieldName == "Cell Cycle State" )
        {
            typeDefinition.insertMember("Cell Cycle State",
                                  HOFFSET(CellSpherical, cellcycleState),
                                  H5::PredType::NATIVE_B8);
        }
        else if ( fieldName == "Cell Cycle Time" )
        {
            typeDefinition.insertMember( "Cell Cycle Time",
                                         HOFFSET(CellSpherical, cycleTime),
                                         H5::PredType::NATIVE_DOUBLE );
        }
        else if ( fieldName == "Dumbbell Division Partner Index" )
        {
            typeDefinition.insertMember( "Dumbbell Division Partner Index",
                                         HOFFSET(CellSpherical, mDivisionPartnerIndex),
                                         H5::PredType::NATIVE_LONG );
        }
        else if ( fieldName == "Cell Subtype" )
        {
            typeDefinition.insertMember( "Cell Subtype",
                                         HOFFSET(CellSpherical, mSubType),
                                         H5::PredType::NATIVE_LONG );
        }
        else if ( fieldName == "Ignore" )
        {}
        else
        {
            warnings << "CellSpherical::ParseHDF5DataFormat:  "
                     << "Unknown/Ignored field in HDF5 data set:\n"
                     << "\t\"" << fieldName << "\"\n";
        }
    }

    return typeDefinition;
}


void
CellSpherical::RestrictionCheckPoint()
{
    // The state transition could be depending on the type or different
    // conditions such as pressure, internal metabolic state, etc.

    switch ( mProliferationControl )
    {
    // pressure dependent transition:
    case ProliferationControlPressure:
        if ( getState(StateQuiescent) )
            setState( Cell::StateQuiescent, (lastPressure>mQuiescentMin) );
        else
            setState( Cell::StateQuiescent, (lastPressure>mQuiescentMax) );
        break;

    default:
    case ProliferationControlNone:
    case ProliferationControlCustom: // only trigger proliferation by external WakeUp() call
        setState( Cell::StateQuiescent );
        break;
    }
}


bool
CellSpherical::getState( Cell::State state )
{
    if (state == Cell::StateSPhase )
    {
        bool myState = false;
        if ( mpModel->time - mDetachTime > mSPhaseStart
             &&
             mpModel->time - mDetachTime < mSPhaseEnd )
            myState = true;

        this->setState(state, myState);

        return myState;
    }
    // else:
    return Cell::getState(state);
}
