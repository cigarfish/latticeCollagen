///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CellSphericalPolar.cpp                                               //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2013-11-05 20:22:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "CellSphericalPolar.h"
#include "CellSpherical.h"
#include "../Model/CSModel.h"
// preliminary:
#include "../Model/ModelCellsSpherical/ModelCellsSpherical.h"

#include <H5Cpp.h>
#include "../../tools/dataIO/hdf5/CSHDF5Helpers.h"

// preliminary:  hardcode max adhesion energy
//  = 2 * M_PI * r^2 * adhesionDensity_default * singleBondEnergy_default / 2
// | with r=1, adhesionDensity_default = 1e13 * length_scale^2 = 1e3
// | singleBondEnergy_default = 1e-19 / energy_scale = 1e-3
// = M_PI
double CellSphericalPolar::adhesionEnergyMax = M_PI;

CellSphericalPolar::CellSphericalPolar(double x, double y, double z)
    : CellSpherical(x,y,z)
{
// Further init is done at AddCell in Model
  mType = ModelElement::TypeCellSphericalPolar;

  mpContactsPointers = new CSListContainer<ModelElement *>( 512 );
}


CellSphericalPolar::~CellSphericalPolar()
{
    delete mpContactsPointers;

    if (mpDivisionPartner)
    {
        mpDivisionPartner->mpDivisionPartner = NULL;
        delete mpDivisionPartner;
    }
}


ModelElement * CellSphericalPolar::Divide( CSModel * parentModel )
{
    // preliminary:  accessing model parameters by casting to ModelCellsSpherical
    // to be replaced by having parameters local to cells.
    ModelCellsSpherical * model = dynamic_cast<ModelCellsSpherical *>(parentModel);

    Vector3f newCellPosition = this->position;

    mDivisionAxis.Set( mPolarDirection );
    if ( model->is2D )
        mDivisionAxis.z = 0.;
    mDivisionAxis.Normalize();

    // preliminary:  random division
    // to be replaced by alignment to closest sinusoid.
    // if ( model->is2D )
    // {
    //     mpRandom->GetRandomUnitVector(&mDivisionAxis.x, &mDivisionAxis.y);
    //     mDivisionAxis.z = 0;
    // }
    // else
    // {

    if ( model->mpGraphBloodVesselNetwork )
    {
        ModelElementVesselSphere * nextVesselElement =
            static_cast<ModelElementVesselSphere *>
            ( model->cells2->NextNeighborByType( this, ModelElement::TypeVesselSphere ) );
        if ( !nextVesselElement )
            mpRandom->GetRandomUnitVector(&mDivisionAxis.x, &mDivisionAxis.y, &mDivisionAxis.z);
        else
        {
            if ( nextVesselElement->mvpSegments.size() == 1 )
            {
                CSGraphEdge *edge = nextVesselElement->mvpSegments[0];
                mDivisionAxis.Set( edge->mpStart->position.x - edge->mpEnd->position.x,
                                   edge->mpStart->position.y - edge->mpEnd->position.y,
                                   edge->mpStart->position.z - edge->mpEnd->position.z );
            }
            else
            {
                double minDistance = std::numeric_limits<double>::max();

                for ( auto edge: nextVesselElement->mvpSegments )
                {
                    Vector3f distVector( position.x - nextVesselElement->position.x,
                                         position.y - nextVesselElement->position.y,
                                         position.z - nextVesselElement->position.z );

                    Vector3f segmentVector( edge->mpStart->position.x - edge->mpEnd->position.x,
                                            edge->mpStart->position.y - edge->mpEnd->position.y,
                                            edge->mpStart->position.z - edge->mpEnd->position.z );

                    double direction = 1;
                    if ( edge->mpStart == nextVesselElement )
                        direction = -1;

                    double distDotSegment = direction * (distVector.x*segmentVector.x + distVector.y*segmentVector.y + distVector.z*segmentVector.z );
                    double distance = 0;
                    if ( distDotSegment <= 0 )
                        distance = distVector.Norm();
                    else
                        distance = std::sqrt(distVector.Norm()*distVector.Norm() - distDotSegment*distDotSegment);

                    if ( distance < minDistance )
                    {
                        mDivisionAxis = segmentVector;
                        minDistance  = distance;
                    }
                    else if ( distance == minDistance )
                    {
                        mDivisionAxis -= segmentVector;
                    }
                }
            }
        } // if ( nextVesselElement )
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

    CellSphericalPolar * newCell = Duplicate( newCellPosition.x,
                                              newCellPosition.y,
                                              newCellPosition.z );

    #pragma region Create daughter cell 2 (by modifying the existing cell)

    // Recaclulate cell cycle time
    SetCycleTimeGaussClamped( model->defaultCellCycleTime,
                              model->defaultCellCycleTimeStandardDeviation );
    newCell->SetCycleTimeGaussClamped( model->defaultCellCycleTime,
                                       model->defaultCellCycleTimeStandardDeviation );

    if ( mUseDumbbell )
    {
        newCell->mDaughterCell = true;

        // At what cell radius to go to dumbbell formation in the next cell cycle.
        mDumbbellInitialRadius = divisionRadius - (divisionRadius - initialRadius) * mDumbbellPeriod/cycleTime;
        newCell->mDumbbellInitialRadius = divisionRadius - (divisionRadius - initialRadius) * mDumbbellPeriod/newCell->cycleTime;

        // calculate shrink rate of the dumbbell partners and growth rate of
        // their distance from each other.
        mElongationFactor   = defaultDivisionDistance  / mDumbbellPeriod;

        deltaRadiusDumbbell = ( initialRadius - mRadius ) / mDumbbellPeriod;

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


CellSphericalPolar *
CellSphericalPolar::Duplicate(double x, double y, double z)
{
    CellSphericalPolar * newCell = new CellSphericalPolar(x,y,z);

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

    newCell->mProliferationControl = this->mProliferationControl;

    if (mUseDumbbell)
    {
        newCell->mDumbbellPeriod = mDumbbellPeriod;
        newCell->mDetachTime = this->mDetachTime;
    }

    newCell->mPolarDirection.Set( this->mPolarDirection );

    return newCell;
}


void
CellSphericalPolar::HDF5DataFormat( H5::CompType & typeDefinition )
{
    // CellSpherical definitions
    H5_DOUBLE_ARRAY( vectorType, 3 );

    typeDefinition.insertMember( "Position", HOFFSET(CellSphericalPolar, position), vectorType );
    typeDefinition.insertMember( "Radius", HOFFSET(CellSphericalPolar, mRadius), H5::PredType::NATIVE_DOUBLE );
    typeDefinition.insertMember( "Cell Cycle State", HOFFSET(CellSphericalPolar, cellcycleState), H5::PredType::NATIVE_B8 );
    typeDefinition.insertMember( "Cell Cycle Time", HOFFSET(CellSphericalPolar, cycleTime), H5::PredType::NATIVE_DOUBLE );
    typeDefinition.insertMember( "Dumbbell Division Partner Index", HOFFSET(CellSphericalPolar, mDivisionPartnerIndex), H5::PredType::NATIVE_LONG );
    typeDefinition.insertMember( "Cell Subtype", HOFFSET(CellSphericalPolar, mSubType), H5::PredType::NATIVE_INT );
    // end of CellSpherical definitions

    // CellSphericalPolar definitions
    typeDefinition.insertMember( "Polar Direction", HOFFSET(CellSphericalPolar, mPolarDirection), vectorType );
    typeDefinition.insertMember( "Pole Region Angle", HOFFSET(CellSphericalPolar, mAngle), H5::PredType::NATIVE_DOUBLE );
}


H5::CompType
CellSphericalPolar::ParseHDF5DataFormat( H5::CompType & inputTypeDefinition,
                                         std::stringstream & /*errors*/,
                                         std::stringstream & warnings )
{
    H5_DOUBLE_ARRAY( vectorType, 3 );

    int numMembers = inputTypeDefinition.getNmembers();

    H5::CompType typeDefinition( sizeof(CellSphericalPolar) );

    // build the data type:
    for ( int i=0; i<numMembers; ++i )
    {
        std::string fieldName = inputTypeDefinition.getMemberName(i);

        if ( fieldName == "Position" )
        {
            H5_DOUBLE_ARRAY(positionType, 3);
            typeDefinition.insertMember( "Position",
                                   HOFFSET(CellSphericalPolar, position),
                                   vectorType );
        }
        else if ( fieldName == "Radius" )
        {
            typeDefinition.insertMember("Radius",
                                  HOFFSET(CellSphericalPolar, mRadius),
                                  H5::PredType::NATIVE_DOUBLE);
        }
        else if ( fieldName == "Cell Cycle State" )
        {
            typeDefinition.insertMember("Cell Cycle State",
                                  HOFFSET(CellSphericalPolar, cellcycleState),
                                  H5::PredType::NATIVE_B8);
        }
        else if ( fieldName == "Cell Cycle Time" )
        {
            typeDefinition.insertMember( "Cell Cycle Time",
                                         HOFFSET(CellSphericalPolar, cycleTime),
                                         H5::PredType::NATIVE_DOUBLE );
        }
        else if ( fieldName == "Dumbbell Division Partner Index" )
        {
            typeDefinition.insertMember( "Dumbbell Division Partner Index",
                                         HOFFSET(CellSphericalPolar, mDivisionPartnerIndex),
                                         H5::PredType::NATIVE_LONG );
        }
        else if ( fieldName == "Cell Subtype" )
        {
            typeDefinition.insertMember( "Cell Subtype",
                                         HOFFSET(CellSphericalPolar, mSubType),
                                         H5::PredType::NATIVE_LONG );
        }
        else if ( fieldName == "Polar Direction" )
        {
            typeDefinition.insertMember( "Polar Direction",
                                         HOFFSET(CellSphericalPolar, mPolarDirection),
                                         vectorType );
        }
        else if ( fieldName == "Pole Region Angle" )
        {
            typeDefinition.insertMember( "Pole Region Angle",
                                         HOFFSET(CellSphericalPolar, mAngle),
                                         H5::PredType::NATIVE_DOUBLE );
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
