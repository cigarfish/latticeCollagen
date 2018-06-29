////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSInteractionFrictionMatrixHertz.cpp                          //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-07-18 23:55:10                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CSInteractionFrictionMatrix.h"

#include "../Cell/CellSpherical.h"
#include "../Elements/ModelElementSphere.h"
#include "../Elements/ModelElementVesselSphere.h"



CSInteractionFrictionMatrix::CSInteractionFrictionMatrix()
    : mpDistance(NULL),
      mpContactArea(NULL),
      mpGammaCellsParallel(NULL),
      mpGammaCellsPerpendicular(NULL),
      mpFrictionMatrices(NULL),
      mFrictionMatrixCounter(0),
      mFrictionMatricesAllocationSize(0)
{}


void
CSInteractionFrictionMatrix::interact( ModelElementSphere *cell_1, ModelElementSphere *cell_2 )
{
    // increment counter for contacts
    ++cell_1->mNumContacts;
    ++cell_2->mNumContacts;

    ++mFrictionMatrixCounter;

    reallocate();

    double distanceUnitVector[3];

    double inverseDistance = 1/(*mpDistance);

    distanceUnitVector[0] = inverseDistance * (cell_2->position.x - cell_1->position.x);
    distanceUnitVector[1] = inverseDistance * (cell_2->position.y - cell_1->position.y);
    distanceUnitVector[2] = inverseDistance * (cell_2->position.z - cell_1->position.z);

    // The surface friction matrix is a matrix that combines the
    // surface friction Term with the projection onto the plane
    // perpendicular to the distance vector.  The projection has to
    // be used to determine the contribution of the parallel motion
    // in the expression
    //
    // Sum_j ( gamma_cell-cell * sigma_ij * (v_i - v_j)_parallel )
    //
    // where
    //   gamma_cell-cell is the friction coefficient for cell-cell
    //                   interaction,
    //   sigma_ij        is the contact area derived by means of the
    //                   Hertz or the JKR model,
    //   v_i, v_j        are the velocity _vectors_ of cell i and j
    //                   respectively.
    //
    mpFrictionMatrix = mpFrictionMatrices + 9 * mFrictionMatrixCounter;

    double surfaceTermPerpendicular = *mpContactArea * (*mpGammaCellsPerpendicular);
    double surfaceTermCombined = *mpContactArea * ( *mpGammaCellsParallel - *mpGammaCellsPerpendicular );


    // diagonal elements:
    *mpFrictionMatrix       = surfaceTermPerpendicular + surfaceTermCombined * distanceUnitVector[0]*distanceUnitVector[0];
    *(mpFrictionMatrix + 4) = surfaceTermPerpendicular + surfaceTermCombined * distanceUnitVector[1]*distanceUnitVector[1];
    *(mpFrictionMatrix + 8) = surfaceTermPerpendicular + surfaceTermCombined * distanceUnitVector[2]*distanceUnitVector[2];

    // off-diagonal (symmetric)
    *(mpFrictionMatrix + 1) = surfaceTermCombined * distanceUnitVector[0]*distanceUnitVector[1];
    *(mpFrictionMatrix + 3) = *(mpFrictionMatrix + 1);

    *(mpFrictionMatrix + 2) = surfaceTermCombined * distanceUnitVector[0]*distanceUnitVector[2];
    *(mpFrictionMatrix + 6) = *(mpFrictionMatrix + 2);

    *(mpFrictionMatrix + 5) = surfaceTermCombined * distanceUnitVector[1]*distanceUnitVector[2];
    *(mpFrictionMatrix + 7) = *(mpFrictionMatrix + 5);

    cell_1->frictionMatrixEntryNumbers.push_back( mFrictionMatrixCounter );
    cell_2->frictionMatrixEntryNumbers.push_back( mFrictionMatrixCounter );

    cell_1->mContacts.push_back( cell_2->mGlobalIndex );
    cell_2->mContacts.push_back( cell_1->mGlobalIndex );
}


void
CSInteractionFrictionMatrix::interact( ModelElementSphere * sphere, ModelElementBarrierTriangle * barrier )
{
    // increment counter for contacts
    ++sphere->mNumContacts;

    ++mFrictionMatrixCounter;

    reallocate();

    double * distanceUnitVector = barrier->getNormalVector();

    // The surface friction matrix is a matrix that combines the
    // surface friction Term with the projection onto the plane
    // perpendicular to the distance vector.  The projection has to
    // be used to determine the contribution of the parallel motion
    // in the expression
    //
    // Sum_j ( gamma_cell-cell * sigma_ij * (v_i - v_j)_parallel )
    //
    // where
    //   gamma_cell-cell is the friction coefficient for cell-cell
    //                   interaction,
    //   sigma_ij        is the contact area derived by means of the
    //                   Hertz or the JKR model,
    //   v_i, v_j        are the velocity _vectors_ of cell i and j
    //                   respectively.
    //
    mpFrictionMatrix = mpFrictionMatrices + 9 * mFrictionMatrixCounter;

    double surfaceTermPerpendicular = *mpContactArea * (*mpGammaCellsPerpendicular);
    double surfaceTermCombined = *mpContactArea * ( *mpGammaCellsParallel - *mpGammaCellsPerpendicular );


    // diagonal elements:
    *mpFrictionMatrix       = surfaceTermPerpendicular + surfaceTermCombined * distanceUnitVector[0]*distanceUnitVector[0];
    *(mpFrictionMatrix + 4) = surfaceTermPerpendicular + surfaceTermCombined * distanceUnitVector[1]*distanceUnitVector[1];
    *(mpFrictionMatrix + 8) = surfaceTermPerpendicular + surfaceTermCombined * distanceUnitVector[2]*distanceUnitVector[2];

    // off-diagonal (symmetric)
    *(mpFrictionMatrix + 1) = surfaceTermCombined * distanceUnitVector[0]*distanceUnitVector[1];
    *(mpFrictionMatrix + 3) = *(mpFrictionMatrix + 1);

    *(mpFrictionMatrix + 2) = surfaceTermCombined * distanceUnitVector[0]*distanceUnitVector[2];
    *(mpFrictionMatrix + 6) = *(mpFrictionMatrix + 2);

    *(mpFrictionMatrix + 5) = surfaceTermCombined * distanceUnitVector[1]*distanceUnitVector[2];
    *(mpFrictionMatrix + 7) = *(mpFrictionMatrix + 5);

    sphere->frictionMatrixEntryNumbers.push_back( mFrictionMatrixCounter );
    sphere->mContacts.push_back( barrier->mGlobalIndex );
}


void
CSInteractionFrictionMatrix::interact( CellSpherical *cell_1, CellSpherical *cell_2 )
{
    ModelElementSphere *sphere1 = static_cast<ModelElementSphere *>(cell_1);
    ModelElementSphere *sphere2 = static_cast<ModelElementSphere *>(cell_2);

    interact( sphere1, sphere2 );
}


void
CSInteractionFrictionMatrix::interact( CellSpherical *cell_1, CellSphericalPolar *cell_2 )
{
    ModelElementSphere *sphere1 = static_cast<ModelElementSphere *>(cell_1);
    ModelElementSphere *sphere2 = static_cast<ModelElementSphere *>(cell_2);

    interact( sphere1, sphere2 );
}


void
CSInteractionFrictionMatrix::interact( CellSpherical * cell, ModelElementBarrierTriangle * barrier )
{
    ModelElementSphere * sphere = static_cast<ModelElementSphere *>(cell);

    interact( sphere, barrier );
}


void
CSInteractionFrictionMatrix::interact( CellSpherical *cell_1, ModelElementVesselSphere *cell_2 )
{
    ModelElementSphere *sphere1 = static_cast<ModelElementSphere *>(cell_1);
    ModelElementSphere *sphere2 = static_cast<ModelElementSphere *>(cell_2);

    interact( sphere1, sphere2 );
}


void
CSInteractionFrictionMatrix::interact( CellSphericalPolar *cell_1, CellSphericalPolar *cell_2 )
{
    ModelElementSphere *sphere1 = static_cast<ModelElementSphere *>(cell_1);
    ModelElementSphere *sphere2 = static_cast<ModelElementSphere *>(cell_2);

    interact( sphere1, sphere2 );
}


void
CSInteractionFrictionMatrix::interact( CellSphericalPolar *cell_1, ModelElementBarrierTriangle *barrier)
{
    ModelElementSphere *sphere = static_cast<ModelElementSphere *>(cell_1);

    interact( sphere, barrier );
}


void
CSInteractionFrictionMatrix::interact( CellSphericalPolar *cell_1, ModelElementVesselSphere *cell_2)
{
    ModelElementSphere *sphere1 = static_cast<ModelElementSphere *>(cell_1);
    ModelElementSphere *sphere2 = static_cast<ModelElementSphere *>(cell_2);

    interact( sphere1, sphere2 );
}


void
CSInteractionFrictionMatrix::interact( ModelElementVesselSphere * vsphere, ModelElementBarrierTriangle * barrier )
{
    ModelElementSphere * sphere = static_cast<ModelElementSphere *>(vsphere);

    interact( sphere, barrier );
}


void
CSInteractionFrictionMatrix::interact( CellSpherical * /*cell*/, ModelElementHollowSphere * /*capsule*/ )
{
    // ++mFrictionMatrixCounter;

    // reallocate();

    // Vector3f normalVector = cell->position;

    // normalVector.Normalize();

    // double surfaceTermPerpendicular = *mpContactArea * (*mpGammaCellsPerpendicular);
    // double surfaceTermCombined      = *mpContactArea * (*mpGammaCellsParallel - *mpGammaCellsPerpendicular);


    // mpFrictionMatrix = mpFrictionMatrices + 9 * mFrictionMatrixCounter;

    // // diagonal elements:
    // *mpFrictionMatrix       = surfaceTermPerpendicular + surfaceTermCombined * normalVector.x*normalVector.x;
    // *(mpFrictionMatrix + 4) = surfaceTermPerpendicular + surfaceTermCombined * normalVector.y*normalVector.y;
    // *(mpFrictionMatrix + 8) = surfaceTermPerpendicular + surfaceTermCombined * normalVector.z*normalVector.z;

    // // off-diagonal (symmetric)
    // *(mpFrictionMatrix + 1) = surfaceTermCombined * normalVector.x*normalVector.y;
    // *(mpFrictionMatrix + 3) = *(mpFrictionMatrix + 1);

    // *(mpFrictionMatrix + 2) = surfaceTermCombined * normalVector.x*normalVector.z;
    // *(mpFrictionMatrix + 6) = *(mpFrictionMatrix + 2);

    // *(mpFrictionMatrix + 5) = surfaceTermCombined * normalVector.y*normalVector.z;
    // *(mpFrictionMatrix + 7) = *(mpFrictionMatrix + 5);

}

// added by Jieling
void
CSInteractionFrictionMatrix::interact( ModelElementLatticeNode * node, CellSpherical * cell)
{
	//ModelElementSphere * sphere1 = static_cast<ModelElementSphere *>(node);
	//ModelElementSphere * sphere2 = static_cast<ModelElementSphere *>(cell);

	//interact( sphere1, sphere2 );
}

void 
CSInteractionFrictionMatrix::interact( ModelElementLatticeNode * node, CellSphericalPolar * cell)
{
	//ModelElementSphere * sphere1 = static_cast<ModelElementSphere *>(node);
	//ModelElementSphere * sphere2 = static_cast<ModelElementSphere *>(cell);

	//interact(sphere1, sphere2);
}

void
CSInteractionFrictionMatrix::interact( ModelElementLatticeNode * node, ModelElementBarrierTriangle * barrier)
{
	//ModelElementSphere * sphere = static_cast<ModelElementSphere *>(node);

	//interact(sphere, barrier);
}

void
CSInteractionFrictionMatrix::interact( ModelElementLatticeNode * node, ModelElementVesselSphere * sphere)
{
	//ModelElementSphere * sphere1 = static_cast<ModelElementSphere *>(node);
	//ModelElementSphere * sphere2 = static_cast<ModelElementSphere *>(sphere);

	//interact(sphere1, sphere2);
}

void
CSInteractionFrictionMatrix::interact( CellSpherical * cell, LatticeSpring * spring)
{

}

void
CSInteractionFrictionMatrix::interact( CellSphericalPolar * cell, LatticeSpring * spring)
{

}