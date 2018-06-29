////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSInteractionHertz.cpp                                        //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-07-18 17:34:24                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CSInteractionHertz.h"

#include "../Cell/CellSpherical.h"
#include "../Cell/CellSphericalPolar.h"
#include "../Elements/ModelElementBarrierTriangle.h"
#include "../Elements/ModelElementHollowSphere.h"

#include "../../tools/model/CSModelTools.h"


#include <cmath>
#include <iostream>

CSInteractionHertz::CSInteractionHertz()
    : CSInteraction(),
      mForceVector(0.,0.,0.),
      mDistance(0.),
      mContactArea(0.),
      mForce(0.),
      mIs2D(false)
{}


void
CSInteractionHertz::interact( ModelElementSphere *cell_1, ModelElementSphere *cell_2 )
{
    mContactArea =0;

    if( cell_1->mStatic && cell_2->mStatic)
        return;

    double force_ij_absolute;

    if (mIs2D)
        mDistance = CSModelTools::GetDistance2D( cell_2->position, cell_1->position );
    else
        mDistance = CSModelTools::GetDistance3D( cell_2->position, cell_1->position );

    if ( ( cell_1->mRadius + cell_2->mRadius ) <= mDistance )
        return;

    // contact area needed for both Hertz force and Cell-Cell friction:
    mContactArea = CSModelTools::GetContactAreaHertz( cell_1->mRadius,
                                                    cell_2->mRadius,
                                                    mDistance);

    cell_1->surfaceContactArea += mContactArea;
    cell_2->surfaceContactArea += mContactArea;

    // Calculate value of (Hertz) force
    mForce = CSModelTools::GetHertzForce( cell_1->poissonRatio, cell_1->youngModulus, cell_1->mRadius,
                                        cell_2->poissonRatio, cell_2->youngModulus, cell_2->mRadius,
                                        mDistance, *mpSingleBondEnergy, *mpAdhesionDensity);


    force_ij_absolute = fabs(mForce);

    // Store cumulated forces
    cell_1->lastForce += mForce;
    cell_2->lastForce += mForce;

    // Store cumulated absolute forces
    cell_1->accumulatedForceAbsolute += force_ij_absolute;
    cell_2->accumulatedForceAbsolute += force_ij_absolute;

    cell_1->accumulatedPressure += force_ij_absolute/mContactArea;
    cell_2->accumulatedPressure += force_ij_absolute/mContactArea;

    if (mIs2D)
    {
        mForceVector.Set(cell_2->position.x - cell_1->position.x,
                         cell_2->position.y - cell_1->position.y,
                         0);
    }
    else
    {
        mForceVector.Set(cell_2->position.x - cell_1->position.x,
                         cell_2->position.y - cell_1->position.y,
                         cell_2->position.z - cell_1->position.z);
    }

    mForceVector.Normalize();
    mForceVector.Multiply(mForce);

    //add directed force
    cell_1->directedForce.Add( -mForceVector.x, -mForceVector.y, -mForceVector.z);
    cell_2->directedForce.Add(  mForceVector.x,  mForceVector.y,  mForceVector.z);
}


void
CSInteractionHertz::interact( ModelElementSphere *cell, ModelElementBarrierTriangle *barrier )
{
  mContactArea =0;

  double force_ij_absolute;

  mDistance = CSModelTools::GetDistance3D( cell->position, barrier->getNormalVector(), barrier->getD() );

//  cell->directedForce.Add( mDistance*barrier->getNormalVector(0) , mDistance*barrier->getNormalVector(1) , mDistance*barrier->getNormalVector(2) );

  if ( cell->mRadius <= mDistance )
    return;

  // if ( useHertz )
  // contact area needed for both Hertz force and Cell-Cell friction:
  mContactArea = CSModelTools::GetContactAreaHertz( cell->mRadius,
                                                  mDistance);
  cell->surfaceContactArea += mContactArea;

  // Calculate value of (Hertz) force
  mForce = CSModelTools::GetHertzForce( cell->poissonRatio, cell->youngModulus, cell->mRadius,
                                      mDistance, *mpSingleBondEnergy, *mpAdhesionDensity);

  //  int j = 0;
  // Get contact area
  // Update remaining (contact free) surface area

  force_ij_absolute = fabs(mForce);

  // Store cumulated forces
  cell->lastForce += mForce;

  // Store cumulated absolute forces
  cell->accumulatedForceAbsolute += force_ij_absolute;

  cell->accumulatedPressure += force_ij_absolute/mContactArea;


  mForceVector.Set( barrier->getNormalVector(0),
                    barrier->getNormalVector(1),
                    barrier->getNormalVector(2));

    // Prelim: These methods could be optimized for 2D
  mForceVector.Multiply(mForce);

  cell->directedForce.Add( mForceVector.x, mForceVector.y, mForceVector.z);

}


void
CSInteractionHertz::interact( CellSphericalPolar *cell_1, CellSphericalPolar *cell_2 )
{
    mContactArea =0;

    if ( cell_1->mpDivisionPartner == cell_2 )
        return;

    if (mIs2D)
        mDistance = CSModelTools::GetDistance2D( cell_2->position, cell_1->position );
    else
        mDistance = CSModelTools::GetDistance3D( cell_2->position, cell_1->position );

    if ( ( cell_1->mRadius + cell_2->mRadius ) <= mDistance )
        return;

    // contact area needed for both Hertz force and Cell-Cell friction:
    mContactArea = CSModelTools::GetContactAreaHertz( cell_1->mRadius,
                                                    cell_2->mRadius,
                                                    mDistance);

    cell_1->surfaceContactArea += mContactArea;
    cell_2->surfaceContactArea += mContactArea;


    double adhesionEnergy = CellSphericalPolar::adhesionEnergyMax * CellSphericalPolar::CalculateOverlap(cell_1,cell_2);
    //In Case Of TEST: CellSphericalPolar::adhesionEnergyMax * fmax( cell_1->mPolarDirection.x * cell_2->mPolarDirection.x + cell_1->mPolarDirection.y * cell_2->mPolarDirection.y + cell_1->mPolarDirection.z * cell_2->mPolarDirection.z ,0);


    cell_1->mAdhesionEnergy += adhesionEnergy;
    cell_2->mAdhesionEnergy += adhesionEnergy;

    // Calculate value of (Hertz) force
    mForce = CSModelTools::GetHertzForce( cell_1->poissonRatio, cell_1->youngModulus, cell_1->mRadius,
                                        cell_2->poissonRatio, cell_2->youngModulus, cell_2->mRadius,
                                        mDistance, 1, adhesionEnergy );

    double force_ij_absolute = fabs(mForce);

    // Store cumulated forces
    cell_1->lastForce += mForce;
    cell_2->lastForce += mForce;

    // Store cumulated absolute forces
    cell_1->accumulatedForceAbsolute += force_ij_absolute;
    cell_2->accumulatedForceAbsolute += force_ij_absolute;

    cell_1->accumulatedPressure += force_ij_absolute/mContactArea;
    cell_2->accumulatedPressure += force_ij_absolute/mContactArea;

    if (mIs2D)
    {
        mForceVector.Set(cell_2->position.x - cell_1->position.x,
                         cell_2->position.y - cell_1->position.y,
                         0);
    }
    else
    {
        mForceVector.Set(cell_2->position.x - cell_1->position.x,
                         cell_2->position.y - cell_1->position.y,
                         cell_2->position.z - cell_1->position.z);
    }

    mForceVector.Normalize();
    mForceVector.Multiply(mForce);

    //add directed force
    if ( cell_1->mpDivisionPartner )
    {
        cell_1->directedForce.Add( -mForceVector.x/2, -mForceVector.y/2, -mForceVector.z/2 );
        cell_1->mpDivisionPartner->directedForce.Add( -mForceVector.x/2, -mForceVector.y/2, -mForceVector.z/2 );
    }
    else
        cell_1->directedForce.Add( -mForceVector.x, -mForceVector.y, -mForceVector.z);

    if ( cell_2->mpDivisionPartner )
    {
        cell_2->directedForce.Add( mForceVector.x/2, mForceVector.y/2, mForceVector.z/2 );
        cell_2->mpDivisionPartner->directedForce.Add( mForceVector.x/2, mForceVector.y/2, mForceVector.z/2 );
    }
    else
        cell_2->directedForce.Add(  mForceVector.x,  mForceVector.y,  mForceVector.z );

    cell_1->mpContactsPointers->push_back(cell_2);
    cell_2->mpContactsPointers->push_back(cell_1);
}


void
CSInteractionHertz::interact( ModelElementSphere *cell_1, ModelElementVesselSphere *cell_2 )
{
    mContactArea =0;

    if( cell_1->mStatic && cell_2->mStatic)
        return;

    double force_ij_absolute;

    if (mIs2D)
        mDistance = CSModelTools::GetDistance2D( cell_2->position, cell_1->position );
    else
        mDistance = CSModelTools::GetDistance3D( cell_2->position, cell_1->position );

    if ( ( cell_1->mRadius + cell_2->mRadius ) <= mDistance )
        return;

    // contact area needed for both Hertz force and Cell-Cell friction:
    mContactArea = CSModelTools::GetContactAreaHertz( cell_1->mRadius,
                                                    cell_2->mRadius,
                                                    mDistance);

    cell_1->surfaceContactArea += mContactArea;
    cell_2->surfaceContactArea += mContactArea;

    // Calculate value of (Hertz) force
    mForce = CSModelTools::GetHertzForcePure( cell_1->poissonRatio, cell_1->youngModulus, cell_1->mRadius,
                                            cell_2->poissonRatio, cell_2->youngModulus, cell_2->mRadius,
                                            mDistance );


    force_ij_absolute = fabs(mForce);

    // Store cumulated forces
    cell_1->lastForce += mForce;
    cell_2->lastForce += mForce;

    // Store cumulated absolute forces
    cell_1->accumulatedForceAbsolute += force_ij_absolute;
    cell_2->accumulatedForceAbsolute += force_ij_absolute;

    cell_1->accumulatedPressure += force_ij_absolute/mContactArea;
    cell_2->accumulatedPressure += force_ij_absolute/mContactArea;

    if (mIs2D)
    {
        mForceVector.Set(cell_2->position.x - cell_1->position.x,
                         cell_2->position.y - cell_1->position.y,
                         0);
    }
    else
    {
        mForceVector.Set(cell_2->position.x - cell_1->position.x,
                         cell_2->position.y - cell_1->position.y,
                         cell_2->position.z - cell_1->position.z);
    }

    mForceVector.Normalize();
    mForceVector.Multiply(mForce);

    //add directed force
    cell_1->directedForce.Add( -mForceVector.x, -mForceVector.y, -mForceVector.z);
    cell_2->directedForce.Add(  mForceVector.x,  mForceVector.y,  mForceVector.z);
}


void
CSInteractionHertz::interact( CellSpherical *cell_1, CellSpherical *cell_2 )
{
    mContactArea =0;

    if( cell_1->mStatic && cell_2->mStatic)
        return;

    if ( cell_1->mpDivisionPartner == cell_2 )
        return;

    double force_ij_absolute;

    if (mIs2D)
        mDistance = CSModelTools::GetDistance2D( cell_2->position, cell_1->position );
    else
        mDistance = CSModelTools::GetDistance3D( cell_2->position, cell_1->position );

    if ( ( cell_1->mRadius + cell_2->mRadius ) <= mDistance )
        return;

    // contact area needed for both Hertz force and Cell-Cell friction:
    mContactArea = CSModelTools::GetContactAreaHertz( cell_1->mRadius,
                                                    cell_2->mRadius,
                                                    mDistance);

    cell_1->surfaceContactArea += mContactArea;
    cell_2->surfaceContactArea += mContactArea;

    // Calculate value of (Hertz) force
    mForce = CSModelTools::GetHertzForce( cell_1->poissonRatio, cell_1->youngModulus, cell_1->mRadius,
                                        cell_2->poissonRatio, cell_2->youngModulus, cell_2->mRadius,
                                        mDistance, *mpSingleBondEnergy, *mpAdhesionDensity);


    force_ij_absolute = fabs(mForce);

    // Store cumulated forces
    cell_1->lastForce += mForce;
    cell_2->lastForce += mForce;

    // Store cumulated absolute forces
    cell_1->accumulatedForceAbsolute += force_ij_absolute;
    cell_2->accumulatedForceAbsolute += force_ij_absolute;

    cell_1->accumulatedPressure += force_ij_absolute/mContactArea;
    cell_2->accumulatedPressure += force_ij_absolute/mContactArea;

    if (mIs2D)
    {
        mForceVector.Set(cell_2->position.x - cell_1->position.x,
                         cell_2->position.y - cell_1->position.y,
                         0);
    }
    else
    {
        mForceVector.Set(cell_2->position.x - cell_1->position.x,
                         cell_2->position.y - cell_1->position.y,
                         cell_2->position.z - cell_1->position.z);
    }

    mForceVector.Normalize();
    mForceVector.Multiply(mForce);

    //add directed force
    if ( cell_1->mpDivisionPartner )
    {
        cell_1->directedForce.Add( -mForceVector.x/2, -mForceVector.y/2, -mForceVector.z/2 );
        cell_1->mpDivisionPartner->directedForce.Add( -mForceVector.x/2, -mForceVector.y/2, -mForceVector.z/2 );
    }
    else
        cell_1->directedForce.Add( -mForceVector.x, -mForceVector.y, -mForceVector.z);

    if ( cell_2->mpDivisionPartner )
    {
        cell_2->directedForce.Add( mForceVector.x/2, mForceVector.y/2, mForceVector.z/2 );
        cell_2->mpDivisionPartner->directedForce.Add( mForceVector.x/2, mForceVector.y/2, mForceVector.z/2 );
    }
    else
        cell_2->directedForce.Add(  mForceVector.x,  mForceVector.y,  mForceVector.z );
}


void
CSInteractionHertz::interact( CellSpherical *cell_1, CellSphericalPolar *cell_2 )
{
    ModelElementSphere *sphere1 = static_cast<ModelElementSphere *>(cell_1);
    ModelElementSphere *sphere2 = static_cast<ModelElementSphere *>(cell_2);

    interact( sphere1, sphere2 );
}


void
CSInteractionHertz::interact( CellSpherical *cell, ModelElementBarrierTriangle *barrier )
{

    ModelElementSphere *sphere = static_cast<ModelElementSphere *>(cell);

    interact( sphere, barrier );
}


void
CSInteractionHertz::interact( CellSphericalPolar *cell_1, ModelElementVesselSphere *cell_2 )
{
    ModelElementSphere * sphere = static_cast<ModelElementSphere *>(cell_1);

    interact( sphere, cell_2 );
}


void
CSInteractionHertz::interact( CellSphericalPolar *cell, ModelElementBarrierTriangle *barrier )
{

    ModelElementSphere *sphere = static_cast<ModelElementSphere *>(cell);

    interact( sphere, barrier );
}


void
CSInteractionHertz::interact( CellSpherical *cell_1, ModelElementVesselSphere *cell_2 )
{
    ModelElementSphere *sphere1 = static_cast<ModelElementSphere *>(cell_1);

    interact( sphere1, cell_2 );
}


void
CSInteractionHertz::interact( ModelElementVesselSphere *cell_1, ModelElementVesselSphere *cell_2 )
{
    // ToDo do not interact with next neighbors in vessel chain!
    ModelElementSphere *sphere1 = static_cast<ModelElementSphere *>(cell_1);

    //exclude neighbors
    for( unsigned int i = 0 ; i < cell_1->mvpBlackList.size() ; i++)
      if( cell_1->mvpBlackList[i] == cell_2 )
      {
          mContactArea =0;
          return;
      }

    interact( sphere1, cell_2 );
}


void
CSInteractionHertz::interact( ModelElementVesselSphere *cell, ModelElementBarrierTriangle *barrier )
{
    ModelElementSphere *sphere = static_cast<ModelElementSphere *>(cell);

    interact( sphere, barrier );
}


void
CSInteractionHertz::interact( CellSpherical * cell, ModelElementHollowSphere * capsule )
{
    mContactArea =0;

    double force_ij_absolut;

    mDistance = capsule->mRadius;

    if (mIs2D)
        mDistance -= std::sqrt( cell->position.x*cell->position.x + cell->position.y*cell->position.y );
    else
        mDistance -= std::sqrt( cell->position.x*cell->position.x + cell->position.y*cell->position.y + cell->position.z*cell->position.z );

    if ( mDistance >= cell->mRadius )
        return;

    // ToDo:  use correct contact area calculation for convex-concave surface contacts:
    // contact area needed for both Hertz force and Cell-Cell friction.
    // Interpreting the capsule as a barrier (BarrierTriangle, i.e.).
    mContactArea = CSModelTools::GetContactAreaHertz( cell->mRadius,
                                                      mDistance);
    if (mContactArea != 0)
        cell->surfaceContactArea += mContactArea;

    // Calculate value of (Hertz) force
    mForce = CSModelTools::GetHertzForce( cell->poissonRatio, cell->youngModulus, cell->mRadius,
                                          mDistance, *mpSingleBondEnergy, *mpAdhesionDensity);


    force_ij_absolut = fabs(mForce);

    // Store cumulated forces
    cell->lastForce += mForce;

    // Store cumulated absolute forces
    cell->accumulatedForceAbsolute   += force_ij_absolut;
    capsule->accumulatedForceAbsolute += force_ij_absolut;

    // radial force vector
    if (mIs2D)
    {
        mForceVector.Set(cell->position.x,
                         cell->position.y,
                         0);
    }
    else
    {
        mForceVector.Set(cell->position.x,
                         cell->position.y,
                         cell->position.z );
    }

    mForceVector.Normalize();
    mForceVector.Multiply(mForce);

    cell->directedForce.Add( -mForceVector.x, -mForceVector.y, -mForceVector.z);
}


void
CSInteractionHertz::interact( CSVoxel * cube, ModelElementSphere *sphere )
{
    // CellSpherical and CellSphericalPolar take precedence
    if ( cube->contains == ModelElement::TypeCellSpherical
         || cube->contains == ModelElement::TypeCellSphericalPolar )
        return;

    BoundingBox * cubeLimits = cube->getBoundingBox();

    // do the two objects actually intersect?
    double dist_squared = sphere->mRadius * sphere->mRadius;

    if ( sphere->position.x < cubeLimits->xmin)
        dist_squared -= (sphere->position.x - cubeLimits->xmin)
            * (sphere->position.x - cubeLimits->xmin);
    else if ( sphere->position.x > cubeLimits->xmax)
        dist_squared -= (sphere->position.x - cubeLimits->xmax)
            * (sphere->position.x - cubeLimits->xmax);

    if (sphere->position.y < cubeLimits->ymin)
        dist_squared -= (sphere->position.y - cubeLimits->ymin)
            * (sphere->position.y - cubeLimits->ymin);
    else if (sphere->position.y > cubeLimits->ymax)
        dist_squared -= (sphere->position.y - cubeLimits->ymax)
            * (sphere->position.y - cubeLimits->ymax);

    if (sphere->position.z < cubeLimits->zmin)
        dist_squared -= (sphere->position.z - cubeLimits->zmin)
            * (sphere->position.z - cubeLimits->zmin);
    else if (sphere->position.z > cubeLimits->zmax)
        dist_squared -= (sphere->position.z - cubeLimits->zmax)
            * (sphere->position.z - cubeLimits->zmax);

    if ( dist_squared <= 0. )
        return;

    cube->contains = sphere->mType;
}

// added by Jieling
void
CSInteractionHertz::interact(ModelElementLatticeNode * node, CellSpherical *cell)
{
	// only add distance and contactArea for cell-cell friction
	// Hertz force is added in cell-fiber interaction part
	/*
	mContactArea = 0;

	if (mIs2D)
		mDistance = CSModelTools::GetDistance2D(node->position, cell->position);
	else
		mDistance = CSModelTools::GetDistance3D(node->position, cell->position);

	if (node->mRadius + cell->mRadius <= mDistance)
		return;

	mContactArea = CSModelTools::GetContactAreaHertz(node->mRadius, cell->mRadius, mDistance);

	node->surfaceContactArea += mContactArea;
	cell->surfaceContactArea += mContactArea;
	*/
}

void 
CSInteractionHertz::interact(ModelElementLatticeNode * node, CellSphericalPolar *cell)
{
	// only add distance and contactArea for cell-cell friction
	// Hertz force is added in cell-fiber interaction part
	/*
	mContactArea = 0;

	if (mIs2D)
		mDistance = CSModelTools::GetDistance2D(node->position, cell->position);
	else
		mDistance = CSModelTools::GetDistance3D(node->position, cell->position);

	if (node->mRadius + cell->mRadius <= mDistance)
		return;

	mContactArea = CSModelTools::GetContactAreaHertz(node->mRadius, cell->mRadius, mDistance);

	node->surfaceContactArea += mContactArea;
	cell->surfaceContactArea += mContactArea;
	*/
}

void
CSInteractionHertz::interact(ModelElementLatticeNode *cell, ModelElementBarrierTriangle *barrier)
{
	//ModelElementSphere *sphere = static_cast<ModelElementSphere *>(cell);

	//interact( sphere, barrier);
}

void 
CSInteractionHertz::interact(ModelElementLatticeNode * node, ModelElementVesselSphere * sphere)
{
	// only add distance and contactArea for cell-cell friction
	// Hertz force is added in cell-fiber interaction part
	/*
	mContactArea = 0;

	if (mIs2D)
		mDistance = CSModelTools::GetDistance2D(node->position, sphere->position);
	else
		mDistance = CSModelTools::GetDistance3D(node->position, sphere->position);

	if (node->mRadius + sphere->mRadius <= mDistance)
		return;

	mContactArea = CSModelTools::GetContactAreaHertz(node->mRadius, sphere->mRadius, mDistance);

	node->surfaceContactArea += mContactArea;
	sphere->surfaceContactArea += mContactArea;
	*/
}

void 
CSInteractionHertz::interact(CellSpherical * cell, LatticeSpring * spring)
{
	mContactArea = 0;

	if (spring->mSpringType == LatticeSpring::TypeLinearSpring)
	{
		LinearSpring *lspring = static_cast<LinearSpring *>(spring);
		double eL = lspring->getLength();
		// 0: Start; 1: End
		double nodex1 = lspring->nodes().at(0)->position.x;
		double nodey1 = lspring->nodes().at(0)->position.y;
		double nodez1 = lspring->nodes().at(0)->position.z;
		double nodex2 = lspring->nodes().at(1)->position.x;
		double nodey2 = lspring->nodes().at(1)->position.y;
		double nodez2 = lspring->nodes().at(1)->position.z;

		double iterx = cell->position.x;
		double itery = cell->position.y;
		double iterz = cell->position.z;

		double iternode1 = sqrt((iterx - nodex1)*(iterx - nodex1) +
			(itery - nodey1)*(itery - nodey1) +
			(iterz - nodez1)*(iterz - nodez1));
		double iternode2 = sqrt((iterx - nodex2)*(iterx - nodex2) +
			(itery - nodey2)*(itery - nodey2) +
			(iterz - nodez2)*(iterz - nodez2));
		bool touchStart = false;
		bool touchEnd = false;
		//
		if (iternode1 <= cell->mRadius) touchStart = true;
		if (iternode2 <= cell->mRadius) touchEnd = true;
		//
		double crosx = (itery - nodey1) * (nodez1 - nodez2) - (iterz - nodez1) * (nodey1 - nodey2);
		double crosy = (iterz - nodez1) * (nodex1 - nodex2) - (iterx - nodex1) * (nodez1 - nodez2);
		double crosz = (iterx - nodex1) * (nodey1 - nodey2) - (itery - nodey1) * (nodex1 - nodex2);
		double distE = sqrt(crosx * crosx + crosy * crosy + crosz * crosz) / eL;
		if (eL < 0.000001) distE = sqrt((iterx - nodex1)*(iterx - nodex1) + (itery - nodey1)*(itery - nodey1) + (iterz - nodez1)*(iterz - nodez1));
		if (distE > cell->mRadius) 
			return;
		// Least squares
		double ac = (nodex2 - nodex1) * (nodex2 - nodex1) +
			(nodey2 - nodey1) * (nodey2 - nodey1) +
			(nodez2 - nodez1) * (nodez2 - nodez1);
		double bc = 2 * (nodex2 - nodex1) * (nodex1 - iterx) +
			2 * (nodey2 - nodey1) * (nodey1 - itery) +
			2 * (nodez2 - nodez1) * (nodez1 - iterz);
		double ts = -bc / 2 / ac;
		// to check if the cell is really in touch with the spring
		if ((ts > 1 || ts < 0) && !touchStart && !touchEnd)
			return;
		// node3: the point with shortest distance to the cell
		double nodex3 = nodex1 + ts * (nodex2 - nodex1);
		double nodey3 = nodey1 + ts * (nodey2 - nodey1);
		double nodez3 = nodez1 + ts * (nodez2 - nodez1);
		// from the center to the edge
		double i3d = sqrt((nodex3 - iterx)*(nodex3 - iterx) +
			(nodey3 - itery)*(nodey3 - itery) +
			(nodez3 - iterz)*(nodez3 - iterz));
		double i3x = (nodex3 - iterx) / i3d;
		double i3y = (nodey3 - itery) / i3d;
		double i3z = (nodez3 - iterz) / i3d; // direction perpendicular to the edge
		// the direction from cell to the node
		double itern1x = nodex1 - iterx;
		double itern1y = nodey1 - itery;
		double itern1z = nodez1 - iterz;
		double itern1d = sqrt(itern1x*itern1x + itern1y*itern1y + itern1z*itern1z);
		itern1x /= itern1d;
		itern1y /= itern1d;
		itern1z /= itern1d; // cell to Start
		double itern2x = nodex2 - iterx;
		double itern2y = nodey2 - itery;
		double itern2z = nodez2 - iterz;
		double itern2d = sqrt(itern2x*itern2x + itern2y*itern2y + itern2z*itern2z);
		itern2x /= itern2d;
		itern2y /= itern2d;
		itern2z /= itern2d; // cell to End

		if (lspring->nodes().at(0)->vesselNeighbor != NULL) // Start
		{
			// the force from the cell upon Start
			double csx = i3x;
			double csy = i3y;
			double csz = i3z;
			double sOverlap = cell->mRadius - distE;
			int contactType = 0; // 0: the edge  1: the end point
			//
			if (touchStart && !touchEnd)
			{
				double productIS = itern1x * i3x + itern1y * i3y + itern1z * i3z;
				if (productIS < sqrt(2) / 2)
				{
					csx = itern1x;
					csy = itern1y;
					csz = itern1z;
					sOverlap = cell->mRadius - itern1d;
					contactType = 1;
				}
			}
			else if (!touchStart && touchEnd)
			{
				double productIE = itern2x * i3x + itern2y * i3y + itern2z * i3z;
				if (productIE < sqrt(2) / 2)
				{
					csx = itern2x;
					csy = itern2y;
					csz = itern2z;
					sOverlap = cell->mRadius - itern2d;
					contactType = 1;
				}
			}
			//
			double vnx = lspring->nodes().at(0)->vesselNeighbor->position.x;
			double vny = lspring->nodes().at(0)->vesselNeighbor->position.y;
			double vnz = lspring->nodes().at(0)->vesselNeighbor->position.z;
			double cosAVmax = 0.;
			int vnMinIndex = -1;
			int vnSize = lspring->nodes().at(0)->vesselNeighbor->mvpNeighbor.size();
			for (int j = 0; j < vnSize; j++)
			{
				if (lspring->nodes().at(0)->vesselNeighbor->mvpNeighbor[j]->mVesselType != 3) continue; // 3: sinusoid
				if (lspring->nodes().at(0)->vesselNeighbor->mvpNeighbor[j]->mRadius > 0.201) continue; // sinusoid radius: 0.2
				double vnjx = lspring->nodes().at(0)->vesselNeighbor->mvpNeighbor[j]->position.x;
				double vnjy = lspring->nodes().at(0)->vesselNeighbor->mvpNeighbor[j]->position.y;
				double vnjz = lspring->nodes().at(0)->vesselNeighbor->mvpNeighbor[j]->position.z;
				double vnnx = vnjx - vnx;
				double vnny = vnjy - vny;
				double vnnz = vnjz - vnz; // vector from bound vessel cell to its neighbor j
				double cosAV = (i3x * vnnx + i3y * vnny + i3z * vnnz) / sqrt(vnnx*vnnx + vnny*vnny + vnnz*vnnz); // normalized vectors
				if (cosAV > 0) // positive projection onto the vessel
				{
					if (cosAV > cosAVmax)
					{
						cosAVmax = cosAV;
						vnMinIndex = j;
					}
				}
			}
			if (vnMinIndex >= 0) // positive projection found
			{
				lspring->nodes().at(0)->touched = true;
				Vector3f vVNMinIndex = lspring->nodes().at(0)->vesselNeighbor->mvpNeighbor[vnMinIndex]->position -
					lspring->nodes().at(0)->vesselNeighbor->position;
				double vNorm = vVNMinIndex.Normalize();
				// calculate value of (Hertz) force
				double mForce = CSModelTools::GetHertzForceSphereCylinder(
					cell->poissonRatio, cell->youngModulus, cell->mRadius,
					lspring->mPoisson, lspring->mYoung, lspring->mRadius,
					sOverlap, *mpSingleBondEnergy, *mpAdhesionDensity, contactType);
				vVNMinIndex.Multiply(mForce * cosAVmax);

				// get contact area
				mContactArea = CSModelTools::GetContactAreaHertzSphereCylinder(cell->mRadius, lspring->mRadius, sOverlap, contactType);

				lspring->nodes().at(0)->surfaceContactArea += mContactArea;
				cell->surfaceContactArea += mContactArea;

				// store cumulated forces
				lspring->nodes().at(0)->lastForce += mForce;
				cell->lastForce += mForce;

				// store cumulated absolute forces
				lspring->nodes().at(0)->accumulatedForceAbsolute += fabs(mForce);
				cell->accumulatedForceAbsolute += fabs(mForce);

				// add directed force
				lspring->nodes().at(0)->directedForce.Add(vVNMinIndex);
				cell->directedForce.Add(vVNMinIndex * -1.);
			}
		}
		if (lspring->nodes().at(1)->vesselNeighbor != NULL) // End
		{
			// the force from the cell upon End
			double csx = i3x;
			double csy = i3y;
			double csz = i3z;
			double sOverlap = cell->mRadius - distE;
			int contactType = 0; // 0: the edge  1: the end point
			//
			if (!touchStart && touchEnd)
			{
				double productIE = itern2x * i3x + itern2y * i3y + itern2z * i3z;
				if (productIE < sqrt(2) / 2)
				{
					csx = itern2x;
					csy = itern2y;
					csz = itern2z;
					sOverlap = cell->mRadius - itern2d;
					contactType = 1;
				}
			}
			else if (touchStart && !touchEnd)
			{
				double productIS = itern1x * i3x + itern1y * i3y + itern1z * i3z;
				if (productIS < sqrt(2) / 2)
				{
					csx = itern1x;
					csy = itern1y;
					csz = itern1z;
					sOverlap = cell->mRadius - itern1d;
					contactType = 1;
				}
			}
			//
			double vnx = lspring->nodes().at(1)->vesselNeighbor->position.x;
			double vny = lspring->nodes().at(1)->vesselNeighbor->position.y;
			double vnz = lspring->nodes().at(1)->vesselNeighbor->position.z;
			double cosAVmax = 0.;
			int vnMinIndex = -1;
			int vnSize = lspring->nodes().at(1)->vesselNeighbor->mvpNeighbor.size();
			for (int j = 0; j < vnSize; j++)
			{
				if (lspring->nodes().at(1)->vesselNeighbor->mvpNeighbor[j]->mVesselType != 3) continue; // 3: sinusoid
				if (lspring->nodes().at(1)->vesselNeighbor->mvpNeighbor[j]->mRadius > 0.201) continue; // sinusoid radius: 0.2
				double vnjx = lspring->nodes().at(1)->vesselNeighbor->mvpNeighbor[j]->position.x;
				double vnjy = lspring->nodes().at(1)->vesselNeighbor->mvpNeighbor[j]->position.y;
				double vnjz = lspring->nodes().at(1)->vesselNeighbor->mvpNeighbor[j]->position.z;
				double vnnx = vnjx - vnx;
				double vnny = vnjy - vny;
				double vnnz = vnjz - vnz; // vector from bound vessel cell to its neighbor j
				double cosAV = (i3x * vnnx + i3y * vnny + i3z * vnnz) / sqrt(vnnx*vnnx + vnny*vnny + vnnz*vnnz); // normalized vectors
				if (cosAV > 0) // positive projection onto the vessel
				{
					if (cosAV > cosAVmax)
					{
						cosAVmax = cosAV;
						vnMinIndex = j;
					}
				}
			}
			if (vnMinIndex >= 0) // positive projection found
			{
				lspring->nodes().at(1)->touched = true;
				Vector3f vVNMinIndex = lspring->nodes().at(1)->vesselNeighbor->mvpNeighbor[vnMinIndex]->position -
					lspring->nodes().at(1)->vesselNeighbor->position;
				double vNorm = vVNMinIndex.Normalize();
				// calculate value of (Hertz) force
				double mForce = CSModelTools::GetHertzForceSphereCylinder(
					cell->poissonRatio, cell->youngModulus, cell->mRadius,
					lspring->mPoisson, lspring->mYoung, lspring->mRadius,
					sOverlap, *mpSingleBondEnergy, *mpAdhesionDensity, contactType);
				vVNMinIndex.Multiply(mForce * cosAVmax);

				// get contact area
				mContactArea = CSModelTools::GetContactAreaHertzSphereCylinder(cell->mRadius, lspring->mRadius, sOverlap, contactType);

				lspring->nodes().at(1)->surfaceContactArea += mContactArea;
				cell->surfaceContactArea += mContactArea;

				// store cumulated forces
				lspring->nodes().at(1)->lastForce += mForce;
				cell->lastForce += mForce;

				// store cumulated absolute forces
				lspring->nodes().at(1)->accumulatedForceAbsolute += fabs(mForce);
				cell->accumulatedForceAbsolute += fabs(mForce);

				// add directed force
				lspring->nodes().at(1)->directedForce.Add(vVNMinIndex);
				cell->directedForce.Add(vVNMinIndex * -1.);
			}
		}
	}
	else
	{
		//RotationalSpring *rspring = static_cast<RotationalSpring *>(spring);
		// 0: Start; 1: Center; 2: End
	}
}

void 
CSInteractionHertz::interact(CellSphericalPolar * cell, LatticeSpring * spring)
{
	CellSpherical * cell_1 = static_cast<CellSpherical *>(cell);

	interact(cell_1, spring);
}