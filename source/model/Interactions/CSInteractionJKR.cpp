////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSInteractionJKR.cpp                                          //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-08-21 16:59:33                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "CSInteractionJKR.h"

#include "../../tools/model/CSModelTools.h"
#include "../Cell/CellSpherical.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>


CSInteractionJKR::CSInteractionJKR()
    : CSInteraction(),
      mContactEstablished(false),
      mForceVector(0.,0.,0.),
      mDistance(0.),
      mContactArea(0.),
      mForce(0.),
      mIs2D(false)
{}


void
CSInteractionJKR::interact( CellSpherical *cell_1, CellSpherical *cell_2 )
{
    mContactArea = 0;

    if ( cell_1->mpDivisionPartner == cell_2 )
        return;

    if (mIs2D)
        mDistance = CSModelTools::GetDistance2D( cell_2->position, cell_1->position );
    else
        mDistance = CSModelTools::GetDistance3D( cell_2->position, cell_1->position );

    double delta = cell_1->mRadius + cell_2->mRadius - mDistance;

    if  (delta < 0.0 && !mContactEstablished)
        return;

    //// PVL
    double E_eff = 1.0
                 / ( (1 - cell_1->poissonRatio * cell_1->poissonRatio)
                     /cell_1->youngModulus
                    +(1-  cell_2->poissonRatio * cell_2->poissonRatio)
                     / cell_2->youngModulus );

    double R_eff = cell_2->mRadius * cell_1->mRadius
                 /(cell_2->mRadius + cell_1->mRadius);

    double adhesion_constant = (*mpSingleBondEnergy) * (*mpAdhesionDensity);

    double dcrit = std::pow(
        std::pow( std::sqrt(3.)* M_PI * adhesion_constant /(8*E_eff), 2) * R_eff,
        0.333333);

    // pulling off distance when two cells are in contact. As long as the overlap
    // is larger than -dcrit then there is still contact !!!
    //Overlap: overlap > 0  => contact.
    if ( delta < -dcrit )
        return;

    //update JKR contact force and area.

    CSModelTools::GetJKRForce(E_eff ,R_eff, delta, adhesion_constant, mForce, mContactArea, dcrit, 1); // 1: CellSpherical - CellSpherical

    cell_1->surfaceContactArea += mContactArea;
    cell_2->surfaceContactArea += mContactArea;

    applyForce( cell_1, cell_2 );
}


void
CSInteractionJKR::interact( CellSpherical *cell, ModelElementBarrierTriangle *barrier )
{
    ModelElementSphere * sphere = static_cast<ModelElementSphere*>(cell);

    interact( sphere, barrier );
}


// non-adhesive interaction -> delegate to
// interact( ModelElementSphere *, ModelElementSphere * );
void
CSInteractionJKR::interact( CellSpherical *cell, ModelElementVesselSphere * sphere )
{
    ModelElementSphere * sphere_1 = static_cast<ModelElementSphere *>(cell);
    ModelElementSphere * sphere_2 = static_cast<ModelElementSphere *>(sphere);

    interact(sphere_1, sphere_2);
}


void
CSInteractionJKR::interact( CellSphericalPolar *cell_1, CellSphericalPolar * cell_2 )
{
    mContactArea = 0;

    if ( cell_1->mpDivisionPartner == cell_2 )
        return;

    if (mIs2D)
        mDistance = CSModelTools::GetDistance2D( cell_2->position, cell_1->position );
    else
        mDistance = CSModelTools::GetDistance3D( cell_2->position, cell_1->position );

    double delta = cell_1->mRadius + cell_2->mRadius - mDistance;

    if  (delta < 0.0 && !mContactEstablished)
        return;

    //// PVL
    double E_eff = 1.0
                 / ( (1 - cell_1->poissonRatio * cell_1->poissonRatio)
                     /cell_1->youngModulus
                    +(1-  cell_2->poissonRatio * cell_2->poissonRatio)
                     / cell_2->youngModulus );

    double R_eff = cell_2->mRadius * cell_1->mRadius
                 /(cell_2->mRadius + cell_1->mRadius);

    double adhesion_constant;
    double dcrit;

    if ( mContactEstablished )
    {
        adhesion_constant = (*mpSingleBondEnergy) * (*mpAdhesionDensity);

        double dcrit = std::pow(
            std::pow( std::sqrt(3.)* M_PI * adhesion_constant /(8*E_eff), 2) * R_eff,
            0.333333);

        // pulling off distance when two cells are in contact. As long as the overlap
        // is larger than -dcrit then there is still contact !!!
        //Overlap: overlap > 0  => contact.
        if ( delta < -dcrit )
            return;

        //update JKR contact force and area.
		//if (delta != delta)
		//{
		//	std::cerr << "	-> delta: " << delta << ", mForce: " << mForce << ", cell1_Radius: " << cell_1->mRadius << ", cell2_Radius: " << cell_2->mRadius;
		//	std::cerr << ", cell1: " << cell_1->position.x << "_" << cell_1->position.y << "_" << cell_1->position.z;
		//	std::cerr << ", cell2: " << cell_2->position.x << "_" << cell_2->position.y << "_" << cell_2->position.z << std::endl;
		//}

        // first calculation assuming adhesive region overlap on the mContactArea:
        CSModelTools::GetJKRForce(E_eff ,R_eff, delta, adhesion_constant, mForce, mContactArea, dcrit, 2); // 2: CellSphericalPolar - CellSphericalPolar

        // does the contact area have an adhesive overlap?
        if ( CellSphericalPolar::CalculateOverlap( cell_1, cell_2, mContactArea / (2*M_PI)) == 0. )
        {
            // recalculate without adhesion -> use Hertz force (without adhesion)
            if ( ( cell_1->mRadius + cell_2->mRadius ) <= mDistance )
            {
                mContactArea = 0;
                return;
            }

            mContactArea = CSModelTools::GetContactAreaHertz( cell_1->mRadius,
                                                              cell_2->mRadius,
                                                              mDistance);
            mForce = CSModelTools::GetHertzForcePure( cell_1->poissonRatio, cell_1->youngModulus, cell_1->mRadius,
                                                      cell_2->poissonRatio, cell_2->youngModulus, cell_2->mRadius,
                                                      mDistance );
        }
    }
    else // mContactEstablished == false
    {
        // new contact, first assuming no adhesion -> use Hertz force
        if ( ( cell_1->mRadius + cell_2->mRadius ) <= mDistance )
        {
            mContactArea = 0;
            return;
        }

        mContactArea = CSModelTools::GetContactAreaHertz( cell_1->mRadius,
                                                          cell_2->mRadius,
                                                          mDistance);

        // does the contact area have an adhesive overlap?
        if ( CellSphericalPolar::CalculateOverlap( cell_1, cell_2, mContactArea / (2*M_PI)) != 0. )
        {
            adhesion_constant = (*mpSingleBondEnergy) * (*mpAdhesionDensity);
            double dcrit = std::pow(
                std::pow( std::sqrt(3.)* M_PI * adhesion_constant /(8*E_eff), 2) * R_eff,0.333333);

            CSModelTools::GetJKRForce(E_eff ,R_eff, delta, adhesion_constant, mForce, mContactArea, dcrit, 2);
        }
        else
        {
            mForce = CSModelTools::GetHertzForcePure( cell_1->poissonRatio, cell_1->youngModulus, cell_1->mRadius,
                                                      cell_2->poissonRatio, cell_2->youngModulus, cell_2->mRadius,
                                                      mDistance );
        }
    }

    cell_1->surfaceContactArea += mContactArea;
    cell_2->surfaceContactArea += mContactArea;

    applyForce( cell_1, cell_2 );
}


void
CSInteractionJKR::interact( CellSphericalPolar *cell, ModelElementBarrierTriangle *barrier)
{
    ModelElementSphere * sphere = static_cast<ModelElementSphere *>(cell);

    interact( sphere, barrier );
}


void
CSInteractionJKR::interact( CellSphericalPolar * cell, ModelElementVesselSphere *sphere )
{
    ModelElementSphere *sphere_1 = static_cast<ModelElementSphere *>(cell);
    ModelElementSphere *sphere_2 = static_cast<ModelElementSphere *>(sphere);

    interact(sphere_1, sphere_2);
}


// interaction with (non-neighbor) vessel spheres are non-adhesive, so use
// interact(ModelElementSphere*,ModelElementSphere*).
void
CSInteractionJKR::interact( ModelElementVesselSphere *sphere_1, ModelElementVesselSphere *sphere_2 )
{
    //exclude neighbors
    for( unsigned int i = 0 ; i < sphere_1->mvpBlackList.size() ; i++)
      if( sphere_1->mvpBlackList[i] == sphere_2 )
      {
          mContactArea =0;
          return;
      }

    ModelElementSphere *sphere1 = static_cast<ModelElementSphere *>(sphere_1);
    ModelElementSphere *sphere2 = static_cast<ModelElementSphere *>(sphere_2);

    interact( sphere1, sphere2 );
}


void
CSInteractionJKR::interact( ModelElementVesselSphere *sphere, ModelElementBarrierTriangle *barrier )
{
    ModelElementSphere * sphere_1 = static_cast<ModelElementSphere*>(sphere);

    interact( sphere_1, barrier );
}


// Method used for non-adhesive interactions between spherical elements,
// therefore it is using CSModelTools::GetHertzForcePure();
void
CSInteractionJKR::interact( ModelElementSphere *sphere_1, ModelElementSphere *sphere_2 )
{
    mContactArea = 0;

    if (mIs2D)
        mDistance = CSModelTools::GetDistance2D( sphere_2->position, sphere_1->position );
    else
        mDistance = CSModelTools::GetDistance3D( sphere_2->position, sphere_1->position );

    if ( ( sphere_1->mRadius + sphere_2->mRadius ) <= mDistance )
        return;

    // contact area needed for both Hertz force and Cell-Cell friction:
    mContactArea = CSModelTools::GetContactAreaHertz( sphere_1->mRadius,
                                                      sphere_2->mRadius,
                                                      mDistance);

    // No adhesion between cells and vessels, use less expensive Hertz
    mForce = CSModelTools::GetHertzForcePure( sphere_1->poissonRatio, sphere_1->youngModulus, sphere_1->mRadius,
                                              sphere_2->poissonRatio, sphere_2->youngModulus, sphere_2->mRadius,
                                              mDistance);

    sphere_1->surfaceContactArea += mContactArea;
    sphere_2->surfaceContactArea += mContactArea;

    applyForce( sphere_1, sphere_2 );
}


// Method used for non-adhesive interactions between spherical elements and
// planar barrier triangles, therefore it is deploying
// CSModelTools::GetHertzForcePure();
void
CSInteractionJKR::interact( ModelElementSphere *sphere, ModelElementBarrierTriangle *barrier )
{
  mContactArea =0;

  double force_ij_absolute;

  mDistance = CSModelTools::GetDistance3D( sphere->position, barrier->getNormalVector(), barrier->getD() );

//  cell->directedForce.Add( mDistance*barrier->getNormalVector(0) , mDistance*barrier->getNormalVector(1) , mDistance*barrier->getNormalVector(2) );

  if ( sphere->mRadius <= mDistance )
    return;

  // if ( useHertz )
  // contact area needed for both Hertz force and Cell-Cell friction:
  mContactArea = CSModelTools::GetContactAreaHertz( sphere->mRadius,
                                                    mDistance);
  sphere->surfaceContactArea += mContactArea;

  // Calculate value of (Hertz) force
  mForce = CSModelTools::GetHertzForce( sphere->poissonRatio, sphere->youngModulus, sphere->mRadius,
                                        mDistance, *mpSingleBondEnergy, *mpAdhesionDensity);

  //  int j = 0;
  // Get contact area
  // Update remaining (contact free) surface area

  force_ij_absolute = fabs(mForce);

  // Store cumulated forces
  sphere->lastForce += mForce;

  // Store cumulated absolute forces
  sphere->accumulatedForceAbsolute += force_ij_absolute;

  sphere->accumulatedPressure += force_ij_absolute/mContactArea;


  mForceVector.Set( barrier->getNormalVector(0),
                    barrier->getNormalVector(1),
                    barrier->getNormalVector(2));

    // Prelim: These methods could be optimized for 2D
  mForceVector.Multiply(mForce);

  sphere->directedForce.Add( mForceVector.x, mForceVector.y, mForceVector.z);
}


void
CSInteractionJKR::interact( CellSpherical * cell, ModelElementHollowSphere *capsule )
{
    ModelElementSphere * sphere = static_cast<ModelElementSphere *>(cell);

    this->interact(sphere,capsule);
}


void
CSInteractionJKR::interact( ModelElementSphere *sphere, ModelElementHollowSphere *capsule )
{
    mContactArea =0;

    mDistance = capsule->mRadius;

    if (mIs2D)
        mDistance -= std::sqrt( (sphere->position.x - capsule->position.x)*(sphere->position.x - capsule->position.x)
                                +(sphere->position.y - capsule->position.y)*(sphere->position.y - capsule->position.y));
    else
        mDistance -= std::sqrt( (sphere->position.x - capsule->position.x)*(sphere->position.x - capsule->position.x)
                                +(sphere->position.y - capsule->position.y)*(sphere->position.y - capsule->position.y)
                                +(sphere->position.z - capsule->position.z)*(sphere->position.z - capsule->position.z) );

    double delta = mDistance - sphere->mRadius;

    if ( delta > 0. && !mContactEstablished )
        return;


    double E_eff = sphere->youngModulus / ( 1 - sphere->poissonRatio * sphere->poissonRatio );
    double adhesion_constant = (*mpSingleBondEnergy) * (*mpAdhesionDensity);

    double dcrit = std::pow( sphere->mRadius
                             * std::pow( std::sqrt(3.) * M_PI * adhesion_constant / (8*E_eff), 2 ),
                             0.333333 );

    if ( delta >= dcrit )
        return;

    CSModelTools::GetJKRForce(E_eff ,sphere->mRadius, -delta, adhesion_constant, mForce, mContactArea, dcrit, 3); // 3: ModelElementSphere - ModelElementHollowSphere

    sphere->surfaceContactArea += mContactArea;

    double force_ij_absolute = fabs(mForce);

    sphere->accumulatedForceAbsolute += force_ij_absolute;

    sphere->accumulatedPressure += force_ij_absolute/mContactArea;

    if ( mIs2D )
    {
        mForceVector.Set( capsule->position.x - sphere->position.x,
                          capsule->position.y - sphere->position.y,
                          0 );
    }
    else
    {
        mForceVector.Set( capsule->position.x - sphere->position.x,
                          capsule->position.y - sphere->position.y,
                          capsule->position.z - sphere->position.z );
    }

    mForceVector.Normalize();
    mForceVector.Multiply( mForce );

    sphere->directedForce.Add( mForceVector.x, mForceVector.y, mForceVector.z );
}

// added by Jieling
void
CSInteractionJKR::interact(ModelElementLatticeNode * node, CellSpherical * cell)
{
	/*
	mContactArea = 0;

	if (mIs2D)
		mDistance = CSModelTools::GetDistance2D(node->position, cell->position);
	else
		mDistance = CSModelTools::GetDistance3D(node->position, cell->position);

	double delta = node->mRadius + cell->mRadius - mDistance;

	if (delta < 0.0 && !mContactEstablished)
		return;

	double E_eff = 1.0
		/ ((1 - node->poissonRatio * node->poissonRatio)
			/ node->youngModulus
			+ (1 - cell->poissonRatio * cell->poissonRatio)
			/ cell->youngModulus);

	double R_eff = cell->mRadius * node->mRadius
		/ (cell->mRadius + node->mRadius);

	double adhesion_constant = (*mpSingleBondEnergy) * (*mpAdhesionDensity);

	double dcrit = std::pow(
		std::pow(std::sqrt(3.)* M_PI * adhesion_constant / (8 * E_eff), 2) * R_eff,
		0.333333);

	if (delta < -dcrit)
		return;

	// only take the contact area
	CSModelTools::GetJKRForce(E_eff, R_eff, delta, adhesion_constant, mForce, mContactArea, dcrit, 1); // 1: CellSpherical - CellSpherical

	node->surfaceContactArea += mContactArea;
	cell->surfaceContactArea += mContactArea;
	*/
}

void
CSInteractionJKR::interact(ModelElementLatticeNode * node, CellSphericalPolar * cell)
{
	/*
	mContactArea = 0;

	if (mIs2D)
		mDistance = CSModelTools::GetDistance2D(node->position, cell->position);
	else
		mDistance = CSModelTools::GetDistance3D(node->position, cell->position);

	double delta = node->mRadius + cell->mRadius - mDistance;

	if (delta < 0.0 && !mContactEstablished)
		return;

	double E_eff = 1.0
		/ ((1 - node->poissonRatio * node->poissonRatio)
			/ node->youngModulus
			+ (1 - cell->poissonRatio * cell->poissonRatio)
			/ cell->youngModulus);

	double R_eff = cell->mRadius * node->mRadius
		/ (cell->mRadius + node->mRadius);

	double adhesion_constant = (*mpSingleBondEnergy) * (*mpAdhesionDensity);

	double dcrit = std::pow(
		std::pow(std::sqrt(3.)* M_PI * adhesion_constant / (8 * E_eff), 2) * R_eff,
		0.333333);

	if (delta < -dcrit)
		return;

	// only take the contact area
	CSModelTools::GetJKRForce(E_eff, R_eff, delta, adhesion_constant, mForce, mContactArea, dcrit, 1); // 1: CellSpherical - CellSpherical

	node->surfaceContactArea += mContactArea;
	cell->surfaceContactArea += mContactArea;
	*/
}

void 
CSInteractionJKR::interact(ModelElementLatticeNode * node, ModelElementBarrierTriangle * barrier)
{
	/*
	mContactArea = 0;

	mDistance = CSModelTools::GetDistance3D( node->position, barrier->getNormalVector(), barrier->getD() );

	if (node->mRadius <= mDistance)
		return;

	mContactArea = CSModelTools::GetContactAreaHertz( node->mRadius, mDistance );

	node->surfaceContactArea += mContactArea;
	*/
}

void 
CSInteractionJKR::interact( ModelElementLatticeNode * node, ModelElementVesselSphere * sphere)
{
	/*
	mContactArea = 0;

	if (mIs2D)
		mDistance = CSModelTools::GetDistance2D(node->position, sphere->position);
	else
		mDistance = CSModelTools::GetDistance3D(node->position, sphere->position);

	double delta = node->mRadius + sphere->mRadius - mDistance;

	if (delta < 0.0 && !mContactEstablished)
		return;

	double E_eff = 1.0
		/ ((1 - node->poissonRatio * node->poissonRatio)
			/ node->youngModulus
			+ (1 - sphere->poissonRatio * sphere->poissonRatio)
			/ sphere->youngModulus);

	double R_eff = sphere->mRadius * node->mRadius
		/ (sphere->mRadius + node->mRadius);

	double adhesion_constant = (*mpSingleBondEnergy) * (*mpAdhesionDensity);

	double dcrit = std::pow(
		std::pow(std::sqrt(3.)* M_PI * adhesion_constant / (8 * E_eff), 2) * R_eff,
		0.333333);

	if (delta < -dcrit)
		return;

	CSModelTools::GetJKRForce(E_eff, R_eff, delta, adhesion_constant, mForce, mContactArea, dcrit, 1);

	node->surfaceContactArea += mContactArea;
	sphere->surfaceContactArea += mContactArea;
	*/
}

void
CSInteractionJKR::interact(CellSpherical * cell, LatticeSpring * spring)
{
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
			double sDistance = distE;
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
					sDistance = itern1d;
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
					sDistance = itern2d;
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
				// calculate value of (JKR) force
				double E_eff = 1.0 / ((1 - cell->poissonRatio * cell->poissonRatio) / cell->youngModulus +
									  (1 - lspring->mPoisson * lspring->mPoisson) / lspring->mYoung);
				double R_eff = (sqrt(cell->mRadius * lspring->mRadius / (cell->mRadius + lspring->mRadius)) + sqrt(cell->mRadius)) * 
							   (sqrt(cell->mRadius * lspring->mRadius / (cell->mRadius + lspring->mRadius)) + sqrt(cell->mRadius)) / 4;

				double adhesion_constant = (*mpSingleBondEnergy) * (*mpAdhesionDensity);

				double dcrit = std::pow( 
					std::pow( std::sqrt(3.) * M_PI * adhesion_constant / (8 * E_eff), 2) * R_eff, 
					0.333333);

				double mForce = 0;

				CSModelTools::GetJKRForceSphereCylinder(E_eff, R_eff, sOverlap, adhesion_constant, mForce, mContactArea, dcrit, contactType);

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
			double sDistance = distE;
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
					sDistance = itern2d;
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
					sDistance = itern1d;
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
				// calculate value of (JKR) force
				double E_eff = 1.0 / ((1 - cell->poissonRatio * cell->poissonRatio) / cell->youngModulus +
					(1 - lspring->mPoisson * lspring->mPoisson) / lspring->mYoung);
				double R_eff = (sqrt(cell->mRadius * lspring->mRadius / (cell->mRadius + lspring->mRadius)) + sqrt(cell->mRadius)) *
					(sqrt(cell->mRadius * lspring->mRadius / (cell->mRadius + lspring->mRadius)) + sqrt(cell->mRadius)) / 4;

				double adhesion_constant = (*mpSingleBondEnergy) * (*mpAdhesionDensity);

				double dcrit = std::pow(
					std::pow(std::sqrt(3.) * M_PI * adhesion_constant / (8 * E_eff), 2) * R_eff,
					0.333333);

				double mForce = 0;

				CSModelTools::GetJKRForceSphereCylinder(E_eff, R_eff, sOverlap, adhesion_constant, mForce, mContactArea, dcrit, contactType);

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
	}
}

void
CSInteractionJKR::interact(CellSphericalPolar * cell, LatticeSpring * spring)
{
	CellSpherical * cell_1 = static_cast<CellSpherical *>(cell);

	interact(cell_1, spring);
}