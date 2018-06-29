////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSInteractionJKR.h                                            //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-08-21 16:41:14                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_INTERACTION_JKR_H
#define CS_INTERACTION_JKR_H

#include "CSInteraction.h"

#include "../BasicDatatypes/Vector.h"

//! \brief Interaction based on the Johnson-Kendall-Roberts (JKR) model

//! This is the implementation of the JKR interaction model.
//! It is closely related to CSInteractionHertz, which calculates
//! the same quantities with a different formula.
//! See CSInteractionHertz.h for more.

//! The JKR force calculation calculates the distance, contact area, and
//! resulting force.
//! These quantities may be needed in subsequent CSInteraction calculations
//! (on the same cell-cell pair) and are therefore stored in member variables.
//! They are meant to be accessed by another CSInteraction via the address of
//! the variable.  I.e. after instantiation of the object of this and the
//! dependent CSInteraction one has to connect the two by passing on the
//! pointer of the variable to the corresponding member variable of the
//! dependent CSInteraction:
//!   CSInteractionHertz hertzForce;
//!   CSInteractionDependent dependentCalculation;
//!   dependentCalculation.mpDistance = &hertzForce.mDistance;
//!
//! The rule is:
//!   pointers for input of already calculated quantities,
//!   variables for output of quantities to be calculated.
//!
//! There will be run-time errors if a dependent CSInteraction is not connected.
//! properly.  Test before committing changes.


class CSInteractionJKR : public CSInteraction
{
public:
    CSInteractionJKR();

    void interact( CellSpherical *,            CellSpherical * );
    void interact( CellSpherical *,            CellSphericalPolar * ) {};
    void interact( CellSpherical *,            CellTriangulated * ) {};
    void interact( CellSpherical *,            ModelElementBarrierTriangle * );
    void interact( CellSpherical *,            ModelElementVesselSphere * );
    void interact( CellSpherical *,            ModelElementHollowSphere * );

    void interact( CellSphericalPolar *,       CellSphericalPolar * );
    void interact( CellSphericalPolar *,       ModelElementVesselSphere *);
    void interact( CellSphericalPolar *,       ModelElementBarrierTriangle *);

    void interact( CellTriangulated *,         CellTriangulated * ) {};

    void interact( ModelElementVesselSphere *, ModelElementVesselSphere * );
    void interact( ModelElementVesselSphere *, ModelElementBarrierTriangle * );

    void interact( ModelElementSphere *, ModelElementSphere *);
    void interact( ModelElementSphere *, ModelElementBarrierTriangle *);
    void interact( ModelElementSphere *, ModelElementHollowSphere * );

    void interact( CSVoxel *, CellSpherical * ) {};
    void interact( CSVoxel *, ModelElementVesselSphere * ) {};

	// added by Jieling
	void interact(ModelElementLatticeNode *, CellSpherical *);
	void interact(ModelElementLatticeNode *, CellSphericalPolar *);
	void interact(ModelElementLatticeNode *, ModelElementBarrierTriangle *);
	void interact(ModelElementLatticeNode *, ModelElementVesselSphere *);
	void interact(CellSpherical *, LatticeSpring *);
	void interact(CellSphericalPolar *, LatticeSpring *);

    // input quantities, therefore pointers
    double * mpSingleBondEnergy;
    double * mpAdhesionDensity;
    bool     mContactEstablished;

    // calculated by this functor
    Vector3f mForceVector;
    double mDistance;
    double mContactArea;
    double mForce;

    bool mIs2D;

private:
    inline void applyForce( ModelElement *, ModelElement * );
};


void
CSInteractionJKR::applyForce( ModelElement *obj1, ModelElement *obj2 )
{
    double force_ij_absolute = fabs(mForce);

    // Store cumulated forces
    obj1->lastForce += mForce;
    obj2->lastForce += mForce;

    // Store cumulated absolute forces
    obj1->accumulatedForceAbsolute += force_ij_absolute;
    obj2->accumulatedForceAbsolute += force_ij_absolute;

    // Store cumulated absolute forces
    obj1->accumulatedPressure += force_ij_absolute/mContactArea;
    obj2->accumulatedPressure += force_ij_absolute/mContactArea;

    if (mIs2D)
    {
        mForceVector.Set(obj2->position.x - obj1->position.x,
                         obj2->position.y - obj1->position.y,
                         0);
    }
    else
    {
        mForceVector.Set(obj2->position.x - obj1->position.x,
                         obj2->position.y - obj1->position.y,
                         obj2->position.z - obj1->position.z);
    }

    // Prelim: These methods could be optimized for 2D
    mForceVector.Normalize();
    mForceVector.Multiply(mForce);

    obj1->directedForce.Add( -mForceVector.x, -mForceVector.y, -mForceVector.z);
    obj2->directedForce.Add(  mForceVector.x,  mForceVector.y,  mForceVector.z);
}


#endif // CS_INTERACTION_JKR_H
