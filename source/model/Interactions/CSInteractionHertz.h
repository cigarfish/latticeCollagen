////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSInteractionHertz.h                                          //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-07-18 17:25:52                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_INTERACTION_HERTZ_H
#define CS_INTERACTION_HERTZ_H

#include "CSInteraction.h"

#include "../BasicDatatypes/Vector.h"

//! \brief Hertz force interaction

//! This is the implementation of the Hertz force calculation as a CSInteraction
//! functor.  This means that this class implements only the operator(), with two
//! arguments.
//! The operators to be implemented are forced by the super-class CSInteraction,
//! which declares all of these as pure virtual methods.

//! The Hertz force calculation calculates the distance, contact area, and
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


class CSInteractionHertz : public CSInteraction
{
public:
    CSInteractionHertz();

    void interact( CellSpherical *,            CellSpherical * );
    void interact( CellSpherical *,            CellSphericalPolar * );
    void interact( CellSpherical *,            CellTriangulated * ) {};
    void interact( CellSpherical *,            ModelElementBarrierTriangle * );
    void interact( CellSpherical *,            ModelElementVesselSphere * );
    void interact( CellSpherical *,            ModelElementHollowSphere * );

    void interact( CellSphericalPolar *,       CellSphericalPolar * );
    void interact( CellSphericalPolar *,       ModelElementVesselSphere * );
    void interact( CellSphericalPolar *,       ModelElementBarrierTriangle * );

    void interact( CellTriangulated *,         CellTriangulated * ) {};

    void interact( ModelElementVesselSphere *, ModelElementVesselSphere * );
    void interact( ModelElementVesselSphere *, ModelElementBarrierTriangle * );

    void interact( ModelElementSphere *, ModelElementSphere * );
    void interact( ModelElementSphere *, ModelElementBarrierTriangle * );
    void interact( ModelElementSphere *, ModelElementVesselSphere * );

    inline void interact( CSVoxel *, CellSpherical * );
    inline void interact( CSVoxel *, ModelElementVesselSphere * );
    void interact( CSVoxel *, ModelElementSphere * );

	// added by Jieling
	void interact( ModelElementLatticeNode *, CellSpherical *);
	void interact( ModelElementLatticeNode *, CellSphericalPolar *);
	void interact( ModelElementLatticeNode *, ModelElementBarrierTriangle *);
	void interact( ModelElementLatticeNode *, ModelElementVesselSphere *);
	void interact( CellSpherical *, LatticeSpring *);
	void interact( CellSphericalPolar *, LatticeSpring *);

    // input quantities, therefore pointers
    double * mpSingleBondEnergy;
    double * mpAdhesionDensity;

    // calculated by this functor
    Vector3f mForceVector;
    double mDistance;
    double mContactArea;
    double mForce;

    bool mIs2D;
};


inline
void
CSInteractionHertz::interact( CSVoxel * voxel, CellSpherical * cell )
{
    ModelElementSphere * sphere = static_cast<ModelElementSphere *>( cell );
    CSInteractionHertz::interact( voxel, sphere );
};

inline
void
CSInteractionHertz::interact( CSVoxel * voxel, ModelElementVesselSphere * element )
{
    ModelElementSphere * sphere = static_cast<ModelElementSphere *>( element );
    CSInteractionHertz::interact( voxel, sphere );
};


#endif // CS_INTERACTION_HERTZ_H

