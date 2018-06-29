////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSInteractionFrictionMatrix.h                                 //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-07-18 18:56:13                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_INTERACTION_FRICTION_MATRIX_H
#define CS_INTERACTION_FRICTION_MATRIX_H

#include "CSInteraction.h"

#include <cstdlib>
#include <cstring>

class ModelElementSphere;

#define FRICTION_MATRIX_ALLOCATION_CHUNK 1024


class CSInteractionFrictionMatrix : public CSInteraction
{
public:
    CSInteractionFrictionMatrix();

    void interact( CellSpherical *,               CellSpherical * );
    void interact( CellSpherical *,               CellSphericalPolar * );
    void interact( CellSpherical *,               CellTriangulated * ) {};
    void interact( CellSpherical *,               ModelElementBarrierTriangle * );
    void interact( CellSpherical *,               ModelElementVesselSphere * );
    void interact( CellSpherical *,               ModelElementHollowSphere * );

    void interact( CellSphericalPolar *,          CellSphericalPolar * );
    void interact( CellSphericalPolar *,          ModelElementBarrierTriangle * );
    void interact( CellSphericalPolar *,          ModelElementVesselSphere * );

    void interact( CellTriangulated *,            CellTriangulated * ) {};

    void interact( ModelElementVesselSphere *,    ModelElementVesselSphere * ) {};
    void interact( ModelElementVesselSphere *,    ModelElementBarrierTriangle * );

    // unified friction for spherical ModelElements
    void interact( ModelElementSphere *, ModelElementSphere * );
    void interact( ModelElementSphere *, ModelElementBarrierTriangle * );

    void interact( CSVoxel *, CellSpherical * ) {};
    void interact( CSVoxel *, ModelElementVesselSphere * ) {};
    void interact( CSVoxel *, ModelElementSphere * ) {};

	// added by Jieling
	void interact( ModelElementLatticeNode *, CellSpherical *);
	void interact( ModelElementLatticeNode *, CellSphericalPolar *);
	void interact( ModelElementLatticeNode *, ModelElementBarrierTriangle *);
	void interact( ModelElementLatticeNode *, ModelElementVesselSphere *);
	void interact( CellSpherical*, LatticeSpring* );
	void interact( CellSphericalPolar*, LatticeSpring* );

    void Reset() { mFrictionMatrixCounter = 0; };

    // input quantities
    // from the Hertz force calculation:
    double * mpDistance;
    double * mpContactArea;
    // parameters from the model:
    double * mpGammaCellsParallel;
    double * mpGammaCellsPerpendicular;

    // output:  The matrix calculated in the last call;
    double * mpFrictionMatrix;

    // the model's friction matrix container:
    double * mpFrictionMatrices;

    bool mIs2D;


private:
    unsigned long mFrictionMatrixCounter;
    unsigned long mFrictionMatricesAllocationSize;
 
    void reallocate()
    {
        if ( mFrictionMatrixCounter >= mFrictionMatricesAllocationSize )
        {
            mFrictionMatricesAllocationSize += FRICTION_MATRIX_ALLOCATION_CHUNK;
            mpFrictionMatrices =
                (double *)realloc( mpFrictionMatrices,
                                   9 * mFrictionMatricesAllocationSize * sizeof(double));
        }
    }
};

#endif // CS_INTERACTION_FRICTION_MATRIX_H
