////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSInteraction.h                                               //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-07-18 17:12:08                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_INTERACTION_H
#define CS_INTERACTION_H


#include <stdio.h>

#include "../Elements/ModelElement.h"

#include "../Cell/CellSpherical.h"
#include "../Cell/CellSphericalPolar.h"
#include "../Cell/CellTriangulated.h"
#include "../Elements/ModelElementBarrierTriangle.h"
#include "../Elements/ModelElementSphere.h"
#include "../Elements/ModelElementVesselSphere.h"
#include "../Elements/ModelElementHollowSphere.h"
#include "../Elements/CSVoxel.h"
// added by Jieling
#include "model/Lattice/ModelElementLatticeNode.h"
#include "model/Lattice/LatticeSpring.h"
#include "model/Lattice/LinearSpring.h"
#include "model/Lattice/RotationalSpring.h"

//! \brief Abstract class for Interaction functors.

//! CSInteraction is the abstract parent class for all interaction classes
//! in CellSys, i.e. for all interaction modules acting on two ModelElements or
//! cells.
//! CSInteraction is designed as a so-called functor, i.e. it implements the
//! void operator() method.
//! Since interactions are between two Cells (in general between two
//! ModelElements), the CSInteraction class declares functions interact() for
//! every permutation of two Cell types for every Cell type defined in CellSys
//! (and later on generalized to all ModelElements).
//! It has to be extended whenever another Cell type has been added to CellSys.
//!
//! The call to the interaction is done with the operator() taking the two
//! arguments of type ModelElement *, the superclass to all Cells and Elements.
//! This operator only needs to be defined here, and decides which interact()
//! method to call depending on the values of the ModelElement::mType members.
//!
//! A derived class has to implement any of these interact() specifications,
//! but can leave all definitions of irrelevant pairs of element types empty {}.
//!
//! Note:
//! A special technique is required to pass on quantities calculated in one
//! CSInteraction and needed in a subsequent CSInteraction acting on the same
//! pair of objects.
//! If a CSInteraction calculates a quantity in the course of computation, it
//! should store it in a member variable.
//! If a CSInteraction needs input from the outside, it should use a pointer
//! stored in a member variable and access the pointer's value during execution.
//! The programmer then has to provide the CSInteraction with the address of the
//! variable containing the value needed.
//! In this way the programmer has to connect output variable of one
//! CSInteraction to the pointer input variable of the subsequent CSInteraction:
//!   CSInteactionA a;  CSInteractionB b;
//!   b.mpInput = &a.mOutput;

class CSInteraction
{
public:
    virtual void interact( CellSpherical *, CellSpherical * ) =0;
    virtual void interact( CellSpherical *, CellSphericalPolar * ) =0;
    virtual void interact( CellSpherical *, CellTriangulated * ) =0;
    virtual void interact( CellSpherical *, ModelElementBarrierTriangle * ) =0;
    virtual void interact( CellSpherical *, ModelElementVesselSphere * ) =0;
    virtual void interact( CellSpherical *, ModelElementHollowSphere * ) =0;

    virtual void interact( CellSphericalPolar *, CellSphericalPolar * ) =0;
    virtual void interact( CellSphericalPolar *, ModelElementVesselSphere *) =0;
    virtual void interact( CellSphericalPolar *, ModelElementBarrierTriangle *) =0;

    virtual void interact( CellTriangulated *, CellTriangulated * ) =0;

    virtual void interact( ModelElementVesselSphere *, ModelElementVesselSphere * ) =0;
    virtual void interact( ModelElementVesselSphere *, ModelElementBarrierTriangle * ) =0;

    virtual void interact( CSVoxel *, CellSpherical *) =0;
    virtual void interact( CSVoxel *, ModelElementVesselSphere * ) =0;

	// added by Jieling
	virtual void interact( ModelElementLatticeNode *, CellSpherical *) = 0;
	virtual void interact( ModelElementLatticeNode *, CellSphericalPolar *) = 0;
	virtual void interact( ModelElementLatticeNode *, ModelElementBarrierTriangle *) = 0;
	virtual void interact( ModelElementLatticeNode *, ModelElementVesselSphere *) = 0;
	virtual void interact( CellSpherical *, LatticeSpring * ) = 0;
	virtual void interact( CellSphericalPolar *, LatticeSpring *) = 0;

    void interact( CellTriangulated *a, CellSpherical *b )
    { interact( b, a ); };
    void interact( CellSphericalPolar *a, CellSpherical *b )
    { interact( b, a ); };
    void interact( ModelElementVesselSphere *a, CellSphericalPolar *b)
    { interact( b, a ); }
    void interact( ModelElementBarrierTriangle *a, CellSpherical *b )
    { interact( b, a ); };
    void interact( ModelElementHollowSphere *a, CellSpherical *b )
    { interact( b, a ); };

    void interact( ModelElementVesselSphere *a , CellSpherical *b )
    { interact( b, a ); };
    void interact( ModelElementBarrierTriangle *a, ModelElementVesselSphere *b )
    { interact( b, a ); };

    void interact( CellSpherical *a, CSVoxel *b )
    { interact( b, a ); };
    void interact( ModelElementVesselSphere *a, CSVoxel *b )
    { interact( b, a ); };

	// added by Jieling
	void interact(CellSpherical *a, ModelElementLatticeNode *b)
	{ interact( b, a ); }
	void interact(CellSphericalPolar *a, ModelElementLatticeNode *b)
	{ interact(b, a); }
	void interact(ModelElementBarrierTriangle *a, ModelElementLatticeNode *b)
	{ interact(b, a); }
	void interact(ModelElementVesselSphere *a, ModelElementLatticeNode *b)
	{ interact(b, a); }

    void operator() ( ModelElement *obj1, ModelElement *obj2 )
    {
        switch (obj1->mType)
        {

        case ModelElement::TypeCellSpherical:
        {
            CellSpherical *cast_object1 = static_cast<CellSpherical *>(obj1);

            switch (obj2->mType)
            {

            case ModelElement::TypeCellSpherical:
            {
                CellSpherical * cast_object2 = static_cast<CellSpherical *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeCellSphericalPolar:
            {
                CellSphericalPolar * cast_object2 = static_cast<CellSphericalPolar *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeCellTriangulated:
            {
                CellTriangulated * cast_object2 = static_cast<CellTriangulated *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeBarrierTriangle:
            {
                ModelElementBarrierTriangle * cast_object2
                    = static_cast<ModelElementBarrierTriangle *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeVesselSphere:
            {
                ModelElementVesselSphere * cast_object2
                    = static_cast<ModelElementVesselSphere *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeHollowSphere:
            {
                ModelElementHollowSphere * cast_object2
                    = static_cast<ModelElementHollowSphere *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

			// added by Jieling
			case ModelElement::TypeCellLatticeNode:
			{
				ModelElementLatticeNode * cast_object2 
					= static_cast<ModelElementLatticeNode *>(obj2);
				this->interact(cast_object1, cast_object2);
				break;
			}

            case ModelElement::TypeSphere:
            case ModelElement::TypeTriangulated:
            case ModelElement::TypeVoxel:
            default:
            case ModelElement::TypeUnspecified:
                break;
            }

            break;
        }


        case ModelElement::TypeCellSphericalPolar:
        {
            CellSphericalPolar *cast_object1
                = static_cast<CellSphericalPolar *>(obj1);

            switch (obj2->mType)
            {

            case ModelElement::TypeCellSphericalPolar:
            {
                CellSphericalPolar * cast_object2
                    = static_cast<CellSphericalPolar *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeCellSpherical:
            {
                CellSpherical * cast_object2
                    = static_cast<CellSpherical *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeVesselSphere:
            {
                ModelElementVesselSphere * cast_object2
                    = static_cast<ModelElementVesselSphere *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeBarrierTriangle:
            {
                ModelElementBarrierTriangle * cast_object2
                    = static_cast<ModelElementBarrierTriangle *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeHollowSphere:
            {
                ModelElementHollowSphere * cast_object2
                    = static_cast<ModelElementHollowSphere *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

			// added by Jieling
			case ModelElement::TypeCellLatticeNode:
			{
				ModelElementLatticeNode * cast_object2
					= static_cast<ModelElementLatticeNode*>(obj2);
				this->interact(cast_object1, cast_object2);
				break;
			}

            case ModelElement::TypeTriangulated:
            case ModelElement::TypeCellTriangulated:
            case ModelElement::TypeVoxel:
            default:
            case ModelElement::TypeUnspecified:
                break;
            }

            break;
        }


        case ModelElement::TypeCellTriangulated:
        {
            CellTriangulated *cast_object1 = static_cast<CellTriangulated *>(obj1);

            switch (obj2->mType)
            {
            case ModelElement::TypeCellSpherical:
            {
                CellSpherical * cast_object2 = static_cast<CellSpherical *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeCellTriangulated:
            {
                CellTriangulated * cast_object2 = static_cast<CellTriangulated *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeSphere:
            case ModelElement::TypeTriangulated:
            case ModelElement::TypeBarrierTriangle:
            case ModelElement::TypeHollowSphere:
            case ModelElement::TypeVesselSphere:
            case ModelElement::TypeVoxel:
            default:
            case ModelElement::TypeUnspecified:
                break;
            }
            break;
        }


        case ModelElement::TypeSphere:
            break;


        case ModelElement::TypeTriangulated:
            break;


        case ModelElement::TypeBarrierTriangle:
        {
            ModelElementBarrierTriangle * cast_object1 =
                static_cast<ModelElementBarrierTriangle *>(obj1);

            switch (obj2->mType)
            {

            case ModelElement::TypeCellSpherical:
            {
                CellSpherical * cast_object2 = static_cast<CellSpherical *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeCellSphericalPolar:
            {
                CellSphericalPolar * cast_object2 = static_cast<CellSphericalPolar *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeVesselSphere:
            {
                ModelElementVesselSphere * cast_object2
                    = static_cast<ModelElementVesselSphere *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

			// added by Jieling
			case ModelElement::TypeCellLatticeNode:
			{
				ModelElementLatticeNode * cast_object2
					= static_cast<ModelElementLatticeNode *>(obj2);
				this->interact( cast_object1, cast_object2);
				break;
			}

            case ModelElement::TypeCellTriangulated:
            case ModelElement::TypeSphere:
            case ModelElement::TypeTriangulated:
            case ModelElement::TypeBarrierTriangle:
            case ModelElement::TypeHollowSphere:
            case ModelElement::TypeVoxel:
            default:
            case ModelElement::TypeUnspecified:
                break;
            }
            break;
        }


        case ModelElement::TypeHollowSphere:
        {
            ModelElementHollowSphere * cast_object1
                = static_cast<ModelElementHollowSphere *>(obj1);

            switch (obj2->mType)
            {

            case ModelElement::TypeCellSpherical:
            {
                CellSpherical * cast_object2
                    = static_cast<CellSpherical *>(obj2);
                this->interact( cast_object1,  cast_object2 );
                break;
            }

            case ModelElement::TypeCellSphericalPolar:
            {
                CellSphericalPolar * cast_object2
                    = static_cast<CellSphericalPolar *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            default:
            case ModelElement::TypeCellTriangulated:
            case ModelElement::TypeSphere:
            case ModelElement::TypeTriangulated:
            case ModelElement::TypeBarrierTriangle:
            case ModelElement::TypeHollowSphere:
            case ModelElement::TypeVoxel:
            case ModelElement::TypeUnspecified:
                break;
            }
            break;
        }

        case ModelElement::TypeVesselSphere:
        {
            ModelElementVesselSphere *cast_object1
                = static_cast<ModelElementVesselSphere *>(obj1);

            switch (obj2->mType)
            {

            case ModelElement::TypeCellSpherical:
            {
                CellSpherical * cast_object2 = static_cast<CellSpherical *>(obj2);
                this->interact( cast_object1, cast_object2);
                break;
            }

            case ModelElement::TypeCellSphericalPolar:
            {
                CellSphericalPolar * cast_object2 = static_cast<CellSphericalPolar *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeBarrierTriangle:
            {
                ModelElementBarrierTriangle * cast_object2
                    = static_cast<ModelElementBarrierTriangle *>(obj2);
                this->interact( cast_object1, cast_object2);
                break;
            }

            case ModelElement::TypeVesselSphere:
            {
                ModelElementVesselSphere * cast_object2
                    = static_cast<ModelElementVesselSphere *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

			// added by Jieling
			case ModelElement::TypeCellLatticeNode:
			{
				ModelElementLatticeNode * cast_object2
					= static_cast<ModelElementLatticeNode *>(obj2);
				this->interact( cast_object1, cast_object2 );
				break;
			}

            case ModelElement::TypeCellTriangulated:
            case ModelElement::TypeSphere:
            case ModelElement::TypeTriangulated:
            case ModelElement::TypeHollowSphere:
            case ModelElement::TypeVoxel:
            default:
            case ModelElement::TypeUnspecified:
                break;
            }
            break;
        }

        case ModelElement::TypeVoxel:
        {
            CSVoxel *cast_object1
                = static_cast<CSVoxel *>(obj1);

            switch ( obj2->mType )
            {
            case ModelElement::TypeVesselSphere:
            {
                ModelElementVesselSphere * cast_object2
                    = static_cast<ModelElementVesselSphere *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeCellSpherical:
            case ModelElement::TypeCellSphericalPolar:
            {
                CellSpherical * cast_object2
                    = static_cast<CellSpherical *>(obj2);
                this->interact( cast_object1, cast_object2 );
                break;
            }

            case ModelElement::TypeBarrierTriangle:
            case ModelElement::TypeCellTriangulated:
            case ModelElement::TypeSphere:
            case ModelElement::TypeTriangulated:
            case ModelElement::TypeHollowSphere:
            default:
            case ModelElement::TypeUnspecified:
                break;
            }
            break;
        }

		// added by Jieling
		case ModelElement::TypeCellLatticeNode:
		{
			ModelElementLatticeNode * cast_object1
				= static_cast<ModelElementLatticeNode *>(obj1);

			switch (obj2->mType)
			{
			case ModelElement::TypeCellSpherical:
			{
				CellSpherical * cast_object2 = static_cast<CellSpherical *>(obj2);
				this->interact(cast_object1, cast_object2);
				break;
			}

			case ModelElement::TypeCellSphericalPolar:
			{
				CellSphericalPolar * cast_object2 = static_cast<CellSphericalPolar *>(obj2);
				this->interact(cast_object1, cast_object2);
				break;
			}
			
			case ModelElement::TypeBarrierTriangle:
			{
				ModelElementBarrierTriangle * cast_object2 = static_cast<ModelElementBarrierTriangle *>(obj2);
				this->interact(cast_object1, cast_object2);
				break;
			}

			case ModelElement::TypeVesselSphere:
			{
				ModelElementVesselSphere * cast_object2 = static_cast<ModelElementVesselSphere *>(obj2);
				this->interact(cast_object1, cast_object2);
				break;
			}

			case ModelElement::TypeCellTriangulated:
			case ModelElement::TypeSphere:
			case ModelElement::TypeTriangulated:
			case ModelElement::TypeHollowSphere:
			case ModelElement::TypeVoxel:
			default:
			case ModelElement::TypeUnspecified:
				break;
			}
			break;
		}

        default:
        case ModelElement::TypeUnspecified:
            break;
        }
    }

	// added by Jieling
	void operator() (CellSpherical *obj1, LatticeSpring *obj2)
	{
		switch (obj1->mType)
		{
		case ModelElement::TypeCellSpherical:
		{
			this->interact(obj1, obj2);
			break;
		}

		case ModelElement::TypeCellSphericalPolar:
		{
			CellSphericalPolar * cast_obj1 = static_cast<CellSphericalPolar*>(obj1);
			this->interact(cast_obj1, obj2);
			break;
		}

		}
	}
};

#endif // CS_INTERACTION_H
