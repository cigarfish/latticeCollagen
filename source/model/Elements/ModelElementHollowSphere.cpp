////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  ModelElementHollowSphere.cpp                                  //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-09-09 12:29:28                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES
#include "ModelElementHollowSphere.h"

#include "../../gui/GLTools/CSGLSphere.h"
#include "../BasicDatatypes/Vector.h"

#include <H5Cpp.h>
#include "../../tools/dataIO/hdf5/CSHDF5Helpers.h"

#include <cmath>

#include <iostream>


ModelElementHollowSphere::ModelElementHollowSphere(double x, double y, double z)
    : ModelElement( x, y, z, ModelElement::TypeHollowSphere ),
      mRadius( 5e-5 ),
      mPressure(0.),
      mStressResponse(1.)
{
    color.alpha = .1;
    mpGLObject = new CSGLSphere( &(this->position), &(this->color), &(this->mRadius) );
}


BoundingBox *
ModelElementHollowSphere::boundingBox()
{
    mBoundingBox.xmin = position.x - mRadius;
    mBoundingBox.xmax = position.x + mRadius;

    mBoundingBox.ymin = position.y - mRadius;
    mBoundingBox.ymax = position.y + mRadius;

    mBoundingBox.zmin = position.z - mRadius;
    mBoundingBox.zmax = position.z + mRadius;

    return &mBoundingBox;
}


void
ModelElementHollowSphere::Update()
{
    if (!mStatic)
    {
        lastForceAbsolute = accumulatedForceAbsolute;

        mPressure = lastForceAbsolute / ( 2*M_PI * mRadius*mRadius );

        mRadius += mStressResponse * mPressure;
    }
}


void
ModelElementHollowSphere::HDF5DataFormat( H5::CompType & typeDefinition )
{
    H5_DOUBLE_ARRAY( positionType, 3 );

    typeDefinition.insertMember( "Position", HOFFSET(ModelElementHollowSphere, position), positionType );
    typeDefinition.insertMember( "Radius", HOFFSET(ModelElementHollowSphere, mRadius), H5::PredType::NATIVE_DOUBLE );

    typeDefinition.insertMember( "Ignore", 0, H5::PredType::NATIVE_B8 );
}


H5::CompType
ModelElementHollowSphere::ParseHDF5DataFormat( H5::CompType & inputTypeDefinition,
                                               std::stringstream & /*errors*/,
                                               std::stringstream & warnings )
{
    int numMembers = inputTypeDefinition.getNmembers();

    H5::CompType typeDefinition( sizeof(ModelElementHollowSphere) );

        // build the data type:
    for ( int i=0; i<numMembers; ++i )
    {
        std::string fieldName = inputTypeDefinition.getMemberName(i);

        if ( fieldName == "Position" )
        {
            H5_DOUBLE_ARRAY(positionType, 3);
            typeDefinition.insertMember( "Position",
                                         HOFFSET(ModelElementHollowSphere, position),
                                         positionType );
        }
        else if ( fieldName == "Radius" )
        {
            typeDefinition.insertMember("Radius",
                                  HOFFSET(ModelElementHollowSphere, mRadius),
                                  H5::PredType::NATIVE_DOUBLE);
        }
        else if ( fieldName == "Ignore" )
        {}
        else
        {
            warnings << "ModelElementHollowsphere::ParseHDF5DataFormat:  "
                     << "Unknown/Ignored field in HDF5 data set:\n"
                     << "\t\"" << fieldName << "\"\n";
        }

    }
    return typeDefinition;
}
