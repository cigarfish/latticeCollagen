///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ModelElementVesselGraph.cpp                                          //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-09-07 20:12:33                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "ModelElementVesselSphere.h"

#include "../../gui/GLTools/CSGLSphere.h"

#include <H5Cpp.h>
#include "../../tools/dataIO/hdf5/CSHDF5Helpers.h"


ModelElementVesselSphere::ModelElementVesselSphere(double x, double y, double z)
    : ModelElementSphere( x, y, z)
{
  this->mType =  ModelElement::TypeVesselSphere;
  this->mRadius = 0.1136;
  this->color.alpha = 1.0;
  this->color.red = 1.0;
  this->color.green = 0.;
  this->color.blue = 0.;

  // added by Jieling
  highlight = 0; // default, none

  mpGLObject = new CSGLSphere( &position, &color, &this->mRadius );
  this->setQuality(10,10);

  this->mVolume = 0;
  this->mConcentration = 0;

}


void ModelElementVesselSphere::setBoundingBox(){
  this->mBoundingBox.xmax = this->position.x + this->mRadius;
  this->mBoundingBox.ymax = this->position.y + this->mRadius;
  this->mBoundingBox.zmax = this->position.z + this->mRadius;
  this->mBoundingBox.xmin = this->position.x - this->mRadius;
  this->mBoundingBox.ymin = this->position.y - this->mRadius;
  this->mBoundingBox.zmin = this->position.z - this->mRadius;
}


BoundingBox *
ModelElementVesselSphere::boundingBox()
{
  setBoundingBox();
  return &mBoundingBox;
}


BoundingBox *
ModelElementVesselSphere::getBoundingBox()
{
  return &mBoundingBox;
}


void
ModelElementVesselSphere::HDF5DataFormat( H5::CompType & typeDefinition)
{
    H5_DOUBLE_ARRAY(positionType, 3);

    typeDefinition.insertMember( "Position", HOFFSET(ModelElementVesselSphere, position), positionType );
    typeDefinition.insertMember( "Radius", HOFFSET(ModelElementVesselSphere, mRadius), H5::PredType::NATIVE_DOUBLE );
    typeDefinition.insertMember( "Vessel Type", HOFFSET(ModelElementVesselSphere, mVesselType), H5::PredType::NATIVE_UINT );
}


H5::CompType
ModelElementVesselSphere::ParseHDF5DataFormat( H5::CompType & inputTypeDefinition,
                                               std::stringstream & /*errors*/,
                                               std::stringstream & warnings )
{
    H5_DOUBLE_ARRAY(positionType, 3);

    int numMembers = inputTypeDefinition.getNmembers();

    H5::CompType typeDefinition( sizeof(ModelElementVesselSphere) );

    for ( int i=0; i<numMembers; ++i )
    {
        std::string fieldName = inputTypeDefinition.getMemberName(i);

        if ( fieldName == "Position" )
        {
            typeDefinition.insertMember( "Position",
                                         HOFFSET(ModelElementVesselSphere, position),
                                         positionType );
        }
        else if ( fieldName == "Radius" )
        {
            typeDefinition.insertMember( "Radius",
                                         HOFFSET(ModelElementVesselSphere, mRadius),
                                         H5::PredType::NATIVE_DOUBLE );
        }
        else if ( fieldName == "Vessel Type" )
        {
            typeDefinition.insertMember( "Vessel Type",
                                         HOFFSET(ModelElementVesselSphere, mVesselType),
                                         H5::PredType::NATIVE_UINT );
        }
        else if ( fieldName == "Ignore" )
        {
        }
        else
        {
            warnings << "ModelElementVesselSphere::ParseHDF5DataFormat:  "
                     << "Unknown/Ignored field in HDF5 data set:\n"
                     << "\t\"" << fieldName << "\"\n";
        }
    }

    return typeDefinition;
}
