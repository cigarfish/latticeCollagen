#include "ModelElementBarrierTriangle.h"
#include "../../../tools/math/mathematics.h"

#include "../../gui/GLTools/CSGLTriangle.h"

#include <sstream>
#include <algorithm>

#include <H5Cpp.h>
#include "../../tools/dataIO/hdf5/CSHDF5Helpers.h"


ModelElementBarrierTriangle::ModelElementBarrierTriangle(double x, double y, double z)
    : ModelElement( x, y, z, ModelElement::TypeBarrierTriangle )
{
  color.alpha = .3;
  mpGLObject = new CSGLTriangle( &position, &color, &mpPoints );

  mD = 0.;

  this->mStatic = 1;
}

ModelElementBarrierTriangle::~ModelElementBarrierTriangle(){}

BoundingBox * ModelElementBarrierTriangle::boundingBox()
{
  return &mBoundingBox;
}

void ModelElementBarrierTriangle::setPoint(double x, double y, double z, int index){

  this->mpPoints[index][0] = x;
  this->mpPoints[index][1] = y;
  this->mpPoints[index][2] = z;

}

double ModelElementBarrierTriangle::getD(){
  return this->mD;
}

double* ModelElementBarrierTriangle::getNormalVector(){
  return this->mpNormalVector;
}

double ModelElementBarrierTriangle::getNormalVector(int i){
  return this->mpNormalVector[i];
}

void ModelElementBarrierTriangle::setNormalVector(){
  crossProduct(mpPoints[0],mpPoints[1],mpPoints[2],mpNormalVector);

  mD = dotmult(mpNormalVector,mpPoints[0]);
}


void ModelElementBarrierTriangle::setBoundingBox(double epsilon){

  double x_min, y_min, z_min;
  double x_max, y_max, z_max;

  x_min = std::min(std::min(this->mpPoints[0][0],this->mpPoints[1][0]),this->mpPoints[2][0]);
  y_min = std::min(std::min(this->mpPoints[0][1],this->mpPoints[1][1]),this->mpPoints[2][1]);
  z_min = std::min(std::min(this->mpPoints[0][2],this->mpPoints[1][2]),this->mpPoints[2][2]);

  x_max = std::max(std::max(this->mpPoints[0][0],this->mpPoints[1][0]),this->mpPoints[2][0]);
  y_max = std::max(std::max(this->mpPoints[0][1],this->mpPoints[1][1]),this->mpPoints[2][1]);
  z_max = std::max(std::max(this->mpPoints[0][2],this->mpPoints[1][2]),this->mpPoints[2][2]);

  this->boundingBox()->xmin = x_min-epsilon;
  this->boundingBox()->ymin = y_min-epsilon;
  this->boundingBox()->zmin = z_min-epsilon;

  this->boundingBox()->xmax = x_max+epsilon;
  this->boundingBox()->ymax = y_max+epsilon;
  this->boundingBox()->zmax = z_max+epsilon;

}


void
ModelElementBarrierTriangle::HDF5DataFormat( H5::CompType & typeDefinition)
{
    const hsize_t dims[] = {3,3};
    H5::ArrayType pointGroup( H5::PredType::NATIVE_DOUBLE, 2, dims );

    typeDefinition.insertMember( "Points", HOFFSET(ModelElementBarrierTriangle, mpPoints), pointGroup );
    typeDefinition.insertMember( "Ignore", 0, H5::PredType::NATIVE_B8 );
}


H5::CompType
ModelElementBarrierTriangle::ParseHDF5DataFormat( H5::CompType & inputTypeDefinition,
                                                  std::stringstream & /*errors*/,
                                                  std::stringstream & warnings )
{
    const hsize_t dims[] = {3,3};
    H5::ArrayType pointGroup( H5::PredType::NATIVE_DOUBLE, 2, dims );

    int numMembers = inputTypeDefinition.getNmembers();

    H5::CompType typeDefinition( sizeof(ModelElementBarrierTriangle) );

    for ( int i=0; i<numMembers; ++i )
    {
        std::string fieldName = inputTypeDefinition.getMemberName(i);

        if ( fieldName == "Points" )
        {
            typeDefinition.insertMember( "Points",
                                         HOFFSET(ModelElementBarrierTriangle, mpPoints),
                                         pointGroup );
        }
        else if ( fieldName == "Ignore" )
        {
        }
        else
        {
            warnings << "ModelElementBarrierTriangle::ParseHDF5DataFormat:  "
                     << "Unknown/Ignored field in HDF5 data set:\n"
                     << "\t\"" << fieldName << "\"\n";
        }
    }

    return typeDefinition;
}
