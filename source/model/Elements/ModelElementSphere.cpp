#include <cstdio>

#include "ModelElementSphere.h"
#include "../../gui/GLTools/CSGLSphere.h"

ModelElementSphere::ModelElementSphere(double x, double y, double z)
    : ModelElement( x, y, z, ModelElement::TypeSphere )
{
	// added by Jieling
	mVelocity.x = 0; mVelocity.y = 0; mVelocity.z = 0;
  mpGLObject = new CSGLSphere( &position, &color, &mRadius );
}

ModelElementSphere::ModelElementSphere(double x, double y, double z, double r,
	double theta, double phi,
	unsigned long index, Type type)
	: ModelElement(x, y, z, index, type),
	mRadius(r),
	//maxoverlap(100.0),
	//criticalCompressionPeriod(0.0),
	//mForceThreshold(1000),
	//minPullDistance(1000),
	//mOverlapThreshold(0.6),
	//density("Density"),
	//deformation("Deformation"),
	//stress(volume),
	//nucleus_hard_force_crit(0.),
	//nucleus_hard_magnitude(0.),
	theta(theta),
	phi(phi)
{
	// HACK XXX
	if (type == ModelElement::TypeSphere)
	{
		mpGLObject = new CSGLSphere(&position, &color, &mRadius);
		boundingBox();
	}

	//setOrientation(theta, phi);
}

ModelElementSphere::ModelElementSphere(double x, double y, double z, double r,
	double theta, double phi, Type type)
	: ModelElementSphere(x, y, z, r, theta, phi, 0, type)
{}

ModelElementSphere::ModelElementSphere(double x, double y, double z, double r,
	double theta, double phi,
	unsigned long index)
	: ModelElementSphere(x, y, z, r, theta, phi, index, ModelElement::TypeSphere)
{}


ModelElementSphere::ModelElementSphere(double x, double y, double z, double r, unsigned long index, Type type)
	: ModelElementSphere(x, y, z, r, 90., 0, index, type)
{}

ModelElementSphere::ModelElementSphere(double x, double y, double z, double r,
	unsigned long index)
	: ModelElementSphere(x, y, z, r, 90., 0, index)
{}

ModelElementSphere::ModelElementSphere(double x, double y, double z, double r,
	Type type)
	: ModelElementSphere(x, y, z, r, 0, type)
{
	// added by Jieling
	mVelocity.x = 0; mVelocity.y = 0; mVelocity.z = 0;
	this->mRadius = r;
	this->mType = type;
	mpGLObject = new CSGLSphere(&position, &color, &mRadius);
}

ModelElementSphere::ModelElementSphere(double x, double y, double z, double r)
	: ModelElementSphere(x, y, z, r, 0)
{}

ModelElementSphere::~ModelElementSphere()
{}

ModelElementSphere::ModelElementSphere(const ModelElementSphere& element)
	: ModelElement(element),
	mRadius(element.mRadius),
	//maxoverlap(element.maxoverlap),
	//criticalCompressionPeriod(element.criticalCompressionPeriod),
	//mForceThreshold(element.mForceThreshold),
	//minPullDistance(element.minPullDistance),
	//mOverlapThreshold(element.mOverlapThreshold),
	//mOverlapThresholdMin(element.mOverlapThresholdMin),
	//density(element.density),
	//deformation(element.deformation),
	//nucleus_hard_force_crit(element.nucleus_hard_force_crit),
	//nucleus_hard_magnitude(element.nucleus_hard_magnitude),
	//stress(element.stress),
	theta(element.theta),
	phi(element.phi)
	//orientation(element.orientation)
{}

ModelElementSphere::ModelElementSphere(ModelElementSphere&& element)
	: ModelElement(element),
	mRadius(element.mRadius),
	//maxoverlap(element.maxoverlap),
	//criticalCompressionPeriod(element.criticalCompressionPeriod),
	//mForceThreshold(element.mForceThreshold),
	//minPullDistance(element.minPullDistance),
	//mOverlapThreshold(element.mOverlapThreshold),
	//mOverlapThresholdMin(element.mOverlapThresholdMin),
	//density(element.density),
	//deformation(element.deformation),
	//nucleus_hard_force_crit(element.nucleus_hard_force_crit),
	//nucleus_hard_magnitude(element.nucleus_hard_magnitude),
	//stress(element.stress),
	theta(element.theta),
	phi(element.phi)
	//orientation(element.orientation)
{}

BoundingBox *
ModelElementSphere::boundingBox()
{
  mBoundingBox.xmin = position.x - mRadius;
  mBoundingBox.xmax = position.x + mRadius;

  mBoundingBox.ymin = position.y - mRadius;
  mBoundingBox.ymax = position.y + mRadius;

  mBoundingBox.zmin = position.z - mRadius;
  mBoundingBox.zmax = position.z + mRadius;

  return &mBoundingBox;
}

void ModelElementSphere::setQuality( int slices, int stacks){
  CSGLSphere *sphere = (CSGLSphere*)this->mpGLObject;
  sphere->setQuality(slices, stacks);
}
