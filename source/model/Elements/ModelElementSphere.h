#ifndef MODEL_ELEMENT_SPHERE_H
#define MODEL_ELEMENT_SPHERE_H

#include "ModelElement.h"

//! Generic class from which all models are derived from
class ModelElementSphere : public ModelElement
{
public:

    double mRadius;
	// added by Jieling
	Vector3f mVelocity;

public:
    // Default constructor
    ModelElementSphere(double x, double y, double z);
	//
	ModelElementSphere(double x, double y, double z, double r, double theta,
		double phi, unsigned long index, Type type);
	ModelElementSphere(double x, double y, double z, double r, double theta,
		double phi, unsigned long index);
	ModelElementSphere(double x, double y, double z, double r, double theta,
		double phi, Type type);
	ModelElementSphere(double x, double y, double z, double r, unsigned long index, Type type);
	ModelElementSphere(double x, double y, double z, double r, unsigned long index);
	ModelElementSphere(double x, double y, double z, double r, Type type);
	ModelElementSphere(double x, double y, double z, double r);
	ModelElementSphere(const ModelElementSphere& element);
	ModelElementSphere(ModelElementSphere&& element);

	virtual ~ModelElementSphere();

    BoundingBox * boundingBox();
    void setQuality( int slices, int stacks);

	virtual void Reset()
	{
		ModelElement::Reset();
		// CAREFUL!!
		//maxoverlap = std::numeric_limits<double>::max();
		//minPullDistance = std::numeric_limits<double>::max();
		//density.Reset();
		//deformation.Reset();
		//stress.Reset();
	}

	// Radius
	//double mRadius;

protected:
	// phi is in (0, 2pi), theta (0, pi) i.e. theta is the inclination and phi
	// is the azimuth
	double theta, phi;
};

#endif
