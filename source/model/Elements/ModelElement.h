#ifndef MODEL_ELEMENT_H
#define MODEL_ELEMENT_H

#include "../BasicDatatypes/Vector.h"
#include "../BasicDatatypes/Color.h"
#include "../BasicDatatypes/BoundingBox.h"
#include "../BasicDatatypes/CSListContainer.h"
#include "../../gui/CSGLObject.h"

#include <vector>
#include <string>
#include <climits>
#include <unordered_map>

#define ID_UNSET ULONG_MAX

//#include "../BasicDatatypes/Lattice.h"

//class Lattice;
//struct latticeContent;

class BoundingBoxList;

//! Generic class from which all models are derived from
class ModelElement
{
public:

    //double currentRadius;

    enum Type {
        TypeUnspecified=0,
        TypeSphere,
        TypeCellSpherical,
        TypeTriangulated,
        TypeCellTriangulated,
        TypeBarrierTriangle,
        TypeHollowSphere,
        TypeVesselSphere,
		TypeECMSphere,
        TypeCellSphericalPolar,
        TypeVoxel,
		// new ones
		TypeCellCylindricalOriented,
		TypeCellLatticeNode
    };

    Type mType;


    ModelElement(double x, double y, double z, Type type=TypeUnspecified);
	//
	ModelElement(double x, double y, double z, unsigned long id,
		Type type = Type::TypeUnspecified);
	ModelElement(const ModelElement& element);
	ModelElement(ModelElement&& element);
	ModelElement& operator=(const ModelElement& element) = delete;
	ModelElement& operator=(ModelElement&& element) = delete;

	virtual ~ModelElement();

    //unsigned long mGlobalIndex;

  /*
   * region Spatial properties
   */

    //! Position of element in model space
    Vector3f position;

	// get type
	std::string getType() const;

    //! Reset quantities that will be calculated anew in each model time step.
    virtual void Reset()
        {
			this->lastForce = 0;
            this->accumulatedForceAbsolute = 0;
            this->lastContactArea = 0;
            this->accumulatedPressure = 0.;
            this->directedForce = 0;
            this->frictionCoefficient = 0;
            this->surfaceContactArea = 0;
            this->mNumContacts = 0;
            this->mContacts.swap(this->mEstablishedContacts);
            this->mContacts.clear();
            this->frictionMatrixEntryNumbers.clear();
        };

	std::unordered_map<std::string, double> attributes;

  /*
   * region Color
   */

    ARGBColor color;
    bool mVisible;

    //! Changes cell color
    void SetColor(float r, float g, float b, float alpha = 1.)
    { color.red = r; color.green = g; color.blue = b; color.alpha = alpha;};

    void SetAlpha( float alpha ) { color.alpha = alpha; };

  /*
   * region Biophysical properties
   */

	//! Young modulus E in dimensionless units
    double youngModulus;

	//! Poisson ratio nu
    double poissonRatio;



    // Mainly for visualisation
    double lastForce;
    double accumulatedForceAbsolute;
    double lastForceAbsolute;
    double lastContactArea;

    //! Variable to sum up pressure from interaction forces.
    //  This member will be reset to zero in Reset().
    double accumulatedPressure;

    //! Pressure of the last timestep.
    //  For visualization.
    double lastPressure;

    // Cumulated Force excerted on cell
    Vector3f directedForce;

    //! Friction coefficient (Accumulated)
    double frictionCoefficient;

	// for solver
	virtual Vector3f getProduct(const Vector3f& velocity) const;
	virtual void getProduct2(double * v_out, const double * v_in) const;
	virtual void setFriction();

    //! Surface contact area (used to calculate friction)
    double surfaceContactArea;

    //! The contact-free area (used to calculate friction with the ECM)
    //  = 4*Pi*currentRadius*currentRadius - surfaceContactArea;
    double freeSurfaceArea;

    //! How many cells overlap with this one.
    int mNumContacts;

	//! Element volume
	double volume;

	// TODO Height TEMP
	double currentHeight;

	// unique id for the simulation objects
	// TODO merge with SimulationObject
	unsigned long mSimObjId;

    //! A list of indices of overlapping cells within ModelCellsSphericals 'cells' vector.
    CSListContainer<unsigned long> mContacts;

    //! The list of contacts established during the last time step.  To be swapped with mContacts.
    CSListContainer<unsigned long> mEstablishedContacts;

    //! The list of elements with bounding boxes intersecting with this element.
    //  Since each intersection should be processed only, once we only save
    //  intersections with elements with an mGlobalIndex larger than
    //  this->mGlobalIndex;
    CSListContainer<unsigned long> mIntersectingList;

	// return the global index
	unsigned long globalIndex() const { return mGlobalIndex; }

    //! Container for the so-called friction matrices
    //! This container will store the entry of a matrix in the Model's array of friction matrices.
    //! The model has to allocate the space for these matrices and an interaction will fill in the matrices
    //! and hand over the entry's offset to this container.
    //!
    //!  These matrices are built up by the coefficient of the friction between this and one overlapping cell
    //!  times the projection matrix for the projection of the velocity vector onto the plane of overlap
    //!  (which is perpendicular to the distance vector between the two cells).
    std::vector<unsigned long> frictionMatrixEntryNumbers;

    virtual BoundingBox * boundingBox() =0;
    // wth does this not work with as const memeber??
    BoundingBox * getBoundingBox() { return & mBoundingBox; };

	// Add forces
	void ApplyForce(const Vector3f& force);
	void ApplyForce(const Vector3f& force,
		const double magnitude);

	void ApplyForce(const Vector3f& force,
		const Vector3f& forceDirection,
		const double magnitude);

	// Stress tensor
	BasicDatatypes::Tensor mVirialStressMatrix;

public:
    CSGLObject * GLObject() const { return mpGLObject; };

public:
	// Contact Detection, do not use this to identify simulation objects use ID!!
	unsigned long mGlobalIndex;

	// make this friend so it can write to mGlobalIndex
	friend BoundingBoxList;

//protected:
public:
    BoundingBox mBoundingBox;
    CSGLObject * mpGLObject;


    bool mStatic;//if Element not change the position, set mStatic to 1

};

#endif
