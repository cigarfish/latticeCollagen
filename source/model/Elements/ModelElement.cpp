
#include "ModelElement.h"
#include "model/BasicDatatypes/Utils.h"

ModelElement::ModelElement(double x, double y, double z, Type type)
    : mType( type ),
      mNumContacts(0)
{
	// Default position: Centered
	position.x = x;
	position.y = y;
	position.z = z;

	// Default color: White
    SetColor( 1., 1., 1., 1. );

    this->mStatic = 0;
    this->mVisible = 1;

}

static
std::unordered_map<std::underlying_type<ModelElement::Type>::type, std::string>
ModelElemTypeMap
{
	{ ModelElement::Type::TypeUnspecified, "TypeUnspecified" },
	{ ModelElement::Type::TypeSphere, "TypeSphere" },
	{ ModelElement::Type::TypeCellSpherical, "TypeCellTriangulated" },
	{ ModelElement::Type::TypeTriangulated, "TypeTriangulated" },
	{ ModelElement::Type::TypeCellTriangulated, "TypeCellTriangulated" },
	{ ModelElement::Type::TypeBarrierTriangle, "TypeBarrierTriangle" },
	{ ModelElement::Type::TypeHollowSphere, "TypeHollowSphere" },
	{ ModelElement::Type::TypeVesselSphere, "TypeVesselSphere" },
	{ ModelElement::Type::TypeCellSphericalPolar, "TypeCellSphericalPolar" },
	{ ModelElement::Type::TypeCellCylindricalOriented, "TypeCellCylindricalOriented" },
	{ ModelElement::Type::TypeCellLatticeNode, "TypeCellLatticeNode" }
};

ModelElement::ModelElement(double x, double y, double z, unsigned long id,
	Type type)
	: mType(type),
	mNumContacts(0),
	//is2D(false),
	mStatic(false),
	currentHeight(0.),
	//mStaticX(false),
	//mStaticY(false),
	//mStaticZ(false),
	//mStaticMask(1., 1., 1.),
	mVisible(true),
	poissonRatio(0.5),
	youngModulus(1000),
	lastForce(0),
	mSimObjId(id),
	//singleBondEnergy(0.),
	//adhesionDensity(0.),
	accumulatedForceAbsolute(0),
	lastForceAbsolute(0),
	lastContactArea(0),
	//mContactSize(0),
	accumulatedPressure(0),
	lastPressure(0),
	frictionCoefficient(0),
	surfaceContactArea(0),
	freeSurfaceArea(0),
	volume(0.),
	//directedForce(0),
	mpGLObject(nullptr),
	position(x, y, z),
	mGlobalIndex(ID_UNSET),
	attributes(),
	mEstablishedContacts(),
	mContacts(),
	mIntersectingList()
	//diagonalFrictionMatrixNumber(-1)
{
	// Default color: White
	SetColor(1., 1., 1., 1.);
}

ModelElement::ModelElement(const ModelElement& other)
	: mType(other.mType),
	mNumContacts(other.mNumContacts),
	//is2D(other.is2D),
	mStatic(other.mStatic),
	//mStaticX(other.mStaticX),
	//mStaticY(other.mStaticY),
	//mStaticZ(other.mStaticZ),
	mVisible(other.mVisible),
	volume(other.volume),
	currentHeight(other.currentHeight),
	poissonRatio(other.poissonRatio),
	youngModulus(other.youngModulus),
	lastForce(other.lastForce),
	mSimObjId(other.mSimObjId),
	//singleBondEnergy(other.singleBondEnergy),
	//adhesionDensity(other.adhesionDensity),
	accumulatedForceAbsolute(other.accumulatedForceAbsolute),
	lastForceAbsolute(other.lastForceAbsolute),
	lastContactArea(other.lastContactArea),
	//mContactSize(other.mContactSize),
	accumulatedPressure(other.accumulatedPressure),
	lastPressure(other.lastPressure),
	frictionCoefficient(other.frictionCoefficient),
	surfaceContactArea(other.surfaceContactArea),
	freeSurfaceArea(other.freeSurfaceArea),
	directedForce(other.directedForce),
	position(other.position),
	color(other.color),
	//mpGLObject(other.mpGLObject),
	mGlobalIndex(other.mGlobalIndex),
	//mStaticMask(other.mStaticMask),
	attributes(other.attributes),
	mEstablishedContacts(other.mEstablishedContacts),
	mContacts(other.mContacts),
	mIntersectingList(other.mIntersectingList)
	//diagonalFrictionMatrixNumber(other.diagonalFrictionMatrixNumber)
{
	// TODO make sure this is copyable!
	// if this is not null it will cause a segfault
	mpGLObject = nullptr;
}

ModelElement::ModelElement(ModelElement&& other)
	: mType(other.mType),
	mNumContacts(other.mNumContacts),
	//is2D(other.is2D),
	mStatic(other.mStatic),
	//mStaticX(other.mStaticX),
	//mStaticY(other.mStaticY),
	//mStaticZ(other.mStaticZ),
	mVisible(other.mVisible),
	volume(other.volume),
	currentHeight(other.currentHeight),
	poissonRatio(other.poissonRatio),
	youngModulus(other.youngModulus),
	lastForce(other.lastForce),
	mSimObjId(other.mSimObjId),
	//singleBondEnergy(other.singleBondEnergy),
	//adhesionDensity(other.adhesionDensity),
	accumulatedForceAbsolute(other.accumulatedForceAbsolute),
	lastForceAbsolute(other.lastForceAbsolute),
	lastContactArea(other.lastContactArea),
	//mContactSize(other.mContactSize),
	accumulatedPressure(other.accumulatedPressure),
	lastPressure(other.lastPressure),
	frictionCoefficient(other.frictionCoefficient),
	surfaceContactArea(other.surfaceContactArea),
	freeSurfaceArea(other.freeSurfaceArea),
	directedForce(other.directedForce),
	position(other.position),
	color(other.color),
	mpGLObject(other.mpGLObject),
	mGlobalIndex(other.mGlobalIndex),
	//mStaticMask(other.mStaticMask),
	attributes(other.attributes),
	mEstablishedContacts(other.mEstablishedContacts),
	mContacts(other.mContacts),
	mIntersectingList(other.mIntersectingList)
	//diagonalFrictionMatrixNumber(other.diagonalFrictionMatrixNumber)
{
	other.mpGLObject = nullptr;
}

ModelElement::~ModelElement()
{
	if (mpGLObject) delete mpGLObject;
}

std::string ModelElement::getType() const
{
	return ModelElemTypeMap[mType];
}

void ModelElement::ApplyForce(const Vector3f& force)
{
	double magnitude = Norm(force);
	directedForce += force;
	accumulatedForceAbsolute += magnitude;
}

void ModelElement::ApplyForce(const Vector3f& force, const double magnitude)
{
	directedForce += force;
	accumulatedForceAbsolute += magnitude;
}

void ModelElement::ApplyForce(const Vector3f& force,
	const Vector3f& forceDirection,
	const double magnitude)
{
	// Apply force
	directedForce += force;
	accumulatedForceAbsolute += magnitude;

	mVirialStressMatrix += outer(force, forceDirection);
}

Vector3f ModelElement::getProduct(const Vector3f& velocity) const
{
	return frictionCoefficient * velocity;
}

void ModelElement::getProduct2(double * v_out, const double * v_in) const
{
	//*v_out = frictionCoefficient.x * (*v_in);
	//*(v_out + 1) = frictionCoefficient.y * *(v_in + 1);
	//*(v_out + 2) = frictionCoefficient.z * *(v_in + 2);
	*v_out = frictionCoefficient * (*v_in);
	*(v_out + 1) = frictionCoefficient * *(v_in + 1);
	*(v_out + 2) = frictionCoefficient * *(v_in + 2);
}

void ModelElement::setFriction()
{}
