// added by Jieling
// For test of ECM as fibers
// 05/02/2018

#ifndef MODEL_ELEMENT_ECM_SPHERE_H
#define MODEL_ELEMENT_ECM_SPHERE_H

#include "ModelElementSphere.h"
#include "../tools/model/BoundingBoxList.h"
#include "../../gui/CSGLArena.h"
#include "../../gui/GLTools/CSGLBar.h"
#include "../tools/math/delaunay.h"
#include "../Elements/ModelElementVesselSphere.h"

#include <vector>
#include <sstream>

#include <Vector.h>

namespace H5 { class CompType; };

struct BoundingBox;
class ModelElementECMSphere;
class ModelElementECMEdge;

class ModelElementECMSphere : public ModelElementSphere
{
public:
	ModelElementECMSphere(double x = 0, double y=0, double z=0);
	~ModelElementECMSphere() {};

	static void HDF5DataFormat(H5::CompType &);
	static H5::CompType ParseHDF5DataFormat(H5::CompType & inputTypeDefinition,
		std::stringstream &,
		std::stringstream & warnings);

	enum ElementType {HSC, ECM};
	ElementType mElementType;

	unsigned int mIndex;

	double mYoungModulus;
	double mPoisonRatio;
	double mStiffness;
	double mPressure;
	double mConcentration;
	double mdConcentration;//concentration difference
	double mVolume;

	bool touched;

	Vector3f cLast;
	Vector3f balongV;

	std::vector<ModelElementECMSphere*> mepNeighbor;
	std::vector<ModelElementECMEdge*> mepEdges;

	BoundingBox * boundingBox();
	BoundingBox * getBoundingBox();

	ModelElementVesselSphere * vesselNeighbor;
	ModelElementVesselSphere * vesselNeighborLast; // neighbor of last step
	//
	ModelElementVesselSphere * alongVesselNeighbor; // for alongVessel
	//
	void setBoundingBox();
	void addNeighbor(ModelElementECMSphere* E);
	void addEdge(ModelElementECMEdge* E);
	void setElementType(ElementType mT);
	void updateVesselNeighbor();
	void alongVessel();
	void updateMigration(double x, double y, double z, double deltaT);

};

class ModelElementECMEdge : public ModelElement
{
public:
	ModelElementECMSphere *mStart;
	ModelElementECMSphere *mEnd;

	short int mIndex;

	double initialLength;
	double mYoungModulus;
	double mPoisonRatio;
	double mStiffness;
	double area; // the cross-section of the edge (assuming it is a cylinder)

	bool touched;

	ModelElementECMEdge(ModelElementECMSphere *p1, ModelElementECMSphere *p2);
	~ModelElementECMEdge() {};

	BoundingBox * boundingBox();
	void setInitialLength();
	void setStiffness();
	double getLength();
	//void setQuality(int slices, int stacks);
};

class ECMGraph
{
public:
	ECMGraph();
	~ECMGraph();
	void createECMNetwork();
	void addNode(ModelElementECMSphere*); // only for initializing the ecm spheres
	void removeNode(ModelElementECMSphere*);
	void setBoundingBoxList(BoundingBoxList *BoundingBoxList);
	void setArena(CSGLArena *Arena);
	void setBoundingBox();
	void updateParameters();

	// list of ECM spheres and edges
	std::vector<ModelElementECMSphere*> meNode;
	std::vector<ModelElementECMEdge*> meEdge;

	BoundingBoxList *mpBoundingBoxList;

	CSGLArena *mpArena;

	// bounding Box
	double mXmin, mXmax, mYmin, mYmax, mZmin, mZmax;

	// physical properties
	double defaultYoungModulus;
	double defaultPoissonRatio;
};

#endif //MODEL_ELEMENT_ECM_SPHERE_H