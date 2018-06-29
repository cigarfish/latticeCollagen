// added by Jieling
// For test of ECM as fibers
// 05/02/2018

#include <fstream>
#include <iostream>
#include <math.h>

#include "ModelElementECMSphere.h"

#include "../../gui/GLTools/CSGLSphere.h"

#include <H5Cpp.h>
#include "../../tools/dataIO/hdf5/CSHDF5Helpers.h"


ModelElementECMSphere::ModelElementECMSphere(double x, double y, double z)
	: ModelElementSphere(x, y, z)
{
	this->mType = ModelElement::TypeECMSphere;
	this->mRadius = 0.0212 * 4; // virtually existed
	this->color.alpha = 1.0;
	this->color.red = 0.5;
	this->color.green = 0.5;
	this->color.blue = 0.5;

	vesselNeighbor = NULL;
	vesselNeighborLast = NULL;

	alongVesselNeighbor = NULL;

	cLast = this->position;
	balongV = this->position;

	touched = false;

	mpGLObject = new CSGLSphere(&position, &color, &this->mRadius);
	this->setQuality(10, 10);

	this->mVolume = 0;
	this->mConcentration = 0;
}

void ModelElementECMSphere::setBoundingBox()
{
	this->mBoundingBox.xmax = this->position.x + this->mRadius;
	this->mBoundingBox.ymax = this->position.y + this->mRadius;
	this->mBoundingBox.zmax = this->position.z + this->mRadius;
	this->mBoundingBox.xmin = this->position.x - this->mRadius;
	this->mBoundingBox.ymin = this->position.y - this->mRadius;
	this->mBoundingBox.zmin = this->position.z - this->mRadius;
}

BoundingBox *
ModelElementECMSphere::boundingBox()
{
	setBoundingBox();
	return &mBoundingBox;
}

BoundingBox *
ModelElementECMSphere::getBoundingBox()
{
	return &mBoundingBox;
}

void
ModelElementECMSphere::HDF5DataFormat(H5::CompType & typeDefinition)
{
	H5_DOUBLE_ARRAY(positionType, 3);

	typeDefinition.insertMember("Position", HOFFSET(ModelElementECMSphere, position), positionType);
	typeDefinition.insertMember("Radius", HOFFSET(ModelElementECMSphere, mRadius), H5::PredType::NATIVE_DOUBLE);
	typeDefinition.insertMember("Element Type", HOFFSET(ModelElementECMSphere, mElementType), H5::PredType::NATIVE_UINT);
}

H5::CompType
ModelElementECMSphere::ParseHDF5DataFormat(H5::CompType & inputTypeDefinition,
	std::stringstream & /*errors*/,
	std::stringstream & warnings)
{
	H5_DOUBLE_ARRAY(positionType, 3);

	int numMembers = inputTypeDefinition.getNmembers();

	H5::CompType typeDefinition(sizeof(ModelElementECMSphere));

	for (int i = 0; i<numMembers; ++i)
	{
		std::string fieldName = inputTypeDefinition.getMemberName(i);

		if (fieldName == "Position")
		{
			typeDefinition.insertMember("Position",
				HOFFSET(ModelElementECMSphere, position),
				positionType);
		}
		else if (fieldName == "Radius")
		{
			typeDefinition.insertMember("Radius",
				HOFFSET(ModelElementECMSphere, mRadius),
				H5::PredType::NATIVE_DOUBLE);
		}
		else if (fieldName == "Element Type")
		{
			typeDefinition.insertMember("Element Type",
				HOFFSET(ModelElementECMSphere, mElementType),
				H5::PredType::NATIVE_UINT);
		}
		else if (fieldName == "Ignore")
		{
		}
		else
		{
			warnings << "ModelElementECMSphere::ParseHDF5DataFormat:  "
				<< "Unknown/Ignored field in HDF5 data set:\n"
				<< "\t\"" << fieldName << "\"\n";
		}
	}

	return typeDefinition;
}

void ModelElementECMSphere::addNeighbor(ModelElementECMSphere *E)
{
	if (mElementType != ElementType::ECM) return;

	std::vector<ModelElementECMSphere*>::iterator findE
		= std::find(mepNeighbor.begin(), mepNeighbor.end(), E);
	if (findE != mepNeighbor.end())
	{
		std::cerr << "ECM sphere already exists in the neighbor sphere vector\n";
		return;
	}
	else
	{
		mepNeighbor.push_back(E);
	}
}

void ModelElementECMSphere::addEdge(ModelElementECMEdge *E)
{
	if (mElementType != ElementType::ECM) return;

	std::vector<ModelElementECMEdge*>::iterator findE
		= std::find(mepEdges.begin(), mepEdges.end(), E);
	if (findE != mepEdges.end())
	{
		std::cerr << "ECM edge already exists in the ECM edge vector\n";
		return;
	}
	else
	{
		mepEdges.push_back(E);
	}
}

void ModelElementECMSphere::setElementType(ElementType mT)
{
	switch (mT)
	{
	case HSC:
		this->color.red = 1.0;
		this->color.green = 0.0;
		this->color.blue = 0.0;
		break;
	case ECM:
		this->color.red = 0.5;
		this->color.green = 0.5;
		this->color.blue = 0.5;
		break;
	}
	this->mElementType = mT;
}

void ModelElementECMSphere::updateVesselNeighbor()
{
	if (vesselNeighbor != NULL)
	{
		int vnNSize = (int)vesselNeighbor->mvpNeighbor.size();
		double disVN = (position.x - vesselNeighbor->position.x) * (position.x - vesselNeighbor->position.x) +
					   (position.y - vesselNeighbor->position.y) * (position.y - vesselNeighbor->position.y) + 
					   (position.z - vesselNeighbor->position.z) * (position.z - vesselNeighbor->position.z);
		int vnIndex = -1;
		int vncIndex = -1; // only for sinusoid connecting to the CV, radius is alwasy > 0.2
		double disVNmin = 1000000.;
		double disVNCmin = 1000000.;
		for (int i = 0; i < vnNSize; i++)
		{
			if (vesselNeighbor->mvpNeighbor[i]->mVesselType != 3) continue; // sinusoid type: 3
			//
			double vnNx = vesselNeighbor->mvpNeighbor[i]->position.x;
			double vnNy = vesselNeighbor->mvpNeighbor[i]->position.y;
			double vnNz = vesselNeighbor->mvpNeighbor[i]->position.z;
			double disVNn = (position.x - vnNx) * (position.x - vnNx) + 
						    (position.y - vnNy) * (position.y - vnNy) + 
						    (position.z - vnNz) * (position.z - vnNz);
			if (vesselNeighbor->mvpNeighbor[i]->mRadius < 0.201) // regular sinusoid node
			{
				if (disVNmin > disVNn)
				{
					disVNmin = disVNn;
					vnIndex = i;
				}
			}
			else
			{
				if (disVNCmin > disVNn)
				{
					disVNCmin = disVNn;
					vncIndex = i;
				}
			}
		}
		if (disVNmin <= disVNCmin)
		{
			if (disVNmin < disVN) // change the bound vessel cell
			{
				cout << "	vesselNeighbor has been changed from " << vesselNeighbor->mIndex << " to ";
				vesselNeighbor->highlight = 0; // none
				vesselNeighbor = vesselNeighbor->mvpNeighbor[vnIndex];
				vesselNeighbor->highlight = 1; // ECM
				cout << vesselNeighbor->mIndex << endl;
				// add alongVessel here
				alongVessel();
			}
		}
		else
		{
			if (disVNCmin < disVN) 
			{
				alongVesselNeighbor = vesselNeighbor->mvpNeighbor[vncIndex];
				// need to keep the vesselNeighbor and move it backward towards the original vesselNeighbor
				double vr = vesselNeighbor->mRadius;
				// node1/2: two end points of the vessel
				double node1x = vesselNeighbor->position.x;
				double node1y = vesselNeighbor->position.y;
				double node1z = vesselNeighbor->position.z;
				double node2x = vesselNeighbor->mvpNeighbor[vncIndex]->position.x;
				double node2y = vesselNeighbor->mvpNeighbor[vncIndex]->position.y;
				double node2z = vesselNeighbor->mvpNeighbor[vncIndex]->position.z;
				// vector for backward
				double bvd = sqrt((node1x - node2x)*(node1x - node2x) + (node1y - node2y)*(node1y - node2y) + (node1z - node2z)*(node1z - node2z));
				double bvx = (node1x - node2x) / bvd;
				double bvy = (node1y - node2y) / bvd;
				double bvz = (node1z - node2z) / bvd;
				// iter: current position
				double iterx = position.x;
				double itery = position.y;
				double iterz = position.z;
				// Least squares
				double a = (node2x - node1x) * (node2x - node1x) +
					(node2y - node1y) * (node2y - node1y) +
					(node2z - node1z) * (node2z - node1z);
				double b = 2 * (node2x - node1x) * (node1x - iterx) +
					2 * (node2y - node1y) * (node1y - itery) +
					2 * (node2z - node1z) * (node1z - iterz);
				double ts = -b / 2 / a;
				// node3: the point with shortest distance to iter
				double node3x = node1x + ts * (node2x - node1x);
				double node3y = node1y + ts * (node2y - node1y);
				double node3z = node1z + ts * (node2z - node1z);
				// distance of node3 to node1
				double d3 = sqrt((node3x - node1x)*(node3x - node1x) + (node3y - node1y)*(node3y - node1y) + (node3z - node1z)*(node3z - node1z));
				// 
				double i3d = sqrt((iterx - node3x)*(iterx - node3x) + (itery - node3y)*(itery - node3y) + (iterz - node3z)*(iterz - node3z));
				double i3x = (iterx - node3x) / i3d;
				double i3y = (itery - node3y) / i3d;
				double i3z = (iterz - node3z) / i3d;
				// np: new position
				double npx = node3x + i3x * vr + bvx * (d3 - bvd / 2);
				double npy = node3y + i3y * vr + bvy * (d3 - bvd / 2);
				double npz = node3z + i3z * vr + bvz * (d3 - bvd / 2);
				//
				position.x = npx;
				position.y = npy;
				position.z = npz;
			}
		}
	}
}

void ModelElementECMSphere::alongVessel()
{
	// make the ECM sphere wrapped around the vessel with fixed distance
	double vr = vesselNeighbor->mRadius;
	int vnNSize = (int)vesselNeighbor->mvpNeighbor.size();
	int vnIndex = -1;
	double disVNmin = 1000000.;
	for (int i = 0; i < vnNSize; i++)
	{
		if (vesselNeighbor->mvpNeighbor[i]->mVesselType != 3) continue; // sinusoid type: 3
		if (vesselNeighbor->mvpNeighbor[i]->mRadius > 0.201) continue; // default sinusoid node radius: 0.2
		double vnNx = vesselNeighbor->mvpNeighbor[i]->position.x;
		double vnNy = vesselNeighbor->mvpNeighbor[i]->position.y;
		double vnNz = vesselNeighbor->mvpNeighbor[i]->position.z;
		double disVNn = (position.x - vnNx) * (position.x - vnNx) +
						(position.y - vnNy) * (position.y - vnNy) +
						(position.z - vnNz) * (position.z - vnNz);
		if (disVNmin > disVNn)
		{
			disVNmin = disVNn;
			vnIndex = i;
		}
	}
	/****************************************
	                o iter
                    |
			o--->---o--->---o
		node1   node3   node2 
	****************************************/
	alongVesselNeighbor = vesselNeighbor->mvpNeighbor[vnIndex];
	// node1/2: two end points of the vessel
	double node1x = vesselNeighbor->position.x;
	double node1y = vesselNeighbor->position.y;
	double node1z = vesselNeighbor->position.z;
	double node2x = vesselNeighbor->mvpNeighbor[vnIndex]->position.x;
	double node2y = vesselNeighbor->mvpNeighbor[vnIndex]->position.y;
	double node2z = vesselNeighbor->mvpNeighbor[vnIndex]->position.z;
	// iter: current position
	double iterx = position.x;
	double itery = position.y;
	double iterz = position.z;
	// Least squares
	double a = (node2x - node1x) * (node2x - node1x) + 
			   (node2y - node1y) * (node2y - node1y) + 
			   (node2z - node1z) * (node2z - node1z);
	double b = 2 * (node2x - node1x) * (node1x - iterx) + 
			   2 * (node2y - node1y) * (node1y - itery) + 
			   2 * (node2z - node1z) * (node1z - iterz);
	double ts = -b / 2 / a;
	// node3: the point with shortest distance to iter
	double node3x = node1x + ts * (node2x - node1x);
	double node3y = node1y + ts * (node2y - node1y);
	double node3z = node1z + ts * (node2z - node1z);
	// 
	double i3d = sqrt((iterx - node3x)*(iterx - node3x) + (itery - node3y)*(itery - node3y) + (iterz - node3z)*(iterz - node3z));
	double i3x = (iterx - node3x) / i3d;
	double i3y = (itery - node3y) / i3d;
	double i3z = (iterz - node3z) / i3d;
	// np: new position
	double npx = node3x + i3x * vr;
	double npy = node3y + i3y * vr;
	double npz = node3z + i3z * vr;
	double checkD = sqrt((position.x - npx)*(position.x - npx) + 
						 (position.y - npy)*(position.y - npy) + 
						 (position.z - npz)*(position.z - npz));
	if (checkD > 0.1)
	{
		cout << "		-> Too large displacement after rotation along the vessel for ECM bound to vessel " << vesselNeighbor->mIndex << ", dist: " << checkD << endl;
		cout << "			-> lastC: " << position.x << ", " << position.y << ", " << position.z << "; cC: " << npx << ", " << npy << ", " << npz << endl;
		cout << "			-> vnC: " << node1x << ", " << node1y << ", " << node1z << "; vNnC: " << node2x << ", " << node2y << ", " << node2z << endl;
	}
	//	cout << "			-> before alongVessel: " << position.x << ", " << position.y << ", " << position.z << endl;
	position.x = npx;
	position.y = npy;
	position.z = npz;
	//	cout << "			-> after alongVessel: " << position.x << ", " << position.y << ", " << position.z << endl;
}

void ModelElementECMSphere::updateMigration(double x, double y, double z, double deltaT)
{
	switch (this->mElementType)
	{
	case HSC:
		// migrate along the vessels 
		// TODO
		break;
	case ECM:
		// migrate along the vessels
		double mDist = sqrt(x*x + y*y + z*z) * deltaT;
		if (mDist < 1e-10)
		{
			position.x += x;
			position.y += y;
			position.z += z;
			return;
		}
		double vnx = vesselNeighbor->position.x;
		double vny = vesselNeighbor->position.y;
		double vnz = vesselNeighbor->position.z;
		int vnSize = vesselNeighbor->mvpNeighbor.size();
		if (vnSize == 0) return;
		double vnDistA = 0.;
		for (int i = 0; i < vnSize; i++)
		{
			double vnix = vesselNeighbor->mvpNeighbor[i]->position.x;
			double vniy = vesselNeighbor->mvpNeighbor[i]->position.y;
			double vniz = vesselNeighbor->mvpNeighbor[i]->position.z;
			double vnDist = sqrt((vnix - vnx)*(vnix - vnx) + (vniy - vny)*(vniy - vny) + (vniz - vnz)*(vniz - vnz));
			vnDistA += vnDist;
		}
		vnDistA /= vnSize;
		//cout << "			-> average vnDist: " << vnDistA << endl;
		double scale = mDist;
		if (mDist > vnDistA)
		{
			// too large migration
			std::cerr << "Warning: the ECM sphere migrates too fast!!" << endl;
			scale = vnDistA;
		}
		int nIndex = -1;
		double cosAmax = 0.;
		for (int i = 0; i < vnSize; i++)
		{
			if (vesselNeighbor->mvpNeighbor[i]->mVesselType != 3) continue; // sinusoid type: 3
			if (vesselNeighbor->mvpNeighbor[i]->mRadius > 0.201) continue; // default sinusoid node radius: 0.2
			double vnix = vesselNeighbor->mvpNeighbor[i]->position.x;
			double vniy = vesselNeighbor->mvpNeighbor[i]->position.y;
			double vniz = vesselNeighbor->mvpNeighbor[i]->position.z;
			double cosA = (vnix - vnx) * x + (vniy - vny) * y + (vniz - vnz) * z;
			cosA /= sqrt(x*x + y*y + z*z) * sqrt((vnix - vnx)*(vnix - vnx) + (vniy - vny)*(vniy - vny) + (vniz - vnz)*(vniz - vnz));
			if (cosA > cosAmax)
			{
				// find the vessel with the best alignment with the migration direction
				cosAmax = cosA;
				nIndex = i;
			}
		}
		if (nIndex >= 0)
		{
			double vnIndexx = vesselNeighbor->mvpNeighbor[nIndex]->position.x;
			double vnIndexy = vesselNeighbor->mvpNeighbor[nIndex]->position.y;
			double vnIndexz = vesselNeighbor->mvpNeighbor[nIndex]->position.z;
			double vNorm = sqrt((vnIndexx-vnx)*(vnIndexx-vnx) + (vnIndexy-vny)*(vnIndexy-vny) + (vnIndexz-vnz)*(vnIndexz-vnz));
			double vnNormx = (vnIndexx - vnx) / vNorm;
			double vnNormy = (vnIndexy - vny) / vNorm;
			double vnNormz = (vnIndexz - vnz) / vNorm;
			position.x += cosAmax * scale * vnNormx;
			position.y += cosAmax * scale * vnNormy;
			position.z += cosAmax * scale * vnNormz;
			//double checkD = sqrt( cosAmax * scale * cosAmax * scale) * sqrt(vnNormx*vnNormx + vnNormy*vnNormy + vnNormz*vnNormz);
			//if (checkD > 0.1) 
			//{ 
			//	cout << "		-> Too large displacement detected for ECM bound to vessel " << vesselNeighbor->mIndex << ", cosAmax: " << cosAmax << ", scale: " << scale << endl; 
			//}
			//cout << "		for ECM sphere_" << position.y << ": " << cosAmax * scale * vnNormx << ", " << cosAmax * scale * vnNormy << ", " << cosAmax * scale * vnNormz << endl;
			// update the bound vessel neighbor after position update
			balongV = position; // record the coordination after mapping the velocity on the the vessel
			updateVesselNeighbor();
			// update the position to wrap around the vessel
			//alongVessel();
		}
		break;
	}
}

ModelElementECMEdge::ModelElementECMEdge(ModelElementECMSphere *p1, ModelElementECMSphere *p2)
	: ModelElement(p1->position.x, p1->position.y, p1->position.z)
{
	mStart = p1;
	mEnd = p2;
	touched = false;
	setInitialLength();
	mpGLObject = new CSGLBar(&p1->position,&p2->position, &color);
}

BoundingBox *
ModelElementECMEdge::boundingBox()
{
	mBoundingBox.xmin = std::min(mStart->position.x, mEnd->position.x);
	mBoundingBox.xmax = std::max(mStart->position.x, mEnd->position.x);
	mBoundingBox.ymin = std::min(mStart->position.y, mEnd->position.y);
	mBoundingBox.ymax = std::max(mStart->position.y, mEnd->position.y);
	mBoundingBox.zmin = std::min(mStart->position.z, mEnd->position.z);
	mBoundingBox.zmax = std::max(mStart->position.z, mEnd->position.z);
	return &mBoundingBox;
}

void ModelElementECMEdge::setInitialLength()
{
	if (mStart != NULL && mEnd != NULL && mStart != mEnd)
	{
		double x1 = mStart->position.x;
		double y1 = mStart->position.y;
		double z1 = mStart->position.z;
		double x2 = mEnd->position.x;
		double y2 = mEnd->position.y;
		double z2 = mEnd->position.z;
		double dist = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
		initialLength = dist;
		// set stiffness
		setStiffness();
	}
}

void ModelElementECMEdge::setStiffness()
{
	if (mStart != NULL && mEnd != NULL)
	{
		double aY = (mStart->mYoungModulus + mEnd->mYoungModulus) / 2;
		double aP = (mStart->mPoisonRatio + mEnd->mPoisonRatio) / 2;
		double aR = (mStart->mRadius + mEnd->mRadius) / 2;
		mStiffness = aY * aR * aR * 3.14159265 / initialLength;
	}
}

double ModelElementECMEdge::getLength()
{
	if (mStart != NULL && mEnd != NULL)
	{
		double x1 = mStart->position.x;
		double y1 = mStart->position.y;
		double z1 = mStart->position.z;
		double x2 = mEnd->position.x;
		double y2 = mEnd->position.y;
		double z2 = mEnd->position.z;
		double dist = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
		return dist;
	}
	return -1.;
}

ECMGraph::ECMGraph()
{
	this->mpBoundingBoxList = NULL;
}

ECMGraph::~ECMGraph()
{
	for (int i = 0; i < (int)meNode.size();i++)
	{
		delete meNode[i];
		meNode[i] = NULL;
	}
	for (int i = 0; i < (int)meEdge.size(); i++)
	{
		delete meEdge[i];
		meEdge[i] = NULL;
	}
}

void ECMGraph::setBoundingBoxList(BoundingBoxList *BoundingBoxList)
{
	this->mpBoundingBoxList = BoundingBoxList;

	if (mpBoundingBoxList)
		for (auto node : meNode)
			mpBoundingBoxList->add(node);
}

void ECMGraph::setBoundingBox()
{
	this->mXmin = 0.;
	this->mXmax = 0.;
	this->mYmin = 0.;
	this->mYmax = 0.;
	this->mZmin = 0.;
	this->mZmax = 0.;

	for (unsigned int i = 0; i < this->meNode.size();i++)
	{
		if (this->mXmin > this->meNode[i]->position.x)
			this->mXmin = this->meNode[i]->position.x;
		if (this->mXmax < this->meNode[i]->position.x)
			this->mXmax = this->meNode[i]->position.x;
		if (this->mYmin > this->meNode[i]->position.y)
			this->mYmin = this->meNode[i]->position.y;
		if (this->mYmax < this->meNode[i]->position.y)
			this->mYmax = this->meNode[i]->position.y;
		if (this->mZmin > this->meNode[i]->position.z)
			this->mZmin = this->meNode[i]->position.z;
		if (this->mZmax < this->meNode[i]->position.z)
			this->mZmax = this->meNode[i]->position.z;
	}
}

void ECMGraph::setArena(CSGLArena *Arena)
{
	mpArena = Arena;
}

void ECMGraph::removeNode(ModelElementECMSphere *E)
{
	if (E->mElementType != ModelElementECMSphere::ECM) return;

	std::vector<ModelElementECMSphere *>::iterator found
		= std::find(meNode.begin(), meNode.end(), E);

	if (found != meNode.end())
		meNode.erase(found);
	else
		std::cerr << "Error in ModelElementECMSphere: Attempted to remove unregistered ECM sphere from ECMSphere vector\n";

	// remove the point from the bounding box
	mpBoundingBoxList->remove(E);

	// remove the point from the arena
	mpArena->removeObject(E->GLObject());

	for (int i = 0;i < (int)E->mepEdges.size();i++)
	{
		// edge to remove: E->mepEdges[i]
		ModelElementECMEdge *edgeI = E->mepEdges[i];
		std::vector<ModelElementECMEdge*>::iterator foundE
			= std::find(meEdge.begin(), meEdge.end(), edgeI);
		if (foundE != meEdge.end())
			meEdge.erase(foundE);
		else
			std::cerr << "Error in ModelElementECMSphere: Attempt to remove unregistered ECM edge from ECMEdge vector\n";
		
		// also to remove the edge from the other end-point
		ModelElementECMSphere *eOther = edgeI->mStart;
		if (eOther == E) eOther = edgeI->mEnd;
		for (int j = 0;j < (int)eOther->mepEdges.size();j++)
		{
			std::vector<ModelElementECMEdge*>::iterator foundOther
				= std::find(eOther->mepEdges.begin(), eOther->mepEdges.end(),edgeI);
			if (foundOther != eOther->mepEdges.end())
				eOther->mepEdges.erase(foundOther);
			else
				std::cerr << "Error in Modcd elElementECMSphere: Attempt to remove unregistered ECM edge from ECMEdge vector of the other endpoint\n";
		}

		// remove the edge from the arena
		mpArena->removeObject(edgeI->GLObject());

		// delete ECM edge
		edgeI->~ModelElementECMEdge();
	}

	// delete ECM sphere
	E->~ModelElementECMSphere();
}

void ECMGraph::addNode(ModelElementECMSphere *E)
{
	if (E->mElementType != ModelElementECMSphere::ECM) return;
	meNode.push_back(E);
}

void ECMGraph::updateParameters()
{
	for (auto sphere : meNode)
	{
		sphere->poissonRatio = defaultPoissonRatio;
		sphere->youngModulus = defaultYoungModulus;
	}
}

void ECMGraph::createECMNetwork()
{
	// delaunay triangulation of the ECM spheres
	CD3DW *Del = new CD3DW();
	int size = (int)meNode.size();
	for (int i = 0; i < size; i++)
	{
		// add ECM spheres into the system one-by-one
		Del->add(meNode[i]->position.x,
				 meNode[i]->position.y,
				 meNode[i]->position.z,
				 meNode[i]->mRadius,
				 i,
				 &meNode[i]);
	}
	CComplex *Complex = new CComplex(Del);
	// add ECM edges
	int Pn = (int)Complex->P.size();
	for (int i = 0; i < Pn; i++)
	{
		CSimplex* p = Complex->P[i];
		int size = p->size();
		if (size == 2) // edges
		{
			C1Simplex *p1 = (C1Simplex*)p;
			int n1 = p1->A->n();
			int n2 = p1->B->n();
			Vector3f v12 = meNode[n1]->position - meNode[n2]->position;
			double dist12 = v12.Norm();
			if (dist12 > 2.5) continue;
			// add ECM edges
			ModelElementECMEdge *edgeN = new ModelElementECMEdge(meNode[n1],meNode[n2]);
			meNode[n1]->addEdge(edgeN);
			meNode[n2]->addEdge(edgeN);
			meNode[n1]->addNeighbor(meNode[n2]);
			meNode[n2]->addNeighbor(meNode[n1]);
			meEdge.push_back(edgeN);
			// add fiber object into the arena
			mpArena->addObject(edgeN->GLObject());
			//cout << " Node " << meNode[n1]->vesselNeighbor->mIndex << " - Node " << meNode[n2]->vesselNeighbor->mIndex << " with distance " << dist12 << endl;
			//cout << "	Node " << meNode[n1]->vesselNeighbor->mIndex << ": " << meNode[n1]->position.x << ", " << meNode[n1]->position.y << ", " << meNode[n1]->position.z << endl;
			//cout << "	Node " << meNode[n2]->vesselNeighbor->mIndex << ": " << meNode[n2]->position.x << ", " << meNode[n2]->position.y << ", " << meNode[n2]->position.z << endl;
		}
	}
	// delete delaunay triangulation
	Del->~CD3DW();
	// delete delaunay simplices
	Complex->~CComplex();
}