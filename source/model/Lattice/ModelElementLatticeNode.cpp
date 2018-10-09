////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  ModelElementLatticeNode.cpp                                   //
//                                                                            //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                //
//    Created:  2016-08-02                                                    //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "ModelElementLatticeNode.h"
#include "tools/model/CSModelTools.h"

// added by Jieling
#include <H5Cpp.h>
#include "../../tools/dataIO/hdf5/CSHDF5Helpers.h"

using namespace CSModelTools;

ModelElementLatticeNode::ModelElementLatticeNode(double x, double y,
                                                 double z, double r)
    : SimulationObject(), ModelElementSphere(x, y, z, r,
                                             ModelElement::TypeCellLatticeNode),
    density(1.),
    hard_magnitude(0.),
    hard_force_crit(0.),
    pulling(false),
    pulling_force(0.)
{
    //mYoung.setPointer(&youngModulus);
    mYoung = 450;
    //mYoung.setUnitVal("N/m^2");
    //mYoung.setXMLPath("YoungsModulus");
    //registerParameter(mYoung);

    // TODO this is really cumbersome! I need to somehow create a structure
    // containing all the required information!
    //mPoisson.setPointer(&poissonRatio);
    //mPoisson.setXMLPath("PoissonRatio");
    mPoisson = 0.4;
    //registerParameter(mPoisson);

    //mForceHardMag.setPointer(&hard_magnitude);
    mForceHardMag = 15;
    //registerParameter(mForceHardMag);

    //mForceHardCrit.setPointer(&hard_force_crit);
    mForceHardCrit = 0.8;
    //registerParameter(mForceHardCrit);

    //mDecayAble.setXMLPath("DecayAble");
    mDecayAble = true;
    //registerParameter(mDecayAble);

	// added by Jieling
	fixed = false;
	touched = false;
	fibreEnd = false;
	top = false;
	bottom = false;
	left = false;
	right = false;
	front = false;
	back = false;
	free = false;

	mMoment = 0.;
	mMomentI = 0.;
	mShearMoment = 0.;
	mShearMomentDirection.x = 0.;
	mShearMomentDirection.y = -1.;
	mShearMomentDirection.z = 0.;
	oCoord.x = x;
	oCoord.y = y;
	oCoord.z = z;

	mStrainTestForce = 0.;

    assert(mType == ModelElement::Type::TypeCellLatticeNode);
}

ModelElementLatticeNode::ModelElementLatticeNode(/*XMLNode& node,*/
                                                 std::stringstream& errors,
                                                 std::stringstream& warnings)
    : ModelElementLatticeNode(0, 0, 0, 0)
{
    load(/*node,*/ errors, warnings);
}

ModelElementLatticeNode::ModelElementLatticeNode()
    : ModelElementLatticeNode(0, 0, 0, 0)
{}

ModelElementLatticeNode::ModelElementLatticeNode(const ModelElementLatticeNode& element)
    : ModelElementSphere(element), SimulationObject(element),
    density(element.density),
    pulling(element.pulling), pulling_force(element.pulling_force)
{}

ModelElementLatticeNode::ModelElementLatticeNode(ModelElementLatticeNode&& element)
    : ModelElementSphere(element), SimulationObject(element),
    density(element.density),
    pulling(element.pulling), pulling_force(element.pulling_force)
{}

ModelElementLatticeNode::~ModelElementLatticeNode()
{}

void
ModelElementLatticeNode::HDF5DataFormat(H5::CompType & typeDefinition)
{
	H5_DOUBLE_ARRAY(positionType, 3);

	typeDefinition.insertMember("Position", HOFFSET(ModelElementLatticeNode, position), positionType);
	typeDefinition.insertMember("Radius", HOFFSET(ModelElementLatticeNode, mRadius), H5::PredType::NATIVE_DOUBLE);
	typeDefinition.insertMember("Element Type", HOFFSET(ModelElementLatticeNode, mType), H5::PredType::NATIVE_UINT);
	typeDefinition.insertMember("Young", HOFFSET(ModelElementLatticeNode, mYoung), H5::PredType::NATIVE_DOUBLE);
	typeDefinition.insertMember("Poisson", HOFFSET(ModelElementLatticeNode, mPoisson), H5::PredType::NATIVE_DOUBLE);
}

H5::CompType
ModelElementLatticeNode::ParseHDF5DataFormat(H5::CompType &inputTypeDefinition,
	std::stringstream & /*errors*/,
	std::stringstream & warnings)
{
	H5_DOUBLE_ARRAY(positionType, 3);

	int numMembers = inputTypeDefinition.getNmembers();

	H5::CompType typeDefinition(sizeof(ModelElementLatticeNode));

	for (int i = 0; i<numMembers; ++i)
	{
		std::string fieldName = inputTypeDefinition.getMemberName(i);

		if (fieldName == "Position")
		{
			typeDefinition.insertMember("Position",
				HOFFSET(ModelElementLatticeNode, position),
				positionType);
		}
		else if (fieldName == "Radius")
		{
			typeDefinition.insertMember("Radius",
				HOFFSET(ModelElementLatticeNode, mRadius),
				H5::PredType::NATIVE_DOUBLE);
		}
		else if (fieldName == "Element Type")
		{
			typeDefinition.insertMember("Element Type",
				HOFFSET(ModelElementLatticeNode, mType),
				H5::PredType::NATIVE_UINT);
		}
		else if (fieldName == "Young")
		{
			typeDefinition.insertMember("Young",
				HOFFSET(ModelElementLatticeNode, mYoung),
				H5::PredType::NATIVE_DOUBLE);
		}
		else if (fieldName == "Poisson")
		{
			typeDefinition.insertMember("Poisson",
				HOFFSET(ModelElementLatticeNode, mPoisson),
				H5::PredType::NATIVE_DOUBLE);
		}
		else if (fieldName == "Ignore")
		{
		}
		else
		{
			warnings << "ModelElementLatticeNode::ParseHDF5DataFormat:  "
				<< "Unknown/Ignored field in HDF5 data set:\n"
				<< "\t\"" << fieldName << "\"\n";
		}
	}

	return typeDefinition;
}

void ModelElementLatticeNode::addNeighbor(ModelElementLatticeNode *E)
{
	std::vector<ModelElementLatticeNode*>::iterator findE
		= std::find(nodes.begin(), nodes.end(), E);
	if (findE != nodes.end())
	{
		std::cerr << "ECM sphere already exists in the neighboring nodes" << std::endl;
		return;
	}
	else
	{
		nodes.push_back(E);
	}
}

void ModelElementLatticeNode::removeNeighbor(ModelElementLatticeNode *E)
{
	std::vector<ModelElementLatticeNode*>::iterator findE
		= std::find(nodes.begin(), nodes.end(), E);
	if (findE == nodes.end())
	{
		std::cerr << "ECM sphere does not exist in the neighboring nodes" << std::endl;
		return;
	}
	else
	{
		nodes.erase(findE);
	}
}

bool ModelElementLatticeNode::checkNeighbor(ModelElementLatticeNode *E)
{
	std::vector<ModelElementLatticeNode*>::iterator findE
		= std::find(nodes.begin(), nodes.end(), E);
	if (findE != nodes.end())
		return true;
	else
		return false;
}

void ModelElementLatticeNode::addSpring(LatticeSpring *E)
{
	std::vector<LatticeSpring*>::iterator findE
		= std::find(springs.begin(), springs.end(), E);
	if (findE != springs.end())
	{
		std::cerr << "ECM spring already exists in the neighboring springs" << std::endl;
		return;
	}
	else
	{
		springs.push_back(E);
	}
}

void ModelElementLatticeNode::removeSpring(LatticeSpring *E)
{
	std::vector<LatticeSpring*>::iterator findE
		= std::find(springs.begin(), springs.end(), E);
	if (findE == springs.end())
	{
		std::cerr << "ECM spring does not exist in the neighboring springs" << std::endl;
		return;
	}
	else
	{
		springs.erase(findE);
	}
}

bool ModelElementLatticeNode::checkSpring(LatticeSpring *E)
{
	std::vector<LatticeSpring*>::iterator findE
		= std::find(springs.begin(), springs.end(), E);
	if (findE != springs.end())
		return true;
	else
		return false;
}

LatticeSpring* ModelElementLatticeNode::locateSpring(ModelElementLatticeNode *E)
{
	LatticeSpring *locateS = NULL;
	for (int i = 0; i < springs.size(); i++)
	{
		if ((springs.at(i)->nodes().at(0) == this && springs.at(i)->nodes().at(1) == E) ||
			(springs.at(i)->nodes().at(1) == this && springs.at(i)->nodes().at(0) == E))
		{
			locateS = springs.at(i);
			break;
		}
	}
	return locateS;
}

void ModelElementLatticeNode::removeRspring(LatticeSpring *E)
{
	std::vector<LatticeSpring*>::iterator findE
		= std::find(Rsprings.begin(), Rsprings.end(), E);
	if (findE == Rsprings.end())
	{
		std::cerr << "ECM spring does not exist in the neighboring springs" << std::endl;
		return;
	}
	else
	{
		Rsprings.erase(findE);
	}
}

bool ModelElementLatticeNode::addFibre(ModelElementFibre *F)
{
	if (F == NULL) // empty pointer
		return false;
		
	std::vector<ModelElementFibre *>::iterator findF
		= std::find(fibres.begin(), fibres.end(), F);
	if (findF != fibres.end())
	{
		std::cerr << "	-> Fibre already bound to the node!" << std::endl;
		return false;
	}
	else
	{
		// set the fibre's end-node to bound state
		if (F->n1 == this)
			F->n1Bound = true;
		else if (F->n2 == this)
			F->n2Bound = true;

		fibres.push_back(F);
	}
	return true;
}

bool ModelElementLatticeNode::removeFibre(ModelElementFibre *F)
{
	// a fibre detach from the node
	if (F != NULL)
		std::cout << "	-> Removing a fibre " << F->n1->mGlobalIndex << "-" << F->n2->mGlobalIndex << " in the node " << mGlobalIndex << std::endl;
	else
		return false;

	std::vector<ModelElementFibre*>::iterator findF
		= std::find(fibres.begin(), fibres.end(), F);
	if (findF != fibres.end())
	{
		// set the fibre's end-node to unbound state
		if (F->n1 == this)
			F->n1Bound = false;
		else if (F->n2 == this)
			F->n2Bound = false;

		fibres.erase(findF);
	}
	else
	{
		std::cerr << "	-> Error in ModelElementLatticeNode: attempt to remove unregistered fibre!" << std::endl;
		return false;
	}

	return true;
}

bool ModelElementLatticeNode::checkFibre(ModelElementFibre *F)
{
	std::vector<ModelElementFibre*>::iterator findF
		= std::find(fibres.begin(), fibres.end(), F);
	if (findF != fibres.end())
		return true;
	else
		return false;
}

void ModelElementLatticeNode::load(/*XMLNode& node,*/ std::stringstream& errors,
                                   std::stringstream& warnings)
{
    // this loads the sim object idn
    SimulationObject::load(/*node,*/ errors, warnings);

    // HACK merge mSimObj with id of SimObject
    mSimObjId = id;

    //getXMLAttribute(node, "X",      position.x,    false);
    //getXMLAttribute(node, "Y",      position.y,    false);
    //getXMLAttribute(node, "Z",      position.z,    false);
    //getXMLAttribute(node, "R",      mRadius   ,    false);
    //getXMLAttribute(node, "Static", mStatic,       false);
    //getXMLAttribute(node, "StaticX", mStaticX,       false);
    //getXMLAttribute(node, "StaticY", mStaticY,       false);
    //getXMLAttribute(node, "StaticZ", mStaticZ,       false);

    attributes["displacement"] = 0.;
    pulling_force = 0.;
    pulling = false;

    // update static status
    //SetStaticMask();

    // process potential tags, they're called tag#
    //int tags = node.nAttribute() - 5;
    //for (int tag = 1; tag < tags; ++tag)
    //{
    //    std::string name = "tag" + std::to_string(tag);

    //    if (!node.isAttributeSet(name.c_str()))
    //        continue;

    //    if (verbose)
    //        std::cout << "Found tag:" << node.getAttribute(name.c_str())
    //                  << std::endl;

    //    attributes[node.getAttribute(name.c_str())] = 1.0;
    //}
}

void ModelElementLatticeNode::print(std::ostream& stream) const
{
    stream << "Node(x="   << position
           << ", f="      << directedForce
           << ", |f|="    << lastForce
           << ", r="      << mRadius
           << ", id="     << id
           << ", E="      << youngModulus
           << ", v="      << poissonRatio
           << ", fh="     << hard_magnitude
           << ", fh_crit="<< hard_force_crit
           //<< ", stMask=" << mStaticMask
           << ", gamma= " << frictionCoefficient
           << ", decayAble=" << std::boolalpha << mDecayAble
           << ", static=" << std::boolalpha << mStatic << ").";
}

bool ModelElementLatticeNode::ready() const
{
    return std::isnormal(frictionCoefficient);
}

void ModelElementLatticeNode::Reset()
{
    ModelElementSphere::Reset();;

    // stress matrix
    mVirialStressMatrix.set(0);

    pulling = 0.;
    pulling_force = 0.;
}

void ModelElementLatticeNode::setFriction(double value)
{
    frictionCoefficient = SurfaceAreaSphere(mRadius) * value;
}

bool ModelElementLatticeNode::isAnchor() const
{
    return mStatic;
}

void ModelElementLatticeNode::decay(const double rate)
{
    density -= rate * density;
}


