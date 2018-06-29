////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  LinarSpring.cpp                                               //
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

#include <vector>

#include "LinearSpring.h"
#include "ModelLattice.h"
#include "BiologyLink.h"
//#include "CSParameterContext.h"
#include "LinearFunction.h"
#include "NonLinearStiffeningFunction.h"

LinearSpring::LinearSpring() // DO NOT USE
    : LatticeSpring(), n1(nullptr), n2(nullptr), length0(0.)
{
    //l0.setXMLPath("l0");
    //registerParameter(l0);

    //Node1_Id.setXMLPath("Node1_Id");
    //registerParameter(Node1_Id);

    //Node2_Id.setXMLPath("Node2_Id");
    //registerParameter(Node2_Id);

    //mDensityDependency.setXMLPath("DensityDependency");
    mDensityDependency = true;
    //registerParameter(mDensityDependency);
	harden = false;
	hardenRatio = 0.5; // added by Jieling
	strainHarden = 0.2;
}

LinearSpring::LinearSpring(ModelElementLatticeNode *N1, ModelElementLatticeNode *N2)
	: LatticeSpring(), n1(N1), n2(N2)
{
	mDensityDependency = true;
	//if (mEquilibrium)
	//	length0 = Norm(n1->position - n2->position);
	//else
	//	length0 = l0;
	initialize();
}

LinearSpring::~LinearSpring()
{
    n1 = nullptr;
    n2 = nullptr;
}

/*LinearSpring::LinearSpring(XMLNode& springNode,
                           std::stringstream& errors,
                           std::stringstream& warnings)
    : LinearSpring()
{
    // can we get rid of this?
    load(springNode, errors, warnings);
}*/

LinearSpring::LinearSpring(const LinearSpring& other)
    : LatticeSpring(other), n1(other.n1), n2(other.n2),
    Node1_Id(other.Node1_Id), Node2_Id(other.Node2_Id),
    l0(other.l0), length0(other.length0)
{
	// added by Jieling
	initialize();
}

LinearSpring::LinearSpring(LinearSpring&& other) // DO NOT USE
    : LatticeSpring(), n1(other.n1), n2(other.n2),
    Node1_Id(other.Node1_Id), Node2_Id(other.Node2_Id),
    l0(other.l0), length0(other.length0)
{
    other.n1 = nullptr;
    other.n2 = nullptr;
}

void LinearSpring::print(std::ostream& stream) const
{
    stream << "LinearSpring(Id="<<id;
    stream << ", l0="<<l0
           <<", Node1_Id="<<Node1_Id<< ", Node2_Id="<<Node2_Id
           <<", position=" << position()
           << ", u=" << displacement()
           << ", " << *msStiffness << ").";
}

void LinearSpring::initialize()
{
	// added by Jieling
	this->mSpringType = LatticeSpring::TypeLinearSpring;

    //n1 = lattice->getNodeById(Node1_Id);
    //n2 = lattice->getNodeById(Node2_Id);

    if (!n1)
    {
        std::cerr << "Node: " /*<< Node1_Id*/ << " does not exist!" << std::endl;
        assert(false);
        return;
    }

    if (!n2)
    {
        std::cerr << "Node: " /*<< Node2_Id*/ << " does not exist!" << std::endl;
        assert(false);
        return;
    }

    if (n1 == n2)
    {
        std::cerr << "Spring( " << id << ") Nodes "
            //<< Node1_Id << " and " << Node2_Id
            << " are the same!" << std::endl;
        assert(false);
        return;
    }

	// added by Jieling
	color.red = 1.; color.green = 1.; color.blue = 1.; color.alpha = 1.; // default color: white
	mpGLObject = new CSGLBar(&n1->position, &n2->position, &color);

    auto OrientationAxis = n2->position - n1->position;

    if (mEquilibrium)
        length0 = Norm(OrientationAxis);
    else
        length0 = l0;

	initialL0 = length0;

	initialLY = std::abs(n2->position.y - n1->position.y);
	if (initialLY == 0) initialLY = 0.000001;

	curStrain = 0.;

	strainHarden = 0.15;

	msStiffness = new LinearStrainFunction();
	double aR = (n1->mRadius + n2->mRadius)/2;
	double aY = (n1->mYoung + n2->mYoung)/2;
	mYoung = aY;
	mPoisson = (n1->mPoisson + n2->mPoisson) / 2;
	mRadius = 0.0212 * 4; // for test

	((LinearStrainFunction*)msStiffness)->setK(aY * aR * aR * 3.14159265 / length0);

	// initiate unbinding rate
	n1_k_off = 0.;
	n2_k_off = 0.;
	n1_unbound = false;
	n2_unbound = false;

    if (!ready() || !isRelaxed())
    {
        std::cout << *this << " is not ready!" << std::endl;
        // Print node info
        std::cout << "Node1=" << n1 << ", " << *n1 << std::endl;
        std::cout << "Node2=" << n2 << ", " << *n2 << std::endl;

        assert(false);
    }
}

void LinearSpring::update(double timeStep)
{
    // compute the orientation of the spring
    auto OrientationAxis = n2->position - n1->position;

    double length = Norm(OrientationAxis);

    double strain = 0.;
	double strainStaining = 0.;
    if (std::isnormal(length))
    {
        // normalize the OrientationAxis
        OrientationAxis /= length;
        strain = length - length0;
		strainStaining = length - initialL0;
    }
    else
    {
        strain = 0;
        std::cerr << "WARNING: " << std::endl;
    }

	curStrain = strain / length0;
	curStrainStaining = strainStaining / initialL0;

    // compute the force
	double strain0 = strain;
	double strain1 = 0.;
	if (strain0 > strainHarden * length0)
	{
		strain1 = strain0 - strainHarden * length0;
		strain0 = strainHarden * length0;
	}
    double forceAbs = ((*msStiffness)(strain0) + (*msStiffness)(strain1) * hardenRatio) / timeStep; // revised by Jieling

    // if this is enabled then we reduce stiffness due to node density
    if (mDensityDependency)
        forceAbs *= (n2->density + n1->density) / 2.;

    Vector3d force = forceAbs * OrientationAxis;

    // save the final forces here
    // check if this is proper!
    //n1->ApplyForce(force,   length * OrientationAxis, forceAbs);
    //n2->ApplyForce(-force, -length * OrientationAxis, forceAbs);

	std::cout << "	-> linearForce on node " << n1->mGlobalIndex << " and node " << n2->mGlobalIndex << ": " << forceAbs << ", strain: " << curStrain << std::endl;

    n1->ApplyForce(force,  abs(forceAbs));
    n2->ApplyForce(-force, abs(forceAbs));

	// added by Jieling
	n1->mLinearForce += force;
	n2->mLinearForce -= force;

    //const double max_force {10000};
    //if (max(abs(force)) > max_force)
    //{
    //    print(std::cout);
    //    std::cout << std::endl;
    //    std::cout << "f=" <<force<<", ";
    //    std::cout << "||f||=" << forceAbs<<std::endl;
    //    std::cout << "OrientAxis:"<<OrientationAxis <<std::endl;
    //    std::cout << "strain:" << strain << std::endl;
    //    std::cout << "n1:" << *n1 << std::endl;
    //    std::cout << "n2:" << *n2<< std::endl;
    //}

    //assert(forceAbs < 1e6);
    //assert(max(abs(force)) < 1e6);
}

// added by Jieling
void LinearSpring::changeN1(ModelElementLatticeNode *N)
{
	n1 = N;
	static_cast<CSGLBar*>(mpGLObject)->changeP1(&N->position);
}

void LinearSpring::changeN2(ModelElementLatticeNode *N)
{
	n2 = N;
	static_cast<CSGLBar*>(mpGLObject)->changeP2(&N->position);
}

void LinearSpring::unboundCheck()
{
	std::string outputFileName1 = "debugFile.txt";
	std::ofstream FO;
	FO.open(outputFileName1.c_str(), std::ios::out | std::ios::app);
	FO << "LinearSpring.cpp / unboundCheck()	starts" << endl;
	if (n1->fibres.size() > 1 && !n1_unbound) // end-node spring and still bound
	{
		double s = getStrain();
		if (s > 1.)
		{
			// unbinding only occurs when stretched
			n1_k_off = 0.025 * (std::exp((s - 1)/0.1) - 1.);
			if (n1_k_off > 1.)
				n1_k_off = 1.;
			else if (n1_k_off < 0.)
				n1_k_off = 0.;
			// random toss here
			double toss = ((double)rand() / (RAND_MAX));
			FO << "	-> s: " << s << ", n1_k_off: " << n1_k_off << ", toss: " << toss << std::endl;
			if (toss < n1_k_off)
			{
				n1_unbound = true;
				std::cout << "	-> unbindg of n1" << std::endl;
				FO << "	-> unbinding of n1" << std::endl;
			}
		}
	}
	if (n2->fibres.size() > 1 && !n2_unbound) // end-node spring and still bound
	{
		double s = getStrain();
		if (s > 1.)
		{
			// unbinding only occurs when stretched
			if (!n1_unbound)
			{
				// if n1 is unbound, then n2 has to bind
				n2_k_off = 0.025 * (std::exp((s - 1)/0.1) - 1.);
				if (n2_k_off > 1.)
					n2_k_off = 1.;
				else if (n2_k_off < 0.)
					n2_k_off = 0.;
				// random toss here
				double toss = ((double)rand() / (RAND_MAX));
				FO << "	-> s: " << s << ", n2_k_off: " << n2_k_off << ", toss: " << toss << std::endl;
				if (toss < n2_k_off)
				{
					n2_unbound = true;
					std::cout << "	-> unbindg of n2" << std::endl;
					FO << "	-> unbinding of n2" << std::endl;
				}
			}
		}
	}
	FO.close();
}

double LinearSpring::getStrainY()
{
	double y = std::abs(n1->position.y - n2->position.y);
	return y / initialLY - 1.;
}

double LinearSpring::getStrain()
{
	double s = (n1->position - n2->position).Norm() / initialL0;
	return s;
}

void LinearSpring::elongation()
{
	auto OrientationAxis = n2->position - n1->position;
	double length = Norm(OrientationAxis);
	double strain = 0.;
	if (std::isnormal(length))
	{
		OrientationAxis /= length;
		strain = length - length0;
	}
	else
	{
		strain = 0;
		std::cerr << "WARNING: " << std::endl;
	}
	strain /= length0;

	double sY = getStrainY();

	if (strain > strainHarden)
	{
		std::cout << " n1: " << n1->latticeIndex << " and n2: " << n2->latticeIndex << " is harden, strain: " << strain << std::endl;
		harden = true;
		// update the default length according to the hardenRatio
		length0 += length0 * (strain - strainHarden) * (1. - hardenRatio);
		//
		double aR = (n1->mRadius + n2->mRadius) / 2;
		double aY = (n1->mYoung + n2->mYoung) / 2;
		((LinearStrainFunction*)msStiffness)->setK(aY * aR * aR * 3.14159265 / length0);
	}
}

bool LinearSpring::ready() const
{
    return ((n1 != nullptr) && (n2 != nullptr)
            && isfinite(n1->position - n2->position)
            && (msStiffness));
}

double LinearSpring::strain() const
{
    return displacement();
}

Vector3d LinearSpring::position(void) const
{
    return  n2->position - n1->position;
}

double LinearSpring::displacement(void) const
{

    return Norm(position()) - length0;
}

double LinearSpring::eq_position() const
{
    return l0;
}

bool LinearSpring::isRelaxed() const
{
    return displacement() < 1e-10;
}

std::vector<unsigned int> LinearSpring::nodeIds() const
{
    return {Node1_Id, Node2_Id};
}

std::vector<ModelElementLatticeNode*> LinearSpring::nodes() const
{
    return {n1, n2};
}

unsigned int LinearSpring::numberOfNodes() const
{
    return 2;
}

BoundingBox *
LinearSpring::boundingBox()
{
	mBoundingBox.xmin = n1->position.x < n2->position.x ? n1->position.x : n2->position.x;
	mBoundingBox.xmax = n1->position.x > n2->position.x ? n1->position.x : n2->position.x;
	
	mBoundingBox.ymin = n1->position.y < n2->position.y ? n1->position.y : n2->position.y;
	mBoundingBox.ymax = n1->position.y > n2->position.y ? n1->position.y : n2->position.y;
	
	mBoundingBox.zmin = n1->position.z < n2->position.z ? n1->position.z : n2->position.z;
	mBoundingBox.zmax = n1->position.z > n2->position.z ? n1->position.z : n2->position.z;
	
	return &mBoundingBox;
}

// added by Jieling
void LinearSpring::strainTestForce()
{
	auto OrientationAxis = n2->position - n1->position;

	double length = Norm(OrientationAxis);

	double strain = 0.;
	if (std::isnormal(length))
	{
		// normalize the OrientationAxis
		OrientationAxis /= length;
		strain = length - length0;
	}
	else
	{
		strain = 0;
		std::cerr << "WARNING: " << std::endl;
	}

	curStrain = strain / length0;

	// compute the force
	double forceAbs = (*msStiffness)(strain);

	// if this is enabled then we reduce stiffness due to node density
	if (mDensityDependency)
		forceAbs *= (n2->density + n1->density) / 2.;

	Vector3d force = forceAbs * OrientationAxis;

	n1->mStrainTestForce += force;
	n2->mStrainTestForce -= force;
}