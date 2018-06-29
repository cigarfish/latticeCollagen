////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  RotationalSpring.cpp                                          //
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

#include "RotationalSpring.h"
#include "ModelLattice.h"
#include "BiologyLink.h"
#include "LinearFunction.h"

#define SIN_LOWER_BOUND 1.e-3

RotationalSpring::RotationalSpring()
    : LatticeSpring(), n1(nullptr), n2(nullptr), center(nullptr), phi0(0.)
{
    //Phi0.setXMLPath("phi0");
    //registerParameter(Phi0);

    //Node1_Id.setXMLPath("Node1_Id");
    //registerParameter(Node1_Id);

    //Node2_Id.setXMLPath("Node2_Id");
    //registerParameter(Node2_Id);

    //NodeC_Id.setXMLPath("NodeC_Id");
    //registerParameter(NodeC_Id);
}

RotationalSpring::RotationalSpring(ModelElementLatticeNode *N1, ModelElementLatticeNode *N2, ModelElementLatticeNode *Center)
	: LatticeSpring(), n1(N1), n2(N2), center(Center), phi0(0.)
{
	initialize();
}

/*RotationalSpring::RotationalSpring(XMLNode& springNode,
                                   std::stringstream& errors,
                                   std::stringstream& warnings)
    : RotationalSpring()
{
    // TODO can we get rid of this?
    load(springNode, errors, warnings);
}*/

RotationalSpring::RotationalSpring(const RotationalSpring& other)
    : LatticeSpring(other),
    n1(other.n1), n2(other.n2), center(other.center),
    phi0(other.phi0),
    Node1_Id(other.Node1_Id), Node2_Id(other.Node2_Id),
    NodeC_Id(other.NodeC_Id), Phi0(other.Phi0)
{}

RotationalSpring::RotationalSpring(RotationalSpring&& other)
    : LatticeSpring(other), n1(other.n1), n2(other.n2),
    center(other.center),
    phi0(other.phi0),
    Node1_Id(other.Node1_Id), Node2_Id(other.Node2_Id),
    NodeC_Id(other.NodeC_Id), Phi0(other.Phi0)
{}

RotationalSpring::~RotationalSpring()
{
    n1     = nullptr;
    n2     = nullptr;
    center = nullptr;
}

void RotationalSpring::print(std::ostream& stream) const
{
    Vector3d nb1 = n1->position - center->position;
    Vector3d nb2 = n2->position - center->position;

    stream << "RotationalSpring(Id="<<id;
    stream << ", phi0="<< phi0
        << ", dphi=" << Angle(nb1, nb2) - phi0
        <<", Node1_Id="<<Node1_Id<<", NodeC_Id="<<NodeC_Id
        <<", Node2_Id="<<Node2_Id<<").";
}

void RotationalSpring::initialize()
{
	// added by Jieling
	this->mSpringType = LatticeSpring::TypeRotationalSpring;

    // get the nodes from the lattice
    //n1     = lattice->getNodeById(Node1_Id);
    //n2     = lattice->getNodeById(Node2_Id);
    //center = lattice->getNodeById(NodeC_Id);

    if (!n1)
    {
        std::cerr << "RotationalSpring(" << id << ") "
            << "Node1: " << Node1_Id << " does not exist!" << std::endl;

        // TODO throw
        assert(false);
        return;
    }

    if (!n2)
    {
        std::cerr << "RotationalSpring(" << id << ") "
            << "Node2: " << Node2_Id << " does not exist!" << std::endl;

        // TODO throw
        assert(false);
        return;
    }

    if (!center)
    {
        std::cerr << "RotationalSpring(" << id << ") "
            << "NodeC: " << NodeC_Id << " does not exist!" << std::endl;

        // TODO throw
        assert(false);
        return;
    }

    if ((n1 == n2) || (n1 == center) || (n2 == center))
    {
        std::cerr << "RotationalSpring(" << id << ") "
                  << "Nodes " << Node1_Id << ", " << Node2_Id << ", "
                  << NodeC_Id << " two of them are the same!" << std::endl;

        throw InternalPluginError("Lattice", "Rotational spring not well defined!");

        return;
    }

    if (mEquilibrium)
    {
        Vector3d n_b  = n1->position - center->position;
        Vector3d n_bp = n2->position - center->position;
        phi0 = Angle(n_b, n_bp);
    }
    else
    {
        phi0 = Phi0;
    }

	msStiffness = new LinearStrainFunction();
	double aR = (n1->mRadius + n2->mRadius + center->mRadius) / 3;
	double aY = (n1->mYoung + n2->mYoung + center->mYoung) / 3;
	mYoung = aY;
	mPoisson = (n1->mPoisson + n2->mPoisson + center->mPoisson) / 3;
	mRadius = 0.0212 * 4; // for test

	((LinearStrainFunction*)msStiffness)->setK(aY * aR * aR * 3.14159265 / phi0 / 1000.); // for test
}

void RotationalSpring::update(double timeStep)
{
    // get new vectors connecting the nodes
    Vector3d nb1 = n1->position - center->position;
    Vector3d nb2 = n2->position - center->position;

    // dot and cross
    double dot = Dot(nb1, nb2);
    double crs = Norm(Cross(nb1, nb2));

    // norms
    double nb1_sq = Norm2Squared(nb1);
    double nb2_sq = Norm2Squared(nb2);

    // norms 2nd
    double nb1_norm = sqrt(nb1_sq);
    double nb2_norm = sqrt(nb2_sq);

    // Compute torque exerted by the torsion spring
    double currentAngle = atan2(crs, dot);
    double TorqueMagnitude = (*msStiffness)(currentAngle - phi0) / timeStep;

    double norm_product = nb1_norm * nb2_norm;
    double c = dot / norm_product;
    // next we divide by s, make sure it is bounded away from zero
    double s = std::max(crs / norm_product, SIN_LOWER_BOUND);

    double a   = - TorqueMagnitude / (2. * s);
    double a11 = a * c / nb1_sq;
    double a12 = - a / norm_product;
    double a22 = a * c / nb2_sq;

    // compute forces
    Vector3d F1 = a11 * nb1 + a12 * nb2;
    Vector3d F3 = a22 * nb2 + a12 * nb1;
    Vector3d F2 = - (F1 + F3);

    // compute force magnitudes
    double Fmag1 = a12 * nb2_norm * s;
    double Fmag2 = Norm(F2);
    double Fmag3 = a12 * nb1_norm * s;

    // apply forces to the nodes
    n1->ApplyForce(F1, Fmag1);
    n2->ApplyForce(F3, Fmag3);

    // center node
    center->ApplyForce(F2, Fmag2);

	// added by Jieling
	n1->mRotationalForce += F1;
	n2->mRotationalForce += F3;
	center->mRotationalForce += F2;

	//std::cout << "	-> rotational force on node " << n1->mGlobalIndex << ", " 
	//		  << n2->mGlobalIndex << " and " << center->mGlobalIndex << ": " 
	//		  << Fmag1 << ", " << Fmag2 << " and " << Fmag3 << ", angle: " << currentAngle - phi0 << std::endl;
}

bool RotationalSpring::ready() const
{
    return ((n1 != nullptr) && (n2 != nullptr)
            && (center != nullptr)
            //&& isfinite(RotationAxis)
            //&& ((RotationAxis != 0) || isRelaxed())
            && (msStiffness));
}

bool RotationalSpring::isRelaxed() const
{
    Vector3d n_b  = n1->position - center->position;
    Vector3d n_bp = n2->position - center->position;

    return (Angle(n_b, n_bp) - phi0) < 1e-10;
}

double RotationalSpring::eq_position() const
{
    return phi0;
}

std::vector<unsigned int> RotationalSpring::nodeIds() const
{
    return {Node1_Id, NodeC_Id, Node2_Id};
}

std::vector<ModelElementLatticeNode*> RotationalSpring::nodes() const
{
    return {n1, center, n2};
}

unsigned int RotationalSpring::numberOfNodes() const
{
    return 3;
}

void RotationalSpring::strainTestForce()
{
	Vector3d nb1 = n1->position - center->position;
	Vector3d nb2 = n2->position - center->position;

	double dot = Dot(nb1, nb2);
	double crs = Norm(Cross(nb1, nb2));

	double nb1_sq = Norm2Squared(nb1);
	double nb2_sq = Norm2Squared(nb2);

	double nb1_norm = sqrt(nb1_sq);
	double nb2_norm = sqrt(nb2_sq);

	double currentAngle = atan2(crs, dot);
	double TorqueMagnitude = (*msStiffness)(currentAngle - phi0);

	double norm_product = nb1_norm * nb2_norm;
	double c = dot / norm_product;

	double s = std::max(crs / norm_product, SIN_LOWER_BOUND);

	double a = -TorqueMagnitude / (2. * s);
	double a11 = a * c / nb1_sq;
	double a12 = -a / norm_product;
	double a22 = a * c / nb2_sq;

	Vector3d F1 = a11 * nb1 + a12 * nb2;
	Vector3d F3 = a22 * nb2 + a12 * nb1;
	Vector3d F2 = -(F1 + F3);

	n1->mStrainTestForce += F1;
	n2->mStrainTestForce += F3;
	center->mStrainTestForce += F2;
}