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

	n10.x = n1->position.x;
	n10.y = n1->position.y;
	n10.z = n1->position.z;
	n20.x = n2->position.x;
	n20.y = n2->position.y;
	n20.z = n2->position.z;
	center0.x = center->position.x;
	center0.y = center->position.y;
	center0.z = center->position.z;

	// calculate the curvature: sqrt(2(1 - cos(theta)))
	double n1cd = (center->position - n1->position).Norm();
	double cn2d = (n2->position - center->position).Norm();
	double n1pn2 = n1->position.x * n2->position.x +
		n1->position.y * n2->position.y +
		n1->position.z * n2->position.z;
	curv0 = sqrt(2 * (1. - n1pn2 / n1cd / cn2d));

	msStiffness = new LinearStrainFunction();
	double aR = (n1->mRadius + n2->mRadius + center->mRadius) / 3;
	double aY = (n1->mYoung + n2->mYoung + center->mYoung) / 3;
	double aS = (n1->mShear + n2->mShear + center->mShear) / 3;
	mYoung = aY;
	mPoisson = (n1->mPoisson + n2->mPoisson + center->mPoisson) / 3;
	double l1 = (n1->position - center->position).Norm();
	double l2 = (n2->position - center->position).Norm();
	fibreLength = l1 + l2;
	double aB = 1. / (
		1. / aY + 12. * 10. / 9. / aS * aR / fibreLength * aR / fibreLength
		);
	mRadius = aR;

	((LinearStrainFunction*)msStiffness)->setK(aB * 3.14159265 * aR * aR * aR * aR / 4.); // EI = Young's modulus * pi * r^4 / 4
}

void RotationalSpring::update(double timeStep)
{
    // get new vectors connecting the nodes
    Vector3d nb1 = n1->position - center->position;
    Vector3d nb2 = n2->position - center->position;

    // dot and cross
    double dot = Dot(nb1, nb2); // cos(theta)*|nb1|*|nb2|
    double crs = Norm(Cross(nb1, nb2)); // sin(theta)*|nb1|*|nb2|

    // norms
    double nb1_sq = Norm2Squared(nb1); // |nb1|^2
    double nb2_sq = Norm2Squared(nb2); // |nb2|^2

    // norms 2nd
    double nb1_norm = sqrt(nb1_sq); // |nb1|
    double nb2_norm = sqrt(nb2_sq); // |nb2|

    // Compute torque exerted by the torsion spring
    double currentAngle = atan2(crs, dot); // theta
    double TorqueMagnitude = (*msStiffness)(currentAngle - phi0); // beta * delta_theta

    double norm_product = nb1_norm * nb2_norm; // |nb1|*|nb2|
    double c = dot / norm_product; // cos(theta)
    // next we divide by s, make sure it is bounded away from zero
    double s = std::max(crs / norm_product, SIN_LOWER_BOUND); // sin(theta)

    double a   = - TorqueMagnitude / (2. * s); // beta * delta_theta / (2 * sin(theta) )
    double a11 = a * c / nb1_sq; // beta * delta_theta / 2 * tan(theta) / |nb1|^2
    double a12 = - a / norm_product; // -beta * delta_theta / (2 * sin(theta)) / (|nb1|*|nb2|)
    double a22 = a * c / nb2_sq;

    // compute forces
    Vector3d F1 = a11 * nb1 + a12 * nb2; // beta * delta_theta / 2 / |nb1| 
    Vector3d F3 = a22 * nb2 + a12 * nb1; // beta * delta_theta / 2 / |nb2|
    Vector3d F2 = - (F1 + F3);

    // compute force magnitudes
    double Fmag1 = a12 * nb2_norm * s; // -beta * delta_theta * |nb2| / 2 / |nb1|
    double Fmag2 = Norm(F2);
    double Fmag3 = a12 * nb1_norm * s; // -beta * delta_theta * |nb1| / 2 / |nb2|

    // apply forces to the nodes
    n1->ApplyForce(F1, Fmag1);
    n2->ApplyForce(F3, Fmag3);

    // center node
    center->ApplyForce(F2, Fmag2);

	// added by Jieling
	n1->mRotationalForce += F1;
	n2->mRotationalForce += F3;
	center->mRotationalForce += F2;
	// update force/pressure for nodes and spring
	n1->accumulatedForceAbsolute += abs(Fmag1);
	n1->accumulatedPressure += abs(Fmag1) / (mRadius * mRadius * 3.14159265);
	n2->accumulatedForceAbsolute += abs(Fmag3);
	n2->accumulatedPressure += abs(Fmag3) / (mRadius * mRadius * 3.14159265);
	center->accumulatedForceAbsolute += abs(Fmag2);
	center->accumulatedPressure += abs(Fmag2) / (mRadius * mRadius * 3.14159265);

	//std::cout << "	-> rotational force on node " << n1->mGlobalIndex << ", " 
	//		  << n2->mGlobalIndex << " and " << center->mGlobalIndex << ": " 
	//		  << Fmag1 << ", " << Fmag2 << " and " << Fmag3 << ", angle: " << currentAngle - phi0 << std::endl;
}

double RotationalSpring::getEI()
{
	return (*msStiffness)(1.);
}

double RotationalSpring::getCurvature()
{
	double curv = 0.;
	/* 
		curvature for the center node

		cx, cy, cz is the coordinate of the curvature circle

		(c - (n1 + center)/2) * (center - n1) = 0
		(c - (n2 + center)/2) * (center - n2) = 0
		| cx  cy  cz  1 | = 0
		| n1x n1y n1z 1 |
		| c_x c_y c_z 1 |
		| n2x n2y n2z 1 |

		| a1 b1 c1 | cx   e1
		| a2 b2 c2 | cy = e2, Ac=e
		| a3 b3 c3 | cz   e3

				     | A D G |
		A^-1 = det(A)| B E H |
					 | C F I |
	*/
	// first check if curvature equals to 0
	double crossx = (n1->position.y - center->position.y) * (n2->position.z - center->position.z) - 
					(n2->position.y - center->position.y) * (n1->position.z - center->position.z);
	double crossy = (n1->position.x - center->position.x) * (n2->position.z - center->position.z) -
					(n2->position.x - center->position.x) * (n1->position.z - center->position.z);
	double crossz = (n1->position.x - center->position.x) * (n2->position.y - center->position.y) -
					(n2->position.x - center->position.x) * (n1->position.y - center->position.y);
	double crossd = crossx * crossx + crossy * crossy + crossz * crossz;
	if (crossd > 0)
	{
		double a1 = n1->position.x - center->position.x;
		double b1 = n1->position.y - center->position.y;
		double c1 = n1->position.z - center->position.z;
		double e1 = (n1->position.x * n1->position.x - center->position.x * center->position.x +
			n1->position.y * n1->position.y - center->position.y * center->position.y +
			n1->position.z * n1->position.z - center->position.z * center->position.z) / 2;
		double a2 = n2->position.x - center->position.x;
		double b2 = n2->position.y - center->position.y;
		double c2 = n2->position.z - center->position.z;
		double e2 = (n2->position.x * n2->position.x - center->position.x * center->position.x +
			n2->position.y * n2->position.y - center->position.y * center->position.y +
			n2->position.z * n2->position.z - center->position.z * center->position.z) / 2;
		double a3 = n1->position.y * (center->position.z - n2->position.z) -
			n1->position.z * (center->position.y - n2->position.y) +
			center->position.y * n2->position.z - center->position.z * n2->position.y;
		double b3 = -n1->position.x * (center->position.z - n2->position.z) +
			n1->position.z * (center->position.x - n2->position.x) -
			center->position.x * n2->position.z - center->position.z * n2->position.x;
		double c3 = n1->position.x * (center->position.y - n2->position.y) -
			n1->position.y * (center->position.x - n2->position.x) +
			center->position.x * n2->position.y - center->position.y * n2->position.x;
		double e3 = n1->position.x * (center->position.y * n2->position.z - center->position.z * n2->position.y) -
			n1->position.y * (center->position.x * n2->position.z - center->position.z * n2->position.x) +
			n1->position.z * (center->position.x * n2->position.y - center->position.y * n2->position.x);
		// the inverse of A
		double detA = a1 * (b2 * c3 - c2 * b3) - 
					  a2 * (b1 * c3 - c1 * b3) +
					  a3 * (b1 * c2 - c1 * b2);
		double A = (b2 * c3 - c2 * b3);
		double D = -(b1 * c3 - c1 * b3);
		double G = (b1 * c2 - c1 * b2);
		double B = -(a2 * c3 - c2 * a3);
		double E = (a1 * c3 - c1 * a3);
		double H = -(a1 * c2 - c1 * a2);
		double C = (a2 * b3 - b2 * a3);
		double F = -(a1 * b3 - b1 * a3);
		double I = (a1 * b2 - b1 * a2);
		double cx = (A * e1 + D * e2 + G * e3) / detA;
		double cy = (B * e1 + E * e2 + H * e3) / detA;
		double cz = (C * e1 + F * e2 + I * e3) / detA;
		double ra = std::sqrt((cx - (n1->position.x + center->position.x) / 2) * (cx - (n1->position.x + center->position.x) / 2) +
			(cy - (n1->position.y + center->position.y) / 2) * (cy - (n1->position.y + center->position.y) / 2) +
			(cz - (n1->position.z + center->position.z) / 2) * (cz - (n1->position.z + center->position.z) / 2));
		curv = 1. / ra;
	}

	return curv;
}

double RotationalSpring::getBendingForce()
{
	double lnx = n2->position.x - center->position.x;
	double lny = n2->position.y - center->position.y;
	double lnz = n2->position.z - center->position.z;
	double lnd = sqrt(lnx * lnx + lny * lny + lnz * lnz);
	//
	double lnm1x = center->position.x - n1->position.x;
	double lnm1y = center->position.y - n1->position.y;
	double lnm1z = center->position.z - n1->position.z;
	double lnm1d = sqrt(lnm1x * lnm1x + lnm1y * lnm1y + lnm1z * lnm1z);
	//
	double lnlnm1 = lnx * lnm1x + lny * lnm1y + lnz * lnm1z;
	//
	double dx = (lnx - lnm1x)/lnd/lnm1d + lnlnm1/lnd/lnm1d * (lnx/lnd/lnd - lnm1x/lnm1d/lnm1d);
	double dy = (lny - lnm1y)/lnd/lnm1d + lnlnm1/lnd/lnm1d * (lny/lnd/lnd - lnm1y/lnm1d/lnm1d);
	double dz = (lnz - lnm1z)/lnd/lnm1d + lnlnm1/lnd/lnm1d * (lnz/lnd/lnd - lnm1z/lnm1d/lnm1d);
	// calculate the moment:
	double curv = sqrt(2 * (1. - lnlnm1 / lnd / lnm1d)) - curv0;
	double theta = acos(lnlnm1 / lnd / lnm1d);
	double moment = (*msStiffness)(curv) / lnd;
	double dcos = 1.;
	if (abs(theta - phi0) < 1e-8)
		dcos = 0.;
	else
	{
		dcos = (theta - phi0) / sin(theta);
	}
	//
	double scaleF = (*msStiffness)(1.) / lnd * dcos;
	// this dx, dy, dz is cos(theta(a+1)) for ra of d(a+1) and d(a), for the center
	// assign it to center
	center->mRotationalForce.x += dx * scaleF;
	center->mRotationalForce.y += dy * scaleF;
	center->mRotationalForce.z += dz * scaleF;
	center->directedForce.x += dx * scaleF;
	center->directedForce.y += dy * scaleF;
	center->directedForce.z += dz * scaleF;
	//
	dx = -1. / lnd / lnm1d * (lnlnm1 * lnx / lnd / lnd - lnm1x);
	dy = -1. / lnd / lnm1d * (lnlnm1 * lny / lnd / lnd - lnm1y);
	dz = -1. / lnd / lnm1d * (lnlnm1 * lnz / lnd / lnd - lnm1z);
	// this dx, dy, dz is cos(theta(a)) for ra of d(a-1) and d(a), for the n2
	// assign it to n2
	n2->mRotationalForce.x += dx * scaleF;
	n2->mRotationalForce.y += dy * scaleF;
	n2->mRotationalForce.z += dz * scaleF;
	n2->directedForce.x += dx * scaleF;
	n2->directedForce.y += dy * scaleF;
	n2->directedForce.z += dz * scaleF;
	//
	dx = 1. / lnd / lnm1d * (lnlnm1 / lnm1d / lnm1d * lnm1x - lnx);
	dy = 1. / lnd / lnm1d * (lnlnm1 / lnm1d / lnm1d * lnm1y - lny);
	dz = 1. / lnd / lnm1d * (lnlnm1 / lnm1d / lnm1d * lnm1z - lnz);
	// this dx, dy, dz is cos(theta(a+2)) for ra of d(a+1) and d(a+2), for the n1
	// assign it to n1
	n1->mRotationalForce.x += dx * scaleF;
	n1->mRotationalForce.y += dy * scaleF;
	n1->mRotationalForce.z += dz * scaleF;
	n1->directedForce.x += dx * scaleF;
	n1->directedForce.y += dy * scaleF;
	n1->directedForce.z += dz * scaleF;

	return moment;
}

void RotationalSpring::reset()
{
	Vector3d n_b = n1->position - center->position;
	Vector3d n_bp = n2->position - center->position;
	// update the angle
	phi0 = Angle(n_b, n_bp);
	// update the curvature
	double n1cd = (center->position - n1->position).Norm();
	double cn2d = (n2->position - center->position).Norm();
	double n1pn2 = n1->position.x * n2->position.x +
		n1->position.y * n2->position.y +
		n1->position.z * n2->position.z;
	curv0 = sqrt(2 * (1. - n1pn2 / n1cd / cn2d));

	double aR = (n1->mRadius + n2->mRadius + center->mRadius) / 3;
	double aY = (n1->mYoung + n2->mYoung + center->mYoung) / 3;
	double aS = (n1->mShear + n2->mShear + center->mShear) / 3;
	mYoung = aY;
	mPoisson = (n1->mPoisson + n2->mPoisson + center->mPoisson) / 3;
	double aB = 1. / (
		1. / aY + 12. * 10. / 9. / aS * aR / fibreLength * aR / fibreLength
		);

	((LinearStrainFunction*)msStiffness)->setK(aY * aR * aR * aR * aR * 3.14159265 / 4.); // for test
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

void RotationalSpring::changeN1(ModelElementLatticeNode* N)
{
	n1 = N;
}

void RotationalSpring::changeN2(ModelElementLatticeNode* N)
{
	n2 = N;
}

void RotationalSpring::changeCenter(ModelElementLatticeNode* N)
{
	center = N;
}

void RotationalSpring::setFibreLengthOnly()
{
	fibreLength = (n1->position - center->position).Norm() + (n2->position - center->position).Norm();
	reset();
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