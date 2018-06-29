#ifndef CS_MODEL_TOOLS_H
#define CS_MODEL_TOOLS_H

#define _USE_MATH_DEFINES
#include "../../model/BasicDatatypes/Vector.h"

#include "../../model/Elements/ModelElementSphere.h"
#include "tools/utils/StringUtils.h"

#include <cmath>
#include <math.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <limits>
#include <tuple>

// Converts degrees to radians.
#define degreesToRadians(angleDegrees) (angleDegrees * M_PI / 180.0)

// Converts radians to degrees.
#define radiansToDegrees(angleRadians) (angleRadians * 180.0 / M_PI)

#define GRAVITATIONAL_ACCELERATION (9.80665)


namespace CSModelTools {

	inline bool inBetween(const double lower, const double upper, const double value)
	{
		return ((value < upper) && (value > lower));
	}

	inline double SurfaceAreaSphere(const double radius)
	{
		return 4. * M_PI * radius * radius;
	}

	inline double VolumeSphere(const double radius)
	{
		return (4 / 3.) * radius * radius * radius * M_PI;
	}

	inline double VolumeRadius(const double volume)
	{
		return pow(3. * volume / (4. * M_PI), 1 / 3.);
	}

	// fix naming!
	inline double SurfaceAreaRod2(const double radius, const double height)
	{
		return 2. * M_PI * radius * height;
	}

	// this includes the spherical caps
	inline double SurfaceAreaRod(const double radius, const double height)
	{
		return 4. * M_PI * radius * radius + 2. * M_PI * radius * height;
	}

	inline double CircumferenceRod(const double radius, const double height)
	{
		return 2. * (M_PI * radius + height);
	}

	inline double VolumeRod(const double radius, const double height)
	{
		return (4 / 3.) * M_PI * radius * radius * radius +
			M_PI * height * radius * radius;
	}

	inline bool notInContactSpheres(const double distance,
		const double r1,
		const double r2)
	{
		return distance > (r1 + r2);
	}

	inline double CircumferenceDisk(const double radius)
	{
		return 2. * radius * M_PI;
	}

	inline double AreaDisk(const double radius)
	{
		return radius * radius * M_PI;
	}

inline
double GetDistance2D(Vector3f p1, Vector3f p2)
{
  return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

inline
double GetDistance3D(Vector3f p1, Vector3f p2)
{
  return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

inline
double GetDistance3DPointToPlane(Vector3f point, double normalVector[3], double D)
{
  return abs(point.x*normalVector[0]+point.y*normalVector[1]+point.z*normalVector[2]-D);
}

inline
double GetDistance3D(Vector3f point, double normalVector[3], double D)
{
  return (point.x*normalVector[0]+point.y*normalVector[1]+point.z*normalVector[2]-D);
}

inline
double GetHertzForce(const double A_poi, const double A_young, const double A_rad, const double B_poi, const double B_young, const double B_rad, const double AB_dist, const double single_bond_energy, const double adhesion_density)
{
    return 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1- B_poi * B_poi) / B_young) *
           sqrt(A_rad * B_rad / (A_rad + B_rad)) * pow((A_rad + B_rad - AB_dist), (double)1.5) - M_PI *
           single_bond_energy * adhesion_density * (A_rad * B_rad / ( A_rad + B_rad));
}

inline
double GetHertzForce(const double A_poi, const double A_young, const double A_rad, const double AB_dist, const double /*single_bond_energy*/, const double /*adhesion_density*/)
{
  return 4. / 3. * A_young / (1- A_poi * A_poi ) * sqrt( A_rad ) * pow( (A_rad - AB_dist) , (double) 1.5);
}

// added by Jieling: between sphere and cylinder
inline
double GetHertzForceSphereCylinder(const double A_poi, const double A_young, const double A_rad, const double B_poi, const double B_young, const double B_rad, const double AB_overlap, const double single_bond_energy, const double adhesion_density, int type)
{
	double F = 0.;
	// A: sphere, B: cylinder
	// type 0: edge, 1: node
	if (type == 0)
		F = 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1 - B_poi * B_poi) / B_young) *
		(sqrt(A_rad * B_rad / (A_rad + B_rad)) + sqrt(A_rad)) / 2 * pow(AB_overlap, (double)1.5) - 
		M_PI * single_bond_energy * adhesion_density
			 * (sqrt(A_rad * B_rad / (A_rad + B_rad)) + sqrt(A_rad)) 
			 * (sqrt(A_rad * B_rad / (A_rad + B_rad)) + sqrt(A_rad));
	else
		F = 2. * B_rad / ((1 - A_poi * A_poi) / A_young + (1 - B_poi * B_poi) / B_young) * AB_overlap -
		M_PI * single_bond_energy * adhesion_density
			 * B_rad * B_rad;
	// cell and cylinder already collide with each other, the force has to be larger than 0
	if (F < 0)
	{
		std::cerr << "	-> Warning: the Hertz force between the cell and ECM fiber is below 0!" << std::endl;
		std::cerr << "		type: " << type << std::endl;
		std::cerr << "		1st term: " << 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1 - B_poi * B_poi) / B_young) *
			(sqrt(A_rad * B_rad / (A_rad + B_rad)) + sqrt(A_rad)) / 2 * pow(AB_overlap, (double)1.5) << std::endl;
		std::cerr << "		2nd term: " << M_PI * single_bond_energy * adhesion_density
			* (sqrt(A_rad * B_rad / (A_rad + B_rad)) + sqrt(A_rad))
			* (sqrt(A_rad * B_rad / (A_rad + B_rad)) + sqrt(A_rad)) << std::endl;
		F = 0.;
	}

	return F;
}

inline 
double GetContactAreaHertzSphereCylinder(const double Ri, const double Rj, const double ijOverlap, int type)
{
	// i: sphere, j: cylinder
	// type 0: edge, 1: node
	// effective radius: sqrt(R) = 1/2 * (sqrt(R1*R2/(R1+R2)) + sqrt(R1))
	if (type == 0)
		return M_PI * (sqrt((Ri * Rj) / (Ri + Rj)) + sqrt(Ri)) 
					* (sqrt((Ri * Rj) / (Ri + Rj)) + sqrt(Ri)) / 4
					* ijOverlap;
	else
		return M_PI * Rj * Rj;
}

inline
double GetJKRForceSphereCylinder(const double E_eff, const double R_eff,
								 const double AB_overlap, const double adhesion_constant, 
								 double &F, double &A, const double dcrit, int type)
{
	// A: sphere, B: cylinder
	A = 1.0 / R_eff;
	double B = sqrt(2 * M_PI * adhesion_constant/E_eff);
	// Initial guess to compute contact radius 
	double aguess = pow( (AB_overlap + 2 * dcrit) * R_eff, 0.25);
	double a, aprev;

	a = aguess;

	for (unsigned int i = 0; i < 10; i++)
	{
		aprev = a;
		// newton iterations: x(n+1) = x(n) - f/f'
		a = a - (A * a * a * a * a - B * a - AB_overlap) / (4 * A * a * a * a + 0.5 * B);
		if (fabs((a - aprev) / aprev) < 1e-6) break; // adjust precision
	}

	a *= a;

	// contact radius
	if (a != a)
	{
		std::cerr << "NaN in JKR!!!" << std::endl;
		std::cerr << "	-> E_eff: " << E_eff << ", R_eff: " << R_eff << ", AB_overlap: " << AB_overlap << std::endl;
		throw;
	}

	A = M_PI * a * a;
	// force
	F = 4. / 3. * (E_eff / R_eff) * a * a * a -
		sqrt(8 * M_PI * E_eff * a * a * a * adhesion_constant);

	// The cell and cylinder are already colliding with each other, the force has to be larger than 0
	if (F < 0)
	{
		std::cerr << "	-> Warning: the JKR force between the cell and ECM fiber is below 0!" << std::endl;
		std::cerr << "		contact radius a: " << a << std::endl;
		std::cerr << "		1st term: " << 4. / 3. * (E_eff / R_eff) * a * a * a << std::endl;
		std::cerr << "		2nd term: " << sqrt(8 * M_PI * E_eff * a * a * a * adhesion_constant) << std::endl;
		F = 0;
	}
}

inline
double GetHertzForce(const double A_poi, const double A_young, const double A_rad, const double B_poi, const double B_young, const double B_rad, const double AB_dist, const double single_bond_energy, const double A_adhesion_density, const double B_adhesion_density)
{
	// Temporary workaround for MSVS (std::min leads to conflicting macros in MSVS)
#if WIN32
	double adhestionDensity;
	if (A_adhesion_density < B_adhesion_density) adhestionDensity = A_adhesion_density;
	else adhestionDensity = B_adhesion_density;

    return 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1- B_poi * B_poi) / B_young) *
    sqrt(A_rad * B_rad / (A_rad + B_rad)) * pow((A_rad + B_rad - AB_dist), (double)1.5) - M_PI *
	single_bond_energy * adhestionDensity * (A_rad * B_rad / ( A_rad + B_rad));
#else
    return 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1- B_poi * B_poi) / B_young) *
    sqrt(A_rad * B_rad / (A_rad + B_rad)) * pow((A_rad + B_rad - AB_dist), (double)1.5) - M_PI *
	single_bond_energy * std::min(A_adhesion_density, B_adhesion_density) * (A_rad * B_rad / ( A_rad + B_rad));
#endif


}


inline
double GetHertzForcePure(const double A_poi, const double A_young, const double A_rad, const double B_poi, const double B_young, const double B_rad, const double AB_dist)
{
    return 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1- B_poi * B_poi) / B_young) *
        sqrt(A_rad * B_rad / (A_rad + B_rad)) * pow((A_rad + B_rad - AB_dist), (double)1.5);
}


inline
void GetJKRForce(const double E_eff,   const double R_eff,
                 const double AB_dist, const double adhesion_constant,
                 double & F, double &A, const double dcrit, int type)
{

  A= 1.0 / R_eff;
  double B = sqrt(2* M_PI*adhesion_constant/E_eff);
  //  Initial guess to compute contact radius  (why 2*d_crit? it works!)
  double aguess = pow( (AB_dist + 2*dcrit)*R_eff, 0.25 );
  double a, aprev;

  a = aguess;

  for (unsigned int i =0; i < 10; i++)
  {
      aprev=a;
      // newton iterations: xn+1 = xn - f/f'
      a = a -(A*a*a*a*a - B*a - AB_dist)/(4 * A * a*a*a + 0.5 * B);
      if (fabs((a-aprev)/aprev) < 1e-6) break; // adjust precision

      // if (i > 8)
      //     std::cerr << "warning :  JKR contact area algorithm"
      //               << "takes a lot of iterations !"
      //               << std::endl;
  }

  // original formula works with 'a', but then there's a sqrt in f,
  // and this is not so effective....
  a *= a;

  // contact radius
  if( a != a )
  {
      std::cerr << "NaN in JKR !!!" << std::endl;
	  std::cerr << "	-> E_eff: " << E_eff << ", R_eff: " << R_eff << ", AB_dist: " << AB_dist << ", adhesion_constant: " << adhesion_constant << ", F: " << F << ", A: " << A << ", dcrit: " << dcrit << ", type: " << type << std::endl;
      throw;
  }

  A = M_PI * a * a;
  //force
  F = 4.0/3.0 * (E_eff/R_eff) * a * a * a
      - sqrt(8 * M_PI * E_eff * a * a * a * adhesion_constant);

} // end JKR


inline
double GetContactAreaHertz(const double Ri, const double Rj, const double dij)
{
  return M_PI *(Ri * Rj / (Ri + Rj)) * ( Ri + Rj - dij);
}

inline
double GetContactAreaHertz(const double R, const double dij)
{
  //A = pi*r²
  return M_PI * R * (R - dij );
}

inline
double ClampDouble(double value, double min, double max)
{
  if (value < min) return min;
  if (value > max) return max;
  return value;
}

inline
bool OverlapSphericalElements(ModelElementSphere *sphere1, ModelElementSphere *sphere2 ){

  //only for 3d
  double dist = GetDistance3D(sphere1->position,sphere2->position);
  dist -= sphere1->mRadius;
  dist -= sphere2->mRadius;

  return ( dist < 0);

}

// gravity + buoyancy force
inline double GetGravity(const double c_rad, const double dens_liq,
	const double dens_cell, const double length_scale,
	const double time_scale)
{
	return GRAVITATIONAL_ACCELERATION * time_scale*time_scale / length_scale * 4. / 3.
		* M_PI * c_rad*c_rad*c_rad * (dens_cell - dens_liq);
}

// increase of volume dV of a cylinder with dr=dh=value returned
inline double GetDeltaRForUniformCylinderVolumeIncrease(const double r, const double h,
	const double dV)
{
	return (sqrt(4 * M_PI*h*h*r*r + 4 * M_PI*h*r*r*r + 4 * h*dV + M_PI*r*r*r*r) -
		2 * sqrt(M_PI)*h*r - sqrt(M_PI)*r*r) / (2 * sqrt(M_PI)*h);
}

// hill function to determine the volume increase according to the cell area
inline double GetDeltaVolumeAreaDependent(const double radius, const double hillArea)
{
	double area = M_PI*radius*radius;
	return pow(area, 5.) / (pow(hillArea, 5.) + pow(area, 5.));
}

// hill function to determine the volume increase according to the ratio between cell area and cell-cell contact surface
inline double GetDeltaVolumeRatioDependent(const double radius, const double height,
	const double hillRatio)
{
	double ratio = radius / height;
	return pow(ratio, 5.) / (pow(hillRatio, 5.) + pow(ratio, 5.));
}

inline double GetOverlapAreaSpheres2D(const double R1, const double R2, const double d)
{
	// http://mathworld.wolfram.com/Circle-CircleIntersection.html
	return R1*R1*acos((d*d + R1*R1 - R2*R2) / (2 * d*R1)) + R2*R2*acos((d*d + R2*R2 - R1*R1) / (2 * d*R2)) -
		.5*sqrt((-d + R1 + R2)*(d + R1 - R2)*(d - R1 + R2)*(d + R1 + R2));
}

inline double GetOverlapAreaSphereRect2D(const double R, const double h)
{
	// http://mathworld.wolfram.com/CircularSegment.html
	return R*R*acos(1 - h / R) - (R - h)*sqrt(2 * R*h - h*h);
}

inline int ConvertTimeUnitToSteps(double val, double timeStep, std::string unit)
{
	lower_case(unit);
	if (unit == "seconds") { return (int)(val / timeStep); }
	else if (unit == "minutes") { return (int)(60 * val / timeStep); }
	else if (unit == "hours") { return (int)(3600 * val / timeStep); }
	else if (unit == "days") { return (int)(86400 * val / timeStep); }
	else { return (int)val; }
}

inline int ConvertTimeUnitToSeconds(double val, std::string unit)
{
	lower_case(unit);
	if (unit == "seconds") { return (int)(val); }
	else if (unit == "minutes") { return (int)(60 * val); }
	else if (unit == "hours") { return (int)(3600 * val); }
	else if (unit == "days") { return (int)(86400 * val); }
	else { return (int)val; }
}

inline double ConvertTimeUnitToDays(double val, std::string unit)
{
	lower_case(unit);
	if (unit == "seconds") { return (val / 86400); }
	else if (unit == "minutes") { return (val / 1440); }
	else if (unit == "hours") { return (val / 24); }
	else if (unit == "days") { return (val); }
	else { return (int)val; }
}

inline double
GetCellDivisionForce(const double A_poi, const double A_young, const double A_rad,
	const double B_poi, const double B_young, const double B_rad,
	const double AB_dist, const double single_bond_energy,
	const double A_adhesion_density,
	const double B_adhesion_density)
{
	const double mThreshold = 0.7;
	if (AB_dist < mThreshold * (A_rad + B_rad))
		return GetHertzForcePure(A_poi, A_young, A_rad, B_poi, B_young,
			B_rad, mThreshold * (A_rad + B_rad));
	else
		return GetHertzForce(A_poi, A_young, A_rad, B_poi, B_young, B_rad,
			AB_dist, single_bond_energy, A_adhesion_density,
			B_adhesion_density);
}

inline double
GetDivisionForce(const double A_poi, const double A_young, const double A_rad,
	const double B_poi, const double B_young, const double B_rad,
	const double AB_dist, const double single_bond_energy,
	const double A_adhesion_density, const double B_adhesion_density,
	const double time_in_div, const double divtime,
	const double divisionRadius)
{
	double fHertz = GetCellDivisionForce(A_poi, A_young, A_rad, B_poi, B_young,
		B_rad, AB_dist, single_bond_energy,
		A_adhesion_density, B_adhesion_density);

	const double K_div = 10.;
	double force = K_div * fHertz * ((2. * divisionRadius * time_in_div / divtime) - AB_dist);
	return force > 50000 ? 50000 : force;
}

inline
double GetContactAreaHertzReff(const double R_eff, const double delta)
{
	return M_PI * R_eff * delta;
}



inline
double GetJKRDeltaCritical(const double E_eff, const double R_eff,
	const double adhesion_constant)
{
	return std::pow(3. * R_eff * std::pow(M_PI * adhesion_constant / (8.*E_eff), 2),
		0.333333333333);
}

inline
void GetJKRContactRadius(const double E_eff, const double R_eff,
	const double AB_dist, const double adhesion_constant,
	double & R, const double dcrit)
{
	double A = 1.0 / R_eff;
	double B = sqrt(2 * M_PI*adhesion_constant / E_eff);
	//  Initial guess to compute contact radius  (why 2*d_crit? it works!)
	double aguess = pow((AB_dist + 2 * dcrit)*R_eff, 0.25);
	double a, aprev;

	a = aguess;

	for (unsigned int i = 0; i < 10; i++)
	{
		aprev = a;
		// Newton's method to solve
		// A a**2 - sqrt(2 pi gamma a / E) - delta = 0
		// 1. Let x = sqrt(a)
		// Then we have to solve
		// f(x) = A x**4 - B x - delta = 0
		//
		// 2. Newton iterations
		// xn+1 = xn - f(xn)/f'(xn)
		//
		// here f'(xn) = 4 A x**3 - B
		//
		a = a - (A*a*a*a*a - B*a - AB_dist) / (4 * A * a*a*a - B);
		if (fabs((a - aprev) / aprev) < 1e-6) break; // adjust precision
	}

	// original formula works with 'a', but then there's a sqrt in f,
	// and this is not so effective....
	a *= a;

	// contact radius
	if (!std::isfinite(a))
	{
		// THROW ACTUAL EXCEPTION
		std::cerr << "NaN in JKR !!!" << std::endl;
		throw;
	}

	R = a;
}

inline
double GetJKRContactArea(const double E_eff, const double R_eff,
	const double AB_dist, const double adhesion_constant,
	const double dcrit)
{
	double a;
	GetJKRContactRadius(E_eff, R_eff, AB_dist, adhesion_constant, a, dcrit);
	return M_PI * a * a;
}

inline
double GetHertzEnergyPure(const double A_poi, const double A_young,
	const double A_rad, const double B_poi,
	const double B_young, const double B_rad,
	const double AB_dist)
{
	return 8. / 15. / ((1 - A_poi * A_poi) / A_young + (1 - B_poi * B_poi) / B_young)
		* std::sqrt(A_rad * B_rad / (A_rad + B_rad))
		* std::pow(A_rad + B_rad - AB_dist, 2.5);
}

inline
double GetHertzEnergyPure(const double E_eff, const double R_eff,
	const double delta)
{
	return 8. / 15. * E_eff * std::sqrt(R_eff) * std::pow(delta, 2.5);
}

inline double
GetHertzForceNoAdhesion(const double A_poi, const double A_young,
	const double A_rad, const double B_poi,
	const double B_young, const double B_rad,
	const double AB_dist)
{
	return 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1 - B_poi * B_poi) / B_young) *
		sqrt(A_rad * B_rad / (A_rad + B_rad)) * pow((A_rad + B_rad - AB_dist), (double) 1.5);
}

inline double
GetHertzForceAdhesion(const double A_rad, const double B_rad, const double AB_dist,
	const double single_bond_energy,
	const double A_adhesion_density, const double B_adhesion_density)
{
	return M_PI * single_bond_energy * std::min(A_adhesion_density, B_adhesion_density)
		* (A_rad * B_rad / (A_rad + B_rad));
}

inline double
GetHertzForceAdhesion(const double A_rad, const double B_rad, const double AB_dist,
	const double single_bond_energy, const double adhesion_density)
{
	return M_PI * single_bond_energy * adhesion_density*(A_rad * B_rad / (A_rad + B_rad));
}

inline double GetHertzForceHollowSphere(const double A_poi, const double A_young, const double A_rad, const double B_poi, const double B_young, const double B_rad, const double AB_dist, const double single_bond_energy, const double adhesion_density) {
	return 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1 - B_poi * B_poi) / B_young) *
		sqrt(A_rad * B_rad / (-A_rad + B_rad)) * pow((-B_rad + AB_dist + A_rad), 1.5) - M_PI *
		single_bond_energy * adhesion_density * (A_rad * B_rad / (-A_rad + B_rad));
}

inline double GetHertzForceHollowSphereNoAdhesion(const double A_poi, const double A_young, const double A_rad, const double B_poi, const double B_young, const double B_rad, const double AB_dist) {
	return 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1 - B_poi * B_poi) / B_young) *
		sqrt(A_rad * B_rad / (-A_rad + B_rad)) * pow((-B_rad + AB_dist + A_rad), 1.5);
}

inline double GetHertzForceHollowCylinder(const double A_poi, const double A_young, const double A_rad, const double B_poi, const double B_young, const double B_rad, const double AB_dist, const double single_bond_energy, const double adhesion_density) {
	return 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1 - B_poi * B_poi) / B_young) *
		sqrt(A_rad * B_rad / (-0.5*A_rad + B_rad)) * pow((-B_rad + AB_dist + A_rad), 1.5) - M_PI *
		single_bond_energy * adhesion_density * (A_rad * B_rad / (-0.5*A_rad + B_rad));
}

inline double GetHertzForceHollowCylinderNoAdhesion(const double A_poi, const double A_young, const double A_rad, const double B_poi, const double B_young, const double B_rad, const double AB_dist) {
	return 4. / 3. / ((1 - A_poi * A_poi) / A_young + (1 - B_poi * B_poi) / B_young) *
		sqrt(A_rad * B_rad / (-0.5*A_rad + B_rad)) * pow((-B_rad + AB_dist + A_rad), 1.5);
}

inline double GetHertzForceLongRange(const double A_rad, const double B_rad, const double AB_dist, const double single_bond_energy, const double adhesion_density, const double sigma_sq) {
	return -M_PI * single_bond_energy * adhesion_density * (A_rad * B_rad) / (A_rad + B_rad) * exp(-0.5 * (A_rad + B_rad - AB_dist)*(A_rad + B_rad - AB_dist) / (sigma_sq));
}

// Explicitly integrating Hertzforce from A_rad+B_rad to AB_dist:
inline
double GetHertzEnergy(const double A_poi, const double A_young, const double A_rad,
	const double B_poi, const double B_young, const double B_rad,
	const double AB_dist, const double single_bond_energy,
	const double adhesion_density)
{
	return 8. / 15. / ((1 - A_poi * A_poi) / A_young + (1 - B_poi * B_poi) / B_young)
		* std::sqrt(A_rad * B_rad / (A_rad + B_rad))
		* std::pow((A_rad + B_rad - AB_dist), 2.5)
		- M_PI * single_bond_energy * adhesion_density *
		(A_rad * B_rad / (A_rad + B_rad)) * (AB_dist - (A_rad + B_rad));
}

inline
double GetHertzEnergy(const double A_poi, const double A_young, const double A_rad,
	const double AB_dist, const double /*single_bond_energy*/,
	const double /*adhesion_density*/)
{
	return 8. / 15. * A_young / (1 - A_poi * A_poi) * sqrt(A_rad)
		* std::pow(A_rad - AB_dist, 2.5);
}

} // END namespace CSModelTools


#endif // CS_MODEL_TOOLS_H
