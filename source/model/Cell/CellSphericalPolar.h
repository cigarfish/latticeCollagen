///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CellSphericalPolar.h                                                 //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2013-11-05 20:22:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////
#ifndef CELLSPHERICALPOLAR_H
#define CELLSPHERICALPOLAR_H

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#include <math.h>
#endif

#include "CellSpherical.h"
#include <iostream>
#include "../../tools/random/Random.h"
#include "../BasicDatatypes/CSListContainer.h"

class CSModel;
namespace H5 {
    class CompType;
}

#define SIGN(x_) (x_ > 0) - (x_ < 0)
using namespace std;

//! Defines a spherical cell with restricted adhesive areas which are located at
//  its poles.  The area is defined by an orientation (i.e. where the poles reside)
//  and a given opening angle relative to the orientation axis.
class CellSphericalPolar : public CellSpherical
{
public:
  CellSphericalPolar( double x=0, double y=0, double z=0 );
  virtual ~CellSphericalPolar();


  virtual void Reset()
  {
      CellSpherical::Reset();

      mpContactsPointers->clear();
     // Test Paul & Noemie
     // mPolarDirection.Set(mNewPolarDirection);
      mAdhesionEnergy = 0.;
  };

  //! Overriding virtual ModelElement * CellSpherical::Divide(Model *)
  virtual ModelElement * Divide( CSModel * parentModel );

  //! Duplicate this cell and postion it at x,y,z.
  //! Make sure all relevant member values are set from this cell
  CellSphericalPolar * Duplicate(double x, double y, double z);

  //! Metropolis algorithm for optimising the adhesion energy
  inline void OrientationMetropolis();

  //! Randomly tilt the polar direction by a given angle.
  inline void probeRandomRotation(double maxAngle);

  //! Calculate overlap of adhesive regions of two CellSphericalPolar cells.
  inline static double CalculateOverlap( CellSphericalPolar *, CellSphericalPolar *, double contactRadius = 0. );

  inline static double CalculateOverlapSimple( CellSphericalPolar *, CellSphericalPolar * );

  //! Set the angle defining the adhesive regions.  Precalculate sin and cos.
  void SetPoleRegionAngle( double angle )
  { mAngle = angle; mSinAngle = std::sin(mAngle); mCosAngle = std::cos(mAngle); };

  //! Get the overall adhesion energy of this cell.
  double GetAdhesionEnergy();

  //!TEST FUNCTION FOR THE METROPOLIS ALGORITHM.
  double GetAdhesionEnergyTest();

  //! Define the hdf5 data format of this class.
  //  Inserts members into the given compound type typeDefinition.
  static void HDF5DataFormat( H5::CompType & typeDefinition );

  //! Builds up the H5::CompType from input data
  static H5::CompType ParseHDF5DataFormat( H5::CompType & inputType,
                                           std::stringstream & errors,
                                           std::stringstream & warnings );


  Vector3f mPolarDirection;    //!< The current orientation of the cell.
  Vector3f mNewPolarDirection; //!< The calculated new orientation of the cell,
                               //   to be updated in Reset() in the next time
                               //   step.

  double mAdhesionEnergy; //!< The adhesion energy.

  CSListContainer<ModelElement *> *mpContactsPointers;  //!< List of all CellSphericalPolar objects in contact with 'this'.

  static double adhesionEnergyMax;

  double maxRotationAngleForMetropolis;

  void setMaxRotationAngleForMetropolis(double maxAngle){
        maxRotationAngleForMetropolis = maxAngle;
    }

protected:
  double mAngle;    //!< The opening angle of the cone defining the adhesion area.
  double mSinAngle; //!< Pre-calculated value for overlap calculations.
  double mCosAngle; //!< Pre-calculated value for overlap calculations.

};
/* UNIT TEST GetAdhesionEnergy function
energy = scalar product of the two axis */
inline
double
CellSphericalPolar::GetAdhesionEnergyTest()
{
    double newAdhesionEnergy =0;
    ModelElement ** contactIterator = mpContactsPointers->begin();
    while ( contactIterator != mpContactsPointers->end() )
    {
        // Preliminary:
        // only take pairs of CellSphericalPolar into account
        if ( (*contactIterator)->mType != ModelElement::TypeCellSphericalPolar )
        continue;

        // calculate adhesion energy between cell and *contactIt
        // first the angle
        // assuming both mPolarDirection vectors are normalized
        CellSphericalPolar* tmpCell = static_cast<CellSphericalPolar *>(*contactIterator);
        newAdhesionEnergy +=  adhesionEnergyMax * fmax( this->mPolarDirection.x * tmpCell->mPolarDirection.x
                            + this->mPolarDirection.y * tmpCell->mPolarDirection.y
                            + this->mPolarDirection.z * tmpCell->mPolarDirection.z ,0);

        ++contactIterator;
    }

    return newAdhesionEnergy;
}

//! Metropolis algorithm for optimising the adhesion energy
/*! Randomly tilt the polar direction by one degree, calculate adhesion energy
 *  and compare it with the adhesion energy calculated before.  Use the new
 *  polar direction if the new adhesion energy is larger, or, if not, randomly
 *  with a probability of the Boltzmann Factor exp(newAdhesionEnergy). */
inline
void
CellSphericalPolar::OrientationMetropolis()
//! From D.Drasdo Chapter III.1: Center-based Single-cells Models ; from book Single-cell-based Models in Biology and Medicine 2007
//! 1. probe a new orientation
//! 2. trial this orientation according to a probability reflecting the relative time scales of orientation changes
//! 3. probability < exp(-∆E/Ft)}
{
    if ( !mNumContacts )
        return;

    double newAdhesionEnergy;

    // t1m-debug
    // double energyChange = 0;
    //!t1m-debug

    Vector3f oldPolarAxis( mPolarDirection );

    probeRandomRotation( maxRotationAngleForMetropolis );

    newAdhesionEnergy = GetAdhesionEnergy();

    // t1m-debug
    // std::cout<<newAdhesionEnergy<<" new"<< std::endl;
    //if (newAdhesionEnergy!=newAdhesionEnergy)
    //    double energyDifference = newAdhesionEnergy -   mAdhesionEnergy;
    //!t1m-debug

    double energyDifference = newAdhesionEnergy - mAdhesionEnergy;

    // t1m-debug
    //cout<<mAdhesionEnergy<<" old"<<endl;
    //cout<<energyDifference<<" Delta"<<endl;
    // if (energyDifference) std::cout << "OrientationMetropolis: cell rotation changed energy by " << energyDifference << " internal units.\n";
    //!t1m-debug

    if ( energyDifference < 0 )
    {
        if ( mpRandom->GetRandomUniform01() < std::exp(1000*energyDifference) )
        {
            // accept newPolarAxis
            mNewPolarDirection.Set( mPolarDirection );
            // t1m-debug
            //std::cout << "accepted energy change:  " << energyDifference << std::endl;
            //energyChange += energyDifference;
            //!t1m-debug
        }
        else
        {
            // use original axis
            mNewPolarDirection.Set( oldPolarAxis );
            //std::cout << "not accepted\n" << std::endl;
        }
    }
    else
    {
        mNewPolarDirection.Set( mPolarDirection );

        // t1m-debug
        // energyChange += energyDifference;
        //!t1m-debug
    }
    mPolarDirection.Set(mNewPolarDirection);

    //TIM:
    // mPolarDirection.Set( oldPolarAxis );

    // t1m-debug
    // if ( energyChange )
    // std::cout << "OrientationMetropolis:  energyChange = " << energyChange << " (positive is better)\n" << std::flush;
    // !t1m-debug
}


//! Randomly tilt the polar direction by a given angle.
inline
void
CellSphericalPolar::probeRandomRotation( double maxAngle )
{
    // the direction of inclination
    double phi;
    // the inclination of the new polar direction wrt the current mPolarDirection
    double theta;

    // the unit vector perpendicular to the mPolarDirection around which to
    // perform the rotation
    Vector3f thetaRotationVector;

    // randomize a rotation within maxAngle (rad) of mPolarDirection

    // phi yields the direction into which the new polar axis is inclined to
    //   with respect to mPolarDirection.
    //   It is relative to one of the two vectors defined by cutting the
    //   plane perpendicular to mPolarDiraction with the x-y plane.

    phi = 2 * M_PI * mpRandom->GetRandomUniform01();

    // Rotation vector for angle theta, i.e. we rotate mPolarDirection on
    //  the plane perpendicular to this vector.  Therefore the vector also
    //  has to be perpendicular to mPolarDirection.
    // vector perpendicular to mPolarDirection, lying in x-y plane (normalized):
    double x = ( mPolarDirection[0] ) ?  mPolarDirection[0] : 1.;
    thetaRotationVector.Set( -mPolarDirection[1]/x, 1., 0. );
    thetaRotationVector.Normalize();

    thetaRotationVector = Rotate( thetaRotationVector, mPolarDirection, phi );

    // Not working - yields angles much too small:
    // Sample theta out of [0, thetaMax] with the distribution of sin²alpha dalpha.
    // Cf. the accepted answer at
    // http://math.stackexchange.com/questions/131336/uniform-random-quaternion-in-a-restricted-angle-range
    // do {
    //     theta = sqrtThree * thetaMaxCubed * m.Random->GetRandomUniform01();
    // } while ( std::sin(theta)*std::sin(theta)/(theta*theta) < mpRandom->GetRandomUniform01() );

    // using the density function 2*M_PI*r*sin(theta) (r=1, 0<=theta<thetaMax) we have a distribution
    // of F(theta) = (1-cos(theta))/(1-cos(thetaMax)).
    // So, we can generate a random number x uniformly distributed over [0,1) and calculate our angle
    // with the inverse function theta = G(x) = acos( 1 - x * (1 - cos(thetaMax)) );
    theta = std::acos( 1 - mpRandom->GetRandomUniform01() * (1-std::cos(maxAngle)) );

    // alternatively one can assume sin(x) ~ x for the angles we sample.
    // then the distribution is F(theta) = theta^2 / thetaMax^2, which gets you G(x) = thetaMax * Sqrt(x);
    // theta = thetaMax * std::sqrt(mRandom.GetRandomUniform01());

    mPolarDirection = Rotate( mPolarDirection, thetaRotationVector, theta );
    mPolarDirection.Normalize();
}


//! Calculate overlap of adhesive regions of two CellSphericalPolar cells.
/*! The area overlap of two spherical caps is approximated by the overlap of the
    circles defined by projecting the caps onto the plane perpendicular to the
    distance vector between the centres of the two cells.  The actual contact
    area defines a third circle.

    Then the overlap of the three circles is calculated according to\n
      MP Fewell:  _Area of Common Overlap of Three Circles_; Oct. 2006,
                  DSTO-TN-0722, Defense Science and Technology Organisation of
                  the Australian Government (Defense Department).
*/
/*! \returns The overlap as a value between 0 and 1, the latter when one cap completely */
inline
double
CellSphericalPolar::CalculateOverlap(CellSphericalPolar *cell1,CellSphericalPolar * cell2, double contactRadius )
{
    // calculate overlap radius and angles
    Vector3f distanceVector( cell2->position.x - cell1->position.x,
                             cell2->position.y - cell1->position.y,
                             cell2->position.z - cell1->position.z  );

    double distanceSquared = distanceVector.x * distanceVector.x
                           + distanceVector.y * distanceVector.y
                           + distanceVector.z * distanceVector.z;

    double distance = distanceVector.Norm();

    distanceVector.Multiply( 1/distance );

    double contactRadiusSquared;

    if ( contactRadius == 0 )
    {
        // calculate geometrical contact radius:
        contactRadiusSquared = (cell1->mRadius + cell2->mRadius + distance)
            * (cell1->mRadius + cell2->mRadius - distance)
            * (cell1->mRadius - cell2->mRadius - distance)
            * (-cell1->mRadius + cell2->mRadius - distance);
        contactRadiusSquared /= 4*distanceSquared;

        if ( contactRadiusSquared == 0. )
            return 0.;

        contactRadius = std::sqrt( contactRadiusSquared );
    }
    else
    {
        contactRadiusSquared = contactRadius * contactRadius;
    }

    double r1DotDistanceVector = cell1->mPolarDirection.x * distanceVector.x
                               + cell1->mPolarDirection.y * distanceVector.y
                               + cell1->mPolarDirection.z * distanceVector.z;

    double r2DotDistanceVector = cell2->mPolarDirection.x * distanceVector.x
                               + cell2->mPolarDirection.y * distanceVector.y
                               + cell2->mPolarDirection.z * distanceVector.z;

    // see if both polar axis vector are inclined to the other cell (or both
    // away).  If not, the second vector direction is flipped.  For symmetric
    // reasons, it does not matter for the latter projection if they incline to
    // the respective other cell or away from it, as long as they both are
    // inclined that way.
    double cellRadius2;
    if ( SIGN(r1DotDistanceVector) == SIGN(r2DotDistanceVector) )
        cellRadius2 = - cell2->mRadius;
    else
        cellRadius2 = cell2->mRadius;


    // project the center of the spherical cap (i.e. the center of the base)
    // spanned by the adhesive region onto the cell-cell contact plane.
    // The formula was taken from tools/math/mathematics.cpp (see reference there)
    // but reduced to our purposes.
    Vector3f projectedCapCenter1( cell1->mPolarDirection.x - distanceVector.x * r1DotDistanceVector,
                                  cell1->mPolarDirection.y - distanceVector.y * r1DotDistanceVector,
                                  cell1->mPolarDirection.z - distanceVector.z * r1DotDistanceVector );
    projectedCapCenter1.Multiply( cell1->mCosAngle * cell1->mRadius );
    double projectedCapRadius1 = cell1->mRadius * cell1->mSinAngle * std::fabs(r1DotDistanceVector);

    double projectedCapCenter1Squared = projectedCapCenter1.x * projectedCapCenter1.x
                                      + projectedCapCenter1.y * projectedCapCenter1.y
                                      + projectedCapCenter1.z * projectedCapCenter1.z;
    double projectedCapCenter1Norm = std::sqrt( projectedCapCenter1Squared );

    if ( projectedCapCenter1Norm - projectedCapRadius1 > contactRadius )
        return 0.;


    Vector3f projectedCapCenter2( cell2->mPolarDirection.x - distanceVector.x * r2DotDistanceVector,
                                  cell2->mPolarDirection.y - distanceVector.y * r2DotDistanceVector,
                                  cell2->mPolarDirection.z - distanceVector.z * r2DotDistanceVector );
    projectedCapCenter2.Multiply( cell2->mCosAngle * cellRadius2 );
    double projectedCapRadius2 = cell2->mRadius * cell2->mSinAngle * std::fabs(r2DotDistanceVector);

    double projectedCapCenter2Squared = projectedCapCenter2.x * projectedCapCenter2.x
                                      + projectedCapCenter2.y * projectedCapCenter2.y
                                      + projectedCapCenter2.z * projectedCapCenter2.z;
    double projectedCapCenter2Norm = std::sqrt( projectedCapCenter2Squared );

    if ( projectedCapCenter2Norm - projectedCapRadius2 > contactRadius )
        return 0.;


    // calculate overlap...
    Vector3f circleCenterDistanceVector( projectedCapCenter2.x - projectedCapCenter1.x,
                                         projectedCapCenter2.y - projectedCapCenter1.y,
                                         projectedCapCenter2.z - projectedCapCenter1.z );

    double circleCenterDistanceSquared = circleCenterDistanceVector.x * circleCenterDistanceVector.x
                                       + circleCenterDistanceVector.y * circleCenterDistanceVector.y
                                       + circleCenterDistanceVector.z * circleCenterDistanceVector.z;

    double circleCenterDistance = std::sqrt( circleCenterDistanceSquared );

    if ( projectedCapRadius2 + projectedCapRadius1 <= circleCenterDistance )
        return 0.;

    circleCenterDistanceVector.Multiply( 1/circleCenterDistance );

    // needed later, so we introduce the intermediates
    double projectedCapRadius1Squared = projectedCapRadius1*projectedCapRadius1;
    double projectedCapRadius2Squared = projectedCapRadius2*projectedCapRadius2;

    // t1m-debug
    //if ( !projectedCapRadius1Squared || !projectedCapRadius2Squared)
    //    projectedCapRadius1Squared=0;
    //!t1m-debug

    double areaOverlap;
    double smallestCircleArea;
    if ( projectedCapRadius1 < projectedCapRadius2 )
    {
        // see if circle 1 is contained in circle 2
        if ( circleCenterDistance < projectedCapRadius2 - projectedCapRadius1 )
        {
            // if circle 1 is also completely contained in the contact circle, we return full overlap.
            if ( projectedCapCenter1Norm + projectedCapRadius1 <= contactRadius )
                return 1.;
            else if (projectedCapCenter1Norm + contactRadius <= projectedCapRadius1 )
            {
                return contactRadiusSquared/projectedCapRadius1Squared;
            }
            else
            {
                // calculate and return the lens of intersection between circle 1 and the contact circle
                areaOverlap = contactRadius * std::acos( (projectedCapCenter1Squared + contactRadiusSquared - projectedCapRadius1Squared) / (2*projectedCapCenter1Norm*contactRadius) )
                            + projectedCapRadius1 * std::acos( (projectedCapCenter1Squared - contactRadiusSquared + projectedCapRadius1Squared) / (2*projectedCapCenter1Norm*projectedCapRadius1) )
                            + .5 * std::sqrt( (-contactRadius + projectedCapRadius1 + projectedCapCenter1Norm)
                                              *(contactRadius - projectedCapRadius1 + projectedCapCenter1Norm)
                                              *(contactRadius + projectedCapRadius1 - projectedCapCenter1Norm)
                                              *(contactRadius + projectedCapRadius1 + projectedCapCenter1Norm) );

                // t1m-debug
                // if (areaOverlap!=areaOverlap)
                //     areaOverlap=0.;
                //!t1m-debug

                return areaOverlap/( 2 * M_PI * projectedCapRadius1Squared );
            }
        }

        smallestCircleArea = 2 * M_PI * projectedCapRadius1Squared;
    }
    else // if projectedCapRadius2 <= projectedCapRadius1
    {
        // see if circle 1 is contained in circle 2
        if ( circleCenterDistance < projectedCapRadius1 - projectedCapRadius2 )
        {
            // if circle 1 is also completely contained in the contact circle, we return full overlap.
            if ( projectedCapCenter2Norm + projectedCapRadius2 <= contactRadius )
                return 1.;
            else if (projectedCapCenter2Norm + contactRadius <= projectedCapRadius2 )
            {
                return contactRadiusSquared/projectedCapRadius2Squared;
            }
            else
            {
                // calculate the lens of intersection between circle 2 and the contact circle
                areaOverlap = contactRadius * std::acos( (projectedCapCenter2Squared + contactRadiusSquared - projectedCapRadius2Squared) / (2*projectedCapCenter2Norm*contactRadius) )
                            + projectedCapRadius2 * std::acos( (projectedCapCenter2Squared - contactRadiusSquared + projectedCapRadius2Squared) / (2*projectedCapCenter2Norm*projectedCapRadius2) )
                            + .5 * std::sqrt( (-contactRadius + projectedCapRadius2 + projectedCapCenter2Norm)
                                              *(contactRadius - projectedCapRadius2 + projectedCapCenter2Norm)
                                              *(contactRadius + projectedCapRadius2 - projectedCapCenter2Norm)
                                              *(contactRadius + projectedCapRadius2 + projectedCapCenter2Norm) );
                // t1m-debug
                // if (areaOverlap!=areaOverlap) // NaN
                //     areaOverlap=0.;
                //!t1m-debug
                return areaOverlap / (2 * M_PI * projectedCapRadius2Squared);
            }
        }

        smallestCircleArea = 2 * M_PI * projectedCapRadius2Squared;
    }


    // construct the two intersection points of the two adhesive area circles
    // wrt the center of the contact area circle (of which the center is used as origin):

    double s =  std::sqrt( (circleCenterDistance + projectedCapRadius1 + projectedCapRadius2)
                           * ( -circleCenterDistance + projectedCapRadius1 + projectedCapRadius2)
                           * ( -circleCenterDistance - projectedCapRadius1 + projectedCapRadius2)
                           * ( -circleCenterDistance + projectedCapRadius1 - projectedCapRadius2));
    double d1 = circleCenterDistanceSquared + projectedCapRadius1Squared - projectedCapRadius2Squared;


    // one of the two vectors in the plane plane (=>perpendicular to distanceVector) that is
    // perpendicular to circleCenterDistanceVector, i.e. distanceVector x circleCenterDistanceVector.
    Vector3f sVector ( distanceVector.y * circleCenterDistanceVector.z - distanceVector.z * circleCenterDistanceVector.y,
                       distanceVector.z * circleCenterDistanceVector.x - distanceVector.x * circleCenterDistanceVector.z,
                       distanceVector.x * circleCenterDistanceVector.y - distanceVector.y * circleCenterDistanceVector.x );

    // scale to length s
    sVector.Normalize();
    sVector.Multiply( s );


    Vector3f intersectionPoint1( (d1 * circleCenterDistanceVector.x + sVector.x)/(2 * circleCenterDistance),
                                 (d1 * circleCenterDistanceVector.y + sVector.y)/(2 * circleCenterDistance),
                                 (d1 * circleCenterDistanceVector.z + sVector.z)/(2 * circleCenterDistance) );
    intersectionPoint1.Add( projectedCapCenter1 );

    double intersectionDistanceSquared1 = intersectionPoint1.x * intersectionPoint1.x
                                        + intersectionPoint1.y * intersectionPoint1.y
                                        + intersectionPoint1.z * intersectionPoint1.z;


    Vector3f intersectionPoint2( (d1 * circleCenterDistanceVector.x - sVector.x)/(2 * circleCenterDistance),
                                 (d1 * circleCenterDistanceVector.y - sVector.y)/(2 * circleCenterDistance),
                                 (d1 * circleCenterDistanceVector.z - sVector.z)/(2 * circleCenterDistance) );
    intersectionPoint2.Add( projectedCapCenter1 );

    double intersectionDistanceSquared2 = intersectionPoint2.x * intersectionPoint2.x
                                        + intersectionPoint2.y * intersectionPoint2.y
                                        + intersectionPoint2.z * intersectionPoint2.z;

    Vector3f * involvedIntersectionPointPtr = &intersectionPoint1;
    if ( intersectionDistanceSquared1 < contactRadiusSquared )
    {
        // we have adhesion area overlap + if the same was true for intersectionDistanceSquared2, the
        // whole circle-circle overlap can be calculated (which is the simpler case!)
        if ( intersectionDistanceSquared2 < contactRadiusSquared )
        {
            d1 /= 2*circleCenterDistance;
            areaOverlap =
                  projectedCapRadius1Squared * std::acos( d1/projectedCapRadius1 )
                + projectedCapRadius2Squared * std::acos( (circleCenterDistance - d1)/projectedCapRadius2 )
                + .5 * std::sqrt(   ( circleCenterDistance + projectedCapRadius1 + projectedCapRadius2 )
                                  * ( circleCenterDistance + projectedCapRadius1 - projectedCapRadius2 )
                                  * ( circleCenterDistance - projectedCapRadius1 + projectedCapRadius2 )
                                  * (-circleCenterDistance + projectedCapRadius1 + projectedCapRadius2 ) );
            // t1m-debug
            // if (areaOverlap!=areaOverlap) // NaN
            //     areaOverlap=0.;
            //!t1m-debug
            return areaOverlap/smallestCircleArea;
        }
    }
    else if ( intersectionDistanceSquared2 < contactRadiusSquared )
        involvedIntersectionPointPtr = &intersectionPoint2;
    else
        return 0.;

    Vector3f & involvedIntersectionPoint = * involvedIntersectionPointPtr;
    // calculate circular triangle area plus chord areas with intersectionPoint1

    // calculate the relevant intersection points with the contact area circle, the intersection
    // of circle 1 with the contact area circle that lies within circle 2 and vice versa.

    // contact area intersecting circle1
    s = std::sqrt(   (  projectedCapCenter1Norm + projectedCapRadius1 + contactRadius )
                     * ( -projectedCapCenter1Norm + projectedCapRadius1 + contactRadius )
                     * ( -projectedCapCenter1Norm - projectedCapRadius1 + contactRadius )
                     * ( -projectedCapCenter1Norm + projectedCapRadius1 - contactRadius ) );
    d1 = (projectedCapCenter1Squared + contactRadiusSquared - projectedCapRadius1Squared)/projectedCapCenter1Norm;

    // the 'height' vector of the intersection lens
    sVector.Set( distanceVector.y * projectedCapCenter1.z - distanceVector.z * projectedCapCenter1.y,
                 distanceVector.z * projectedCapCenter1.x - distanceVector.x * projectedCapCenter1.z,
                 distanceVector.x * projectedCapCenter1.y - distanceVector.y * projectedCapCenter1.x );
    sVector.Normalize();
    sVector.Multiply( s );

    // the two intersection points between circle 1 and the contact area circle
    Vector3f x01a( (d1 * projectedCapCenter1.x - sVector.x)/(2*projectedCapCenter1Norm),
                   (d1 * projectedCapCenter1.y - sVector.y)/(2*projectedCapCenter1Norm),
                   (d1 * projectedCapCenter1.z - sVector.z)/(2*projectedCapCenter1Norm) );

    Vector3f x01b( (d1 * projectedCapCenter1.x + sVector.x)/(2*projectedCapCenter1Norm),
                   (d1 * projectedCapCenter1.y + sVector.y)/(2*projectedCapCenter1Norm),
                   (d1 * projectedCapCenter1.z + sVector.z)/(2*projectedCapCenter1Norm) );


    // contact area intersecting circle2

    // s is the length of the lens created by the contact circle and
    // circle 2 (times 2*projectedCapCenter2Squared, which will be
    // diveded later)
    s = std::sqrt( (  projectedCapCenter2Norm + projectedCapRadius2 + contactRadius )
                   * ( -projectedCapCenter2Norm + projectedCapRadius2 + contactRadius )
                   * ( -projectedCapCenter2Norm - projectedCapRadius2 + contactRadius )
                   * ( -projectedCapCenter2Norm + projectedCapRadius2 - contactRadius ) );
    d1 = (projectedCapCenter2Squared + contactRadiusSquared - projectedCapRadius2Squared)/projectedCapCenter2Norm;

    // the 'height' vector of the intersection lens
    sVector.Set( distanceVector.y * projectedCapCenter2.z - distanceVector.z * projectedCapCenter2.y,
                 distanceVector.z * projectedCapCenter2.x - distanceVector.x * projectedCapCenter2.z,
                 distanceVector.x * projectedCapCenter2.y - distanceVector.y * projectedCapCenter2.x );
    sVector.Normalize();
    sVector.Multiply( s );

    // the two intersection points between circle 2 and the contact area circle
    Vector3f x02a( (d1 * projectedCapCenter2.x - sVector.x)/(2*projectedCapCenter2Norm),
                   (d1 * projectedCapCenter2.y - sVector.y)/(2*projectedCapCenter2Norm),
                   (d1 * projectedCapCenter2.z - sVector.z)/(2*projectedCapCenter2Norm) );

    Vector3f x02b( (d1 * projectedCapCenter2.x + sVector.x)/(2*projectedCapCenter2Norm),
                   (d1 * projectedCapCenter2.y + sVector.y)/(2*projectedCapCenter2Norm),
                   (d1 * projectedCapCenter2.z + sVector.z)/(2*projectedCapCenter2Norm) );

    Vector3f * x01ptr, * x02ptr;

    Vector3f intersection01aDistanceFromCenter2( x01a.x - projectedCapCenter2.x,
                                                 x01a.y - projectedCapCenter2.y,
                                                 x01a.z - projectedCapCenter2.z );
    Vector3f intersection01bDistanceFromCenter2( x01b.x - projectedCapCenter2.x,
                                                 x01b.y - projectedCapCenter2.y,
                                                 x01b.z - projectedCapCenter2.z );
    if ( intersection01aDistanceFromCenter2.Norm() < intersection01bDistanceFromCenter2.Norm() )
        x01ptr = &x01a;
    else
        x01ptr = &x01b;

    Vector3f intersection02aDistanceFromCenter1( x02a.x - projectedCapCenter1.x,
                                                 x02a.y - projectedCapCenter1.y,
                                                 x02a.z - projectedCapCenter1.z );
    Vector3f intersection02bDistanceFromCenter1( x02b.x - projectedCapCenter1.x,
                                                 x02b.y - projectedCapCenter1.y,
                                                 x02b.z - projectedCapCenter1.z );

    if ( intersection02aDistanceFromCenter1.Norm() < intersection02bDistanceFromCenter1.Norm() )
        x02ptr = &x02a;
    else
        x02ptr = &x02b;

    Vector3f & x01 = *x01ptr;
    Vector3f & x02 = *x02ptr;

    // the chords of the circular triangle, c0 is opposite to
    // involvedIntersectionPoint (x12), c1 opp. to x02 and c2 opp to x01
    double c0 = std::sqrt( (  x01.x - x02.x)*(x01.x - x02.x)
                           + (x01.y - x02.y)*(x01.y - x02.y)
                           + (x01.z - x02.z)*(x01.z - x02.z) );
    double c1 = std::sqrt(  (  involvedIntersectionPoint.x - x01.x)*(involvedIntersectionPoint.x - x01.x)
                            + (involvedIntersectionPoint.y - x01.y)*(involvedIntersectionPoint.y - x01.y)
                            + (involvedIntersectionPoint.z - x01.z)*(involvedIntersectionPoint.z - x01.z) );
    double c2 = std::sqrt(  (  involvedIntersectionPoint.x - x02.x)*(involvedIntersectionPoint.x - x02.x)
                            + (involvedIntersectionPoint.y - x02.y)*(involvedIntersectionPoint.y - x02.y)
                            + (involvedIntersectionPoint.z - x02.z)*(involvedIntersectionPoint.z - x02.z) );

    // ToDo:  Case in which more than half a circle is involved.

    // spherical triangle area
    areaOverlap = .5 * sqrt( (c0 + c1 + c2)*(c0 + c1 - c2)*(c0 - c1 + c2)*(-c0 + c1 + c2) )
                + contactRadiusSquared        * std::asin( c0 / (2*contactRadius) )
                + projectedCapRadius1Squared  * std::asin( c1 / (2*projectedCapRadius1) )
                + projectedCapRadius2Squared  * std::asin( c2 / (2*projectedCapRadius2) )
                - .25 * ( (  c0 * std::sqrt(4 * contactRadiusSquared - c0*c0 ) )
                          + (c1 * std::sqrt(4 * projectedCapRadius1Squared - c1*c1) )
                          + (c2 * std::sqrt(4 * projectedCapRadius2Squared - c2*c2) ));
    // t1m-debug
    // if (areaOverlap!=areaOverlap) // NaN
    //     areaOverlap=0.;
    //!t1m-debug

    return areaOverlap/smallestCircleArea;
}


inline
double
CellSphericalPolar::CalculateOverlapSimple(CellSphericalPolar *cell1, CellSphericalPolar *cell2)
{
    Vector3f distanceVector = cell1->position - cell2->position;

    distanceVector.Normalize();

    double polarAxis1DotDistanceVector = std::fabs( cell1->mPolarDirection[0] * distanceVector[0] +
                                                    cell1->mPolarDirection[1] * distanceVector[1] +
                                                    cell1->mPolarDirection[2] * distanceVector[2] );

    if ( polarAxis1DotDistanceVector < cell1->mCosAngle )
        return 0.;

    double polarAxis2DotDistanceVector = std::fabs( cell2->mPolarDirection[0] * distanceVector[0] +
                                                    cell2->mPolarDirection[1] * distanceVector[1] +
                                                    cell2->mPolarDirection[2] * distanceVector[2] );

    if ( polarAxis2DotDistanceVector < cell2->mCosAngle )
        return 0.;

    return 1.;
}


//! Get the overall adhesion energy of this cell.
/*! Each CellSphericalPolar-CellSphericalPolar interaction contributes to the
    adhesion energy by the value
      adhesionEnergyMax * CalculateOverlap( this, otherCell );
*/
inline
double
CellSphericalPolar::GetAdhesionEnergy()
{
    double newAdhesionEnergy =0;
    ModelElement ** contactIterator = mpContactsPointers->begin();
    while ( contactIterator != mpContactsPointers->end() )
    {
        // Preliminary:
        // only take pairs of CellSphericalPolar into account
        if ( (*contactIterator)->mType != ModelElement::TypeCellSphericalPolar )
            continue;

        // calculate adhesion energy between cell and *contactIt
        // first the angle
        // assuming both mPolarDirection vectors are normalized
        newAdhesionEnergy += adhesionEnergyMax
                           * CalculateOverlap( this, static_cast<CellSphericalPolar *>(*contactIterator) );

        ++contactIterator;
    }

    return newAdhesionEnergy;
}


#endif
