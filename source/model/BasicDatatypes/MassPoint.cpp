///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  MassPoint.cpp                                                        //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-05-31 22:00:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/**
 * @file	MassPoint.cpp
 * @date	31.05.2012
 * @author	Johannes Neitsch (MassPoint.cpp)
 * @brief	mass point
 *
 * extended description
 *
 */

/** short funtion discription
 * @param	paramter_name	parameter_description
 * @return	return_description
 *
 * extended description
 *
 */


#include "MassPoint.h"

#include "../../tools/math/mathematics.h"

class Triangle;

MassPoint::MassPoint()
{
	// velocity
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;

	// force
	F[0] = 0.0;
	F[1] = 0.0;
	F[2] = 0.0;

	// acceleration
	a[0] = 0.0;
	a[1] = 0.0;
	a[2] = 0.0;

	// old acceleration
	a_old[0] = 0.0;
	a_old[1] = 0.0;
	a_old[2] = 0.0;

	// mass of this point and its inverse
	point_mass = 1.0;
	point_mass_invers = 1.0;

	// voronou area
	mAreaVoronoi = 0.;

	// visualisation
	r = 0;
	g = 0;
	b = 255;
	transparency = 1.;

	// cell division
	divide_number = 0;

  radius = 0.007;

  this->F_repulse[0] = 0.;
  this->F_repulse[1] = 0.;
  this->F_repulse[2] = 0.;


}


MassPoint::~MassPoint()
{
	neighbourSprings.clear();

}

void MassPoint::SetColor(unsigned char r_, unsigned char g_, unsigned char b_)
{
	r = r_;
	g = g_;
	b = b_;
}


void MassPoint::copyfrom(MassPoint *mass_basic){

  this->mSpringsOnSurface = mass_basic->mSpringsOnSurface;

  this->mIndex = mass_basic->mIndex;

  position[0] = mass_basic->position[0];
  position[1] = mass_basic->position[1];
  position[2] = mass_basic->position[2];

  point_mass = mass_basic->point_mass;
  point_mass_invers = mass_basic->point_mass_invers;

  r = mass_basic->r;
  g = mass_basic->g;
  b = mass_basic->b;

  transparency =  mass_basic->transparency;
  mAreaVoronoi = mass_basic->mAreaVoronoi;
/*
  for( unsigned int i = 0 ; i < mass_basic->mvpTriangles.size() ; i++){

    this->mvpTriangles.push_back( mass_basic->mvpTriangles[i] );
    this->mvpPointsL.push_back( mass_basic->mvpPointsL[i] );
    this->mvpPointsR.push_back( mass_basic->mvpPointsR[i] );

  }
*/




}


void MassPoint::changeLength(double m_position[3], double l_c){
	
	// moves a new cell to a new point and rescales (all) springs with the same value
	double v[3] = { position[0]-m_position[0] , position[1]-m_position[1] , position[2]-m_position[2]};
	double d = norm(v);

	v[0] = v[0]*l_c/d;
	v[1] = v[1]*l_c/d;
	v[2] = v[2]*l_c/d;

	
	position[0] = m_position[0]+v[0];
	position[1] = m_position[1]+v[1];
	position[2] = m_position[2]+v[2];

}
