///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Spring.h                                                             //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-05-31 22:00:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef SPRING_H
#define SPRING_H

class MassPoint;
class Triangle;

class Spring
{
public:

  // Connect mass points to spring
  MassPoint *start;
  MassPoint *end;

  short mIndex;

  // Length and direction
  double l0;            // Length of spring
  double direction[3];
  double dlength;       // change of length
  double length;
  double l0_init;       // stores initial l0

  double F_debug_secondDamper;

	// Constants of viscoelatic elemtent
	double k;             //spring constant
	double nu;
	
	// Additional damper
//	double l_damper;
	double eta_damper;
//	double second_damper;

	//second spring
	double k2;

	double force_KV[3];

	bool interacts;
	Triangle *interacts_with;


	// Force
	double F;
	double force[3];

	double store_k;
  double store_old_k;
	double store_l0;
  double store_old_l0;

	// Visualisation
	// Color
    unsigned char r,g,b;
	double transparency;
	void SetColor(unsigned char r, unsigned char g, unsigned char b);// Changes MassPoint color
	double radius;    // Radius of the cylinder

	// Cell-division
	int divideCircle; // Diecide witch side during cell division
		
	Spring();         // Initialise
	void copyfrom(Spring * spring_basic); // Copy a spring
	
  unsigned int index;

};

#endif
