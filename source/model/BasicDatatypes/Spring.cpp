///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Spring.cpp                                                           //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-05-31 22:00:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "Spring.h"

Spring::Spring()
{
	// Initialise sata 
	dlength = 0.;     // change of length
	
	direction[0]=0.;  // direction
	direction[1]=0.;
	direction[2]=0.;

	// Visualisation of cylinder
	// Color
	r = 255;
	g = 0;
	b = 0;
	radius = 0.005;
	transparency = 1.;

	// Cell division
	divideCircle = 1;
	
	// Additional damper
	//l_damper = 0.;
  length = 0.;
}


void Spring::copyfrom(Spring *spring_basic){


  this->mIndex = spring_basic->mIndex;

	start = spring_basic->start;
	end = spring_basic->end;

	l0 = spring_basic->l0;

	l0_init = spring_basic->l0_init;

	k = spring_basic->k;
	k2 = spring_basic->k2;
	nu = spring_basic->nu;
	
	direction[0] = spring_basic->direction[0];
	direction[1] = spring_basic->direction[1];
	direction[2] = spring_basic->direction[2];
	
	dlength = spring_basic->dlength;

//	l_damper = spring_basic->l_damper;
	eta_damper = spring_basic->eta_damper;
//	second_damper = spring_basic->second_damper;
	
    r = spring_basic->r;
    g = spring_basic->g;
    b = spring_basic->b;
	transparency = spring_basic->transparency;

	radius = spring_basic->radius;

	divideCircle = spring_basic->divideCircle;

}

void Spring::SetColor(unsigned char r_, unsigned char g_, unsigned char b_)
{
	r = r_;
	g = g_;
	b = b_;
}
