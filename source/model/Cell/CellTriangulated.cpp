///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CellTriangulated.cpp                                                 //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-05-31 22:00:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#include "../../Core.h"


#include "CellTriangulated.h"

#include "../../tools/math/mathematics.h"
#include "../../tools/sort/sort.h"
#include "../../tools/random/Random.h"
#include "../../tools/Tools.h"

#include "../BasicDatatypes/MassPoint.h"
#include "../BasicDatatypes/Spring.h"
#include "../BasicDatatypes/Triangle.h"
#include "../BasicDatatypes/Triangle.cpp"

#include "../tools/new_triangulation/Geometry/Triangulation.hpp"
#include "../../gui/GLTools/TriangulatedCell.h"

#include "../../tools/new_triangulation/LinearAlgebra/SparseMatrix.hpp"
#include "../../tools/new_triangulation/LinearAlgebra/SparseVector.hpp"
#include "../../tools/new_triangulation/LinearAlgebra/Solver.hpp"

using namespace std;
#include <iostream>
#include <fstream>
#include <cmath>

#include <queue>

#define PI 3.141592653589793238462643383279502

#define mod_3( a) (a<3 ? a : (a-3))

CellTriangulated::CellTriangulated()
    : Cell(),
      ModelElementTriangulatedCell( 0,0,0 )
{
  mType = ModelElement::TypeCellTriangulated;

  can_divide = 0;
  time_relax=0.;
  mpGLObject = new TriangulatedCell( this );

  pressure_high = 0;
  pressure_low_since = pressure_last_number;
  able_to_move = 1;

  count_stretch = 0;
 A = NULL;

 eta = 10;
 gamma = 0.8;


}

//color
void CellTriangulated::color_triangle_deform(){

	double min = 20;//cell->radius_cell;
	double max = 0;

	for( unsigned int i = 0 ; i < mass_huell.size() ;i++){


		double tmp = norm(mass_huell[i]->position) ;//- cell->radius_cell;

//		cout << tmp << "\t" << min << "\t" << max << "\t" << (tmp<min) << "\t" << (tmp>max)<<  endl << endl;


		if ( tmp > max ){
			max = tmp;
		}
		if ( tmp < min ){
			min = tmp;
		}
	}

	double diff = max -min;


	for( unsigned int i = 0 ; i < mass_huell.size(); i++){


		double tmp = norm(mass_huell[i]->position);

		if (tmp < (min + 0.25 * diff)){
	    	  mass_huell[i]->r = 0.;
	    	  mass_huell[i]->g = 4 * (tmp - min) / diff *255;
	    	  mass_huell[i]->b = 1. *255;
		}
	    else  if (tmp < (min + 0.5 * diff)){
	    	mass_huell[i]->r = 0. ;
	    	mass_huell[i]->g = 1. *255;
	    	mass_huell[i]->b = (1 + 4 * (min + 0.25 * diff - tmp) / diff)*255;
	    }
        else  if (tmp < (min + 0.75 * diff)){
	    	mass_huell[i]->r = 4 * (tmp- min - 0.5 * diff) / diff*255;
	    	mass_huell[i]->g = 1.*255;
	    	mass_huell[i]->b = 0.;
        }
	    else {
	    	mass_huell[i]->r = 1.*255;
	    	mass_huell[i]->g = (1 + 4 * (min + 0.75 * diff - tmp) / diff)*255;
	    	mass_huell[i]->b = 0.;
	    }

	}

	for( unsigned int i = 0 ; i < triangle.size() ; i++){
		triangle[i]->transparency = 255;

		triangle[i]->r = (triangle[i]->points[0]->r+triangle[i]->points[1]->r+triangle[i]->points[2]->r)/3.;
		triangle[i]->g = (triangle[i]->points[0]->g+triangle[i]->points[1]->g+triangle[i]->points[2]->g)/3.;
		triangle[i]->b = (triangle[i]->points[0]->b+triangle[i]->points[1]->b+triangle[i]->points[2]->b)/3.;


	}


}

//force
void CellTriangulated::addInnerStress(){

	for( unsigned int i = 0 ; i < triangle.size() ; i++){


			triangle[i]->setArea(p);


			triangle[i]->points[0]->F[0] += triangle[i]->mNormalvector[0];
			triangle[i]->points[0]->F[1] += triangle[i]->mNormalvector[1];
			triangle[i]->points[0]->F[2] += triangle[i]->mNormalvector[2];

			triangle[i]->points[1]->F[0] += triangle[i]->mNormalvector[0];
			triangle[i]->points[1]->F[1] += triangle[i]->mNormalvector[1];
			triangle[i]->points[1]->F[2] += triangle[i]->mNormalvector[2];

			triangle[i]->points[2]->F[0] += triangle[i]->mNormalvector[0];
			triangle[i]->points[2]->F[1] += triangle[i]->mNormalvector[1];
			triangle[i]->points[2]->F[2] += triangle[i]->mNormalvector[2];



		}

		


}
void CellTriangulated::addConstantForce(double x, double y, double z){
	
	for( unsigned int i = 0 ; i < mass.size() ; i++){
		mass[i]->F[0] += x;
		mass[i]->F[1] += y;
		mass[i]->F[2] += z;
	}


}

/*
void CellTriangulated::updateDamper( double timeStep ){

	for(unsigned int i = 0 ; i < spring_huell.size() ; i++)
		spring_huell[i]->l_damper += 1/spring_huell[i]->eta_damper*spring_huell[i]->F*timeStep;

	for(unsigned int i = 0 ; i < spring_cytoscelett.size() ; i++)
		spring_cytoscelett[i]->l_damper += 1/spring_cytoscelett[i]->eta_damper*spring_cytoscelett[i]->F*timeStep;


}
*/

//calc
void CellTriangulated::calcAcceleration(int number){

	for(unsigned int i = 0 ; i < mass.size() ; i++){

		this->mass[i]->a_old[0] = this->mass[i]->a[0];
		this->mass[i]->a_old[1] = this->mass[i]->a[1];
		this->mass[i]->a_old[2] = this->mass[i]->a[2];

		this->mass[i]->a[0] = this->mass[i]->F[0] * this->mass[i]->point_mass_invers;
		this->mass[i]->a[1] = this->mass[i]->F[1] * this->mass[i]->point_mass_invers;
		this->mass[i]->a[2] = this->mass[i]->F[2] * this->mass[i]->point_mass_invers;


	}
}
void CellTriangulated::calcForce_froce3d_viscous(double timeStep){
	//force from springs
	for( unsigned int i = 0 ; i < spring.size() ; i++ ){
		spring[i]->force[0] +=	spring[i]->nu * (spring[i]->end->v[0] - spring[i]->start->v[0]) * spring[i]->direction[0] ;
		spring[i]->force[1] +=	spring[i]->nu * (spring[i]->end->v[1] - spring[i]->start->v[1]) * spring[i]->direction[1] ;
		spring[i]->force[2] +=	spring[i]->nu * (spring[i]->end->v[2] - spring[i]->start->v[2]) * spring[i]->direction[2] ;
		
		spring[i]->force_KV[0] = spring[i]->force[0];
		spring[i]->force_KV[1] = spring[i]->force[1];
		spring[i]->force_KV[2] = spring[i]->force[2];

	}
}
void CellTriangulated::calcForce_froce3d_spring(double timeStep){
	//force from liquid in Kelvin-Voigt
	for( unsigned int i = 0 ; i < spring.size() ; i++ ){
		spring[i]->force[0] = spring[i]->k * spring[i]->dlength;
		spring[i]->force[1] = spring[i]->k * spring[i]->dlength;
		spring[i]->force[2] = spring[i]->k * spring[i]->dlength;


	}
}

void CellTriangulated::updateSecondDamper(double timeStep){
	//change length of springs
	
/*	for( unsigned int i = 0 ; i < mass.size() ; i++){

		mass[i]->norm_force = norm(mass[i]->F);

	}
*/
	/*
	for( unsigned int i = 0 ; i < spring_huell.size() ; i++ ){
		
		//start
		double tmp[3];
		if( spring_huell[0]->start->norm_force != 0.){
			tmp[0] = spring_huell[0]->start->F[0] * spring_huell[0]->start->F[0] * spring_huell[i]->direction[0] / spring_huell[0]->start->norm_force;
			tmp[1] = spring_huell[1]->start->F[1] * spring_huell[1]->start->F[1] * spring_huell[i]->direction[1] / spring_huell[0]->start->norm_force;
			tmp[2] = spring_huell[2]->start->F[2] * spring_huell[2]->start->F[2] * spring_huell[i]->direction[2] / spring_huell[0]->start->norm_force;
		}

		double norm_tmp = norm(tmp);

		//end
		double direction_coefficient =  (dotmult(spring_huell[0]->start->F,spring_huell[i]->direction) > 0 )*2-1;//
		//spring_huell[i]->l0 += direction_coefficient/spring_huell[i]->eta_damper*spring[i]->F*timeStep;
		spring_huell[i]->l0 += direction_coefficient/spring_huell[i]->eta_damper*norm_tmp*timeStep;

		if( spring_huell[0]->end->norm_force != 0.){
			tmp[0] = spring_huell[0]->end->F[0] * spring_huell[0]->end->F[0] * spring_huell[i]->direction[0] / spring_huell[0]->end->norm_force;
			tmp[1] = spring_huell[1]->end->F[1] * spring_huell[1]->end->F[1] * spring_huell[i]->direction[1] / spring_huell[0]->end->norm_force;
			tmp[2] = spring_huell[2]->end->F[2] * spring_huell[2]->end->F[2] * spring_huell[i]->direction[2] / spring_huell[0]->end->norm_force;
		}

		norm_tmp = norm(tmp);

		direction_coefficient =  (dotmult(spring_huell[0]->end->F,spring_huell[i]->direction) > 0 )*2-1;//
		//spring_huell[i]->l0 += direction_coefficient/spring_huell[i]->eta_damper*spring[i]->F*timeStep;
		spring_huell[i]->l0 += direction_coefficient/spring_huell[i]->eta_damper*norm_tmp*timeStep;

		


	}
		for( unsigned int i = 0 ; i < spring_cytoscelett.size() ; i++ ){
		
		//start
		double tmp[3];
		

		if( spring_cytoscelett[0]->end->norm_force != 0.){
			tmp[0] = spring_cytoscelett[0]->end->F[0] * spring_cytoscelett[0]->end->F[0] * spring_cytoscelett[i]->direction[0] / spring_cytoscelett[0]->end->norm_force;
			tmp[1] = spring_cytoscelett[1]->end->F[1] * spring_cytoscelett[1]->end->F[1] * spring_cytoscelett[i]->direction[1] / spring_cytoscelett[0]->end->norm_force;
			tmp[2] = spring_cytoscelett[2]->end->F[2] * spring_cytoscelett[2]->end->F[2] * spring_cytoscelett[i]->direction[2] / spring_cytoscelett[0]->end->norm_force;
		}

		double norm_tmp = norm(tmp);

		double direction_coefficient =  (dotmult(spring_cytoscelett[0]->end->F,spring_cytoscelett[i]->direction) > 0 )*2-1;//
		//spring_huell[i]->l0 += direction_coefficient/spring_huell[i]->eta_damper*spring[i]->F*timeStep;
		spring_huell[i]->l0 += direction_coefficient/spring_cytoscelett[i]->eta_damper*norm_tmp*timeStep;

		


	}
	*/
	
	/*
	for( unsigned int i = 0 ; i < spring_huell.size() ; i++ ){

		spring_huell[i]->F = norm(spring_huell[i]->force);
		double tmp[3];
		if( spring_huell[i]->F != 0.){
			tmp[0] = spring_huell[i]->force[0] * spring_huell[i]->force[0] * spring_huell[i]->direction[0] / spring_huell[i]->F;
			tmp[1] = spring_huell[i]->force[1] * spring_huell[i]->force[1] * spring_huell[i]->direction[1] / spring_huell[i]->F;
			tmp[2] = spring_huell[i]->force[2] * spring_huell[i]->force[2] * spring_huell[i]->direction[2] / spring_huell[i]->F;
		}

		double norm_tmp = norm(tmp);

		double direction_coefficient =  (dotmult(spring_huell[i]->force,spring_huell[i]->direction) > 0 )*2-1;//
		//spring_huell[i]->l0 += direction_coefficient/spring_huell[i]->eta_damper*spring[i]->F*timeStep;
		spring_huell[i]->l0 += direction_coefficient/spring_huell[i]->eta_damper*norm_tmp*timeStep;
	}




	for( unsigned int i = 0 ; i < spring_cytoscelett.size() ; i++ ){
		
		spring_cytoscelett[i]->F = norm(spring_cytoscelett[i]->force);
		double tmp[3];
		if( spring_huell[i]->F != 0.){
			tmp[0] = spring_cytoscelett[i]->force[0] * spring_cytoscelett[i]->force[0] * spring_cytoscelett[i]->direction[0] / spring_cytoscelett[i]->F;
			tmp[1] = spring_cytoscelett[i]->force[1] * spring_cytoscelett[i]->force[1] * spring_cytoscelett[i]->direction[1] / spring_cytoscelett[i]->F;
			tmp[2] = spring_cytoscelett[i]->force[2] * spring_cytoscelett[i]->force[2] * spring_cytoscelett[i]->direction[2] / spring_cytoscelett[i]->F;
		}

		double norm_tmp = norm(tmp);

		spring_cytoscelett[i]->F = norm(spring_cytoscelett[i]->force);

		double direction_coefficient =  (dotmult(spring_cytoscelett[i]->force,spring_huell[i]->direction) > 0 )*2-1;//
		//spring_cytoscelett[i]->l0 += direction_coefficient/spring_cytoscelett[i]->eta_damper*spring[i]->F*timeStep;
		spring_cytoscelett[i]->l0 += direction_coefficient/spring_cytoscelett[i]->eta_damper*norm_tmp*timeStep;

	}
	*/

	
	for( unsigned int i = 0 ; i < spring_huell.size() ; i++ ){
//		double force[3];
//		double tmp_d_init = (spring_huell[i]->l0_init + spring_huell[i]->dlength - spring_huell[i]->l0)*spring_huell[i]->k2;
//		tmp_d_init *= tmp_d_init*tmp_d_init;


//		force[0] = spring_huell[i]->force[0]-tmp_d_init;
//		force[1] = spring_huell[i]->force[1]-tmp_d_init;
//		force[2] = spring_huell[i]->force[2]-tmp_d_init;
		
		spring_huell[i]->force_KV[0] *= spring_huell[i]->direction[0];
		spring_huell[i]->force_KV[1] *= spring_huell[i]->direction[1];
		spring_huell[i]->force_KV[2] *= spring_huell[i]->direction[2];

		spring_huell[i]->F = norm(spring_huell[i]->force_KV);
//		spring_huell[i]->F = norm(force);


		
		//if( spring_huell[i]->F*0.001 > spring_huell[i]->l0 ){
		double direction_coefficient =  (dotmult(spring_huell[i]->force_KV,spring_huell[i]->direction) > 0 )*2-1;//

		//spring_huell[i]->l0 += direction_coefficient/spring_huell[i]->eta_damper*spring[i]->F*timeStep;
		spring_huell[i]->l0 += direction_coefficient/spring_huell[i]->eta_damper*spring_huell[i]->F*timeStep;
	//	}
	}




	for( unsigned int i = 0 ; i < spring_cytoscelett.size() ; i++ ){

//		double force[3];
//		double tmp_d_init = (spring_cytoscelett[i]->l0_init + spring_cytoscelett[i]->dlength - spring_cytoscelett[i]->l0)*spring_cytoscelett[i]->k2;
		//tmp_d_init *= tmp_d_init*tmp_d_init;

//		force[0] = spring_cytoscelett[i]->force[0]-tmp_d_init;
//		force[1] = spring_cytoscelett[i]->force[1]-tmp_d_init;
//		force[2] = spring_cytoscelett[i]->force[2]-tmp_d_init;

		spring_cytoscelett[i]->force_KV[0] *= spring_cytoscelett[i]->direction[0];
		spring_cytoscelett[i]->force_KV[1] *= spring_cytoscelett[i]->direction[1];
		spring_cytoscelett[i]->force_KV[2] *= spring_cytoscelett[i]->direction[2];


		

		spring_cytoscelett[i]->F = norm(spring_cytoscelett[i]->force_KV);
//		spring_cytoscelett[i]->F = norm(force);

		//if( spring_cytoscelett[i]->F*0.000001 > spring_cytoscelett[i]->l0 ){
		


		double direction_coefficient =  (dotmult(spring_cytoscelett[i]->force_KV,spring_cytoscelett[i]->direction) > 0 )*2-1;//
		//spring_cytoscelett[i]->l0 += direction_coefficient/spring_cytoscelett[i]->eta_damper*spring[i]->F*timeStep;
		spring_cytoscelett[i]->l0 += direction_coefficient/spring_cytoscelett[i]->eta_damper*spring_cytoscelett[i]->F*timeStep;
	//	}

	}
	
	
}
void CellTriangulated::projection_to_spring(){
	for( unsigned int i = 0 ; i < spring.size() ; i++ ){
		spring[i]->force[0] *= spring[i]->direction[0];
		spring[i]->force[1] *= spring[i]->direction[1];
		spring[i]->force[2] *= spring[i]->direction[2];
	}
}
void CellTriangulated::force_add_to_points(){

  for( unsigned int i = 0 ; i < mass_huell.size() ; i++){
    mass[i]->F_debug_spring[0] = 0.;
    mass[i]->F_debug_spring[1] = 0.;
    mass[i]->F_debug_spring[2] = 0.;
  }

	for( unsigned int i = 0 ; i < spring.size() ; i++){
		
		spring[i]->end->F[0] -= spring[i]->force[0];
		spring[i]->end->F[1] -= spring[i]->force[1];
		spring[i]->end->F[2] -= spring[i]->force[2];

    spring[i]->end->F_debug_spring[0] -= spring[i]->force[0];
    spring[i]->end->F_debug_spring[1] -= spring[i]->force[1];
    spring[i]->end->F_debug_spring[2] -= spring[i]->force[2];

		spring[i]->start->F[0] += spring[i]->force[0];
		spring[i]->start->F[1] += spring[i]->force[1];
		spring[i]->start->F[2] += spring[i]->force[2];

    spring[i]->start->F_debug_spring[0] += spring[i]->force[0];
    spring[i]->start->F_debug_spring[1] += spring[i]->force[1];
    spring[i]->start->F_debug_spring[2] += spring[i]->force[2];
	}
}
void CellTriangulated::calcForce_force3d_second_damper(){
	for( unsigned int i = 0 ; i < spring_huell.size() ; i++ ){
		
//		double tmp_d_init = -(spring_huell[i]->dlength + spring_huell[i]->l0-spring_huell[i]->l0_init);
		double tmp_d_init = -(spring_huell[i]->l0_init - spring_huell[i]->dlength - spring_huell[i]->l0);
	//	tmp_d_init *= tmp_d_init*tmp_d_init;
		spring_huell[i]->force[0] += spring_huell[i]->k2 * tmp_d_init;
		spring_huell[i]->force[1] += spring_huell[i]->k2 * tmp_d_init;
		spring_huell[i]->force[2] += spring_huell[i]->k2 * tmp_d_init;

    spring_huell[i]->F_debug_secondDamper = spring_huell[i]->k2 * tmp_d_init;

		
	}
	for( unsigned int i = 0 ; i < spring_cytoscelett.size() ; i++ ){
		double tmp_d_init = -(spring_cytoscelett[i]->l0_init - spring_cytoscelett[i]->dlength - spring_cytoscelett[i]->l0);
		//	double tmp_d_init = spring_cytoscelett[i]->dlength;
//		tmp_d_init *= tmp_d_init*tmp_d_init;
		
		spring_cytoscelett[i]->force[0] += spring_cytoscelett[i]->k2 * tmp_d_init ;
		spring_cytoscelett[i]->force[1] += spring_cytoscelett[i]->k2 * tmp_d_init ;
		spring_cytoscelett[i]->force[2] += spring_cytoscelett[i]->k2 * tmp_d_init ;

    spring_cytoscelett[i]->F_debug_secondDamper = spring_cytoscelett[i]->k2 * tmp_d_init;
	
	}
}
void CellTriangulated::calcForce_force3d(double timeStep){


	
	

}
void CellTriangulated::calcForce_medium(){

  //force from medium
  for( unsigned int i = 0 ; i < mass.size(); i++){
    mass[i]->F[0] -= mass[i]->v[0] * nu_medium ;
    mass[i]->F[1] -= mass[i]->v[1] * nu_medium ;
    mass[i]->F[2] -= mass[i]->v[2] * nu_medium ;


    mass[i]->F_debug_medium[0] = mass[i]->v[0] * nu_medium;
    mass[i]->F_debug_medium[1] = mass[i]->v[1] * nu_medium;
    mass[i]->F_debug_medium[2] = mass[i]->v[2] * nu_medium;

  }

}
void CellTriangulated::calcForce(double timeStep){

	//must be first, set F to zero
//	if( this->p != 0.)
//		this->addInnerStress();
  this->calcForce_medium(); 
	this->calcForce_froce3d_spring(timeStep);
	//force from liquid in Kelvin-Voigt
	this->calcForce_froce3d_viscous(timeStep);





//	this->calcForce_force3d_second_damper();
	this->projection_to_spring();
	this->force_add_to_points();

  for( unsigned int i = 0 ; i < this->mass.size() ; i++){
    this->mass[i]->F_springs_k[0] = this->mass[i]->F[0];
    this->mass[i]->F_springs_k[1] = this->mass[i]->F[1];
    this->mass[i]->F_springs_k[2] = this->mass[i]->F[2];
  }

   

//	this->updateDamper(timeStep);
//	this->updateDamper_new(timeStep);

}
void CellTriangulated::calcLength(){

  for( unsigned int i = 0 ; i < spring.size() ; i++){

    // calc direction of a spring
    spring[i]->direction[0] = spring[i]->end->position[0] - spring[i]->start->position[0];
    spring[i]->direction[1] = spring[i]->end->position[1] - spring[i]->start->position[1];
    spring[i]->direction[2] = spring[i]->end->position[2] - spring[i]->start->position[2];

    // calc difference of length
    spring[i]->length = normvector(spring[i]->direction);
    spring[i]->dlength = spring[i]->length - spring[i]->l0;// - spring[i]->l_damper;
  }// for all springs
}
void CellTriangulated::calcLengthWithNormDirection(){

  for( unsigned int i = 0 ; i < spring.size() ; i++){

    // calc direction of a spring
    spring[i]->direction[0] = spring[i]->end->position[0] - spring[i]->start->position[0];
    spring[i]->direction[1] = spring[i]->end->position[1] - spring[i]->start->position[1];
    spring[i]->direction[2] = spring[i]->end->position[2] - spring[i]->start->position[2];

    //length of connection
    double l = normvector(spring[i]->direction);

    //normalize direction
    spring[i]->direction[0] /= l;
    spring[i]->direction[1] /= l;
    spring[i]->direction[2] /= l;

    // calc difference of length
    spring[i]->dlength = l - spring[i]->l0;// - spring[i]->l_damper;

  }// for all springs
}
double CellTriangulated::calcVolume(){
  volume = 0.;
  for( unsigned int i = 0 ; i < triangle.size() ;i++ ){
    triangle[i]->calcVolume(mass[0]->position);
    volume += triangle[i]->Volume;
   }
  return volume;
}
void CellTriangulated::clearForce(){
  for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
    this->mass[i]->F[0] = 0.;
    this->mass[i]->F[1] = 0.;
    this->mass[i]->F[2] = 0.;
  }
}

void CellTriangulated::setAndSolveSystem(){

	//Debug Johanes
	double left = 0 ; int left_index = 0;
  for( unsigned int i = 0 ; i < this->mass.size() ; i++){
	  if( this->mass[i]->position[0] < left ){
		  left = this->mass[i]->position[0];
		  left_index = i;
	  }
  }
  //end Debug


  for( int j = 0 ; j < 3 ; j++){//for each dimension

    double *b=(double*) malloc( sizeof(double) * this->mass.size());

    for( unsigned int i = 0 ; i < this->mass.size() ; i++){

      b[i] = 0.;

      for( unsigned int k = 0 ; k < this->mass[i]->neighbourSprings.size(); k++){
        if( this->mass[i]->neighbourSprings[k]->start == this->mass[i] )
          b[i] += this->mass[i]->neighbourSprings[k]->k * this->mass[i]->neighbourSprings[k]->direction[j] *this->mass[i]->neighbourSprings[k]->dlength;
        else
          b[i] -= this->mass[i]->neighbourSprings[k]->k * this->mass[i]->neighbourSprings[k]->direction[j] *this->mass[i]->neighbourSprings[k]->dlength;
      }


      b[i] += this->mass[i]->F[j];

    }//over all mass points

    //multiply inverse of connection matrix with b ;Ax=b--> x = A^-1*b
    for( unsigned int i = 0 ; i < this->mass.size() ; i++){

      this->mass[i]->v[j] = 0.;
      for( unsigned int k = 0 ; k < this->mass.size(); k++)
        this->mass[i]->v[j] += this->A[i][k] * b[k];


    }

    free(b);


  }//for all dimensions

}

void CellTriangulated::removeLastStep(){

//  this->status = this->mStoreOldStatus;
//  this->can_divide = this->store_can_divide;


  for( unsigned int i = 0 ; i < mass.size() ; i++){

    this->mass[i]->F_max = this->mass[i]->store_old_F_max;

    this->mass[i]->dist_move[0] = this->mass[i]->store_old_dist_move[0];
    this->mass[i]->dist_move[1] = this->mass[i]->store_old_dist_move[1];
    this->mass[i]->dist_move[2] = this->mass[i]->store_old_dist_move[2];

    mass[i]->position[0] = mass[i]->store_position[0];
    mass[i]->position[1] = mass[i]->store_position[1];
    mass[i]->position[2] = mass[i]->store_position[2];

    mass[i]->a[0] = mass[i]->store_a[0];
    mass[i]->a[1] = mass[i]->store_a[1];
    mass[i]->a[2] = mass[i]->store_a[2];

    mass[i]->a_old[0] = mass[i]->store_a_old[0];
		mass[i]->a_old[1] = mass[i]->store_a_old[1];
		mass[i]->a_old[2] = mass[i]->store_a_old[2];

		mass[i]->a_new[0] = mass[i]->store_a_new[0];
		mass[i]->a_new[1] = mass[i]->store_a_new[1];
		mass[i]->a_new[2] = mass[i]->store_a_new[2];
    
    mass[i]->v[0] = mass[i]->store_v[0];
		mass[i]->v[1] = mass[i]->store_v[1];
		mass[i]->v[2] = mass[i]->store_v[2];

    /*
    mass[i]->position[0] = mass[i]->store_old_position[0];
    mass[i]->position[1] = mass[i]->store_old_position[1];
    mass[i]->position[2] = mass[i]->store_old_position[2];

    mass[i]->store_position[0] = mass[i]->store_old_position[0];
    mass[i]->store_position[1] = mass[i]->store_old_position[1];
    mass[i]->store_position[2] = mass[i]->store_old_position[2];
    

    mass[i]->a[0] = mass[i]->store_old_a[0];
    mass[i]->a[1] = mass[i]->store_old_a[1];
    mass[i]->a[2] = mass[i]->store_old_a[2];

    mass[i]->store_a[0] = mass[i]->store_old_a[0];
    mass[i]->store_a[1] = mass[i]->store_old_a[1];
    mass[i]->store_a[2] = mass[i]->store_old_a[2];
    


		mass[i]->a_old[0] = mass[i]->store_old_a_old[0];
		mass[i]->a_old[1] = mass[i]->store_old_a_old[1];
		mass[i]->a_old[2] = mass[i]->store_old_a_old[2];

		mass[i]->store_a_old[0] = mass[i]->store_old_a_old[0];
		mass[i]->store_a_old[1] = mass[i]->store_old_a_old[1];
		mass[i]->store_a_old[2] = mass[i]->store_old_a_old[2];
    



		mass[i]->a_new[0] = mass[i]->store_old_a_new[0];
		mass[i]->a_new[1] = mass[i]->store_old_a_new[1];
		mass[i]->a_new[2] = mass[i]->store_old_a_new[2];
		
    
    mass[i]->store_v[0] = mass[i]->store_old_v[0];
		mass[i]->store_v[1] = mass[i]->store_old_v[1];
		mass[i]->store_v[2] = mass[i]->store_old_v[2];
  
		mass[i]->v[0] = mass[i]->store_old_v[0];
		mass[i]->v[1] = mass[i]->store_old_v[1];
		mass[i]->v[2] = mass[i]->store_old_v[2];
      */
    mass[i]->F_max = this->mass[i]->store_F_repuslive;

//    mass[i]->store_F_max = this->mass[i]->store_old_F_repulsive;

    mass[i]->F[0] = this->mass[i]->store_F[0];
    mass[i]->F[1] = this->mass[i]->store_F[1];
    mass[i]->F[2] = this->mass[i]->store_F[2];
    /*
    mass[i]->store_F[0] = this->mass[i]->store_old_F[0];
    mass[i]->store_F[1] = this->mass[i]->store_old_F[1];
    mass[i]->store_F[2] = this->mass[i]->store_old_F[2];
    */
   /*
    this->mass[i]->neighbourSprings.clear();
    for( unsigned int j = 0 ; j < this->mass[i]->mStoreNeighbourSprings.size() ; j++ ){
      this->mass[i]->neighbourSprings.push_back( this->mass[i]->mStoreNeighbourSprings[j] );
    }
    */
  }//for all mass points
  

  /*
  this->spring_divide.clear();
  for( unsigned int j = 0 ; j < this->mStoreOldSpring_divide.size() ; j++ ){
    mStoreOldSpring_divide[j]->k = mStoreOldSpring_divide[j]->store_k;
    this->spring_divide.push_back(this->mStoreOldSpring_divide[j]);
  }
  
  this->spring_divide_cytoskeleton.clear();
  for( unsigned int j = 0 ; j < this->mStoreOldSpring_divide_cytoskeleton.size() ; j++ ){
    mStoreOldSpring_divide_cytoskeleton[j]->k = mStoreOldSpring_divide_cytoskeleton[j]->store_k;
    this->spring_divide_cytoskeleton.push_back(this->mStoreOldSpring_divide_cytoskeleton[j]);
  }
  
  this->spring.clear();
  for( unsigned int j = 0 ; j < this->mStoreOldSpring.size() ; j++ ){
    mStoreOldSpring[j]->l0 = mStoreOldSpring[j]->store_l0;
    this->spring.push_back(this->mStoreOldSpring[j]);
  }
  
  */

  for( unsigned int i = 0 ; i < this->triangle.size() ; i++){
    this->triangle[i]->mNormalvector[0] = this->triangle[i]->store_normalvetor[0];
    this->triangle[i]->mNormalvector[1] = this->triangle[i]->store_normalvetor[1];
    this->triangle[i]->mNormalvector[2] = this->triangle[i]->store_normalvetor[2];
  }
  
    for( unsigned int i = 0 ; i < this->spring.size() ; i++){

    this->spring[i]->l0 = this->spring[i]->store_l0;
  }


}
void CellTriangulated::storeLastStep(){
  
//  this->mStoreOldStatus = this->status;
//  this->store_can_divide = this->can_divide;

  for( unsigned int i = 0 ; i < mass.size() ; i++){

 //   this->mass[i]->store_old_F_max = this->mass[i]->store_F_max;
    this->mass[i]->store_F_max = this->mass[i]->F_max;
    /*
    this->mass[i]->store_old_dist_move[0] = this->mass[i]->store_dist_move[0];
    this->mass[i]->store_old_dist_move[1] = this->mass[i]->store_dist_move[1];
    this->mass[i]->store_old_dist_move[2] = this->mass[i]->store_dist_move[2];
    */
    this->mass[i]->store_dist_move[0] = this->mass[i]->dist_move[0];
    this->mass[i]->store_dist_move[1] = this->mass[i]->dist_move[1];
    this->mass[i]->store_dist_move[2] = this->mass[i]->dist_move[2];
    /*
    mass[i]->store_old_position[0] = mass[i]->store_position[0];
    mass[i]->store_old_position[1] = mass[i]->store_position[1];
    mass[i]->store_old_position[2] = mass[i]->store_position[2];
    */
    mass[i]->store_position[0] = mass[i]->position[0];
    mass[i]->store_position[1] = mass[i]->position[1];
    mass[i]->store_position[2] = mass[i]->position[2];


    /*
    mass[i]->store_old_a[0] = mass[i]->store_a[0];
    mass[i]->store_old_a[1] = mass[i]->store_a[1];
    mass[i]->store_old_a[2] = mass[i]->store_a[2];
    */
    mass[i]->store_a[0] = mass[i]->a[0];
    mass[i]->store_a[1] = mass[i]->a[1];
    mass[i]->store_a[2] = mass[i]->a[2];


    /*
    mass[i]->store_old_a_old[0] = mass[i]->a_old[0];
    mass[i]->store_old_a_old[1] = mass[i]->a_old[1];
    mass[i]->store_old_a_old[2] = mass[i]->a_old[2];
    */
    mass[i]->store_a_old[0] = mass[i]->a_old[0];
    mass[i]->store_a_old[1] = mass[i]->a_old[1];
    mass[i]->store_a_old[2] = mass[i]->a_old[2];

    /*

    mass[i]->store_old_a_new[0] = mass[i]->store_a_new[0];
    mass[i]->store_old_a_new[1] = mass[i]->store_a_new[1];
    mass[i]->store_old_a_new[2] = mass[i]->store_a_new[2];
    */
    mass[i]->store_a_new[0] = mass[i]->a_new[0];
    mass[i]->store_a_new[1] = mass[i]->a_new[1];
    mass[i]->store_a_new[2] = mass[i]->a_new[2];


/*
    mass[i]->store_old_v[0] = mass[i]->store_v[0];
    mass[i]->store_old_v[1] = mass[i]->store_v[1];
    mass[i]->store_old_v[2] = mass[i]->store_v[2];
*/
    mass[i]->store_v[0] = mass[i]->v[0];
    mass[i]->store_v[1] = mass[i]->v[1];
    mass[i]->store_v[2] = mass[i]->v[2];

//    mass[i]->store_old_F_repulsive = this->mass[i]->store_F_repuslive;

    mass[i]->store_F_repuslive = this->mass[i]->F_max;
/*
    mass[i]->store_old_F[0] = this->mass[i]->store_F[0];
    mass[i]->store_old_F[1] = this->mass[i]->store_F[1];
    mass[i]->store_old_F[2] = this->mass[i]->store_F[2];
  */
    mass[i]->store_F[0] = this->mass[i]->F[0];
    mass[i]->store_F[1] = this->mass[i]->F[1];
    mass[i]->store_F[2] = this->mass[i]->F[2];



  /*  this->mass[i]->mStoreNeighbourSprings.clear();
    for( unsigned int j = 0 ; j < this->mass[i]->neighbourSprings.size() ; j++ ){
      this->mass[i]->mStoreNeighbourSprings.push_back( this->mass[i]->neighbourSprings[j] );
    }
    */
  }//for all mass points


  for( unsigned int i = 0 ; i < this->triangle.size() ; i++){
/*
    this->triangle[i]->store_old_normalvector[0] = this->triangle[i]->store_normalvetor[0];
    this->triangle[i]->store_old_normalvector[1] = this->triangle[i]->store_normalvetor[1];
    this->triangle[i]->store_old_normalvector[2] = this->triangle[i]->store_normalvetor[2];
    */
    this->triangle[i]->store_normalvetor[0] = this->triangle[i]->mNormalvector[0];
    this->triangle[i]->store_normalvetor[1] = this->triangle[i]->mNormalvector[1];
    this->triangle[i]->store_normalvetor[2] = this->triangle[i]->mNormalvector[2];
  }

  /*

  this->mStoreOldSpring_divide.clear();
  for( unsigned int j = 0 ; j < this->spring_divide.size() ; j++ ){
    spring_divide[j]->store_k = spring_divide[j]->k;
    this->mStoreOldSpring_divide.push_back(this->spring_divide[j]);
  }

  this->mStoreOldSpring_divide_cytoskeleton.clear();
  for( unsigned int j = 0 ; j < this->spring_divide_cytoskeleton.size() ; j++ ){
  spring_divide_cytoskeleton[j]->store_k = spring_divide_cytoskeleton[j]->k;
    this->mStoreOldSpring_divide_cytoskeleton.push_back(this->spring_divide_cytoskeleton[j]);
  }

  this->mStoreOldSpring.clear();
  for( unsigned int j = 0 ; j < this->spring.size() ; j++ ){
    spring[j]->store_l0 = spring[j]->l0;
    this->mStoreOldSpring.push_back(this->spring[j]);
  }
  
  */
  
  for( unsigned int i = 0 ; i < this->spring.size() ; i++){

  //  this->spring[i]->store_old_l0 = this->spring[i]->store_l0;

    this->spring[i]->store_l0 = this->spring[i]->l0;
  }

}


void CellTriangulated::updateDamper_new( double timeStep){
	
	for( unsigned int i = 0 ; i < spring_huell.size() ; i++ ){
		double force[3];
		force[0] = spring_huell[i]->end->F[0] - spring_huell[i]->start->F[0];
		force[1] = spring_huell[i]->end->F[1] - spring_huell[i]->start->F[1];
		force[2] = spring_huell[i]->end->F[2] - spring_huell[i]->start->F[2];

		double force_damper = dotmult(spring_huell[i]->force,spring_huell[i]->direction);

		double tmp_length =  spring_huell[i]->dlength - spring_huell[i]->l0 - spring_huell[i]->l0_init;

		force_damper -= tmp_length*spring_huell[i]->k2;

		spring_huell[i]->l0 += force_damper/spring_huell[i]->eta_damper*timeStep;
	}




	for( unsigned int i = 0 ; i < spring_cytoscelett.size() ; i++ ){
		double force[3];
		force[0] = spring_cytoscelett[i]->end->F[0] - spring_cytoscelett[i]->start->F[0];
		force[1] = spring_cytoscelett[i]->end->F[1] - spring_cytoscelett[i]->start->F[1];
		force[2] = spring_cytoscelett[i]->end->F[2] - spring_cytoscelett[i]->start->F[2];

		double force_damper = dotmult(spring_cytoscelett[i]->force,spring_cytoscelett[i]->direction);

		double tmp_length =  spring_cytoscelett[i]->dlength - spring_cytoscelett[i]->l0 - spring_cytoscelett[i]->l0_init;

		force_damper -= tmp_length*spring_cytoscelett[i]->k2;

		spring_cytoscelett[i]->l0 += force_damper/spring_cytoscelett[i]->eta_damper*timeStep;
	
	}
}

//move
double CellTriangulated::updateVelocity(int number, double delta_t){

if( type != 1){	
	if( number == 1){
		for(unsigned int i = 0 ; i < mass.size() ; i++){

			mass[i]->v[0] += mass[i]->a[0] *delta_t;
			mass[i]->v[1] += mass[i]->a[1] *delta_t;
			mass[i]->v[2] += mass[i]->a[2] *delta_t;

                 this->mass[i]->dist_v[0] =  mass[i]->a[0] *delta_t;
      this->mass[i]->dist_v[1] =  mass[i]->a[1] *delta_t;
      this->mass[i]->dist_v[2] =  mass[i]->a[2] *delta_t;



		}
	}
	else if( number ==2){
		for(unsigned int i = 0 ; i < mass.size() ; i++){
			mass[i]->v[0] = mass[i]->v[0] + 0.5*( mass[i]->a[0]+mass[i]->a_old[0]) *delta_t;
			mass[i]->v[1] = mass[i]->v[1] + 0.5*( mass[i]->a[1]+mass[i]->a_old[1]) *delta_t;
			mass[i]->v[2] = mass[i]->v[2] + 0.5*( mass[i]->a[2]+mass[i]->a_old[2]) *delta_t;
/*
           this->mass[i]->dist_v[0] =  0.5*( mass[i]->a[0]+mass[i]->a_old[0]) *delta_t;
      this->mass[i]->dist_v[1] =  0.5*( mass[i]->a[1]+mass[i]->a_old[1]) *delta_t;
      this->mass[i]->dist_v[2] =  0.5*( mass[i]->a[2]+mass[i]->a_old[2]) *delta_t;
      */

		//	double tmp_dist = norm(mass[i]->v);
		//	if( tmp_dist > 1. && abs(dist(mass[i]->v,mass[i]->v_stored) / tmp_dist)> 1.0 )
		//		this->able_to_move = 0;

		}
	}else if( number == 4){
		for(unsigned int i = 0 ; i < mass.size() ; i++){
			//beeman predictor
      mass[i]->store_v[0] = mass[i]->v[0];
			mass[i]->store_v[1] = mass[i]->v[1];
			mass[i]->store_v[2] = mass[i]->v[2];
			
			mass[i]->v[0] += 3./2.*mass[i]->a[0] *delta_t - 0.5*mass[i]->a_old[0]*delta_t;
			mass[i]->v[1] += 3./2.*mass[i]->a[1] *delta_t - 0.5*mass[i]->a_old[1]*delta_t;
			mass[i]->v[2] += 3./2.*mass[i]->a[2] *delta_t - 0.5*mass[i]->a_old[2]*delta_t;

		}
	}else{
		for(unsigned int i = 0 ; i < mass.size() ; i++){
			//beeman corrector
			mass[i]->v[0] = mass[i]->v[0] + 1./3.*mass[i]->a_new[0] *delta_t + 5./6.*mass[i]->a[0]*delta_t - 1./6.*mass[i]->a_old[0]*delta_t;
			mass[i]->v[1] = mass[i]->v[0] + 1./3.*mass[i]->a_new[1] *delta_t + 5./6.*mass[i]->a[1]*delta_t - 1./6.*mass[i]->a_old[1]*delta_t;
			mass[i]->v[2] = mass[i]->v[0] + 1./3.*mass[i]->a_new[2] *delta_t + 5./6.*mass[i]->a[2]*delta_t - 1./6.*mass[i]->a_old[2]*delta_t;

		}

	}
}	
	return 0.;
}
double CellTriangulated::updatePosition(int number, double delta_t){

// if( this->able_to_move == 0)
//    return 0;

  double max = 0.;
  if( this->type != 1){

  if( number ==1){
    mBoundingBox.xmax = mass[0]->position[0]+ mass[0]->v[0] *delta_t;
    mBoundingBox.ymax = mass[0]->position[1]+ mass[0]->v[1] *delta_t;
    mBoundingBox.zmax = mass[0]->position[2]+ mass[0]->v[2] *delta_t;
    mBoundingBox.xmin = mBoundingBox.xmax;
    mBoundingBox.ymin = mBoundingBox.ymax;
    mBoundingBox.zmin = mBoundingBox.zmax;


    for(unsigned int i = 0 ; i < mass.size() ; i++){

      double dp[3];

      dp[0] = (mass[i]->v[0] *delta_t);
      dp[1] = (mass[i]->v[1] *delta_t);
      dp[2] = (mass[i]->v[2] *delta_t);

      mass[i]->position[0] += dp[0];
      mass[i]->position[1] += dp[1];
      mass[i]->position[2] += dp[2];

      this->mass[i]->dist_move[0] = dp[0];
      this->mass[i]->dist_move[1] = dp[1];
      this->mass[i]->dist_move[2] = dp[2];



      double norm_dp = norm(dp);
      if( norm_dp > max)
        max = norm_dp;

      if( mBoundingBox.xmin >  mass[i]->position[0] )
        mBoundingBox.xmin = mass[i]->position[0];
      if ( mBoundingBox.xmax <  mass[i]->position[0] )
        mBoundingBox.xmax = mass[i]->position[0];
      if ( mBoundingBox.ymin >  mass[i]->position[1] )
        mBoundingBox.ymin = mass[i]->position[1];
      if ( mBoundingBox.ymax <  mass[i]->position[1] )
        mBoundingBox.ymax = mass[i]->position[1];
      if ( mBoundingBox.zmin >  mass[i]->position[2] )
        mBoundingBox.zmin = mass[i]->position[2];
      if ( mBoundingBox.zmax <  mass[i]->position[2] )
        mBoundingBox.zmax = mass[i]->position[2];

 //     mass[i]->mAreaVoronoi = 0.;

    }
    if( dim == 2)
      mass[0]->position[2] = 0.;
    if( dim == 1){
      mass[0]->position[1] = 0.;
      mass[0]->position[2] = 0.;
    }
  }else if( number == 2){
/*
    mBoundingBox.xmax = (mass[0]->v[0] *delta_t + 0.5*mass[0]->a_old[0]*delta_t*delta_t);
    mBoundingBox.ymax = (mass[0]->v[1] *delta_t + 0.5*mass[0]->a_old[1]*delta_t*delta_t);
    mBoundingBox.zmax = (mass[0]->v[2] *delta_t + 0.5*mass[0]->a_old[2]*delta_t*delta_t);
*/
    mBoundingBox.xmax = (mass[0]->v[0] *delta_t + 0.5*mass[0]->a[0]*delta_t*delta_t);
    mBoundingBox.ymax = (mass[0]->v[1] *delta_t + 0.5*mass[0]->a[1]*delta_t*delta_t);
    mBoundingBox.zmax = (mass[0]->v[2] *delta_t + 0.5*mass[0]->a[2]*delta_t*delta_t);

    mBoundingBox.xmin = mBoundingBox.xmax;
    mBoundingBox.ymin = mBoundingBox.ymax;
    mBoundingBox.zmin = mBoundingBox.zmax;

    for(unsigned int i = 0 ; i < mass.size() ; i++){
/*
      mass[i]->position[0] += (mass[i]->v[0] *delta_t + 0.5*mass[i]->a_old[0]*delta_t*delta_t);
      mass[i]->position[1] += (mass[i]->v[1] *delta_t + 0.5*mass[i]->a_old[1]*delta_t*delta_t);
      mass[i]->position[2] += (mass[i]->v[2] *delta_t + 0.5*mass[i]->a_old[2]*delta_t*delta_t);
*/
      double dp[3];
      dp[0] = (mass[i]->v[0] *delta_t + 0.5*mass[i]->a[0]*delta_t*delta_t);
      dp[1] = (mass[i]->v[1] *delta_t + 0.5*mass[i]->a[1]*delta_t*delta_t);
      dp[2] = (mass[i]->v[2] *delta_t + 0.5*mass[i]->a[2]*delta_t*delta_t);

      mass[i]->position[0] += dp[0];
      mass[i]->position[1] += dp[1];
      mass[i]->position[2] += dp[2];

      this->mass[i]->dist_move[0] = dp[0];
      this->mass[i]->dist_move[1] = dp[1];
      this->mass[i]->dist_move[2] = dp[2];

      double norm_dp = norm(dp);

      this->mass[i]->log_move_max = norm_dp;

      if( norm_dp > max)
        max = norm_dp;

      if ( mBoundingBox.xmin >  mass[i]->position[0] )
        mBoundingBox.xmin = mass[i]->position[0];
      if ( mBoundingBox.xmax <  mass[i]->position[0] )
        mBoundingBox.xmax = mass[i]->position[0];
      if ( mBoundingBox.ymin >  mass[i]->position[1] )
        mBoundingBox.ymin = mass[i]->position[1];
      if ( mBoundingBox.ymax <  mass[i]->position[1] )
        mBoundingBox.ymax = mass[i]->position[1];
      if ( mBoundingBox.zmin >  mass[i]->position[2] )
        mBoundingBox.zmin = mass[i]->position[2];
      if ( mBoundingBox.zmax <  mass[i]->position[2] )
        mBoundingBox.zmax = mass[i]->position[2];

      if( dim == 2)
        mass[0]->position[2] = 0.;
      if( dim == 1){
        mass[0]->position[1] = 0.;
        mass[0]->position[2] = 0.;
      }
    }
  }else if( number == 11 ){


  }

  else{

		mBoundingBox.xmax = mass[0]->position[0]+ mass[0]->v[0] *delta_t +(2./3.*mass[0]->a[0]-1./6.*mass[0]->a_old[0])*delta_t*delta_t;
		mBoundingBox.ymax = mass[0]->position[1]+ mass[0]->v[1] *delta_t +(2./3.*mass[0]->a[1]-1./6.*mass[0]->a_old[1])*delta_t*delta_t;
		mBoundingBox.zmax = mass[0]->position[2]+ mass[0]->v[2] *delta_t +(2./3.*mass[0]->a[2]-1./6.*mass[0]->a_old[2])*delta_t*delta_t;
		mBoundingBox.xmin = mBoundingBox.xmax;
		mBoundingBox.ymin = mBoundingBox.ymax;
		mBoundingBox.zmin = mBoundingBox.zmax;


		for(unsigned int i = 0 ; i < mass.size() ; i++){

			mass[i]->position[0] += (mass[i]->v[0] *delta_t) +(2./3.*mass[i]->a[0]-1./6.*mass[i]->a_old[0])*delta_t*delta_t;
			mass[i]->position[1] += (mass[i]->v[1] *delta_t) +(2./3.*mass[i]->a[1]-1./6.*mass[i]->a_old[1])*delta_t*delta_t;
			mass[i]->position[2] += (mass[i]->v[2] *delta_t) +(2./3.*mass[i]->a[2]-1./6.*mass[i]->a_old[2])*delta_t*delta_t;

			if ( mBoundingBox.xmin >  mass[i]->position[0] )
				mBoundingBox.xmin = mass[i]->position[0];
			if ( mBoundingBox.xmax <  mass[i]->position[0] )
				mBoundingBox.xmax = mass[i]->position[0];
			if ( mBoundingBox.ymin >  mass[i]->position[1] )
				mBoundingBox.ymin = mass[i]->position[1];
			if ( mBoundingBox.ymax <  mass[i]->position[1] )
				mBoundingBox.ymax = mass[i]->position[1];
			if ( mBoundingBox.zmin >  mass[i]->position[2] )
				mBoundingBox.zmin = mass[i]->position[2];
			if ( mBoundingBox.zmax <  mass[i]->position[2] )
				mBoundingBox.zmax = mass[i]->position[2];



	//		mass[i]->mAreaVoronoi = 0.;

		}
		if( dim == 2)
			mass[0]->position[2] = 0.;
		if( dim == 1){
			mass[0]->position[1] = 0.;
			mass[0]->position[2] = 0.;
		}



	}

    for( unsigned int i = 0 ; i < triangle.size() ; i++)
      triangle[i]->setNormalVector();
  }

	double tmp = 1./8.;
	mBoundingBox.xmin -= tmp;
	mBoundingBox.xmax += tmp;
	mBoundingBox.ymin -= tmp;
	mBoundingBox.ymax += tmp;
	mBoundingBox.zmin -= tmp;
	mBoundingBox.zmax += tmp;


  return max;
}

//Cell division
void   CellTriangulated::addDivideCircle(double v_x, double v_y, double v_z, int number , double eta){

	//direction of divide circle
	double normalvector[3] = {v_x,v_y,v_z};
	if( dim<3)
		normalvector[2] = 0.;
	if( dim==1)
		normalvector[1] = 0.;
	normvector(normalvector);


	map<double,int> mymap_1;
	map<double,int>::iterator it;


	//select points of divideCircle
	for( unsigned int i = 0 ; i < mass_huell.size() ; i++){
		
		//normalvector of point
		double point[3] = {	mass_huell[i]->position[0] - mass[0]->position[0],
							mass_huell[i]->position[1] - mass[0]->position[1],
							mass_huell[i]->position[2] - mass[0]->position[2]};
		normvector(point);

		//degree between divide circle and point
		double mass_dist2 = acos(dotmult(normalvector, point))*180/PI-90;

		mymap_1.insert ( pair<double,int>(abs(mass_dist2),i ));
		
		//deside witch side
		if( mass_dist2 < 0 )
			mass_huell[i]->divide_number = 1;
		else
			mass_huell[i]->divide_number = 2;
	

	}//end select points



	int *mass_reverenz = new int[number];

	int count =0;
	for ( it = mymap_1.begin() ; it != mymap_1.end() && count < number; it++ ){

		mass_reverenz[count] = (*it).second;
		
		//declare divide circle	
		mass_huell[mass_reverenz[count]]->divide_number = 0;

		count++;

	}




	//extract order of points of divideCircle
	
	//reference point of divideCircle( first point)
	double  direction[3] = {	mass_huell[mass_reverenz[0]]->position[0] - mass[0]->position[0],
								mass_huell[mass_reverenz[0]]->position[1] - mass[0]->position[1],
								mass_huell[mass_reverenz[0]]->position[2] - mass[0]->position[2]};
	normvector(direction);

	map<double,int> mymap_2;

	mymap_2.insert ( pair<double,int>(0,mass_reverenz[0] ));


	

	//set a order of points of divideCircle
	for( int i = 1 ; i < number ; i++){

		//direction of point (norm)
		double point[3] = {	mass_huell[mass_reverenz[i]]->position[0] - mass[0]->position[0],
							mass_huell[mass_reverenz[i]]->position[1] - mass[0]->position[1],
							mass_huell[mass_reverenz[i]]->position[2] - mass[0]->position[2]};
		normvector(point);

		double divideCircle_tmp;

		double v_tmp[3];
		crossProduct(direction,point,v_tmp);
			normvector(v_tmp);

		divideCircle_tmp = acos(dotmult(direction, point))*180/PI;

		if( (acos(dotmult(normalvector, v_tmp))*180/PI) < 90)
			divideCircle_tmp = 360-divideCircle_tmp;

		mymap_2.insert ( pair<double,int>(divideCircle_tmp,mass_reverenz[i] ));

	}


	count = 0;
	for ( it = mymap_2.begin() ; it != mymap_2.end(); it++ ){
		mass_reverenz[count] = (*it).second;
		count++;
	}










	//change color of selected mass points
	for(  int i = 0 ; i < number ; i++){
		mass_huell[mass_reverenz[i]]->r = 102;//1.*255;
		mass_huell[mass_reverenz[i]]->g = 205;//0.65*255;
		mass_huell[mass_reverenz[i]]->b = 0.;
		mass_huell[mass_reverenz[i]]->divide_number = 4;//4 is for circle
	}


	//detect for each spring which side
	for( unsigned int i = 0 ; i < spring_huell.size() ; i++){
		
		if( spring_huell[i]->start->divide_number == spring_huell[i]->end->divide_number){
//			if( spring[i]->start->divide_number != 0){
				spring_huell[i]->divideCircle = spring_huell[i]->start->divide_number;

	//		}
		}
		else{
			if( spring_huell[i]->start->divide_number == 4 ){
				spring_huell[i]->divideCircle = spring_huell[i]->end->divide_number;

			}
			if( spring_huell[i]->end->divide_number == 4 ){
				spring_huell[i]->divideCircle = spring_huell[i]->start->divide_number;
			}
		}
		

				
				if( spring_huell[i]->divideCircle == 0)
					spring_huell[i]->divideCircle = 4;
	}


//triangle
	for( unsigned int i = 0 ; i < triangle.size() ; i++) {

		if( triangle[i]->points[0]->divide_number != 4 ){
			triangle[i]->divideCircle = triangle[i]->points[0]->divide_number;
		}else{
			if( triangle[i]->points[1]->divide_number != 4){
				triangle[i]->divideCircle = triangle[i]->points[1]->divide_number;
			}else
				triangle[i]->divideCircle = triangle[i]->points[2]->divide_number;
		}

	}



	
/*

vector<int> index;

	for( unsigned int i = 0 ; i < spring_huell.size() ; i++ ){
		if( spring_huell[i]->divideCircle == 4 ){
		vector<int> tri;
			for( unsigned int j = 0 ; j < triangle.size() ; j++ ){
				if( triangle[j]->springs[0] == spring_huell[i]   || 
					 triangle[j]->springs[1] == spring_huell[i]   ||
					 triangle[j]->springs[2] == spring_huell[i] ) {//&&
					 //triangle[j]->divideCircle != 4 ){
				
							//	change=2;
						

						 tri.push_back(j);					
								spring_huell[i]->divideCircle = triangle[j]->divideCircle;
								//spring_huell[i]->r = 0;
								//triangle[j]->r = 0;
				}	
				


		}

			int a = 9;
			if(((triangle[tri[0]]->points[0]->divide_number == 4 ||
				triangle[tri[0]]->points[1]->divide_number == 4 ||
				triangle[tri[0]]->points[2]->divide_number == 4 )&&(
				triangle[tri[1]]->points[0]->divide_number != 4 ||
				triangle[tri[1]]->points[1]->divide_number != 4 ||
				triangle[tri[1]]->points[2]->divide_number != 4 ))||

				((triangle[tri[0]]->points[0]->divide_number != 4 ||
				triangle[tri[0]]->points[1]->divide_number != 4 ||
				triangle[tri[0]]->points[2]->divide_number != 4 )&&(
				triangle[tri[1]]->points[0]->divide_number == 4 ||
				triangle[tri[1]]->points[1]->divide_number == 4 ||
				triangle[tri[1]]->points[2]->divide_number == 4 ))
				&&
				(triangle[tri[0]]->divideCircle !=4 || triangle[tri[1]]->divideCircle!=4))
				index.push_back(i); 
			tri.clear();
		}
		
	}




for( unsigned int i = 0 ; i < index.size() ; i++ ){
			for( unsigned int j = 0 ; j < triangle.size() ; j++ ){
				if( (triangle[j]->springs[0] == spring_huell[index[i]]   || 
					 triangle[j]->springs[1] == spring_huell[index[i]]   ||
					 triangle[j]->springs[2] == spring_huell[index[i]] )  ){
							
							//	change=2;
						 triangle[j]->divideCircle = spring_huell[index[i]]->divideCircle;
	//							spring_huell[i]->divideCircle = triangle[j]->divideCircle;
								//spring_huell[i]->r = 0;
								//triangle[j]->r = 0;
					

				}

			}

	}
*/


/*
	//triangle
	for( unsigned int i = 0 ; i < triangle.size() ; i++) {
		if( triangle[i]->divideCircle == 4 ){
	
		if( triangle[i]->springs[0]->divideCircle != 4 ){
			triangle[i]->divideCircle = triangle[i]->springs[0]->divideCircle;
		}else{
			if( triangle[i]->springs[1]->divideCircle != 4){
				triangle[i]->divideCircle = triangle[i]->springs[1]->divideCircle;
			}else{
				if( triangle[i]->springs[2]->divideCircle != 4){
					triangle[i]->divideCircle = triangle[i]->springs[2]->divideCircle;
					}
			}
		}


		}

		}
*/	
	//change = 0;
//}








//color change

/*
for( unsigned int i = 0 ; i < spring_huell.size() ; i++ ){
				if( spring_huell[i]->divideCircle == 1 )
					spring_huell[i]->b = 255;
				else{
				if( spring_huell[i]->divideCircle == 2 ){
					spring_huell[i]->b = 100;
					spring_huell[i]->g = 100;
				}
				else{
					if( spring_huell[i]->divideCircle == 4 )
						spring_huell[i]->r = 0;
				}
				}



for( unsigned int i = 0 ; i < triangle.size() ; i++ ){
		if( triangle[i]->divideCircle == 1){
			triangle[i]->g = 255;
		}
		else{
		if( triangle[i]->divideCircle == 2){
			triangle[i]->b = 255;
		}
		else{

triangle[i]->r = 0;

		}

		}
}
*/
		

 //contracting ring
	//add new springs
	for( int i = 0 ; i < (number-1) ; i++){
		Spring *spring = new Spring();
		spring->k = 0;
		spring->nu = eta;
		spring->start = mass_huell[mass_reverenz[i]];
		spring->end = mass_huell[mass_reverenz[i+1]];

		spring->r = 102.;//1.0*255;
		spring->g = 205.;//165;//1.0*255;
		spring->radius = 0.006;
		spring->l0 = 0.;
		spring->divideCircle = 0;
		
		spring->mIndex = this->spring.size();

		spring->eta_damper = eta_damper; 

		this->spring.push_back(spring);
		this->spring_divide.push_back(spring);

		mass_huell[mass_reverenz[i]]->neighbourSprings.push_back(spring);
		mass_huell[mass_reverenz[i+1]]->neighbourSprings.push_back(spring);


		
		Triangle *triangle_tmp = new Triangle();
		triangle_tmp->points[0] = spring->start;
		triangle_tmp->points[1] = spring->end;
		triangle_tmp->points[2] = mass[0];
		triangle_tmp->divideCircle = 3;

		triangle_divide_ring.push_back(triangle_tmp);
		


	}

	//add last spring
	Spring *spring = new Spring();

	spring->k = 0;
	spring->nu = eta;
	spring->start =mass_huell[ mass_reverenz[0]];
	spring->end = mass_huell[mass_reverenz[number-1]];

		spring->r = 102;//1.0*255;
		spring->g = 205;//1.0*255;
	spring->radius = 0.006;
	spring->l0 = 0.;
	spring->divideCircle = 0;

	this->spring.push_back(spring);
	this->spring_divide.push_back(spring);

	mass_huell[mass_reverenz[0]]->neighbourSprings.push_back(spring);
	mass_huell[mass_reverenz[number-1]]->neighbourSprings.push_back(spring);

	
		Triangle *triangle_tmp = new Triangle();
		triangle_tmp->points[0] = spring->start;
		triangle_tmp->points[1] = spring->end;
		triangle_tmp->points[2] = mass[0];
		triangle_tmp->divideCircle = 3;

		triangle_divide_ring.push_back(triangle_tmp);
	

//contracting cytoskeleton
	//add new springs
	for( int i = 0 ; i < number ; i++){
		Spring *spring = new Spring();
		spring->k = 0;
		spring->nu = eta;
		spring->start = mass[0];
		spring->end = mass_huell[mass_reverenz[i]];


		spring->r = 102;//1.0*255;
		spring->g = 205;//1.0*255;


		spring->radius = 0.006;
		spring->l0 = 0.;
		spring->divideCircle = 0;

		this->spring.push_back(spring);
		this->spring_divide_cytoskeleton.push_back(spring);

		mass[0]->neighbourSprings.push_back(spring);
		mass_huell[mass_reverenz[i]]->neighbourSprings.push_back(spring);

	}


	



}
double CellTriangulated::divideCircleLength(){

	double tmp = 0.;
	for( unsigned int i = 0 ; i < spring_divide.size() ; i++)
		tmp += spring_divide[i]->dlength;

	return tmp;

}
double CellTriangulated::calcRadius(){
	
	double tmp = 0.0;

	for( unsigned int i = 0; i < mass_huell.size() ; i++){
		double v[3];
		v[0] = mass_huell[i]->position[0] - position.x;
		v[1] = mass_huell[i]->position[1] - position.y;
		v[2] = mass_huell[i]->position[2] - position.z;
		tmp += norm(v);
	}

	tmp /=  mass_huell.size();

	return tmp;

}

void CellTriangulated::statusDepend( double timeStep, double time){

  double r;
  switch( status){
    case 2000:{
			for( unsigned int i = 1 ; i < mass.size() ; i++ ){			

			}
		break;}
      case 204://growth in cube
      if( this->presure_inside < pressure_threshold ){
        double stretch = 1+0.0256*timeStep;

   
        radius_cell*= stretch;


        v_reference = 4./3.*PI*radius_cell*radius_cell*radius_cell*volume_correction;

        for( unsigned int i = 0 ; i < this->spring.size(); i++){
          this->spring[i]->l0 *= stretch;
          this->spring[i]->l0_init*=stretch;	
        }
      }


      for( unsigned int i = 0 ; i < this->mass_huell.size() ; i++){
        double tmp_dist;
        
        double F[3] = {0.,0.,0.};

        //x
        tmp_dist = (1.0-this->mass_huell[i]->position[0]);
        F[0] -= pow(0.0455/tmp_dist,12)-2*pow(0.0455/tmp_dist,6);
        tmp_dist = abs(-1.0-this->mass_huell[i]->position[0]);
        F[0] += pow(0.0455/tmp_dist,12)-2*pow(0.0455/tmp_dist,6);
        //y
        tmp_dist = (1.0-this->mass_huell[i]->position[1]);
        F[1] -= pow(0.0455/tmp_dist,12)-2*pow(0.0455/tmp_dist,6);
        tmp_dist = abs(-1.0-this->mass_huell[i]->position[1]);
        F[1] += pow(0.0455/tmp_dist,12)-2*pow(0.0455/tmp_dist,6);
        //z
        tmp_dist = (1.0-this->mass_huell[i]->position[2]);
        F[2] -= pow(0.0455/tmp_dist,12)-2*pow(0.0455/tmp_dist,6);
        tmp_dist = abs(-1.0-this->mass_huell[i]->position[2]);
        F[2] += pow(0.0455/tmp_dist,12)-2*pow(0.0455/tmp_dist,6);

        if( norm(F) > 1000){
          this->mass_huell[i]->radius = 0.05;
          this->mass_huell[i]->r = 1;
        }

        this->mass_huell[i]->F[0] += F[0];
        this->mass_huell[i]->F[1] += F[1];
        this->mass_huell[i]->F[2] += F[2];

      }



      break;





		#pragma region pressure inside
		case 202:
			
			this->addInnerStress();


			r = this->calcRadius();
			p = radius_cell*radius_cell*radius_cell/r/r/r*p;

			radius_cell = r;

		break;
		#pragma endregion
		#pragma region push two cells together
	/*	case 201:
			if( mass[0]->position[0] < 0 ){
				mass[0]->F[0]+= 5;

			}
			if( mass[0]->position[0] > 0 ){
				mass[0]->F[0]-= 5;

			}
			break;*/
		#pragma endregion
		#pragma region strech cell with zwo points
		case 200:{
			if( time < 4. ){
			for( unsigned i = 0 ; i < mass_huell.size() ; i++){

				mass_huell[i]->F[0] += mass_huell[i]->F_force_profile[0];
				mass_huell[i]->F[1] += mass_huell[i]->F_force_profile[1];
				mass_huell[i]->F[2] += mass_huell[i]->F_force_profile[2];

			}
			}
			//set Force of cell nucleus to zero
			mass[0]->F[0] = 0.;
			mass[0]->F[1] = 0.;
			mass[0]->F[2] = 0.;
			
			break;}   
		#pragma endregion

		#pragma region old cell cycle
		case 1:
			for( unsigned int j = 0 ; j < spring_divide.size() ; j++)
				spring_divide[j]->k += 0.05;
			break;
		case 10:
//			if( cells->t_tmp < 15 ){
//				cells->x_left +=(cells->delta_t/8);
//				cells->x_right -= (cells->delta_t/8);
//			}
			break;
		case 11:
			addConstantForce(0.25,0,0);
			break;
		case 12:
			for( unsigned int i = 0 ; i < spring.size(); i++){
				double stretch = 1.000001;
				spring[i]->l0*= stretch;
				radius_cell*= stretch;

			}
			break;
		case 13:
			for( unsigned int i = 0 ; i < triangle.size(); i++){

				triangle[i]->setArea(p);

				triangle[i]->points[0]->F[0] += triangle[i]->mNormalvector[0];
				triangle[i]->points[0]->F[1] += triangle[i]->mNormalvector[1];
				triangle[i]->points[0]->F[2] += triangle[i]->mNormalvector[2];

				triangle[i]->points[1]->F[0] += triangle[i]->mNormalvector[0];
				triangle[i]->points[1]->F[1] += triangle[i]->mNormalvector[1];
				triangle[i]->points[1]->F[2] += triangle[i]->mNormalvector[2];

				triangle[i]->points[2]->F[0] += triangle[i]->mNormalvector[0];
				triangle[i]->points[2]->F[1] += triangle[i]->mNormalvector[1];
				triangle[i]->points[2]->F[2] += triangle[i]->mNormalvector[2];




			/*	double l = norm(mass_huell[i]->position);

				double nf = p/l;

				mass_huell[i]->F[0] -=    mass_huell[i]->position[0]*nf;
				mass_huell[i]->F[1] -=    mass_huell[i]->position[1]*nf;
				mass_huell[i]->F[2] -=    mass_huell[i]->position[2]*nf;

*/




			}
			break;
		#pragma endregion


			}//switch
		#pragma region general
			if( dim == 1){
				for( unsigned int i = 0 ; i < this->mass.size() ; i++){

					
				double l = sqrt(mass[i]->position[2]*mass[i]->position[2]+mass[i]->position[1]*mass[i]->position[1]);
				
				if ( l > 0.5){
				//	mass[i]->position[1] = mass[i]->position[1]/(1+l-0.5)*500;
				//	mass[i]->position[2] = mass[i]->position[2]/(1+l-0.5)*500;
				//	mass[i]->F[1] -= 	mass[i]->position[1]*500;
				//	mass[i]->F[2] -= 	mass[i]->position[2]*500;
					mass[i]->F[1] -= 	mass[i]->position[1]/l*(l-0.5)*500;
					mass[i]->F[2] -= 	mass[i]->position[2]/l*(l-0.5)*500;
				}
			}
			
			}
			if( dim == 2){
				for( unsigned int i = 0 ; i < this->mass.size() ; i++){
				if ( abs(mass[i]->position[2]) > 0.5){
					if ( mass[i]->position[2] > 0.5)
						//mass[i]->position[2] = 0.5;
						mass[i]->F[2] -= (mass[i]->position[2]-0.5)*500;
					else
						//mass[i]->position[2] = -0.5;
						mass[i]->F[2] -= (mass[i]->position[2]+0.5)*500;

				}
			}
			
			}//D2
			
		#pragma endregion

	



}
void CellTriangulated::updateCellCycle( double timeStep, double time){

  double r;
  switch( status){

    case 202://pressure inside
      {
      this->addInnerStress();
      r = this->calcRadius();
      p = radius_cell*radius_cell*radius_cell/r/r/r*p;
      radius_cell = r;
      break;
      }
  //CellCycle
   case 103://growth to original size (Volume)
     {
     this->v_reference += 200. * timeStep;

     if( this->v_reference <= this->init_volume ) 
       status = 102;

     break;
     }
    case 102://relax
      {

      if( presure_inside < pressure_threshold ){
        pressure_low_since++;
        if( pressure_low_since > pressure_last_number)
          pressure_high = 0;
      }
      else{
        pressure_low_since = 0;
        pressure_high = 1;
      }

      if(time_relax>relax_time)
        if( !pressure_high)
          status = 100;
      break;
      }
    case 101://shrink divide circle
      {
      {//if( presure_inside < pressure_threshold ){
        double stretch = 125.;
        for( unsigned int j = 0 ; j < spring_divide.size() ; j++)
          spring_divide[j]->k += stretch*timeStep;//0.5;
        for( unsigned int j = 0 ; j < spring_divide_cytoskeleton.size() ; j++)
          spring_divide_cytoskeleton[j]->k += stretch*timeStep;//0.5;
          //if ( divideCircleLength() < 6*meanSpringLength)
          //if ( divideCircleLength() < (sqrt(12*PI*meanSpringLength)) )//why length instead of area?
        double mean_r = 0.;
        for( unsigned int i = 0 ; i < spring_divide_cytoskeleton.size() ; i++){

        mean_r += spring_divide_cytoskeleton[i]->dlength;

        }

        mean_r /= spring_divide_cytoskeleton.size();

      if( mean_r < 0.33)
        //if ( divideCircleLength() < (sqrt(6*PI*sqrt(3.))*meanSpringLength) )
        can_divide = 1;
    }
    break;
    }
    case 100://stretch cell
      {
    	   if( presure_inside < pressure_threshold ){
        double stretch = 1+0.256*timeStep;//1.001 with 2^-8; 0.256//slow 0.03
//			double stretch = 1+0.256*timeStep;//1.001 with 2^-8; 0.256//slow 0.03
        radius_cell*= stretch;
        v_reference = 4./3.*PI*radius_cell*radius_cell*radius_cell*volume_correction;

        for( unsigned int i = 0 ; i < spring.size(); i++){
         spring[i]->l0*= stretch;
         spring[i]->l0_init*=stretch;
      }


      if( v_reference >= 1.0472){
        status = 101;
        //random direction!!!!
        double x = core->tools->random->GetRandomUniform01()*2-1;
        double y = core->tools->random->GetRandomUniform01()*2-1;
        double z = core->tools->random->GetRandomUniform01()*2-1;

        double u = radius_cell*PI*2;
        u *= 1.3; //correction amount of mass points
        double anz = (int) u/meanSpringLength;
        //this->addDivideCircle(x,y,z,(int)mass.size()/4,0.);//hier anderen Anzahl w�hlen!!!!
        this->addDivideCircle(x,y,z,anz,0.);
        this->setMatrixA();
      }
    }
    break;
    }
  }//switch

}

void CellTriangulated::checkTriangle(){
	
	double force = 5.;

	for( unsigned int i = 0 ; i < this->triangle.size() ; i++){
	
/*		double A[3][3];
		A[0] = triangle[i]->points[0]->position;
		A[1] = triangle[i]->points[1]->position;
		A[2] = triangle[i]->points[2]->position;
*/



		if( det(triangle[i]->points[0]->position,triangle[i]->points[1]->position,triangle[i]->points[2]->position) < 0. ){


					double v0[3],v1[3],v2[3];

		v0[0] = this->triangle[i]->points[0]->position[0] - this->triangle[i]->points[1]->position[0];
		v0[1] = this->triangle[i]->points[0]->position[1] - this->triangle[i]->points[1]->position[1];
		v0[2] = this->triangle[i]->points[0]->position[2] - this->triangle[i]->points[1]->position[2];

		v1[0] = this->triangle[i]->points[0]->position[0] - this->triangle[i]->points[2]->position[0];
		v1[1] = this->triangle[i]->points[0]->position[1] - this->triangle[i]->points[2]->position[1];
		v1[2] = this->triangle[i]->points[0]->position[2] - this->triangle[i]->points[2]->position[2];

		v2[0] = this->triangle[i]->points[1]->position[0] - this->triangle[i]->points[2]->position[0];
		v2[1] = this->triangle[i]->points[1]->position[1] - this->triangle[i]->points[2]->position[1];
		v2[2] = this->triangle[i]->points[1]->position[2] - this->triangle[i]->points[2]->position[2];

			double vn0[3],vn1[3],vn2[3];

		vn0[0] = v0[0]+v1[0];
		vn0[1] = v0[1]+v1[1];
		vn0[2] = v0[2]+v1[2];
		normvector(vn0);

		vn1[0] = -v0[0]+v2[0];
		vn1[1] = -v0[1]+v2[1];
		vn1[2] = -v0[2]+v2[2];
		normvector(vn1);

		vn2[0] = -v1[0]-v2[0];
		vn2[1] = -v1[1]-v2[1];
		vn2[2] = -v1[2]-v2[2];
		normvector(vn2);

		this->triangle[i]->points[0]->F[0] += vn0[0]*force;
		this->triangle[i]->points[0]->F[1] += vn0[1]*force;
		this->triangle[i]->points[0]->F[2] += vn0[2]*force;

		this->triangle[i]->points[1]->F[0] += vn1[0]*force;
		this->triangle[i]->points[1]->F[1] += vn1[1]*force;
		this->triangle[i]->points[1]->F[2] += vn1[2]*force;

		this->triangle[i]->points[2]->F[0] += vn2[0]*force;
		this->triangle[i]->points[2]->F[1] += vn2[1]*force;
		this->triangle[i]->points[2]->F[2] += vn2[2]*force;

		}


	}


}

void CellTriangulated::calcInteractionTriangle(CellTriangulated *cell_static){
	
	double force_factor = 300.;//50.;


for( unsigned int i = 0 ; i < spring_huell.size(); i++ ){//check each mass point
		for( unsigned int j = 0 ; j < cell_static->triangle.size(); j++ ){
		
			//check parallel of triangle and line
			if( dotmult(spring_huell[i]->direction,cell_static->triangle[j]->mNormalvector) != 0) {


		double u[3],v[3],w[3];

			v[0] = cell_static->triangle[j]->points[1]->position[0] - cell_static->triangle[j]->points[0]->position[0];
			v[1] = cell_static->triangle[j]->points[1]->position[1] - cell_static->triangle[j]->points[0]->position[1];
			v[2] = cell_static->triangle[j]->points[1]->position[2] - cell_static->triangle[j]->points[0]->position[2];

			u[0] = cell_static->triangle[j]->points[2]->position[0] - cell_static->triangle[j]->points[0]->position[0];
			u[1] = cell_static->triangle[j]->points[2]->position[1] - cell_static->triangle[j]->points[0]->position[1];
			u[2] = cell_static->triangle[j]->points[2]->position[2] - cell_static->triangle[j]->points[0]->position[2];

			w[0] = spring_huell[i]->start->position[0] - cell_static->triangle[j]->points[0]->position[0];
			w[1] = spring_huell[i]->start->position[1] - cell_static->triangle[j]->points[0]->position[1];
			w[2] = spring_huell[i]->start->position[2] - cell_static->triangle[j]->points[0]->position[2];


			
			double cross_d_v[3];
			crossProduct(spring_huell[i]->direction,v,cross_d_v);
			double cross_w_u[3];
			crossProduct(w,u,cross_w_u);

			double invers_d_v_u = 1./dotmult(cross_d_v,u);

			double trs[3];
			
			trs[0] = invers_d_v_u * dotmult(cross_w_u,v);
			if( trs[0] >= 0 )	{

			trs[1] = invers_d_v_u * dotmult(cross_d_v,w);
			if( trs[1] > 0 && trs[1] < 1 ){
			
			
			trs[2] = invers_d_v_u * dotmult(cross_w_u,spring_huell[i]->direction);


			

		
			if ( trs[2] > 0 && trs[2] < 1)	{
		
			double sum_r_s = trs[1]+trs[2];
			if( sum_r_s <= 1 )	{

	
			/*
			double s[3];
			s[0] = spring_huell[i]->start->position[0] + trs[0] * spring_huell[i]->direction[0];
			s[1] = spring_huell[i]->start->position[1] + trs[0] * spring_huell[i]->direction[1];
			s[2] = spring_huell[i]->start->position[2] + trs[0] * spring_huell[i]->direction[2];
			*/
		//	if( dist(s,spring_huell[i]->start->position) <= dist(spring_huell[i]->start->position,spring_huell[i]->end->position) ){
				if( trs[0] <= spring_huell[i]->l0){
		

										double normv[3];
					cell_static->triangle[j]->calcMean(normv);
					normv[0] -= cell_static->mass[0]->position[0];
					normv[1] -= cell_static->mass[1]->position[1];
					normv[2] -= cell_static->mass[2]->position[2];
					normvector(normv);






			//interaction!

			double nv_tmp[3] = {normv[0]*force_factor,
								normv[1]*force_factor,
								normv[2]*force_factor};

			//mass point inside tetrohedron


			nv_tmp[0] /=2.;
			nv_tmp[1] /=2.;
			nv_tmp[2] /=2.;

			spring_huell[i]->start->F[0] += nv_tmp[0];
			spring_huell[i]->start->F[1] += nv_tmp[1];
			spring_huell[i]->start->F[2] += nv_tmp[2];
			spring_huell[i]->end->F[0] += nv_tmp[0];
			spring_huell[i]->end->F[1] += nv_tmp[1];
			spring_huell[i]->end->F[2] += nv_tmp[2];


			nv_tmp[0] *=2./3.;
			nv_tmp[1] *=2./3.;
			nv_tmp[2] *=2./3.;

			cell_static->triangle[j]->points[0]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[0]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[0]->F[2] -= nv_tmp[2];

			cell_static->triangle[j]->points[1]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[1]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[1]->F[2] -= nv_tmp[2];

			cell_static->triangle[j]->points[2]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[2]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[2]->F[2] -= nv_tmp[2];


			}}}
			}//check length
			}//intersection
			}//check angle normalvectors


		}//cell_static->triangle
	}//spring_huell








}//calcInteractionTriangle
void CellTriangulated::calcInteractionTriangle_basic(CellTriangulated *cell_static){
	
	double force_factor = 300.;//50.;


for( unsigned int i = 0 ; i < spring_huell.size(); i++ ){//check each mass point
	for( unsigned int j = 0 ; j < cell_static->triangle.size(); j++ ){
		
			//check parallel of triangle and line
			if( dotmult(spring_huell[i]->direction,cell_static->triangle[j]->mNormalvector) != 0) {
				double intersection[3];
//				if( cell_static->triangle[j]->check_intersection_triangle_section(spring_huell[i]->direction,spring_huell[i]->start->position,spring_huell[i]->end->position,spring_huell[i]->l0,intersection) == true){
				
						double u[3],v[3];
		
			v[0] = cell_static->triangle[j]->points[1]->position[0] - cell_static->triangle[j]->points[0]->position[0];
			v[1] = cell_static->triangle[j]->points[1]->position[1] - cell_static->triangle[j]->points[0]->position[1];
			v[2] = cell_static->triangle[j]->points[1]->position[2] - cell_static->triangle[j]->points[0]->position[2];

			u[0] = cell_static->triangle[j]->points[2]->position[0] - cell_static->triangle[j]->points[0]->position[0];
			u[1] = cell_static->triangle[j]->points[2]->position[1] - cell_static->triangle[j]->points[0]->position[1];
			u[2] = cell_static->triangle[j]->points[2]->position[2] - cell_static->triangle[j]->points[0]->position[2];

				
				if( cell_static->triangle[j]->check_intersection_triangle_section(spring_huell[i]->direction,spring_huell[i]->start->position,spring_huell[i]->end->position,spring_huell[i]->l0,intersection) == true){
		//		if(	check_intersection_triangle_section(cell_static->triangle[j]->points[0]->position,v,u,cell_static->triangle[j]->normalvector,spring_huell[i]->start->position,spring_huell[i]->end->position,spring_huell[i]->l0,intersection)){

					double normv[3];
					cell_static->triangle[j]->calcMean(normv);
					normv[0] -= cell_static->mass[0]->position[0];
					normv[1] -= cell_static->mass[1]->position[1];
					normv[2] -= cell_static->mass[2]->position[2];
					normvector(normv);


					//interaction!

					double nv_tmp[3] = {normv[0]*force_factor,
					normv[1]*force_factor,
					normv[2]*force_factor};
					
					//mass point inside tetrohedron

					nv_tmp[0] /=2.;
					nv_tmp[1] /=2.;
					nv_tmp[2] /=2.;

					spring_huell[i]->start->F[0] += nv_tmp[0];
					spring_huell[i]->start->F[1] += nv_tmp[1];
					spring_huell[i]->start->F[2] += nv_tmp[2];
					spring_huell[i]->end->F[0] += nv_tmp[0];
					spring_huell[i]->end->F[1] += nv_tmp[1];
					spring_huell[i]->end->F[2] += nv_tmp[2];


					nv_tmp[0] *=2./3.;
					nv_tmp[1] *=2./3.;
					nv_tmp[2] *=2./3.;

					cell_static->triangle[j]->points[0]->F[0] -= nv_tmp[0];
					cell_static->triangle[j]->points[0]->F[1] -= nv_tmp[1];
					cell_static->triangle[j]->points[0]->F[2] -= nv_tmp[2];

					cell_static->triangle[j]->points[1]->F[0] -= nv_tmp[0];
					cell_static->triangle[j]->points[1]->F[1] -= nv_tmp[1];
					cell_static->triangle[j]->points[1]->F[2] -= nv_tmp[2];

					cell_static->triangle[j]->points[2]->F[0] -= nv_tmp[0];
					cell_static->triangle[j]->points[2]->F[1] -= nv_tmp[1];
					cell_static->triangle[j]->points[2]->F[2] -= nv_tmp[2];


			
			}//intersection
			}//check angle normalvectors


		}//cell_static->triangle
	}//spring_huell








}//calcInteractionTriangle

void CellTriangulated::calcInteractionTetrahedral_basic(CellTriangulated *cell_static){
				double F_faktor = 110.;
	
	for( unsigned int i = 0 ; i < mass.size(); i++ ){//check each mass point

			double p0i[3];
			
			p0i[0] = mass[i]->position[0] - cell_static->mass[0]->position[0];
			p0i[1] = mass[i]->position[1] - cell_static->mass[0]->position[1];
			p0i[2] = mass[i]->position[2] - cell_static->mass[0]->position[2];


		for( unsigned int j = 0 ; j < cell_static->triangle.size()   ; j++){//in each tetraeder
		//	if( cell_static->triangle[j]->check_point_inside_tetrahedron(this->mass[i]->position,p0i) == true ){
		//	if( cell_static->triangle[j]->check_point_inside_tetrahedron(this->mass[i]->position) == true ){

				if( check_intersection_tetrahedral_naive(cell_static->triangle[j],this->triangle[i])==true){

				double nv_tmp[3] = {cell_static->triangle[j]->mNormalvector[0]*F_faktor,
									cell_static->triangle[j]->mNormalvector[1]*F_faktor,
									cell_static->triangle[j]->mNormalvector[2]*F_faktor};

			//mass point inside tetrohedron

			mass[i]->F[0] += nv_tmp[0];
			mass[i]->F[1] += nv_tmp[1];
			mass[i]->F[2] += nv_tmp[2];

			nv_tmp[0] /=3.;
			nv_tmp[1] /=3.;
			nv_tmp[2] /=3.;

			cell_static->triangle[j]->points[0]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[0]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[0]->F[2] -= nv_tmp[2];

			cell_static->triangle[j]->points[1]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[1]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[1]->F[2] -= nv_tmp[2];

			cell_static->triangle[j]->points[2]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[2]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[2]->F[2] -= nv_tmp[2];

			break;
			}
		}
	}
	
}


void CellTriangulated::calcInteractionTetrahedral(CellTriangulated *cell_static){
				double F_faktor = 110.;
	
	for( unsigned int i = 0 ; i < mass.size(); i++ ){//check each mass point

			double p0i[3];
			
			p0i[0] = mass[i]->position[0] - cell_static->mass[0]->position[0];
			p0i[1] = mass[i]->position[1] - cell_static->mass[0]->position[1];
			p0i[2] = mass[i]->position[2] - cell_static->mass[0]->position[2];


		for( unsigned int j = 0 ; j < cell_static->triangle.size()   ; j++){//in each tetraeder
			double norm_invador = norm(mass[i]->position,cell_static->mass[0]->position);

			double m00[3] = {cell_static->mass[0]->position[0] - cell_static->triangle[j]->points[0]->position[0],
							 cell_static->mass[0]->position[1] - cell_static->triangle[j]->points[0]->position[1],
							 cell_static->mass[0]->position[2] - cell_static->triangle[j]->points[0]->position[2]};
	
			double m01[3] = {cell_static->mass[0]->position[0] - cell_static->triangle[j]->points[1]->position[0],
							 cell_static->mass[0]->position[1] - cell_static->triangle[j]->points[1]->position[1],
							 cell_static->mass[0]->position[2] - cell_static->triangle[j]->points[1]->position[2]};

			double m02[3] = {cell_static->mass[0]->position[0] - cell_static->triangle[j]->points[2]->position[0],
							 cell_static->mass[0]->position[1] - cell_static->triangle[j]->points[2]->position[1],
							 cell_static->mass[0]->position[2] - cell_static->triangle[j]->points[2]->position[2]};
	



			if( !(norm(m02) < norm_invador && norm(m01) < norm_invador && norm(m00) < norm_invador)){
		
			double pm0i[3];
			
			pm0i[0] = mass[i]->position[0] - cell_static->triangle[j]->points[0]->position[0];
			pm0i[1] = mass[i]->position[1] - cell_static->triangle[j]->points[0]->position[1];
			pm0i[2] = mass[i]->position[2] - cell_static->triangle[j]->points[0]->position[2];
			
			//check if point is on right side of triangle
			//first triangle at surface of cell
			if( (dotmult(cell_static->triangle[j]->mNormalvector, pm0i) < 0 )){

			//all 3 other sides


			double normvector_m0_2_1[3];
			double normvector_m0_1_0[3];
			double normvector_m0_0_2[3];

			crossProduct(cell_static->mass[0]->position,cell_static->triangle[j]->points[0]->position,cell_static->triangle[j]->points[2]->position,normvector_m0_0_2);
			if( dotmult( normvector_m0_0_2, p0i) < 0 ){
				//break;
		


			crossProduct(cell_static->mass[0]->position,cell_static->triangle[j]->points[2]->position,cell_static->triangle[j]->points[1]->position,normvector_m0_2_1);

			if( (dotmult( normvector_m0_2_1, p0i) <= 0 )){
				//break;
		



			crossProduct(cell_static->mass[0]->position,cell_static->triangle[j]->points[1]->position,cell_static->triangle[j]->points[0]->position,normvector_m0_1_0);

			if( (dotmult( normvector_m0_1_0, p0i) <= 0 )){
				//break;
		/*
				double R1 = 1./3.*(norm(m00)+norm(m01)+norm(m02));
				double tmp[3];
					tmp[0] = mass[i]->position[0] - mass[0]->position[0];
					tmp[1] = mass[i]->position[1] - mass[0]->position[1];
					tmp[2] = mass[i]->position[2] - mass[0]->position[2];
				double R2 = norm(tmp);
				double E_star = 2.*450./(1.-0.4*0.4);
				double R_star = 1./(1./R1+1./R2);

				double n[4];
					n[0] = cell_static->triangle[j]->normalvector[0];
					n[1] = cell_static->triangle[j]->normalvector[1];
					n[2] = cell_static->triangle[j]->normalvector[2];
					n[3] = n[0]* cell_static->triangle[j]->points[0]->position[0]
						  +n[1]* cell_static->triangle[j]->points[0]->position[1]
						  +n[2]* cell_static->triangle[j]->points[0]->position[2];
				
				double d = n[0]*mass[i]->position[0]
						  +n[1]*mass[i]->position[1]
						  +n[2]*mass[i]->position[2]
						  -n[4];
				d = abs(d);
				double F_faktor = 4./3.*E_star*R_star*pow(d,1.5)*1.;//10^12
				*/

				double nv_tmp[3] = {cell_static->triangle[j]->mNormalvector[0]*F_faktor,
									cell_static->triangle[j]->mNormalvector[1]*F_faktor,
									cell_static->triangle[j]->mNormalvector[2]*F_faktor};

			//mass point inside tetrohedron

			mass[i]->F[0] += nv_tmp[0];
			mass[i]->F[1] += nv_tmp[1];
			mass[i]->F[2] += nv_tmp[2];

			nv_tmp[0] /=3.;
			nv_tmp[1] /=3.;
			nv_tmp[2] /=3.;

			cell_static->triangle[j]->points[0]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[0]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[0]->F[2] -= nv_tmp[2];

			cell_static->triangle[j]->points[1]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[1]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[1]->F[2] -= nv_tmp[2];

			cell_static->triangle[j]->points[2]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[2]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[2]->F[2] -= nv_tmp[2];

			break;
			}}}}}
		}
	}
	
}
void CellTriangulated::calcInteractionTetrahedral_2(CellTriangulated *cell_static){
	
	double force_factor = 1.;
	/*
	for( unsigned int i = 0 ; i < spring_cytoscelett.size(); i++ ){//check each mass point
		for( unsigned int j = 0 ; j < cell_static->triangle.size(); j++ ){
		
			//check parallel of triangle and line
			if( dotmult(spring_cytoscelett[i]->direction,cell_static->triangle[j]->normalvector) != 0) {


		double u[3],v[3],w[3];

			v[0] = cell_static->triangle[j]->points[1]->position[0] - cell_static->triangle[j]->points[0]->position[0];
			v[1] = cell_static->triangle[j]->points[1]->position[1] - cell_static->triangle[j]->points[0]->position[1];
			v[2] = cell_static->triangle[j]->points[1]->position[2] - cell_static->triangle[j]->points[0]->position[2];

			u[0] = cell_static->triangle[j]->points[2]->position[0] - cell_static->triangle[j]->points[0]->position[0];
			u[1] = cell_static->triangle[j]->points[2]->position[1] - cell_static->triangle[j]->points[0]->position[1];
			u[2] = cell_static->triangle[j]->points[2]->position[2] - cell_static->triangle[j]->points[0]->position[2];

			w[0] = spring_cytoscelett[i]->start->position[0] - cell_static->triangle[j]->points[0]->position[0];
			w[1] = spring_cytoscelett[i]->start->position[1] - cell_static->triangle[j]->points[0]->position[1];
			w[2] = spring_cytoscelett[i]->start->position[2] - cell_static->triangle[j]->points[0]->position[2];


			
			double cross_d_v[3];
			crossProduct(spring_cytoscelett[i]->direction,v,cross_d_v);
			double cross_w_u[3];
			crossProduct(w,u,cross_w_u);

			double invers_d_v_u = 1./dotmult(cross_d_v,u);

			double trs[3];
			trs[0] = invers_d_v_u * dotmult(cross_w_u,v);
			trs[1] = invers_d_v_u * dotmult(cross_d_v,w);
			trs[2] = invers_d_v_u * dotmult(cross_w_u,spring_cytoscelett[i]->direction);


			if( trs[0] >= 0 )	{

		
			if (!( trs[1] < 0 || trs[1] > 1 || trs[2] < 0 || trs[2] > 1))	{
		
			double sum_r_s = trs[1]+trs[2];
			if( sum_r_s <= 1 )	{

	
			
			double s[3];
			s[0] = spring_cytoscelett[i]->start->position[0] + trs[0] * spring_cytoscelett[i]->direction[0];
			s[1] = spring_cytoscelett[i]->start->position[1] + trs[0] * spring_cytoscelett[i]->direction[1];
			s[2] = spring_cytoscelett[i]->start->position[2] + trs[0] * spring_cytoscelett[i]->direction[2];

			if( dist(s,spring_cytoscelett[i]->start->position) <= dist(spring_cytoscelett[i]->start->position,spring_cytoscelett[i]->end->position) ){

		








			//interaction!

			double nv_tmp[3] = {cell_static->triangle[j]->normalvector[0]*force_factor,
								cell_static->triangle[j]->normalvector[1]*force_factor,
								cell_static->triangle[j]->normalvector[2]*force_factor};

			//mass point inside tetrohedron

			spring_cytoscelett[i]->end->F[0] += nv_tmp[0];
			spring_cytoscelett[i]->end->F[1] += nv_tmp[1];
			spring_cytoscelett[i]->end->F[2] += nv_tmp[2];

			nv_tmp[0] /=3.;
			nv_tmp[1] /=3.;
			nv_tmp[2] /=3.;

			cell_static->triangle[j]->points[0]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[0]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[0]->F[2] -= nv_tmp[2];

			cell_static->triangle[j]->points[1]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[1]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[1]->F[2] -= nv_tmp[2];

			cell_static->triangle[j]->points[2]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[2]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[2]->F[2] -= nv_tmp[2];


			}}
			}//check length
			}//intersection
			}//check angle normalvectors


		}//cell_static->triangle
	}//spring_huell

	*/

for( unsigned int i = 0 ; i < spring_cytoscelett.size(); i++ ){//check each mass point
		for( unsigned int j = 0 ; j < cell_static->triangle.size(); j++ ){
		
			//check parallel of triangle and line
			if( dotmult(spring_cytoscelett[i]->direction,cell_static->triangle[j]->mNormalvector) != 0) {


		double u[3],v[3],w[3];

			v[0] = cell_static->triangle[j]->points[1]->position[0] - cell_static->triangle[j]->points[0]->position[0];
			v[1] = cell_static->triangle[j]->points[1]->position[1] - cell_static->triangle[j]->points[0]->position[1];
			v[2] = cell_static->triangle[j]->points[1]->position[2] - cell_static->triangle[j]->points[0]->position[2];

			u[0] = cell_static->triangle[j]->points[2]->position[0] - cell_static->triangle[j]->points[0]->position[0];
			u[1] = cell_static->triangle[j]->points[2]->position[1] - cell_static->triangle[j]->points[0]->position[1];
			u[2] = cell_static->triangle[j]->points[2]->position[2] - cell_static->triangle[j]->points[0]->position[2];

			w[0] = spring_cytoscelett[i]->start->position[0] - cell_static->triangle[j]->points[0]->position[0];
			w[1] = spring_cytoscelett[i]->start->position[1] - cell_static->triangle[j]->points[0]->position[1];
			w[2] = spring_cytoscelett[i]->start->position[2] - cell_static->triangle[j]->points[0]->position[2];


			
			double cross_d_v[3];
			crossProduct(spring_cytoscelett[i]->direction,v,cross_d_v);
			double cross_w_u[3];
			crossProduct(w,u,cross_w_u);

			double invers_d_v_u = 1./dotmult(cross_d_v,u);

			double trs[3];
			
			trs[0] = invers_d_v_u * dotmult(cross_w_u,v);
			if( trs[0] >= 0 )	{

			trs[1] = invers_d_v_u * dotmult(cross_d_v,w);
			if( trs[1] > 0 && trs[1] < 1 ){
			
			
			trs[2] = invers_d_v_u * dotmult(cross_w_u,spring_cytoscelett[i]->direction);


			

		
			if ( trs[2] > 0 && trs[2] < 1)	{
		
			double sum_r_s = trs[1]+trs[2];
			if( sum_r_s <= 1 )	{

	
			

		//	if( dist(s,spring_huell[i]->start->position) <= dist(spring_huell[i]->start->position,spring_huell[i]->end->position) ){
				if( trs[0] <= spring_cytoscelett[i]->l0){
		
								double s[3];
			s[0] = spring_cytoscelett[i]->start->position[0] + trs[0] * spring_cytoscelett[i]->direction[0];
			s[1] = spring_cytoscelett[i]->start->position[1] + trs[0] * spring_cytoscelett[i]->direction[1];
			s[2] = spring_cytoscelett[i]->start->position[2] + trs[0] * spring_cytoscelett[i]->direction[2];
			

					double deep =  dist(s,spring_cytoscelett[i]->end->position);//spring_huell[i]->l0 - trs[0];

					deep *= 22.;
				
					
					double correction_fator = force_factor* pow(1.+deep,3);
					/*
						ofstream F;
						F.open( "deep.dat", ios::out|ios::app );
						F << deep << "\t" << dist(s,spring_cytoscelett[i]->end->position) << "\t" <<dist(spring_cytoscelett[i]->end->position,spring_cytoscelett[i]->start->position) << "\t"  << (spring_cytoscelett[i]->l0 - trs[0]) <<  endl;
						F.close();
					*/

			//interaction!

					double normv[3];
					cell_static->triangle[j]->calcMean(normv);
					normv[0] -= cell_static->mass[0]->position[0];
					normv[1] -= cell_static->mass[1]->position[1];
					normv[2] -= cell_static->mass[2]->position[2];
					normvector(normv);

	/*		double nv_tmp[3] = {cell_static->triangle[j]->normalvector[0]*correction_fator,
								cell_static->triangle[j]->normalvector[1]*correction_fator,
								cell_static->triangle[j]->normalvector[2]*correction_fator};
	*/

			double nv_tmp[3] = {normv[0]*correction_fator,
								normv[1]*correction_fator,
								normv[2]*correction_fator};

			//mass point inside tetrohedron


			nv_tmp[0] /=2.;
			nv_tmp[1] /=2.;
			nv_tmp[2] /=2.;

	
			spring_cytoscelett[i]->end->F[0] += nv_tmp[0];
			spring_cytoscelett[i]->end->F[1] += nv_tmp[1];
			spring_cytoscelett[i]->end->F[2] += nv_tmp[2];


			nv_tmp[0] *=2./3.;
			nv_tmp[1] *=2./3.;
			nv_tmp[2] *=2./3.;

			cell_static->triangle[j]->points[0]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[0]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[0]->F[2] -= nv_tmp[2];

			cell_static->triangle[j]->points[1]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[1]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[1]->F[2] -= nv_tmp[2];

			cell_static->triangle[j]->points[2]->F[0] -= nv_tmp[0];
			cell_static->triangle[j]->points[2]->F[1] -= nv_tmp[1];
			cell_static->triangle[j]->points[2]->F[2] -= nv_tmp[2];


			}}}
			}//check length
			}//intersection
			}//check angle normalvectors


		}//cell_static->triangle
	}//spring_huell






	
}
bool CellTriangulated::calcInteraction_repulsive(CellTriangulated *cell_static){

  //check all mass points
  for( unsigned int i = 0 ; i < this->mass_huell.size() ; i++){

     int source = 0; //0 point,1 spring, 2 triangle
     int i0;
     int i1;
     int i2;

      double s0 = 10;
      double s1 = 10;
      double s2 = 10;

    //save F_repulsive from timestep bevor
    this->mass_huell[i]->F_repulsive_old = this->mass_huell[i]->F_max;
    
    //calc dist
    double tmp_dist = dist(this->mass_huell[i]->position,cell_static->mass_huell[0]->position);

    int tmp_int = 0;

    //check all other points (static)
    for( unsigned int j = 1 ; j < cell_static->mass_huell.size() ; j++){

      double tmp_for_dist = dist(this->mass_huell[i]->position,cell_static->mass_huell[j]->position);

      if( tmp_for_dist < tmp_dist){
        tmp_dist = tmp_for_dist;
        tmp_int = j;
        s0 = tmp_for_dist;
        i0 = j;
      }

    }//all static mass points


    double dir[3];
    dir[0] = cell_static->mass_huell[tmp_int]->position[0] - cell_static->mass[0]->position[0];
    dir[1] = cell_static->mass_huell[tmp_int]->position[1] - cell_static->mass[0]->position[1];
    dir[2] = cell_static->mass_huell[tmp_int]->position[2] - cell_static->mass[0]->position[2];
    normvector(dir);

    this->mass_huell[i]->dist_point = tmp_dist;
    
   
    //check distance to all triangles
    for( unsigned int j = 0 ; j < cell_static->triangle.size(); j++){
      double intersection[3];

      double tmp_normalv[3];
      tmp_normalv[0] = -cell_static->triangle[j]->mNormalvector[0];
      tmp_normalv[1] = -cell_static->triangle[j]->mNormalvector[1];
      tmp_normalv[2] = -cell_static->triangle[j]->mNormalvector[2];


      bool bt = intersect_triangle_section(//_area
          cell_static->triangle[j]->points[0]->position,
          cell_static->triangle[j]->points[1]->position,
          cell_static->triangle[j]->points[2]->position,
          this->mass_huell[i]->position,
          tmp_normalv,
          intersection);


     if( bt  == 1){

    	 double tmp_tmp_dist = dist(intersection,this->mass_huell[i]->position);

    

    	 if( tmp_tmp_dist < tmp_dist ){

         source = 2;

         tmp_int = j;
          tmp_dist = tmp_tmp_dist;
          dir[0] = cell_static->triangle[j]->mNormalvector[0];
          dir[1] = cell_static->triangle[j]->mNormalvector[1];
          dir[2] = cell_static->triangle[j]->mNormalvector[2];

          this->mass_huell[i]->dist_triangle = tmp_dist;


        }

       if( s2 > tmp_tmp_dist){
          s2 = tmp_tmp_dist;
          i2 = j;
       }
          

      }//if interact


    }//all triangles
    

    //check all springs from hull
    for( unsigned int j = 0 ; j < cell_static->spring_huell.size() ; j++ ){
  
      double intersection[3];


      double r =  dist_point_section( cell_static->spring_huell[j]->start->position,
                                      cell_static->spring_huell[j]->end->position,
                                      this->mass_huell[i]->position,
                                      intersection);
      if( r != -1 ){

        double d = dist(intersection,this->mass_huell[i]->position);

        if( s1 > d){
          s1 = d;
          i1 = j;
        }

        if( d < tmp_dist ){

          tmp_dist = d;

          dir[0] = this->mass_huell[i]->position[0] - intersection[0];
          dir[1] = this->mass_huell[i]->position[1] - intersection[1];
          dir[2] = this->mass_huell[i]->position[2] - intersection[2];
          
          tmp_int = j;

          norm(dir);

          source = 1;

        }


      }


    }





		if( tmp_dist < 1 ){


 
//      double F_LJ = pow(0.0455/tmp_dist,12)-2*pow(0.0455/tmp_dist,6);
      double F_LJ = pow(0.0455/tmp_dist,12);//24*(2*pow(0.0455/tmp_dist,13));//
//this->mass_huell[i]->F_repulsive_old = F_LJ;

      //only debug (visualization)
      if( F_LJ > 1000.){
        this->mass_huell[i]->radius = 0.03;
        this->mass_huell[i]->b = 1;
        //	F_LJ = 1000.;
        //		this->able_move = 0;
        this->able_to_move = 0;
/*
        this->mass_huell[i]->v[0] =0.;
        this->mass_huell[i]->v[1] =0.;
        this->mass_huell[i]->v[2] =0.;

        this->mass_huell[i]->store_old_v[0] =0.;
        this->mass_huell[i]->store_old_v[1] =0.;
        this->mass_huell[i]->store_old_v[2] =0.;
        
*/

      }







      if( source == 0 ){
        cell_static->mass_huell[tmp_int]->F[0] -= dir[0] * F_LJ * 0.3333;
        cell_static->mass_huell[tmp_int]->F[1] -= dir[1] * F_LJ * 0.3333;
        cell_static->mass_huell[tmp_int]->F[2] -= dir[2] * F_LJ * 0.3333;

        cell_static->mass_huell[tmp_int]->store_old_v[0] = 0;
        cell_static->mass_huell[tmp_int]->store_old_v[1] = 0;
        cell_static->mass_huell[tmp_int]->store_old_v[2] = 0;


      }

      if( source == 1){


        cell_static->spring_huell[tmp_int]->start->F[0] -= dir[0] * F_LJ * 0.5;
        cell_static->spring_huell[tmp_int]->start->F[1] -= dir[1] * F_LJ * 0.5;
        cell_static->spring_huell[tmp_int]->start->F[2] -= dir[2] * F_LJ * 0.5;

        cell_static->spring_huell[tmp_int]->end->F[0] -= dir[0] * F_LJ * 0.5;
        cell_static->spring_huell[tmp_int]->end->F[1] -= dir[1] * F_LJ * 0.5;
        cell_static->spring_huell[tmp_int]->end->F[2] -= dir[2] * F_LJ * 0.5;


        cell_static->spring_huell[tmp_int]->start->store_old_v[0] = 0;
        cell_static->spring_huell[tmp_int]->start->store_old_v[1] = 0;
        cell_static->spring_huell[tmp_int]->start->store_old_v[2] = 0;

        cell_static->spring_huell[tmp_int]->end->store_old_v[0] = 0;
        cell_static->spring_huell[tmp_int]->end->store_old_v[1] = 0;
        cell_static->spring_huell[tmp_int]->end->store_old_v[2] = 0;


      }

      if( source == 2){
        for( int j = 0 ; j < 3 ; j++){
          cell_static->triangle[tmp_int]->points[j]->F[0] -= dir[0] * F_LJ * 0.3333;
          cell_static->triangle[tmp_int]->points[j]->F[1] -= dir[1] * F_LJ * 0.3333;
          cell_static->triangle[tmp_int]->points[j]->F[2] -= dir[2] * F_LJ * 0.3333;
        
          cell_static->triangle[tmp_int]->points[j]->store_old_v[0] = 0;
          cell_static->triangle[tmp_int]->points[j]->store_old_v[1] = 0;
          cell_static->triangle[tmp_int]->points[j]->store_old_v[2] = 0;
        }
      }





      this->mass_huell[i]->F[0] += dir[0] * F_LJ;
      this->mass_huell[i]->F[1] += dir[1] * F_LJ;
      this->mass_huell[i]->F[2] += dir[2] * F_LJ;




      
/*
      this->mass_huell[i]->F_repulse[0] = dir[0] * F_LJ;
      this->mass_huell[i]->F_repulse[1] = dir[1] * F_LJ;
      this->mass_huell[i]->F_repulse[2] = dir[2] * F_LJ;

      this->mass_huell[i]->dist_max = tmp_dist;
      this->mass_huell[i]->F_max = F_LJ;
*/
    }
  }
  return this->able_to_move;
}
void CellTriangulated::calcInteraction(CellTriangulated *cell_static){
				double F_faktor = 30.;
				for( unsigned int i = 0 ; i < (this->triangle.size()-1); i++ ){//check each mass point




		for( unsigned int j = i+1 ; j < cell_static->triangle.size()   ; j++){//in each tetraeder
		//	if( cell_static->triangle[j]->check_point_inside_tetrahedron(this->mass[i]->position,p0i) == true ){
			//if( cell_static->triangle[j]->check_point_inside_tetrahedron(this->mass[i]->position) == true ){

				double V =check_intersection_tetrahedral_naive(cell_static->triangle[j],this->triangle[i]);
	
				
		//	double V=5.;
				if( V > 0 ){

					this->triangle[i]->interaction = 1;
					cell_static->triangle[j]->interaction = 1;

			//		V += 1.;
			//		V *= V*V;	
				/*	
					double p0[3] = {1.,0.,0.};
					double p1[3] = {0.,1.,0.};
					double p2[3] = {0.,0.,1.};
					double p[3] = {0.,0.,0.};
					double point[3] = {0.0,0.0,0.0}; 
					V = point_inside_tetrahedron(p0,p1,p2,p,point);

					ofstream F;
					F.open("inter_V.dat", ios::out|ios::app );
				
					F << V<< endl;
				F.close();
		
		*/
	/*
				double nv_tmp[3] = {0.,0.,0.};
		
				
				
				double tmp_this[3] = {0.,0.,0.};
				for( int l = 0 ; l < 3 ; l++){
				tmp_this[0] += this->triangle[i]->points[l]->position[0]*0.3333;
				tmp_this[1] += this->triangle[i]->points[l]->position[1]*0.3333;
				tmp_this[2] += this->triangle[i]->points[l]->position[2]*0.3333;
				}
				double tmp_static[3] = {0.,0.,0.};
				for( int l = 0 ; l < 3 ; l++){
				tmp_static[0] += cell_static->triangle[j]->points[l]->position[0]*0.3333;
				tmp_static[1] += cell_static->triangle[j]->points[l]->position[1]*0.3333;
				tmp_static[2] += cell_static->triangle[j]->points[l]->position[2]*0.3333;
				}


				nv_tmp[0] = tmp_this[0] - tmp_static[0];
				nv_tmp[1] = tmp_this[1] - tmp_static[1];
				nv_tmp[2] = tmp_this[2] - tmp_static[2];
				normvector(nv_tmp);
				*/
		/*
				double nv_tmp[3] = {0.,0.,0.};
				double tmp[3];
				for( int l = 0 ; l < 3 ; l++){
				tmp[0] = this->triangle[i]->points[l]->position[0] - cell_static->triangle[j]->points[l]->position[0];
				tmp[1] = this->triangle[i]->points[l]->position[1] - cell_static->triangle[j]->points[l]->position[1];
				tmp[2] = this->triangle[i]->points[l]->position[2] - cell_static->triangle[j]->points[l]->position[2];
				normvector(tmp);

				nv_tmp[0] += tmp[0]*0.333;
				nv_tmp[1] += tmp[1]*0.333;
				nv_tmp[2] += tmp[2]*0.333;
				}
				*/

			double nv_tmp[3] = {0.,0.,0.};
			for( int l = 0 ; l < 3 ; l++){
				nv_tmp[0] += this->triangle[i]->points[l]->position[0]*0.3333;
				nv_tmp[1] += this->triangle[i]->points[l]->position[1]*0.3333;
				nv_tmp[2] += this->triangle[i]->points[l]->position[2]*0.3333;
			}

			nv_tmp[0] -= this->mass[0]->position[0];
			nv_tmp[1] -= this->mass[0]->position[1];
			nv_tmp[2] -= this->mass[0]->position[2];

			normvector(nv_tmp);
			nv_tmp[0] *= F_faktor;
			nv_tmp[1] *= F_faktor;
			nv_tmp[2] *= F_faktor;
			/*
			nv_tmp[0] = this->triangle[i]->normalvector[0]*F_faktor,
			nv_tmp[1] = this->triangle[i]->normalvector[1]*F_faktor,
			nv_tmp[2] = this->triangle[i]->normalvector[2]*F_faktor;
		*/
			//mass point inside tetrohedron

			
			cell_static->triangle[j]->points[0]->F[0] += nv_tmp[0];
			cell_static->triangle[j]->points[0]->F[1] += nv_tmp[1];
			cell_static->triangle[j]->points[0]->F[2] += nv_tmp[2];

			cell_static->triangle[j]->points[1]->F[0] += nv_tmp[0];
			cell_static->triangle[j]->points[1]->F[1] += nv_tmp[1];
			cell_static->triangle[j]->points[1]->F[2] += nv_tmp[2];

			cell_static->triangle[j]->points[2]->F[0] += nv_tmp[0];
			cell_static->triangle[j]->points[2]->F[1] += nv_tmp[1];
			cell_static->triangle[j]->points[2]->F[2] += nv_tmp[2];
			

			cell_static->triangle[j]->next_triangle_points[0]->F[0] += nv_tmp[0];
			cell_static->triangle[j]->next_triangle_points[0]->F[1] += nv_tmp[1];
			cell_static->triangle[j]->next_triangle_points[0]->F[2] += nv_tmp[2];

			cell_static->triangle[j]->next_triangle_points[1]->F[0] += nv_tmp[0];
			cell_static->triangle[j]->next_triangle_points[1]->F[1] += nv_tmp[1];
			cell_static->triangle[j]->next_triangle_points[1]->F[2] += nv_tmp[2];

			cell_static->triangle[j]->next_triangle_points[2]->F[0] += nv_tmp[0];
			cell_static->triangle[j]->next_triangle_points[2]->F[1] += nv_tmp[1];
      cell_static->triangle[j]->next_triangle_points[2]->F[2] += nv_tmp[2];

			nv_tmp[0] = 0.;
			nv_tmp[1] = 0.;
			nv_tmp[2] = 0.;

			for( int l = 0 ; l < 3 ; l++){
				nv_tmp[0] += cell_static->triangle[j]->points[l]->position[0]*0.3333;
				nv_tmp[1] += cell_static->triangle[j]->points[l]->position[1]*0.3333;
				nv_tmp[2] += cell_static->triangle[j]->points[l]->position[2]*0.3333;
			}

			nv_tmp[0] -= cell_static->mass[0]->position[0];
			nv_tmp[1] -= cell_static->mass[0]->position[1];
			nv_tmp[2] -= cell_static->mass[0]->position[2];
			normvector(nv_tmp);
			nv_tmp[0] *= -F_faktor;
			nv_tmp[1] *= -F_faktor;
			nv_tmp[2] *= -F_faktor;
	
/*
			nv_tmp[0] = -cell_static->triangle[j]->normalvector[0]*F_faktor,
			nv_tmp[1] = -cell_static->triangle[j]->normalvector[1]*F_faktor,
			nv_tmp[2] = -cell_static->triangle[j]->normalvector[2]*F_faktor;
*/		

			this->triangle[i]->points[0]->F[0] -= nv_tmp[0];
			this->triangle[i]->points[0]->F[1] -= nv_tmp[1];
			this->triangle[i]->points[0]->F[2] -= nv_tmp[2];

			this->triangle[i]->points[1]->F[0] -= nv_tmp[0];
			this->triangle[i]->points[1]->F[1] -= nv_tmp[1];
			this->triangle[i]->points[1]->F[2] -= nv_tmp[2];

			this->triangle[i]->points[2]->F[0] -= nv_tmp[0];
			this->triangle[i]->points[2]->F[1] -= nv_tmp[1];
			this->triangle[i]->points[2]->F[2] -= nv_tmp[2];


			this->triangle[i]->next_triangle_points[0]->F[0] -= nv_tmp[0];
			this->triangle[i]->next_triangle_points[0]->F[1] -= nv_tmp[1];
			this->triangle[i]->next_triangle_points[0]->F[2] -= nv_tmp[2];

			this->triangle[i]->next_triangle_points[1]->F[0] -= nv_tmp[0];
			this->triangle[i]->next_triangle_points[1]->F[1] -= nv_tmp[1];
			this->triangle[i]->next_triangle_points[1]->F[2] -= nv_tmp[2];

			this->triangle[i]->next_triangle_points[2]->F[0] -= nv_tmp[0];
			this->triangle[i]->next_triangle_points[2]->F[1] -= nv_tmp[1];
			this->triangle[i]->next_triangle_points[2]->F[2] -= nv_tmp[2];

			break;
			}
		}
	}


}
void CellTriangulated::calcSelfInteraction(){
	double F_faktor = 2.;
	for( unsigned int i = 0 ; i < (this->triangle.size()); i++ ){//check each triangle
		for( int j = 0 ; j < 3 ; j++){//check each mass point of triangle
			bool interact = 0;
			//check point inside --> interact = 1;
			interact = point_inside_tetrahedron(
					this->triangle[i]->points[0]->position,
					this->triangle[i]->points[1]->position,
					this->triangle[i]->points[2]->position,
					this->mass[0]->position,
					this->triangle[i]->next_triangle_points[j]->position);
			if( interact == 0){
				//check intersect spring_cytoskeleton with triangle(side)
				double intersection[3];
				interact = intersect_triangle_section(
								this->triangle[i]->points[0]->position,
								this->triangle[i]->points[1]->position,
								this->triangle[i]->points[2]->position,
								this->mass[0]->position,
								this->triangle[i]->next_triangle_points[j]->position,
								intersection);
			}

			//check if intersects
			if( interact == 1){
				//add force
				double force[3];// = {0.,0.,0.};

				double f_tmp = norm(this->triangle[i]->next_triangle_points[j]->F);

				int a = mod_3(j);
				int b = mod_3(j+1);
				int c = mod_3(j+2);

				force[0] = (this->triangle[i]->points[mod_3(j)]->position[0] + this->triangle[i]->points[mod_3(j+1)]->position[0] ) *0.5;
				force[1] = (this->triangle[i]->points[mod_3(j)]->position[1] + this->triangle[i]->points[mod_3(j+1)]->position[1] ) *0.5;
				force[2] = (this->triangle[i]->points[mod_3(j)]->position[2] + this->triangle[i]->points[mod_3(j+1)]->position[2] ) *0.5;
				force[0] -= this->triangle[i]->points[mod_3(j+2)]->position[0];
				force[1] -= this->triangle[i]->points[mod_3(j+2)]->position[1];
				force[2] -= this->triangle[i]->points[mod_3(j+2)]->position[2];

				this->triangle[i]->next_triangle_points[j]->F[0] += F_faktor *force[0]*f_tmp;
				this->triangle[i]->next_triangle_points[j]->F[1] += F_faktor *force[1]*f_tmp;
				this->triangle[i]->next_triangle_points[j]->F[2] += F_faktor *force[2]*f_tmp;

			}
		}
	}
	}
void CellTriangulated::calcSelfInteraction_angle(){
  for( unsigned int i = 0 ; i < this->triangle.size() ; i++ ){
    for( int j = 0 ; j < 3 ; j++){
      double tmp_v[3];
      tmp_v[0] = (this->triangle[i]->mNormalvector[0]+this->triangle[i]->triangles[j]->mNormalvector[0])*0.5;


    }
  }
}

void CellTriangulated::calcSelfInteraction_all(){

  this->able_to_move = 1;
  double intersection[3];
  
  for( unsigned int j = 0 ; j < this->mass_huell.size() ; j++){
    for( unsigned int i = 0 ; i < this->triangle.size() ; i++){
        
      if( this->triangle[i]->points[0] != mass_huell[j] &&  this->triangle[i]->points[1] != mass_huell[j] &&  this->triangle[i]->points[2] != mass_huell[j]){
        if( point_inside_tetrahedron(
                    this->triangle[i]->points[0]->position,
                    this->triangle[i]->points[1]->position,
                    this->triangle[i]->points[2]->position,
                    this->mass[0]->position,
                    this->mass_huell[j]->position)
          ||
          intersect_triangle_section(
                    this->triangle[i]->points[0]->position,
                    this->triangle[i]->points[1]->position,
                    this->triangle[i]->points[2]->position,
                    this->mass[0]->position,
                    this->mass_huell[j]->position,
                    intersection)
         ){
          this->able_to_move = 0;
          return;
        }//if
      }//if
    }//for
  }//for





  /*
  for( unsigned int i = 0 ; i < this->mass.size() ; i++ )
    this->mass[i]->self_interaction = 0;

  for( unsigned int i = 0 ; i < (this->triangle.size()-1); i++ ){//check each triangle
    for( int j = i+1 ; j < (this->triangle.size()) ; j++){//check with each other triangle

      for( int k = 0 ; k < 3 ; k++){
      bool interact = 0;
      //check point inside --> interact = 1;
      interact = point_inside_tetrahedron(
                  this->triangle[i]->points[0]->position,
                  this->triangle[i]->points[1]->position,
                  this->triangle[i]->points[2]->position,
                  this->mass[0]->position,
                  this->triangle[i]->next_triangle[j]->position);

      if( interact == 0){
        //check intersect spring_cytoskeleton with triangle(side)
        double intersection[3];
        interact = intersect_triangle_section(
                    this->triangle[i]->points[0]->position,
                    this->triangle[i]->points[1]->position,
                    this->triangle[i]->points[2]->position,
                    this->mass[0]->position,
                    this->triangle[i]->next_triangle[j]->position,
                    intersection);
      }

			//check if intersects
			if( interact == 1){
				//add force
				double force[3];// = {0.,0.,0.};

				double f_tmp = norm(this->triangle[i]->next_triangle[j]->F);

				int a = mod_3(j);
				int b = mod_3(j+1);
				int c = mod_3(j+2);

				force[0] = (this->triangle[i]->points[mod_3(j)]->position[0] + this->triangle[i]->points[mod_3(j+1)]->position[0] ) *0.5;
				force[1] = (this->triangle[i]->points[mod_3(j)]->position[1] + this->triangle[i]->points[mod_3(j+1)]->position[1] ) *0.5;
				force[2] = (this->triangle[i]->points[mod_3(j)]->position[2] + this->triangle[i]->points[mod_3(j+1)]->position[2] ) *0.5;
				force[0] -= this->triangle[i]->points[mod_3(j+2)]->position[0];
				force[1] -= this->triangle[i]->points[mod_3(j+2)]->position[1];
				force[2] -= this->triangle[i]->points[mod_3(j+2)]->position[2];

				this->triangle[i]->next_triangle[j]->F[0] += F_faktor *force[0]*f_tmp;
				this->triangle[i]->next_triangle[j]->F[1] += F_faktor *force[1]*f_tmp;
				this->triangle[i]->next_triangle[j]->F[2] += F_faktor *force[2]*f_tmp;

			}
		}
	}
  */
	}
void CellTriangulated::calcSelfInteraction_all2(){

//  for( unsigned int i = 0 ; i < this->triangle.size() ; i++ )
//   this->triangle[i]->transparency = 0.;


  this->able_to_move = 1;
  double intersection[3];
  
  for( unsigned int j = 0 ; j < this->mass_huell.size() ; j++){
    if ( this->mass_huell[j]->divide_number != 4){
      for( unsigned int i = 0 ; i < this->triangle.size() ; i++){
        if(    this->triangle[i]->points[0]->divide_number != 4 
            && this->triangle[i]->points[1]->divide_number != 4 
            && this->triangle[i]->points[2]->divide_number != 4){
          if( this->triangle[i]->points[0] != mass_huell[j] &&  this->triangle[i]->points[1] != mass_huell[j] &&  this->triangle[i]->points[2] != mass_huell[j]){
            if(/* point_inside_tetrahedron(
                  this->triangle[i]->points[0]->position,
                  this->triangle[i]->points[1]->position,
                  this->triangle[i]->points[2]->position,
                  this->mass[0]->position,
                  this->mass_huell[j]->position)
                  ||*/
                intersect_triangle_section(
                  this->triangle[i]->points[0]->position,
                  this->triangle[i]->points[1]->position,
                  this->triangle[i]->points[2]->position,
                  this->mass[0]->position,
                  this->mass_huell[j]->position,
                  intersection)
            ){
              this->triangle[j]->transparency = 0.85;
              this->able_to_move = 0;
              return;
            }//if inside somehow
          }//if all triangle points huell points
        }//if triangle divide_number
      }//all triangles
    }//if divide_number
  }//all mass points

}

bool CellTriangulated::calcSelfInteraction_reverseStep(){

  return 1;

	}

bool CellTriangulated::check_overlap(CellTriangulated *cell_static){

  for( unsigned int i = 0 ; i < this->mass_huell.size() ; i++){
    for( unsigned int j = 0 ; j < cell_static->triangle.size() ; j++){
      if( point_inside_tetrahedron(
              cell_static->triangle[j]->points[0]->position,
              cell_static->triangle[j]->points[1]->position,
              cell_static->triangle[j]->points[2]->position,
              cell_static->mass[0]->position,
              this->mass_huell[i]->position)){
          this->mass_huell[i]->g = 255;
      //    this->mass_huell[i]->radius = 0.0222;

          double intersection[3];

          double tmp_norm[3];
          tmp_norm[0] = -cell_static->triangle[j]->mNormalvector[0];
          tmp_norm[1] = -cell_static->triangle[j]->mNormalvector[1];
          tmp_norm[2] = -cell_static->triangle[j]->mNormalvector[2];

          intersect_triangle_section(//_area
            cell_static->triangle[j]->points[0]->position,
            cell_static->triangle[j]->points[1]->position,
            cell_static->triangle[j]->points[2]->position,
            this->mass_huell[i]->position,
            tmp_norm,
            intersection);

          this->mass_huell[i]->length_inside = dist(intersection,this->mass_huell[i]->position);

                    tmp_norm[0] = -cell_static->triangle[j]->store_normalvetor[0];
          tmp_norm[1] = -cell_static->triangle[j]->store_normalvetor[1];
          tmp_norm[2] = -cell_static->triangle[j]->store_normalvetor[2];


                    intersect_triangle_section(//_area
            cell_static->triangle[j]->points[0]->store_position,
            cell_static->triangle[j]->points[1]->store_position,
            cell_static->triangle[j]->points[2]->store_position,
            this->mass_huell[i]->store_position,
            tmp_norm,
            intersection);

this->mass_huell[i]->length_inside_old = dist(intersection,this->mass_huell[i]->store_position);


          cell_static->triangle[j]->points[0]->radius = 0.0211;
          cell_static->triangle[j]->points[1]->radius = 0.0211;
          cell_static->triangle[j]->points[2]->radius = 0.0211;
          
         cell_static->triangle[j]->points[0]->b = 1;
         cell_static->triangle[j]->points[1]->b = 1;
         cell_static->triangle[j]->points[2]->b = 1;

         this->mass_huell[i]->radius = 0.03;
         this->mass_huell[i]->b = 1;

          return 1;
      }
    }
  }
  return 0;

}

bool CellTriangulated::check_overlap_past(CellTriangulated *cell_static){

  for( unsigned int i = 0 ; i < this->mass_huell.size() ; i++){
    for( unsigned int j = 0 ; j < cell_static->triangle.size() ; j++){
      if( point_inside_tetrahedron(
              cell_static->triangle[j]->points[0]->store_position,
              cell_static->triangle[j]->points[1]->store_position,
              cell_static->triangle[j]->points[2]->store_position,
              cell_static->mass[0]->store_position,
              this->mass_huell[i]->store_position)){
       //   this->mass_huell[i]->b = 1;
        //  this->mass_huell[i]->radius = 0.01;
    	  return 1;
      }
    }
  }
  return 0;

}


//create new cells
void CellTriangulated::triangulate(){

  point_mass_invers = 1./point_mass;

  // Create Triangulation
  Triangulation<3> *myTriangulation = new Triangulation<3>();

  // Add N Points to Spherical Surface
  double center[3] = {position.x,position.y,position.z};//{position[0],y,z};

  myTriangulation->setPointsOnSphericalSurface( radius_cell, mass_number, 1, center);

	//	center[0] += 13;
	//	myTriangulation->setPointsOnSphericalSurface( radius, m_number, 0.0, center);
	//	m_number*=2;

	// Set Domain Size (automatically)
	myTriangulation->setDomain();

	// Set Frame Points
	myTriangulation->setFramePoints();


	// Triangulate Vertices
	myTriangulation->triangulate();
	
	/*insertVoronoiCell( myTriangulation, new Vertex(-7., -7., 7.));
	insertVoronoiCell( myTriangulation, new Vertex(-7., -6., 8.));
	insertVoronoiCell( myTriangulation, new Vertex(-7., -8., 6.));*/
	//	myTriangulation->printToPovray( "test.pov", true, false, true, false, false);

	// Creates the Convex Hull of
	double maximalVertexDistanceOnSurface = 2*radius_cell;
	myTriangulation->getConvexHull(maximalVertexDistanceOnSurface);
	//myTriangulation->printToPovray( "convexHull.pov", true, false, false, false, true);


	std::cout << std::endl << std::endl << std::endl << "Triangulierung abgeschlossen" << std::endl;



	//my use

	MassPoint *nucleo = new MassPoint();
	mass.push_back( nucleo );
//	mass_nucleo.push_back( nucleo );
	nucleo->position[0] = position.x;
	nucleo->position[1] = position.y;
	nucleo->position[2] = position.z;

	//initialize mass vector
	//initialize and store positions in mass points
	for( unsigned int i = 0 ; i < mass_number ; i++){

		MassPoint *temp_huell = new MassPoint();
		mass_huell.push_back( temp_huell );

		for( int j = 0 ; j < myTriangulation->getCountConvexFaces() ; j++){
			int tmp = -1;
			if( (myTriangulation->getConvexFace(j)[0]->getIndex() == (int)i) )
				tmp = 0;
			if( (myTriangulation->getConvexFace(j)[1]->getIndex() == (int)i) )
				tmp = 1;
			if( (myTriangulation->getConvexFace(j)[2]->getIndex() == (int)i) )
				tmp = 2;

			if( tmp >= 0){

				mass_huell[i]->position[0] = myTriangulation->getConvexFace(j)[tmp]->x();
				mass_huell[i]->position[1] = myTriangulation->getConvexFace(j)[tmp]->y();
				mass_huell[i]->position[2] = myTriangulation->getConvexFace(j)[tmp]->z();

				//connect points
				Spring *tmp_spring = new Spring();

				double v[3];


				v[0] = mass_huell[i]->position[0] - nucleo->position[0];
				v[1] = mass_huell[i]->position[1] - nucleo->position[1];
				v[2] = mass_huell[i]->position[2] - nucleo->position[2];

				tmp_spring->l0 = norm(v);
				tmp_spring->l0_init = tmp_spring->l0;
				tmp_spring->k = k_cytoskelett;
				tmp_spring->k2 = k2;
				tmp_spring->nu = nu_cytoskelett;
				tmp_spring->start = nucleo;
				tmp_spring->end = temp_huell;

//				tmp_spring->l_damper = 0.;
				tmp_spring->eta_damper = eta_damper;
		

				
				temp_huell->neighbourSprings.push_back( tmp_spring );
				nucleo->neighbourSprings.push_back( tmp_spring );

				spring_cytoscelett.push_back( tmp_spring );

				break;
			}
		}
	}


		for( unsigned int i = 0 ; i < this->mass_huell.size(); i++)
			this->mass.push_back( this->mass_huell[i]);



		std::vector< std::vector<int> > count_springs_v(mass_number);

		for( unsigned int i = 0 ; i < mass_number ; i++){

			for( int j = 0 ; j < myTriangulation->getCountConvexFaces() ; j++){

				//search for each mass point
				int found = -1;
				if( (myTriangulation->getConvexFace(j)[0]->getIndex() == (int)i) )
					found = 0;
				if( (myTriangulation->getConvexFace(j)[1]->getIndex() == (int)i) )
					found = 1;
				if( (myTriangulation->getConvexFace(j)[2]->getIndex() == (int)i) )
					found = 2;

				//found mass point in one convex face
				if( found >= 0){//maybe this and d-loop together
					for( int d = 0 ; d < 3 ; d++){
						if(myTriangulation->getConvexFace(j)[d]->getIndex() > (int)i){
							int tmp10 = 0;
							for(unsigned int l = 0 ; l < count_springs_v[i].size() ; l++){
								if( count_springs_v[i][l] == myTriangulation->getConvexFace(j)[d]->getIndex() ){
									tmp10 = 1;
									break;
								}
							}//alle bisherigen durchgehen
							if( tmp10 == 0 ){

								Spring *tmp_spring_huell = new Spring();


								tmp_spring_huell->k = k_hull;
								tmp_spring_huell->nu= nu_hull;

								
								tmp_spring_huell->start = mass_huell[i];
								tmp_spring_huell->end = mass_huell[myTriangulation->getConvexFace(j)[d]->getIndex()];
								
//								tmp_spring_huell->l_damper = 0.;
							

								tmp_spring_huell->eta_damper = eta_damper;

								this->mass_huell[i]->neighbourSprings.push_back( tmp_spring_huell );
								this->mass_huell[myTriangulation->getConvexFace(j)[d]->getIndex()]->neighbourSprings.push_back( tmp_spring_huell );

								count_springs_v[i].push_back(myTriangulation->getConvexFace(j)[d]->getIndex());
								count_springs_v[myTriangulation->getConvexFace(j)[d]->getIndex() ].push_back(i);

								double v[3];
								v[0] = this->mass_huell[i]->position[0] - this->mass_huell[myTriangulation->getConvexFace(j)[d]->getIndex()]->position[0];
								v[1] = this->mass_huell[i]->position[1] - this->mass_huell[myTriangulation->getConvexFace(j)[d]->getIndex()]->position[1];
								v[2] = this->mass_huell[i]->position[2] - this->mass_huell[myTriangulation->getConvexFace(j)[d]->getIndex()]->position[2];

								tmp_spring_huell->l0 = norm(v);
								tmp_spring_huell->l0_init = tmp_spring_huell->l0;
								tmp_spring_huell->k2 = k2;



								this->spring_huell.push_back(tmp_spring_huell);

							}//neue Feder gefunden und hinzugefuegt
						}
					}//for all 3 points in triangle
				}//mass point found in one triangle
			}//search mass points in all triangles
		}//looking for all mass points


		for( unsigned int i = 0 ; i < this->spring_huell.size(); i++){
			spring.push_back( this->spring_huell[i]);

		}

//TODO
		//delete count_spring_v

	//cytoskelett to spring
		for( unsigned int i = 0 ; i < this->spring_cytoscelett.size() ; i ++)
			spring.push_back(this->spring_cytoscelett[i]);




        for( unsigned int i = 0 ; i < mass.size() ; i++ )
      mass[i]->index = i;
    for( unsigned int i = 0 ; i < spring.size() ; i++ )
      spring[i]->index = i;



    for( unsigned int i = 0 ; i < this->mass.size() ; i++ )
      mass[i]->mSpringsOnSurface = this->mass[i]->neighbourSprings.size() -1;


//l0

		//add triangles to cell
		for( int i = 0 ; i < myTriangulation->getCountConvexFaces() ; i++){

				Triangle *tmp_triangle_huell = new Triangle();
				triangle.push_back( tmp_triangle_huell );


				tmp_triangle_huell->points[0] = mass_huell[myTriangulation->getConvexFace(i)[0]->getIndex()];
				tmp_triangle_huell->points[1] = mass_huell[myTriangulation->getConvexFace(i)[1]->getIndex()];
				tmp_triangle_huell->points[2] = mass_huell[myTriangulation->getConvexFace(i)[2]->getIndex()];

				for( unsigned int j = 0 ; j < tmp_triangle_huell->points[0]->neighbourSprings.size() ; j++){
					if( tmp_triangle_huell->points[0]->neighbourSprings[j]->start == tmp_triangle_huell->points[1]  || 
						tmp_triangle_huell->points[0]->neighbourSprings[j]->end == tmp_triangle_huell->points[1] ){
							tmp_triangle_huell->springs[0] = tmp_triangle_huell->points[0]->neighbourSprings[j];
						break;
					}

				}
				for( unsigned int j = 0 ; j < tmp_triangle_huell->points[1]->neighbourSprings.size() ; j++){
					if( tmp_triangle_huell->points[1]->neighbourSprings[j]->start == tmp_triangle_huell->points[2]  || 
						tmp_triangle_huell->points[1]->neighbourSprings[j]->end == tmp_triangle_huell->points[2] ){
							tmp_triangle_huell->springs[1] = tmp_triangle_huell->points[1]->neighbourSprings[j];
						break;
					}

				}
				for( unsigned int j = 0 ; j < tmp_triangle_huell->points[2]->neighbourSprings.size() ; j++){
					if( tmp_triangle_huell->points[2]->neighbourSprings[j]->start == tmp_triangle_huell->points[0]  || 
						tmp_triangle_huell->points[2]->neighbourSprings[j]->end == tmp_triangle_huell->points[0] ){
							tmp_triangle_huell->springs[2] = tmp_triangle_huell->points[2]->neighbourSprings[j];
						break;
					}

				}
				for( unsigned int j = 0 ; j < tmp_triangle_huell->points[0]->neighbourSprings.size() ; j++){
					if( tmp_triangle_huell->points[0]->neighbourSprings[j]->start == mass[0]  || 
						tmp_triangle_huell->points[0]->neighbourSprings[j]->end == mass[0] ){
							tmp_triangle_huell->springs_cytoskeleton[0] = tmp_triangle_huell->points[0]->neighbourSprings[j];
						break;
					}

				}
				for( unsigned int j = 0 ; j < tmp_triangle_huell->points[1]->neighbourSprings.size() ; j++){
					if( tmp_triangle_huell->points[1]->neighbourSprings[j]->start == mass[0]  || 
						tmp_triangle_huell->points[1]->neighbourSprings[j]->end == mass[0] ){
							tmp_triangle_huell->springs_cytoskeleton[1] = tmp_triangle_huell->points[1]->neighbourSprings[j];
						break;
					}

				}
				for( unsigned int j = 0 ; j < tmp_triangle_huell->points[2]->neighbourSprings.size() ; j++){
					if( tmp_triangle_huell->points[2]->neighbourSprings[j]->start == mass[0]  || 
						tmp_triangle_huell->points[2]->neighbourSprings[j]->end == mass[0] ){
							tmp_triangle_huell->springs_cytoskeleton[2] = tmp_triangle_huell->points[2]->neighbourSprings[j];
						break;
					}

				}

				

			
				tmp_triangle_huell->nuc = mass[0];
				
				(*tmp_triangle_huell).checkAndChangeOrientation();
				tmp_triangle_huell->setNormalVector();

		}


        for( unsigned int i = 0 ; i < triangle.size() ; i++ )
      triangle[i]->index = i;
	

//hier noch springs heraussuchen
		for( unsigned int i = 0 ;  i < triangle.size() ; i++ ){
			
			//spring cytoskeleton
			for( unsigned j = 0 ; j < triangle[i]->points[0]->neighbourSprings.size() ; j++){
				if( (triangle[i]->points[0]->neighbourSprings[j]->start == triangle[i]->points[0] &&
					 triangle[i]->points[0]->neighbourSprings[j]->end   == triangle[i]->points[1])|| 
					(triangle[i]->points[0]->neighbourSprings[j]->start == triangle[i]->points[1] &&
					 triangle[i]->points[0]->neighbourSprings[j]->end   == triangle[i]->points[0]) ){
					
					triangle[i]->springs[0]=triangle[i]->points[0]->neighbourSprings[j];
					break;
				}//if
			}//for
			for( unsigned j = 0 ; j < triangle[i]->points[1]->neighbourSprings.size() ; j++){
				if( (triangle[i]->points[1]->neighbourSprings[j]->start == triangle[i]->points[1] &&
					 triangle[i]->points[1]->neighbourSprings[j]->end   == triangle[i]->points[2])|| 
					(triangle[i]->points[1]->neighbourSprings[j]->start == triangle[i]->points[2] &&
					 triangle[i]->points[1]->neighbourSprings[j]->end   == triangle[i]->points[1]) ){
					
					triangle[i]->springs[1]=triangle[i]->points[1]->neighbourSprings[j];
					break;
				}//if

			}//for
			for( unsigned j = 0 ; j < triangle[i]->points[2]->neighbourSprings.size() ; j++){
				if( (triangle[i]->points[2]->neighbourSprings[j]->start == triangle[i]->points[0] &&
					 triangle[i]->points[2]->neighbourSprings[j]->end   == triangle[i]->points[2])|| 
					(triangle[i]->points[2]->neighbourSprings[j]->start == triangle[i]->points[2] &&
					 triangle[i]->points[2]->neighbourSprings[j]->end   == triangle[i]->points[0]) ){
					
					triangle[i]->springs[2]=triangle[i]->points[2]->neighbourSprings[j];
					break;
				}//if
			}//for

		}//for all triangled



    //search third point of neigbourung triangle
		for( unsigned int i = 0 ;  i < triangle.size() ; i++ ){
			triangle[i]->next_triangle_points[0] = 0;
			triangle[i]->next_triangle_points[1] = 0;
			triangle[i]->next_triangle_points[2] = 0;
		}

    
    //search next triangle and set third not equal point
    for( unsigned int i = 0 ;  i < triangle.size() ; i++ ){ //check all triangles
      for( int k = 0 ; k < 3 ; k++){                        //all sides of triangles
       for( unsigned j = 0 ; j < triangle.size() ; j++){    //search in all other triangles
//          if( triangle[i]->equal(triangle[j]) == 2 ){       //check equal
         if( triangle[i] != triangle[j]){
         if( triangle[j]->equal(this->triangle[i]->points[k], this->triangle[i]->points[mod_3(k+1)]) == 1 ){       //check equal
            //triangle[i]->difference(triangle[j]);           //set neigbouring triangle
            triangle[i]->triangles[k] = triangle[j];
            triangle[i]->next_triangle_points[k] = triangle[j]->points[ triangle[j]->missing_point( triangle[i] ) ];
          }
         }}
      }
    }
    

		// Destroy Triangulation
		delete myTriangulation;

		meanSpringLength = 0.;
		for( unsigned int i = 0; i < spring_huell.size() ; i++)
			meanSpringLength += spring_huell[i]->l0;
		
		meanSpringLength /= (double)spring_huell.size();

	v_reference  = 0.;
	for( unsigned int i = 0 ; i < triangle.size() ;i++ ){
		triangle[i]->calcVolume(mass[0]->position);
		v_reference  += triangle[i]->Volume;
	}

  double v_tmp = 4./3.*PI*0.5*0.5*0.5;
  volume_correction = v_reference/v_tmp;


  //set index of elements
  //of mass points
  for( unsigned int i = 0 ; i < this->mass.size() ; i++ )
    this->mass[i]->mIndex = i;
  //of springs
  for( unsigned int i = 0 ; i < this->spring.size() ; i++ )
    this->spring[i]->mIndex = i;
  //of triangles
  for( unsigned int i = 0 ; i < this->triangle.size() ; i++ )
    this->triangle[i]->mIndex = i;


//TODO move more stuff to createInternalStructureOfCell()

  this->createInternalStructureOfCell();

}

void CellTriangulated::setMatrixA(){



  //initialize connection matrix
  int m = this->mass.size();

  double **tmp = (double**) malloc( sizeof(double*)*m);
  double **Ao = (double**) malloc( sizeof(double*)*m);
  for( int i = 0 ; i < m ; i++){
    Ao[i] = (double*) malloc( sizeof(double)*m);
    tmp[i] = (double*) malloc( sizeof(double)*m);
  }

  if( A == NULL){
    A = (double**) malloc( sizeof(double*)*m);
    for( int i = 0 ; i < m ; i++)
      A[i] = (double*) malloc( sizeof(double)*m);
  }



//  int* debug = (int*) malloc( sizeof(int)*m);

  //set connection matrix (count connections)
  for( int i = 0 ; i < m ; i++ )
    for( int j = 0 ; j < m ; j++ )
      tmp[i][j] = 0.;

  for( unsigned int i = 0 ; i < this->spring.size() ; i++ ){
    tmp[this->spring[i]->end->mIndex][this->spring[i]->start->mIndex] += 1;
    tmp[this->spring[i]->start->mIndex][this->spring[i]->end->mIndex] += 1;
    tmp[this->spring[i]->start->mIndex][this->spring[i]->start->mIndex] += 1;
    tmp[this->spring[i]->end->mIndex][this->spring[i]->end->mIndex] += 1;
  }



//  double gamma = 0.8;
//  double eta = 10.;

  SparseMatrix<double> *sm = new SparseMatrix<double>(m,m);
/*
  for( int i = 0 ; i < m ; i++){
    for( int j = 0 ; j < m ; j++){
      if( tmp[i][j] == 0)
        continue;
      if( tmp[i][j] == 1){
        sm->set(i,j,this->eta);
        Ao[i][j] = this->eta;
        }
      else{
        sm->set(i,j,-tmp[i][j]*this->eta+this->gamma);
        Ao[i][j] = -tmp[i][j]*this->eta+this->gamma;
        if( i != j)
          fprintf(stderr, "found problems in matrix\n");
      }
    }
  }*/

    //Debug Johannes : Test new Matrix
    for( int i = 0 ; i < m ; i++){
        for( int j = 0 ; j < m ; j++){
          if( tmp[i][j] == 1){
            sm->set(i,j,-this->eta);
    //        Ao[i][j] = eta;
            }
          }
    //      eye->set(i,j,0);
      sm->set(i,i,tmp[i][i]*this->eta + this->gamma);
      }






//solve system for each column
  for( int i = 0 ; i < m ; i++){

    double *b = (double*) malloc( sizeof(double)*m);
    double *x = (double*) malloc( sizeof(double)*m);

    for( int j = 0 ; j < m ; j++){
      if( j == i )
        b[j] = 1;
      else
        b[j] = 0;
    }

    Solver<double> *S = new Solver<double>(m,Solver<double>::BiCGSTAB);//BiCGSTAB);CG
     S->setIterations( 1000);
     S->setError( 1e-20);
     S->solve(sm,b,x);// solve( sm, b, sv);//debug[i] =
    delete S;
    free(b);

    for( int j = 0 ; j < m ; j++){
      A[j][i] = x[j];
    }
    free(x);


  }

   /*
    //debug output of connections
    ofstream F;
    F.open("../../output/A.txt", ios::out|ios::app );
    F << "connection matrix - start" << std::endl;
           for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
              for( unsigned int j = 0 ; j < this->mass.size() ; j++ ){
                F << tmp[i][j] << "\t";
              }
              F << std::endl;
            }
    F << "connection matrix - end" << std::endl << std::endl;

    F << "phy. connection matrix - start" << std::endl;
           for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
             for( unsigned int j = 0 ; j < this->mass.size() ; j++ ){
               F << sm->get(i,j) << "\t";
             }
             F << std::endl;
           }
    F << "phy. connection matrix - end" << std::endl << std::endl;

    F << "invert phy. connection matrix - start" << std::endl;
           for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
             for( unsigned int j = 0 ; j < this->mass.size() ; j++ ){
               F << A[i][j] << "\t";
             }
             F << std::endl;
           }
    F << "invert phy. connection matrix - end" << std::endl << std::endl;

  //mutliply A*inv(A) --> eye??
    for( int i = 0 ; i < m ; i++){
           for( int j = 0 ; j< m ; j++){
                   Ao[i][j] = 0;
                   for( int k = 0 ; k < m ; k++){
                           Ao[i][j] += sm->get(i,k)*A[k][j];
                   }
           }
    }


    F << "control matrix - start" << std::endl;
           for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
             for( unsigned int j = 0 ; j < this->mass.size() ; j++ ){
               F << Ao[i][j] << "\t";
             }
             F << std::endl;
           }
    F << "control matrix - end" << std::endl << std::endl;


  */


/*
  //debug output of connections
  ofstream F;
  F.open("../../output/tmp.txt", ios::out|ios::app );

  F << "connection matrix (count) " << std::endl;
  for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
	    for( unsigned int j = 0 ; j < this->mass.size() ; j++ ){
	      F << tmp[i][j] << "\t";
	    }
	    F << std::endl;
	  }
  F << "end connection matrix " << std::endl << std::endl;


  F << "connection matrix (phys) " << std::endl;
  for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
	    for( unsigned int j = 0 ; j < this->mass.size() ; j++ ){
	      F <<  sm->get(i,j) << "\t";
	    }
	    F << std::endl;
	  }
  F << "end connection matrix " << std::endl << std::endl;


  F << "inverse connection matrix (phys) " << std::endl;
  for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
	    for( unsigned int j = 0 ; j < this->mass.size() ; j++ ){
	      F <<  A[i][j] << "\t";
	    }
	    F << std::endl;
	  }
  F << "end connection matrix " << std::endl << std::endl;

	  double **tmp2 = (double**) malloc( sizeof(double*)*m);
	  for( int i = 0 ; i < m ; i++){
	    tmp2[i] = (double*) malloc( sizeof(double)*m);
	  }

	  for( int i = 0 ; i < m ; i++){
		  for( int j = 0 ; j< m ; j++){
			  tmp2[i][j] = 0.;
			  for( int k = 0 ; k < m ; k++){
				  tmp2[i][j] += Ao[i][k]*A[k][j];
			  }
		  }
	  }
	  F << std::endl;
	  F << std::endl;
	  F << std::endl;
	  F << "product of matrices" << std::endl;
	  for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
	      for( unsigned int j = 0 ; j < this->mass.size() ; j++ ){
	        F << tmp2[i][j] << "\t";
	      }
	      F << std::endl;
	    }
	  F << "end of matricies product" << std::endl;

	  */
/*
	  F << std::endl;
	  F << std::endl;
	  F << std::endl;
	  for( unsigned int j = 0 ; j < this->mass.size() ; j++ ){
	      F << debug[j] << "\t";
	    }
	    F << std::endl;


	  F << std::endl;
	  F << std::endl;
	  F << std::endl;
*/

	  for( int i = 0 ; i < m ; i++){
		  free(tmp[i]);
	  }
	  free(tmp);






/*

	  for( int i = 0 ; i < m ; i++){
	  	  free(tmp2[i]);
	    }
	    free(tmp2);






	  for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
	    for( unsigned int j = 0 ; j < this->mass.size() ; j++ ){
	      F << A[i][j] << "\t";
	    }
	    F << std::endl;

	  }


	  F.close();
*/

}

void CellTriangulated::copyfrom( CellTriangulated *default_cell, double x, double y, double z, double r){

  mpRandom = default_cell->mpRandom;

//volume
  this->v_reference = default_cell->v_reference;
  this->init_volume = default_cell->v_reference;
  volume_correction = default_cell->volume_correction;

//cell type
  type = default_cell->type;
  status = default_cell->status;//0 normal, 1 divide circle//10 walls//11 gegen Wand

//center and radii
  position.x = default_cell->position.x + x;
  position.y = default_cell->position.y + y;
  position.z = default_cell->position.z + z;
  radius_cell = r;

//physical properties
  this->eta = default_cell->eta;
  this->gamma = default_cell->gamma;
  this->pressure_threshold = default_cell->pressure_threshold;
  nu_medium = default_cell->nu_medium;//for enviroment
  p = default_cell->p;//inner stress
  k_hull = default_cell->k_hull;
  nu_hull = default_cell->nu_hull;
  k_cytoskelett = default_cell->k_cytoskelett;
  nu_cytoskelett = default_cell->nu_cytoskelett;
  k2 = default_cell->k2;
  eta_damper = default_cell->eta_damper;

  point_mass = default_cell->point_mass;
  point_mass_invers = default_cell->point_mass_invers;
	
  dim = default_cell->dim;
  mass_number = default_cell->mass_number;

  meanSpringLength = default_cell->meanSpringLength;

//cell cycle time
  SetCycleTimeGaussClamped(1.00,0.4);
  relax_time = this->cycleTime;
  //relax_time = default_cell->relax_time;
  relax_time = 2.;




//copy Elements
  for( unsigned int i = 0 ; i < default_cell->mass.size() ; i++ ){
    default_cell->mass[i]->divide_number = i;

    MassPoint *mass_tmp = new MassPoint();

    mass_tmp->copyfrom(default_cell->mass[i]);

    mass_tmp->position[0] *= r*2;
    mass_tmp->position[1] *= r*2;
    mass_tmp->position[2] *= r*2;

    mass_tmp->position[0] += x;
    mass_tmp->position[1] += y;
    mass_tmp->position[2] += z;

    mass.push_back(mass_tmp);

    if( i != 0 )
      mass_huell.push_back(mass_tmp);

  }

	bool tmp = 0;
	for( unsigned int i = 0 ; i < default_cell->spring.size() ; i++){
		
		default_cell->spring[i]->divideCircle = i;

		if( default_cell->spring[i] == default_cell->spring_cytoscelett[0] )
			tmp = 1;
	

		Spring *spring_tmp = new Spring();
		spring_tmp->copyfrom(default_cell->spring[i]);
		spring.push_back(spring_tmp);

		if( tmp == 0)
			spring_huell.push_back(spring_tmp);
		else
			spring_cytoscelett.push_back(spring_tmp);

	//connection springs ans mass points
		spring_tmp->start = mass[ spring_tmp->start->divide_number ];
		spring_tmp->end = mass[ spring_tmp->end->divide_number ];

		spring_tmp->start->neighbourSprings.push_back(spring_tmp);
		spring_tmp->end->neighbourSprings.push_back(spring_tmp);
		


	}	
	
	
	//vector <Triangle*> triangle;
	for( unsigned int l = 0 ; l < default_cell->triangle.size() ; l++ ){
		Triangle *triangle_tmp = new Triangle();
		triangle_tmp->copyfrom(default_cell->triangle[l]);
		
		triangle_tmp->points[0] = mass[ default_cell->triangle[l]->points[0]->index ];
		triangle_tmp->points[1] = mass[ default_cell->triangle[l]->points[1]->index ];
		triangle_tmp->points[2] = mass[ default_cell->triangle[l]->points[2]->index ];

		triangle_tmp->springs[0] = spring[ default_cell->triangle[l]->springs[0]->index ];
		triangle_tmp->springs[1] = spring[ default_cell->triangle[l]->springs[1]->index ];
		triangle_tmp->springs[2] = spring[ default_cell->triangle[l]->springs[2]->index ];

		triangle_tmp->springs_cytoskeleton[0] = spring[ default_cell->triangle[l]->springs_cytoskeleton[0]->index ];
		triangle_tmp->springs_cytoskeleton[1] = spring[ default_cell->triangle[l]->springs_cytoskeleton[1]->index ];
		triangle_tmp->springs_cytoskeleton[2] = spring[ default_cell->triangle[l]->springs_cytoskeleton[2]->index ];

    
		triangle_tmp->next_triangle_points[0] = mass[ default_cell->triangle[l]->next_triangle_points[0]->index ];
		triangle_tmp->next_triangle_points[1] = mass[ default_cell->triangle[l]->next_triangle_points[1]->index ];
		triangle_tmp->next_triangle_points[2] = mass[ default_cell->triangle[l]->next_triangle_points[2]->index ];
		



		triangle_tmp->nuc = this->mass[0];

		triangle.push_back(triangle_tmp);
	}

  //copy all neigborships ralation
  	for( unsigned int l = 0 ; l < default_cell->triangle.size() ; l++ ){
     // int a =  default_cell->triangle[l]->triangles[0]->index;
    //triangle_tmp = triangle[0];
    triangle[l]->triangles[0] = triangle[ default_cell->triangle[l]->triangles[0]->index];
    triangle[l]->triangles[1] = triangle[ default_cell->triangle[l]->triangles[1]->index];
    triangle[l]->triangles[2] = triangle[ default_cell->triangle[l]->triangles[2]->index];
    
    }

//copy system specifiy matrix
  this->A = (double**) malloc( sizeof(double*)*default_cell->mass.size());
  for( unsigned int i = 0 ; i < default_cell->mass.size() ; i++)
    this->A[i] = (double*) malloc( sizeof(double)*default_cell->mass.size());

  for( unsigned int i = 0 ; i < this->mass.size() ; i++){
    for( unsigned int j = 0 ; j < this->mass.size() ; j++){
      this->A[i][j] = default_cell->A[i][j];
    }
  }


  //copy internal structure (calculate voronoi area)
  for( unsigned int i = 0 ; i < default_cell->mass.size() ; i++ ){
    for( unsigned int j = 0 ; j < default_cell->mass[i]->mvpTriangles.size() ; j++){
      this->mass[i]->mvpPointsL.push_back( default_cell->mass[i]->mvpPointsL[j]);
      this->mass[i]->mvpPointsR.push_back( default_cell->mass[i]->mvpPointsR[j]);
      this->mass[i]->mvpTriangles.push_back( this->triangle[default_cell->mass[i]->mvpTriangles[j]->mIndex] );
      this->mass[i]->mvpSprings.push_back( this->spring[ default_cell->mass[i]->mvpSprings[j]->mIndex ] );

    }
  }


}
void CellTriangulated::projectToCell(CellTriangulated *originalCell, int side){

	
	//calc default length of projected cell
	double l_c = 0.;
	int l_c_count = 0;
	for( unsigned int i = 0 ; i < originalCell->mass_huell.size() ; i++ ){
		if( originalCell->mass_huell[i]->divide_number == side || originalCell->mass_huell[i]->divide_number == 4)
			l_c += dist(originalCell->mass[0]->position,originalCell->mass_huell[i]->position);
		l_c_count++;

	}
	l_c /= l_c_count;




for( unsigned int i = 0 ; i < mass_huell.size() ; i++ ){

	bool found = 0;

		for( unsigned int j = 0 ; j < originalCell->triangle.size() ; j++ ){
			if( (originalCell->triangle[j]->divideCircle == side || originalCell->triangle[j]->divideCircle == 3|| originalCell->triangle[j]->divideCircle == 0|| originalCell->triangle[j]->divideCircle == 4   )
				&& originalCell->triangle[j]->intersectionWithLine(mass[0],mass_huell[i])  ){
				
					found = 1;
				originalCell->triangle[j]->CalcIntersectionWithLine(mass[0],mass_huell[i]);
				break;
			}



		}



		if( found == 0){
			
			bool found_next = 0;

				for( unsigned int j = 0 ; j < originalCell->triangle_divide_ring.size() ; j++ ){
					if( originalCell->triangle_divide_ring[j]->intersectionWithLine(mass[0],mass_huell[i])  ){
					found_next = 1;	
					originalCell->triangle_divide_ring[j]->CalcIntersectionWithLine(mass[0],mass_huell[i]);
					break;
				}
			}
			if( found_next == 0 )
				mass_huell[i]->changeLength(mass[0]->position,l_c);
			
			
			
			


		}

	
	}

        /*
  for( unsigned int i = 0 ; i < this->spring.size() ; i++){

    this->spring[i]->l0 = dist(this->spring[i]->start->position,this->spring[i]->end->position);

  }
  */


/*
  double shrink = 0.90;



  for( unsigned int i = 0 ; i < this->mass_huell.size() ; i++){

    double dir[3];

    dir[0] = this->mass_huell[i]->position[0] - this->mass[0]->position[0];
    dir[1] = this->mass_huell[i]->position[1] - this->mass[0]->position[1];
    dir[2] = this->mass_huell[i]->position[2] - this->mass[0]->position[2];

    this->mass_huell[i]->position[0] = this->mass[0]->position[0]+dir[0]*shrink;
    this->mass_huell[i]->position[1] = this->mass[0]->position[1]+dir[1]*shrink;
    this->mass_huell[i]->position[2] = this->mass[0]->position[2]+dir[2]*shrink;
  
  }

 
    this->volume = this->calcVolume();
    this->v_reference = this->volume;

    */

	this->calcLength();



}

//old routines
void CellTriangulated::insertPoint(Triangle *triangle){

	MassPoint *mass = inCircle(triangle->points[0],triangle->points[1],triangle->points[2]);

	double r = triangle->meanLength(this->mass[0]->position);

	double v[3];
	v[0] = mass->position[0] - position.x;
	v[1] = mass->position[1] - position.y;
	v[2] = mass->position[2] - position.z;

	double l = norm(v);

	mass->position[0] = mass->position[0]/l*r;
	mass->position[1] = mass->position[1]/l*r;
	mass->position[2] = mass->position[2]/l*r;

	this->mass.push_back(mass);

	//insert 3 springs
	for( int i = 0 ; i < 3 ; i++){

		Spring *spring = new Spring();
		spring->start = mass;
		spring->end = triangle->points[i];
		this->spring.push_back(spring);
	}



	//remove triangle
	//add 3 triangles






	//checking triangulation







}
void CellTriangulated::removeSpring(int i){

	for( unsigned int j = 0 ; j < spring[i]->start->neighbourSprings.size(); j++){
		if(  spring[i]->start->neighbourSprings[j] == spring[i] ){

			std::vector<Spring*>::iterator it = spring[i]->start->neighbourSprings.begin()+j;
			spring[i]->start->neighbourSprings.erase(it);
			break;
		}
	}
	for( unsigned int j = 0 ; j < spring[i]->end->neighbourSprings.size(); j++){
		if(  spring[i]->end->neighbourSprings[j] == spring[i] ){

			std::vector<Spring*>::iterator it = spring[i]->end->neighbourSprings.begin()+j;
			spring[i]->end->neighbourSprings.erase(it);
			break;
		}
	}

	std::vector<Spring*>::iterator it = spring.begin()+i;
	spring.erase(it);
//TODO noch alle Komponenten der feder l�schen
}



void CellTriangulated::createInternalStructureOfCell(){

  //add triangles to node
  for( unsigned int i = 0 ; i < this->triangle.size() ; i++ ){

    //add the triangles to the adjacent points
    this->triangle[i]->points[0]->mvpTriangles.push_back(this->triangle[i]);
    this->triangle[i]->points[1]->mvpTriangles.push_back(this->triangle[i]);
    this->triangle[i]->points[2]->mvpTriangles.push_back(this->triangle[i]);

    //add the remaining points of the triangle
    this->triangle[i]->points[0]->mvpPointsL.push_back(1);
    this->triangle[i]->points[0]->mvpPointsR.push_back(2);
    this->triangle[i]->points[1]->mvpPointsL.push_back(2);
    this->triangle[i]->points[1]->mvpPointsR.push_back(0);
    this->triangle[i]->points[2]->mvpPointsL.push_back(0);
    this->triangle[i]->points[2]->mvpPointsR.push_back(1);

  }
  //adds the first triangle again
  for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
    if( this->mass[i]->mvpTriangles.size() > 0 ){
      this->mass[i]->mvpTriangles.push_back( this->mass[i]->mvpTriangles[0] );
      this->mass[i]->mvpPointsL.push_back( this->mass[i]->mvpPointsL[0] );
      this->mass[i]->mvpPointsR.push_back( this->mass[i]->mvpPointsR[0] );
    }
  }

  //order the triangles belongs to a node
  for( unsigned int i = 0 ; i < this->mass.size() ; i++){
    for( int j = 0 ; j < ((int)this->mass[i]->mvpTriangles.size()-1) ; j++ ){//first element is stored twice (second at last position)
      bool found = false;
      for( int k = (j+1) ; ( k < ((int)this->mass[i]->mvpTriangles.size()-1 ) && (found == false)) ; k++ ){
        if( this->mass[i]->mvpTriangles[j]->points[ this->mass[i]->mvpPointsR[j] ] == this->mass[i]->mvpTriangles[k]->points[ this->mass[i]->mvpPointsL[k] ] ){
          Triangle *m = this->mass[i]->mvpTriangles[j+1];
          int left = this->mass[i]->mvpPointsL[j+1];
          int right = this->mass[i]->mvpPointsR[j+1];
          this->mass[i]->mvpTriangles[j+1] = this->mass[i]->mvpTriangles[k];
          this->mass[i]->mvpPointsL[j+1] = this->mass[i]->mvpPointsL[k];
          this->mass[i]->mvpPointsR[j+1] = this->mass[i]->mvpPointsR[k];
          this->mass[i]->mvpTriangles[k] = m;
          this->mass[i]->mvpPointsL[k] = left;
          this->mass[i]->mvpPointsR[k] = right;
          found = true;
        }//if
      }//all next triangles
      if( found == false && (j < this->mass[i]->mvpTriangles.size()-2) )
        fprintf( stderr , "error while sorting triangles %i %i %i \n", i,j, this->mass[i]->mvpTriangles.size());
    }//all triangles of node
  }//all nodes

  //add ordered all neigbouring springs (at surface)
  for( unsigned int i = 0 ; i < this->mass.size() ; i++ ){
    for( unsigned int j = 0 ; j < this->mass[i]->mvpTriangles.size() ; j++ ){
      MassPoint *m0 = this->mass[i]->mvpTriangles[j]->points[ this->mass[i]->mvpPointsR[j] ];
      MassPoint *m1 = this->mass[i];
      //search spring
      bool found = false;
      for( unsigned int k = 0 ; k < this->spring_huell.size() ; k++ ){
        if( (this->spring_huell[k]->start == m0 && this->spring_huell[k]->end == m1) || (this->spring_huell[k]->start == m1 && this->spring_huell[k]->end == m0) ){
          this->mass[i]->mvpSprings.push_back( this->spring_huell[k] );
          found = true;
          break;
        }
      }
      if( found == false ){
        fprintf(stderr,"did not found corresponding spring\n");
      }

    }
  }


}
