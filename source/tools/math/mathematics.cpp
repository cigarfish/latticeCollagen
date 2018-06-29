/**
 * @file	mathematics.cpp
 * @date	11.11.2011
 * @author	Johannes Neitsch
 * @brief	collection of mathematic algorithms
 *
 * extended discription
 *
 */

/** short funtion discription
 * @param	paramter_name	parameter_discription
 * @return	return_discription
 *
 * extendet discription
 *
 */

#include <cmath>

#include <iostream>

#include <vector>

#include "mathematics.h"

#include "../../model/BasicDatatypes/MassPoint.h"
#include "../../model/BasicDatatypes/Spring.h"

using namespace std;

void add_columns(double A[][3], double D[]){
	D[0] = A[0][0] + A[0][1] + A[0][2];
	D[1] = A[1][0] + A[1][1] + A[1][2];
	D[2] = A[2][0] + A[2][1] + A[2][2];
}

void mult_matrix_3times3_component(double A[][3], double B[][3]){
	A[0][0] *= B[0][0];	A[0][1] *= B[0][1];	A[0][2] *= B[0][2];
	A[1][0] *= B[1][0];	A[1][1] *= B[1][1];	A[1][2] *= B[1][2];
	A[2][0] *= B[2][0];	A[2][1] *= B[2][1];	A[2][2] *= B[2][2];
}

void mult_matrix_3times3_vector(double A[][3], double D[], double back[]){
	back[0] = A[0][0] *D[0] + A[0][1] *D[1] + A[0][2] *D[2];
	back[1] = A[1][0] *D[0] + A[1][1] *D[1] + A[1][2] *D[2];
	back[2] = A[2][0] *D[0] + A[2][1] *D[1] + A[2][2] *D[2];
}


void half_matrix_3times3(double A[][3], double B[][3]){
	B[0][0] = A[0][0]*0.5;	B[0][1] = A[0][1]*0.5;	B[0][2] = A[0][2]*0.5;
	B[1][0] = A[1][0]*0.5;	B[1][1] = A[1][1]*0.5;	B[1][2] = A[1][2]*0.5;
	B[2][0] = A[2][0]*0.5;	B[2][1] = A[2][1]*0.5;	B[2][2] = A[2][2]*0.5;
}


double det(double A[][3]){
return ( A[0][0]*(A[2][2]*A[1][1]-A[2][1]*A[1][2]) - A[1][0]*(A[2][2]*A[0][1]-A[2][1]*A[0][2]) + A[2][0]*(A[1][2]*A[0][1]-A[1][1]*A[0][2]) );
}

double det(double A0[3],double A1[3],double A2[3]){
return ( A0[0]*(A2[2]*A1[1]-A2[1]*A1[2]) - A1[0]*(A2[2]*A0[1]-A2[1]*A0[2]) + A2[0]*(A1[2]*A0[1]-A1[1]*A0[2]) );
}

void inv_matrix_3times3(double A[][3], double A_inv[][3]){



	double det_A_inv = 1./( A[0][0]*(A[2][2]*A[1][1]-A[2][1]*A[1][2]) - A[1][0]*(A[2][2]*A[0][1]-A[2][1]*A[0][2]) + A[2][0]*(A[1][2]*A[0][1]-A[1][1]*A[0][2]) );

	A_inv[0][0] = det_A_inv *(  A[2][2]*A[1][1]-A[2][1]*A[1][2] );
	A_inv[0][1] = det_A_inv *( -A[2][2]*A[0][1]+A[2][1]*A[0][2] );
	A_inv[0][2] = det_A_inv *(  A[1][2]*A[0][1]-A[1][1]*A[0][2] );

	A_inv[1][0] = det_A_inv *( -A[2][2]*A[1][0]+A[2][0]*A[1][2] );
	A_inv[1][1] = det_A_inv *(  A[2][2]*A[0][0]-A[2][0]*A[0][2] );
	A_inv[1][2] = det_A_inv *( -A[1][2]*A[0][0]+A[1][0]*A[0][2] );

	A_inv[2][0] = det_A_inv *(  A[2][1]*A[1][0]-A[2][0]*A[1][1] );
	A_inv[2][1] = det_A_inv *( -A[2][1]*A[0][0]+A[2][0]*A[0][1] );
	A_inv[2][2] = det_A_inv *(  A[1][1]*A[0][0]-A[1][0]*A[0][1] );

}

double outCircle(MassPoint* m0, MassPoint*m1, MassPoint* m2, MassPoint* m_return){

	double A[3][3];

	A[0][0] = m1->position[0] - m0->position[0];
	A[0][1] = m1->position[1] - m0->position[1];
	A[0][2] = m1->position[2] - m0->position[2];

	A[1][0] = m2->position[0] - m0->position[0];
	A[1][1] = m2->position[1] - m0->position[1];
	A[1][2] = m2->position[2] - m0->position[2];

	A[2][0] = m2->position[0] - m1->position[0];
	A[2][1] = m2->position[1] - m1->position[1];
	A[2][2] = m2->position[2] - m1->position[2];


	double A_inv[3][3];
	inv_matrix_3times3( A, A_inv);


	double D[3];

	double B[3][3];
	half_matrix_3times3( A,  B);

	B[0][0] += m0->position[0];
	B[0][1] += m0->position[1];
	B[0][2] += m0->position[2];

	B[1][0] += m1->position[0];
	B[1][1] += m1->position[1];
	B[1][2] += m1->position[2];

	B[2][0] += m2->position[0];
	B[2][1] += m2->position[1];
	B[2][2] += m2->position[2];

	mult_matrix_3times3_component( A, B);//store in A
	add_columns( A ,D);

	double m[3];
	mult_matrix_3times3_vector(A_inv, D,m);

	cout << "Mitpkt" << m[0] << "\t" << m[1] << "\t" << m[2] << endl;
	double tmp[3] = {m[0]-m0->position[0],m[1]-m0->position[1],m[2]-m0->position[2]};


	return norm(tmp);


}




MassPoint* inCircle(MassPoint* m0, MassPoint*m1, MassPoint* m2){

	MassPoint *incirclePoint = new MassPoint();


	double v01[3];
	double v12[3];

	v01[0] = m1->position[0] - m0->position[0];
	v01[1] = m1->position[1] - m0->position[1];
	v01[2] = m1->position[2] - m0->position[2];

	v12[0] = m2->position[0] - m1->position[0];
	v12[1] = m2->position[1] - m1->position[1];
	v12[2] = m2->position[2] - m1->position[2];


	double v1[3];
	double v2[2];

	v1[0] = - v12[0] - v01[0]/2.;
	v1[1] = - v12[1] - v01[1]/2.;
	v1[2] = - v12[2] - v01[2]/2.;

	v2[0] = v01[0] + v12[0]/2.;
	v2[1] = v01[1] + v12[1]/2.;

	//intersection between v1 and v2

	double k = (v1[1]*(m2->position[0]-m0->position[0])-v1[0]*(m2->position[1]-m0->position[1]))/(v2[0]*v1[1]-v1[0]*v2[1]);

	incirclePoint->position[0] = m2->position[0] + k*v1[0];
	incirclePoint->position[1] = m2->position[1] + k*v1[1];
	incirclePoint->position[2] = m2->position[2] + k*v1[2];



	return incirclePoint;
}


double norm(double v[3]){
	return std::sqrt( v[0]*v[0]+v[1]*v[1]+v[2]*v[2] );


}
double norm(double v[3],double m0[3]){

	double x = v[0]-m0[0];
	double y = v[1]-m0[1];
	double z = v[2]-m0[2];
	return std::sqrt( x*x+y*y+z*z );

}



double normvector(double v[3]){
	double tmp = norm(v);
	if (tmp == 0) return 0;
	v[0] /= tmp;
	v[1] /= tmp;
	v[2] /= tmp;
	return tmp;
}



double dotmult(double v0[3], double v1[3]){
	return v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2];
}



void crossProduct( double v0[3], double v1[3], double v_r[3]){
	v_r[0] = v0[1]*v1[2] - v0[2]*v1[1];
	v_r[1] = v0[2]*v1[0] - v0[0]*v1[2];
	v_r[2] = v0[0]*v1[1] - v0[1]*v1[0];
}

void crossProduct( double v0[3], double v1[3], double v2[3], double v_r[3]){

  double v_0[3];
  double v_1[3];

  v_0[0] = v1[0] - v0[0];
  v_0[1] = v1[1] - v0[1];
  v_0[2] = v1[2] - v0[2];

  v_1[0] = v2[0] - v0[0];
  v_1[1] = v2[1] - v0[1];
  v_1[2] = v2[2] - v0[2];

  v_r[0]  = v_0[1]*v_1[2] - v_0[2]*v_1[1];
  v_r[1]  = v_0[2]*v_1[0] - v_0[0]*v_1[2];
  v_r[2]  = v_0[0]*v_1[1] - v_0[1]*v_1[0];

  normvector(v_r);

}

/*
double calcRadius(CellTriangulated *zelle){

	double tmp = 0.;
	for( unsigned int i = 0 ; i < zelle->mass.size() ; i++ ){
		tmp +=norm(zelle->mass[i]->position);

	}


	return tmp/zelle->mass.size();
}
*/
double dist(double v0[3], double v1[3]){
  
  double v[3] = { v0[0]-v1[0] , v0[1]-v1[1] , v0[2]-v1[2]};
  
  return norm(v);

}
void dist(double v0[3], double v1[3], double v_return[3]){

	v_return[0] = v1[0]-v0[0];
	v_return[1] = v1[1]-v0[1];
	v_return[2] = v1[2]-v0[2];

}

/*
double dist(Vector3f v0, Vector3f v1){
	double v[3] = { v0.x - v1.x , v0.y-v1.y , v0.z-v1.z};
	return norm(v);
}
*/

double max_points(std::vector <MassPoint*> mass){

	double max = 0;
	for( unsigned int i = 0 ; i < mass.size() ; i++ ){
		for( unsigned int j = i+1 ; j < mass.size() ; j++){
			double tmp = dist(mass[i]->position,mass[j]->position);

			if( tmp > max)
				max = tmp;


		}

	}
return max;
}

void mean(double a[], double b[], double c[], double r[]){
  
  for( int i = 0 ; i < 3 ; i++){
    r[i] = a[i] + b[i] + c[i];
    r[i] /= 3.;
  }

}


double rotateCylinder(double tmp_rot[3], Spring *s){
	double vx = s->end->position[0]-s->start->position[0];// x2-x1;
	double vy = s->end->position[1]-s->start->position[1];// y2-y1;
	double vz = s->end->position[2]-s->start->position[2];// z2-z1;

//handle the degenerate case of z1 == z2 with an approximation
if(vz == 0)
    vz = .0001;

double v = sqrt( vx*vx + vy*vy + vz*vz );
double ax = 57.2957795*acos( vz/v );
if ( vz < 0.0 )
    ax = -ax;
double rx = -vy*vz;
double ry = vx*vz;

//draw the cylinder body
	//glTranslatef( x1,y1,z1 );
	
	//glRotatef(ax, rx, ry, 0.0);
	tmp_rot[0] = ax;
	tmp_rot[1] = rx;
	tmp_rot[2] = ry;
	return v;

}
double rotateCylinder(double tmp_rot[3],double a[3],double b[3], double scale){
	double vx = (b[0]-a[0])/scale;// x2-x1;
	double vy = (b[1]-a[1])/scale;// y2-y1;
	double vz = (b[2]-a[2])/scale;// z2-z1;

//handle the degenerate case of z1 == z2 with an approximation
if(vz == 0)
    vz = .0001;

double v = sqrt( vx*vx + vy*vy + vz*vz );
double ax = 57.2957795*acos( vz/v );
if ( vz < 0.0 )
    ax = -ax;
double rx = -vy*vz;
double ry = vx*vz;

//draw the cylinder body
	//glTranslatef( x1,y1,z1 );
	
	//glRotatef(ax, rx, ry, 0.0);
	tmp_rot[0] = ax;
	tmp_rot[1] = rx;
	tmp_rot[2] = ry;
	return v;

}

double rotateCylinder(double tmp_rot[3],float a[3],float b[3], double scale){
	double vx = (b[0]-a[0])/scale;// x2-x1;
	double vy = (b[1]-a[1])/scale;// y2-y1;
	double vz = (b[2]-a[2])/scale;// z2-z1;

//handle the degenerate case of z1 == z2 with an approximation
if(vz == 0)
    vz = .0001;

double v = sqrt( vx*vx + vy*vy + vz*vz );
double ax = 57.2957795*acos( vz/v );
if ( vz < 0.0 )
    ax = -ax;
double rx = -vy*vz;
double ry = vx*vz;

//draw the cylinder body
	//glTranslatef( x1,y1,z1 );
	
	//glRotatef(ax, rx, ry, 0.0);
	tmp_rot[0] = ax;
	tmp_rot[1] = rx;
	tmp_rot[2] = ry;
	return v;

}



double rotateCylinder(double tmp_rot[3],double a[3],double b[3]){
	double vx = b[0]-a[0];// x2-x1;
	double vy = b[1]-a[1];// y2-y1;
	double vz = b[2]-a[2];// z2-z1;

//handle the degenerate case of z1 == z2 with an approximation
if(vz == 0)
    vz = .0001;

double v = sqrt( vx*vx + vy*vy + vz*vz );
double ax = 57.2957795*acos( vz/v );
if ( vz < 0.0 )
    ax = -ax;
double rx = -vy*vz;
double ry = vx*vz;

//draw the cylinder body
	//glTranslatef( x1,y1,z1 );
	
	//glRotatef(ax, rx, ry, 0.0);
	tmp_rot[0] = ax;
	tmp_rot[1] = rx;
	tmp_rot[2] = ry;
	return v;

}

bool intersect_triangle_section_area(double p0[3],double p1[3],double p2[3],double start[3],double end[3],double intersection[3]){

		
		
		double u[3],v[3];
		u[0] = p2[0] - p0[0];
		u[1] = p2[1] - p0[1];
		u[2] = p2[2] - p0[2];
		v[0] = p1[0] - p0[0];
		v[1] = p1[1] - p0[1];
		v[2] = p1[2] - p0[2];

		double direction[3];
		direction[0] = end[0] - start[0];
		direction[1] = end[1] - start[1];
		direction[2] = end[2] - start[2];


		double normalvector[3];
		crossProduct(v,u,normalvector);
		//check parallel of triangle and line
		if( dotmult(direction,normalvector) == 0.) return 0;//check angle normalvectors

    double tmp_v[3];
    tmp_v[0] = p0[0]-start[0];
    tmp_v[1] = p0[1]-start[1];
    tmp_v[2] = p0[2]-start[2];

    double t = dotmult(tmp_v,normalvector) / dotmult(direction,normalvector);

    if( t < 0. ) return 0;

    double tmp_interaction[3];
    tmp_interaction[0] = start[0] +t*direction[0];
    tmp_interaction[1] = start[1] +t*direction[1];
    tmp_interaction[2] = start[2] +t*direction[2];

    double area_all = area(p0,p1,p2);

    double a0 = area(tmp_interaction,p0,p1)/area_all;
    double a1 = area(tmp_interaction,p1,p2)/area_all;
    double a2 = area(tmp_interaction,p2,p0)/area_all;
    
    if( (a0+a1+a2-0.001) < 1){
      intersection[0] = tmp_interaction[0];
      intersection[1] = tmp_interaction[1];
      intersection[2] = tmp_interaction[2];
    }



	return false;

}

double area(double p0[3],double p1[2], double p2[3]){

  double tmp_ba[3];
  tmp_ba[0] = p1[0] - p0[0];
  tmp_ba[1] = p1[1] - p0[1];
  tmp_ba[2] = p1[2] - p0[2];

  double tmp_ca[3];
  tmp_ca[0] = p2[0]-p0[0];
  tmp_ca[1] = p2[1]-p0[1];
  tmp_ca[2] = p2[2]-p0[2];
  
  double tmp_cross[3];

  crossProduct(tmp_ba,tmp_ca,tmp_cross);
  
  return 0.5* norm(tmp_cross);

}

bool intersect_triangle_section(double p0[3],double p1[3],double p2[3],double start[3],double end[3],double intersection[3]){

		
		
		double u[3],v[3];
		u[0] = p2[0] - p0[0];
		u[1] = p2[1] - p0[1];
		u[2] = p2[2] - p0[2];
		v[0] = p1[0] - p0[0];
		v[1] = p1[1] - p0[1];
		v[2] = p1[2] - p0[2];

		double direction[3];
		direction[0] = end[0];// - start[0];
		direction[1] = end[1];// - start[1];
		direction[2] = end[2];// - start[2];
    
    

		double length = norm(direction);



		double normalvector[3];
		crossProduct(u,v,normalvector);
		//check parallel of triangle and line
		if( dotmult(direction,normalvector) == 0.) return 0;//check angle normalvectors
    
		double w[3];
		w[0] = start[0] - p0[0];
		w[1] = start[1] - p0[1];
		w[2] = start[2] - p0[2];

		double cross_d_v[3];
		crossProduct(direction,v,cross_d_v);
		
		double cross_w_u[3];
		crossProduct(w,u,cross_w_u);

		double invers_d_v_u = 1./dotmult(cross_d_v,u);

		double trs[3];
			
		trs[0] = invers_d_v_u * dotmult(cross_w_u,v);
		if( trs[0] >= 0. )	{

		trs[1] = invers_d_v_u * dotmult(cross_d_v,w);
		if( trs[1] >= 0. && trs[1] <= 1. ){
			
			
		trs[2] = invers_d_v_u * dotmult(cross_w_u,direction);


		if ( trs[2] >= 0. && trs[2] <= 1.)	{
		
		double sum_r_s = trs[1]+trs[2];
		if( sum_r_s <= 1. )	{

	
		
		if( trs[0] <= length){
			
			intersection[0] = start[0] + trs[0] * direction[0];
			intersection[1] = start[1] + trs[0] * direction[1];
			intersection[2] = start[2] + trs[0] * direction[2];
			return true;


		}}}
		}//check length
	}//intersection

	return false;

}

double det_4( double *A[4]){
	
	//calculates a det of a special formed 4x4 matrix
	
	// A00 A01 A02 1
	// A10 A11 A12 1
	// A20 A21 A22 1
	// A30 A31 A32 1
	
	return	A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + 
			A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0] - 
			A[0][0]*A[1][1]*A[3][2] + A[0][0]*A[1][2]*A[3][1] + A[0][1]*A[1][0]*A[3][2] - 
			A[0][1]*A[1][2]*A[3][0] - A[0][2]*A[1][0]*A[3][1] + A[0][2]*A[1][1]*A[3][0] + 
			A[0][0]*A[2][1]*A[3][2] - A[0][0]*A[2][2]*A[3][1] - A[0][1]*A[2][0]*A[3][2] + 
			A[0][1]*A[2][2]*A[3][0] + A[0][2]*A[2][0]*A[3][1] - A[0][2]*A[2][1]*A[3][0] - 
			A[1][0]*A[2][1]*A[3][2] + A[1][0]*A[2][2]*A[3][1] + A[1][1]*A[2][0]*A[3][2] - 
			A[1][1]*A[2][2]*A[3][0] - A[1][2]*A[2][0]*A[3][1] + A[1][2]*A[2][1]*A[3][0];
	 
}



bool point_inside_tetrahedron(double p0[3],double p1[3],double p2[3],double p[3],double point[3]){
	
//	double A[4][3];
//	A[0] = p0;
	
	double *A[4];
	
	
	A[0] = p0;
	A[1] = p1;
	A[2] = p2;
	A[3] = p;
	
	double vz = det_4(A);
	
	if( vz < 0 ){
		vz = -1.;
	}else{
		if( vz > 0 )
			vz = 1.;
	}
	
	A[0] = point;
	
	double vz_tmp = det_4(A);

	if( vz_tmp < 0 ){
		vz_tmp = -1.;
	}else{
		if( vz_tmp > 0 )
			vz_tmp = 1.;
	}

	if( abs(vz - vz_tmp) > 1.) return 0;

	A[0] = p0;
	A[1] = point;

	vz_tmp = det_4(A);
	if( vz_tmp < 0 ){
		vz_tmp = -1.;
	}else{
		if( vz_tmp > 0 )
			vz_tmp = 1.;
	}

	if( abs(vz - vz_tmp) > 1.) return 0;

	A[1] = p1;
	A[2] = point;

	vz_tmp = det_4(A);
	if( vz_tmp < 0 ){
		vz_tmp = -1.;
	}else{
		if( vz_tmp > 0 )
			vz_tmp = 1.;
	}

	if( abs(vz - vz_tmp) > 1.) return 0;

	A[2] = p2;
	A[3] = point;

	
	vz_tmp = det_4(A);
	if( vz_tmp < 0 ){
		vz_tmp = -1.;
	}else{
		if( vz_tmp > 0 )
			vz_tmp = 1.;
	}

	if( abs(vz - vz_tmp) > 1.) return 0;



	return 1;
}

//return -1 if edge point of section is closer, else distance
//dir = end- start
double dist_point_section(double start[3],double end[3], double p[3], double intersection[3]){

  double dir[3];
  dir[0] = end[0] - start[0];
  dir[1] = end[1] - start[1];
  dir[2] = end[2] - start[2];

  double d = dir[0]*p[0] + dir[1]*p[1] + dir[2]*p[2];
   
//  if( d < -0.001 || d > 1.001)
 //   return -1;

  double rs = dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2];

  double numbers = d - dir[0]*start[0] - dir[1]*start[1] - dir[2]*start[2];

 
  
  double r = numbers/rs;

  intersection[0] = start[0] + r * dir[0];
  intersection[1] = start[1] + r * dir[1];
  intersection[2] = start[2] + r * dir[2];

  double l = dist(start,end);

  double l1 = dist(start,intersection);

  if( l1 > l)
    return -1;

  double l2 = dist( end, intersection);
    if( l2 > l)
      return -1;

  return r;
}


//gives the angle of two normalized vectors
double angle(double v1[3], double v2[3]){
  return acos(dotmult(v1,v2));
}
