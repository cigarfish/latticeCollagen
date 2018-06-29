///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Triangle.cpp                                                         //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-05-31 22:00:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "MassPoint.h"
#include "Spring.h"
#include "Triangle.h"

#include "Vector.h"
#include "../../tools/math/mathematics.h"
#include <cmath>

Triangle::Triangle() {

	divideCircle = 0;

	r = 255;
	g = 255;
	b = 0;
	transparency = 0.0;
	interaction = 0;
	mArea = 0.;
	Volume = 0.;

}
void Triangle::setNormalVector(){
	
	// set the normalvector of a triangle
	crossProduct(points[0]->position,points[1]->position,points[2]->position,mNormalvector);
	normvector(mNormalvector);

}


void Triangle::calcMean(double mean[3]){
	mean[0] = (points[0]->position[0] + points[1]->position[0] + points[2]->position[0])/3.;
	mean[1] = (points[1]->position[1] + points[1]->position[1] + points[2]->position[1])/3.;
	mean[2] = (points[2]->position[2] + points[1]->position[2] + points[2]->position[2])/3.;

}


bool Triangle::intersectionWithLine(MassPoint *m0, MassPoint *m1){

	double d[3];
	dist( m0->position,m1->position,d );
	double e2[3];
	dist( points[0]->position,points[2]->position,e2 );
	double e1[3];
	dist( points[0]->position,points[1]->position,e1 );
	double t[3];
	dist( points[0]->position,m0->position,t );
	
	double p[3] = {0,0,0};
	crossProduct( d, e2,p);
	double D = dotmult(p,e1);

	if( D == 0)
		return 0;
	double alpha = dotmult(p,t)/D;

	if( alpha < 0 || alpha > 1) 
		return 0;

	double q[3] = {0,0,0};
	crossProduct( t, e1,q);
	double beta = dotmult(q,d)/D;

	if( beta < 0 || (alpha+beta)>1 )
		return 0;
	
	double lambda = dotmult(q,e2)/D;

	if( lambda <= 0)
		return 0;

	return 1;

}

void Triangle::CalcIntersectionWithLine(MassPoint *m, MassPoint *m_new){
	
	
	double d[3];
	dist(m->position,m_new->position,d);
	
	double e2[3];
	dist(points[0]->position,points[2]->position,e2);
	double e1[3];
	dist(points[0]->position,points[1]->position,e1);
	double t[3];
	dist(points[0]->position,m->position,t);
	
	double p[3] = {0,0,0};
	crossProduct( d,e2,p );
	double D = dotmult(p,e1);

	double q[3] = {0,0,0};
	crossProduct( t,e1,q );
	
	double lambda = dotmult(q,e2)/D;

	m_new->position[0] = m->position[0]+lambda*d[0];
	m_new->position[1] = m->position[1]+lambda*d[1];
	m_new->position[2] = m->position[2]+lambda*d[2];
	
}


void Triangle::checkAndChangeOrientation(){

	double v0[3];
	double v1[3];
	
	dist(points[0]->position, points[1]->position, v0);
	dist(points[0]->position, points[2]->position, v1);

	double v[3];
	crossProduct(v0,v1,v);


	v0[0] = points[1]->position[0] + v[0];
	v0[1] = points[1]->position[1] + v[1];
	v0[2] = points[1]->position[2] + v[2];

	v1[0] = points[1]->position[0] - v[0];
	v1[1] = points[1]->position[1] - v[1];
	v1[2] = points[1]->position[2] - v[2];


	double nv0 = norm(v0);
	double nv1 = norm(v1);


	if( nv0 < nv1){
		MassPoint *temp = this->points[2];
		this->points[2] = this->points[1];
		this->points[1] = temp;
		
		//TODO
		//~temp;
	}

}

void Triangle::setArea(double p){

	//calculate side-length of the triangle
	double v01[3], v02[3], v12[3];

	dist(points[0]->position,points[1]->position,v01);
	dist(points[0]->position,points[2]->position,v01);
	dist(points[1]->position,points[2]->position,v01);

	double nv01 = norm(v01);
	double nv02 = norm(v02);
	double nv12 = norm(v12);

	double s = (nv01+nv02+nv12)*0.5;

	//caluclate area
	mArea  = p*sqrt(s*(s-nv01)*(s-nv02)*(s-nv12))/3;

	crossProduct(v01,v02,mNormalvector);

	double tmp = mArea/norm(mNormalvector);

	//1/3 of normalvector
	mNormalvector[0] *= tmp;
	mNormalvector[1] *= tmp;
	mNormalvector[2] *= tmp;
	
}

void Triangle::calcVolume(double top[3]){

	double sp0[3];
	double sp1[3];
	double sp2[3];

	dist( top, points[0]->position, sp0);
	dist( top, points[1]->position, sp1);
	dist( top, points[2]->position, sp2);

	double crossprod[3];
	crossProduct( sp0, sp1, crossprod);
	double dotm = fabs(dotmult(sp2,crossprod));

	Volume = dotm/6.;
}

double Triangle::meanLength(double center[3]){
		
	double v[3];
	v[0] = points[0]->position[0] - center[0];
	v[1] = points[0]->position[1] - center[1];
	v[2] = points[0]->position[2] - center[2];

	double l1 = norm(v);

	v[0] = points[1]->position[0] - center[0];
	v[1] = points[1]->position[1] - center[1];
	v[2] = points[1]->position[2] - center[2];

	double l2 = norm(v);

	v[0] = points[2]->position[0] - center[0];
	v[1] = points[2]->position[1] - center[1];
	v[2] = points[2]->position[2] - center[2];

	double l3 = norm(v);
		
	return (l1+l2+l3)/3.;
}

void Triangle::copyfrom( Triangle *triangle_basic){

  this->mIndex = triangle_basic->mIndex;

	//Data
	mNormalvector[0] = triangle_basic->mNormalvector[0];
	mNormalvector[1] = triangle_basic->mNormalvector[1];
	mNormalvector[2] = triangle_basic->mNormalvector[2];

	mArea = triangle_basic->mArea;
  this->index = triangle_basic->index;

	divideCircle = triangle_basic->divideCircle;

	//visualisation
    r = triangle_basic->r;
    g = triangle_basic->g;
	b = triangle_basic->b;
	transparency = triangle_basic->transparency;

}

bool Triangle::check_intersection_triangle_section2(double p0[3],double v[3], double u[3], double direction[3], double start[3], double end[3], double length, double intersection[3]){

		double normalvector[3];
		crossProduct(u,v,normalvector);
		//check parallel of triangle and line
		if( dotmult(direction,normalvector) == 0) return 0;//check angle normalvectors


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
		if( trs[0] >= 0 )	{

		trs[1] = invers_d_v_u * dotmult(cross_d_v,w);
		if( trs[1] > 0 && trs[1] < 1 ){
			
			
		trs[2] = invers_d_v_u * dotmult(cross_w_u,direction);


		if ( trs[2] > 0 && trs[2] < 1)	{
		
		double sum_r_s = trs[1]+trs[2];
		if( sum_r_s <= 1 )	{

	
		
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


bool Triangle::check_intersection_triangle_section( double direction[3], double start[3], double end[3], double length, double intersection[3]){
		//check parallel of triangle and line
		if( dotmult(direction,this->mNormalvector) == 0) return 0;//check angle normalvectors


		double u[3],v[3],w[3];
		
			v[0] = this->points[1]->position[0] - this->points[0]->position[0];
			v[1] = this->points[1]->position[1] - this->points[0]->position[1];
			v[2] = this->points[1]->position[2] - this->points[0]->position[2];

			u[0] = this->points[2]->position[0] - this->points[0]->position[0];
			u[1] = this->points[2]->position[1] - this->points[0]->position[1];
			u[2] = this->points[2]->position[2] - this->points[0]->position[2];



			w[0] = start[0] - this->points[0]->position[0];
			w[1] = start[1] - this->points[0]->position[1];
			w[2] = start[2] - this->points[0]->position[2];


			
			double cross_d_v[3];
			crossProduct(direction,v,cross_d_v);
			double cross_w_u[3];
			crossProduct(w,u,cross_w_u);

			double invers_d_v_u = 1./dotmult(cross_d_v,u);

			double trs[3];
			
			trs[0] = invers_d_v_u * dotmult(cross_w_u,v);
			if( trs[0] >= 0 )	{

			trs[1] = invers_d_v_u * dotmult(cross_d_v,w);
			if( trs[1] > 0 && trs[1] < 1 ){
			
			
			trs[2] = invers_d_v_u * dotmult(cross_w_u,direction);


			

		
			if ( trs[2] > 0 && trs[2] < 1)	{
		
			double sum_r_s = trs[1]+trs[2];
			if( sum_r_s <= 1 )	{

	
		
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

bool Triangle::check_point_inside_tetrahedron(double point[3]){

	double p0i[3];
			
	p0i[0] = point[0] - this->nuc->position[0];
	p0i[1] = point[1] - this->nuc->position[1];
	p0i[2] = point[2] - this->nuc->position[2];

	return this->check_point_inside_tetrahedron(point, p0i);

}
bool Triangle::check_point_inside_tetrahedron(double point[3], double p0i[3]){
	
	double norm_invador = norm(point,this->nuc->position);

	double m00[3] = {this->nuc->position[0] - this->points[0]->position[0],
					 this->nuc->position[1] - this->points[0]->position[1],
					 this->nuc->position[2] - this->points[0]->position[2]};
	
	double m01[3] = {this->nuc->position[0] - this->points[1]->position[0],
					 this->nuc->position[1] - this->points[1]->position[1],
					 this->nuc->position[2] - this->points[1]->position[2]};

	double m02[3] = {this->nuc->position[0] - this->points[2]->position[0],
					 this->nuc->position[1] - this->points[2]->position[1],
					 this->nuc->position[2] - this->points[2]->position[2]};
	
	if( !(norm(m02) < norm_invador && norm(m01) < norm_invador && norm(m00) < norm_invador)){
		
		double pm0i[3];
			
		pm0i[0] = point[0] - this->points[0]->position[0];
		pm0i[1] = point[1] - this->points[0]->position[1];
		pm0i[2] = point[2] - this->points[0]->position[2];
			
		//check if point is on right side of triangle
		//first triangle at surface of cell
		if( (dotmult(this->mNormalvector, pm0i) < 0 )){

		//all 3 other sides


		double normvector_m0_2_1[3];
		double normvector_m0_1_0[3];
		double normvector_m0_0_2[3];

		crossProduct(this->nuc->position,this->points[0]->position,this->points[2]->position,normvector_m0_0_2);
		if( dotmult( normvector_m0_0_2, p0i) < 0 ){
				//break;
	
		crossProduct(this->nuc->position,this->points[2]->position,this->points[1]->position,normvector_m0_2_1);

		if( (dotmult( normvector_m0_2_1, p0i) <= 0 )){
			//break;

		crossProduct(this->nuc->position,this->points[1]->position,this->points[0]->position,normvector_m0_1_0);

		if( (dotmult( normvector_m0_1_0, p0i) <= 0 )){
	
			return true;

		}}}}}
		
	return false;
}

int Triangle::missing_point( Triangle *t){
  
  for( int i = 0 ; i < 3 ; i++){
    int score = 0;
    for( int j = 0 ; j < 3 ; j++ ){
      if( this->points[i] == t->points[j])
        score++;
    }
    if( score == 0)
      return i;
  }
  return -1;
}


int Triangle::equal(Triangle *triangle){
	//gives the amount of same vertex of two triangles
	int score = 0;
	
	for( int i = 0 ; i < 3 ; i++ ){
		for( int j = 0 ; j < 3 ; j++ ){
			if( this->points[i] == triangle->points[j] )
				score++;
		}
	}
	return score;
}

bool Triangle::equal(MassPoint *p1, MassPoint *p2){
	//gives the amount of same vertex of two triangles
	int score = 0;
	
	for( int i = 0 ; i < 3 ; i++ ){
    if( this->points[i] == p1 || this->points[i] == p2){
      score++;
      if( score == 2)
        return 1;
    }
 }
    return 0;
}



void Triangle::difference(Triangle *triangle){
	//two triangles has one side together, gives the index of the missing vertex
	
	int pos = 0;
	int v[2];
	int next_m;

	for( int i = 0 ; i < 3 ; i++ ){
		bool found = 0;
		for(int j = 0 ; j < 3 ; j++ ){
			if( this->points[j] == triangle->points[i] ){
				v[pos] = j;
				pos++;
				found = 1;
			}
		}
		if( found == 0){
			next_m = i;
		}
	}

	if( (v[0] == 0 && v[1] == 1) || (v[1] == 0 && v[0] == 1) ){
		this->next_triangle_points[0] = triangle->points[next_m];
		return;
	}
	if( (v[0] == 1 && v[1] == 2) || (v[1] == 1 && v[0] == 2) ){
		this->next_triangle_points[1] = triangle->points[next_m];
		return;
	}
	if( (v[0] == 0 && v[1] == 2) || (v[1] == 0 && v[0] == 2) ){
		this->next_triangle_points[2] = triangle->points[next_m];
		return;
	}

	
	/*
	this->next_triangle[0] = triangle->points[0];
	this->next_triangle[1] = triangle->points[1];
	this->next_triangle[2] = triangle->points[2];
	*/


}



double check_intersection_tetrahedral_naive(Triangle *triangle1,Triangle *triangle2){


	vector <double> pos_x;
	vector <double> pos_y;
	vector <double> pos_z;

	
	double intersection[3];
	
	//check edges of hull
	if( intersect_triangle_section(
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->points[0]->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->points[0]->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	
	//ceck intersection of sides triangles and springs

	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[2]->position,
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle1->points[0]->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[2]->position,
		triangle1->points[0]->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->points[0]->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[2]->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	//second cell
		if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[2]->position,
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle2->points[0]->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[2]->position,
		triangle2->points[0]->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->points[0]->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[2]->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}

	

	//ceck intersection fo spring_cytskeleton

	if( intersect_triangle_section(
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->nuc->position,
		triangle2->points[0]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->nuc->position,
		triangle2->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->nuc->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->nuc->position,
		triangle1->points[0]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->nuc->position,
		triangle1->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->nuc->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}



	
	
	//ceck intersection spring_cyto with side sides of tetraheron
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle1->nuc->position,
		triangle1->points[0]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle1->nuc->position,
		triangle1->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[1]->position,
		triangle1->nuc->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[2]->position,
		triangle1->nuc->position,
		triangle1->points[0]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[2]->position,
		triangle1->nuc->position,
		triangle1->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[0]->position,
		triangle2->points[2]->position,
		triangle1->nuc->position,
		triangle1->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->nuc->position,
		triangle1->points[0]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->nuc->position,
		triangle1->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle2->nuc->position,
		triangle2->points[1]->position,
		triangle2->points[2]->position,
		triangle1->nuc->position,
		triangle1->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	//other cell
		if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle2->nuc->position,
		triangle2->points[0]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle2->nuc->position,
		triangle2->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[1]->position,
		triangle2->nuc->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[2]->position,
		triangle2->nuc->position,
		triangle2->points[0]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[2]->position,
		triangle2->nuc->position,
		triangle2->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[0]->position,
		triangle1->points[2]->position,
		triangle2->nuc->position,
		triangle2->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->nuc->position,
		triangle2->points[0]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->nuc->position,
		triangle2->points[1]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}
	if( intersect_triangle_section(
		triangle1->nuc->position,
		triangle1->points[1]->position,
		triangle1->points[2]->position,
		triangle2->nuc->position,
		triangle2->points[2]->position,
		intersection))
	{
		pos_x.push_back(intersection[0]);
		pos_y.push_back(intersection[1]);
		pos_z.push_back(intersection[2]);
	}

	



	

	if( point_inside_tetrahedron(
			triangle1->points[0]->position,
			triangle1->points[1]->position,
			triangle1->points[2]->position,
			triangle1->nuc->position,
			triangle2->points[0]->position)){
		pos_x.push_back(triangle2->points[0]->position[0]);
		pos_y.push_back(triangle2->points[0]->position[1]);
		pos_z.push_back(triangle2->points[0]->position[2]);
	}
	
	if( point_inside_tetrahedron(
			triangle1->points[0]->position,
			triangle1->points[1]->position,
			triangle1->points[2]->position,
			triangle1->nuc->position,
			triangle2->points[1]->position)){
		pos_x.push_back(triangle2->points[1]->position[0]);
		pos_y.push_back(triangle2->points[1]->position[1]);
		pos_z.push_back(triangle2->points[1]->position[2]);
	}
	
	if( point_inside_tetrahedron(
			triangle1->points[0]->position,
			triangle1->points[1]->position,
			triangle1->points[2]->position,
			triangle1->nuc->position,
			triangle2->points[2]->position)){
		pos_x.push_back(triangle2->points[2]->position[0]);
		pos_y.push_back(triangle2->points[2]->position[1]);
		pos_z.push_back(triangle2->points[2]->position[2]);
	}
		if( point_inside_tetrahedron(
			triangle1->points[0]->position,
			triangle1->points[1]->position,
			triangle1->points[2]->position,
			triangle1->nuc->position,
			triangle2->nuc->position)){
		pos_x.push_back(triangle2->nuc->position[0]);
		pos_y.push_back(triangle2->nuc->position[1]);
		pos_z.push_back(triangle2->nuc->position[2]);
	}
	

	if( point_inside_tetrahedron(
			triangle2->points[0]->position,
			triangle2->points[1]->position,
			triangle2->points[2]->position,
			triangle2->nuc->position,
			triangle1->points[0]->position)){
		pos_x.push_back(triangle1->points[0]->position[0]);
		pos_y.push_back(triangle1->points[0]->position[1]);
		pos_z.push_back(triangle1->points[0]->position[2]);
	}

	if( point_inside_tetrahedron(
			triangle2->points[0]->position,
			triangle2->points[1]->position,
			triangle2->points[2]->position,
			triangle2->nuc->position,
			triangle1->points[1]->position)){
		pos_x.push_back(triangle1->points[1]->position[0]);
		pos_y.push_back(triangle1->points[1]->position[1]);
		pos_z.push_back(triangle1->points[1]->position[2]);
	}

	if( point_inside_tetrahedron(
			triangle2->points[0]->position,
			triangle2->points[1]->position,
			triangle2->points[2]->position,
			triangle2->nuc->position,
			triangle1->points[2]->position)){
		pos_x.push_back(triangle1->points[2]->position[0]);
		pos_y.push_back(triangle1->points[2]->position[1]);
		pos_z.push_back(triangle1->points[2]->position[2]);
	}

	if( point_inside_tetrahedron(
			triangle2->points[0]->position,
			triangle2->points[1]->position,
			triangle2->points[2]->position,
			triangle2->nuc->position,
			triangle1->nuc->position)){
		pos_x.push_back(triangle1->nuc->position[0]);
		pos_y.push_back(triangle1->nuc->position[1]);
		pos_z.push_back(triangle1->nuc->position[2]);
	}

	
	



	//calc volume von point list

	if( pos_x.size() == 0)	{
		pos_x.clear();
		pos_y.clear();
		pos_z.clear();
		return 0;
	}else{
				double V = 0.;
				//calc volume
				//Triangulation<3> *myTriangulation = new Triangulation<3>();
				//VoronoiDiagram *myVoronoiDiagram = new VoronoiDiagram();
			/*
				for( unsigned int i = 0 ; i< pos_x.size() ; i++){
					myTriangulation->add(new Vertex<3>(pos_x[i],pos_y[i],pos_z[i]));
					

				}
				

					myTriangulation->setDomain();

					// Set Frame Points
					myTriangulation->setFramePoints();
					myTriangulation->triangulate();
					myTriangulation->getConvexHull(10.);
				
		
					
					for( int j = 0 ; j < myTriangulation->getCountConvexFaces() ; j++){
					
						Triangle *tr = new Triangle();
						for( int l = 0 ; l < 3 ; l++){
						tr->points[0]->position[l] = myTriangulation->getConvexFace(j)[l]->x(0);
						tr->points[0]->position[l] = myTriangulation->getConvexFace(j)[l]->x(1);
						tr->points[0]->position[l] = myTriangulation->getConvexFace(j)[l]->x(2);
						}
						double top[3];
						top[0] = myTriangulation->getCentralCell()->x(0);
						top[1] = myTriangulation->getCentralCell()->x(1);
						top[2] = myTriangulation->getCentralCell()->x(2);
						 tr->calcVolume(top);
						 V += tr->Volume;
					}*/
					//my
				//	delete myTriangulation;
					
			/*		
					double mean[3] = {0.,0.,0.};
					for( unsigned int i = 0 ; i < pos_x.size(); i++){
						mean[0] += pos_x[i];
						mean[1] += pos_y[i];
						mean[2] += pos_z[i];
					}

				mean[0] /= pos_x.size();
				mean[1] /= pos_y.size();
				mean[2] /= pos_z.size();

				double mean_r  = 0.;


				for( unsigned int i = 0 ; i < pos_x.size(); i++){
					double v[3] = { pos_x[i]-mean[0] , pos_y[1]-mean[1] , pos_z[2]-mean[2]};
					mean_r += norm(v);
				}

				mean_r /= pos_x.size();

			
				V  = 4./3.*PI*mean_r*mean_r*mean_r;
			
			*/
			

				double min_x = pos_x[0];
				double max_x = pos_x[0];
				double min_y = pos_y[0];
				double max_y = pos_y[0];
				double min_z = pos_z[0];
				double max_z = pos_z[0];


				for( unsigned int i = 1 ; i < pos_x.size() ; i++){
					if( pos_x[i] < min_x )
						min_x = pos_x[i];
					if( pos_x[i] > max_x )
						max_x = pos_x[i];
					if( pos_y[i] < min_y )
						min_y = pos_y[i];
					if( pos_y[i] > max_y )
						max_y = pos_y[i];
					if( pos_z[i] < min_z)
						min_z = pos_z[i];
					if( pos_z[i] > max_z )
						max_z = pos_z[i];
				}

				V = (max_x-min_x) * (max_y-min_y) * (max_z-min_z);



				pos_x.clear();
				pos_y.clear();
				pos_z.clear();

				return  V;
			
			}			
}
