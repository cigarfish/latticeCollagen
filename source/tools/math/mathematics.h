/**
 * @file	mathematics.h
 * @date	11.11.2011
 * @author	Johannes Neitsch
 * @brief	collection of mathematic algorithm
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

#ifndef MATHEMATICS_H_
#define MATHEMATICS_H_


//#include "../BasicDatatypes/MassPoint.h"
//#include "../BasicDatatypes/Spring.h"
//#include "../Cell/CellTriangulated.h"
#include <vector>
class MassPoint;
class Spring;


//#define min( a, b)	(a<b?a:b)
#define MIN( a, b)	(a<b?a:b)
//#define max( a, b)	(a>b?a:b)
#define MAX( a, b)	(a>b?a:b)



void add_columns(double A[][3], double D[]);

void mult_matrix_3times3_vector(double A[][3], double D[], double back[]);

void mult_matrix_3times3_component(double A[][3], double B[][3]);

void half_matrix_3times3(double A[][3], double B[][3]);

void inv_matrix_3times3(double A[][3], double A_inv[][3]);

void mean(double a[], double b[], double c[], double r[]);

double det(double A[][3]);
double det(double A0[3],double A1[3],double A2[3]);

double outCircle(MassPoint* m0, MassPoint*m1, MassPoint* m2, MassPoint* m_return);

MassPoint* inCircle(MassPoint* m1, MassPoint*m2, MassPoint* m3);

double norm(double v[3]);
double norm(double v[3],double m0[3]);
double normvector(double v[3]);


void crossProduct( double v1[3], double v2[3], double v_r[3]);
void crossProduct( double v0[3], double v1[3], double v2[3], double v_r[3]);
double dotmult(double v0[3], double v1[3]);

double area(double p0[3],double p1[2], double p2[3]);

/*
double calcRadius(CellTriangulated *zelle);
*/
double max_points(std::vector <MassPoint*> mass);
double dist(double v0[3], double v1[3]);
void dist(double v0[3], double v1[3], double V_return[3]);
//double dist(Vector3f v0, Vector3f v1);

double rotateCylinder(double tmp_rot[3], Spring *s);
double rotateCylinder(double tmp_rot[3],double a[3],double b[3]);
double rotateCylinder(double tmp_rot[3],double a[3],double b[3], double scale);
double rotateCylinder(double tmp_rot[3],float a[3],float b[3], double scale);

bool intersect_triangle_section(double p0[3],double p1[3],double p2[3],double start[3],double end[3],double intersection[3]);
bool intersect_triangle_section_area(double p0[3],double p1[3],double p2[3],double start[3],double end[3],double intersection[3]);
bool point_inside_tetrahedron(double p0[3],double p1[3],double p2[3],double p[3],double point[3]);

double det_4(double *A[4]);


double dist_point_section(double start[3],double end[3], double p[3], double intersection[3]);


double angle(double v1[3], double v2[3]);

#endif /* MATHEMATICS_H_ */


