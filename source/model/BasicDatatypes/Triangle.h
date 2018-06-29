///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Triangle.h                                                           //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-05-31 22:00:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef TRIANGLE_H
#define TRIANGLE_H


class MassPoint;
class Triangle
{
public:
/*
      p2
      |\
      | \
      |  \
      |   \
      |    \
      |     \
      |      \
    s2|       \ s1
      |        \
      |         \
      |          \
      |           \
      |------------\
     p0    s0       p1
*/

  short mIndex;

  // connectiong mass points
  MassPoint* points[3];
  // connection springs
  Spring* springs[3];
  Spring* springs_cytoskeleton[3];
  Triangle* triangles[3];

  MassPoint* next_triangle_points[3];

  MassPoint *nuc;

  double store_normalvetor[3];
  double store_old_normalvector[3];
  
  double mAngels[3];//corresponding angels in triangle

  double mNormalvector[3];//normalvector
  double mArea;//area

  void setNormalVector();
  void setArea(double  p);

	// volume of tetraeder
	double Volume;
	void calcVolume(double top[3]);

	// cell division
	int divideCircle;

  unsigned int index;


	//visualisation
    unsigned char r;
    unsigned char g;
    unsigned char b;
	double transparency;
    void SetColor(unsigned char r, unsigned char g, unsigned char b);
	bool interaction;


	Triangle();

	int equal(Triangle *triangle);
  bool equal(MassPoint *p1, MassPoint *p2);
	void difference(Triangle *triangle);

	void checkAndChangeOrientation();
	void copyfrom( Triangle *triangle_basic);
	
	void calcMean(double center[3]);
	
	double meanLength(double mean[3]);// calc the mean length of all points (3) to a given point

	bool check_intersection_triangle_section( double direction[3], double start[3], double end[3], double length, double intersection[3]);
	bool check_intersection_triangle_section2(double p0[3],double v[3], double u[3], double direction[3], double start[3], double end[3], double length, double intersection[3]);
	bool check_point_inside_tetrahedron(double point[3], double p0i[3]);
	bool check_point_inside_tetrahedron(double point[3]);

	bool intersectionWithLine(MassPoint *m0, MassPoint *m1);       // calc if there is intersection of a line between two points and this triangle
	void CalcIntersectionWithLine(MassPoint *m, MassPoint *m_new); // calcs the intersection point
	
	double check_intersection_tetrahedral_naive(Triangle *triangle1,Triangle *triangle2);

  int missing_point( Triangle *t);
	};

#endif
