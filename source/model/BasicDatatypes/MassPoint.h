///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  MassPoint.h                                                          //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-05-31 22:00:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef MASS_Point_H
#define MASS_Point_H

#include <vector>

using namespace::std;

class Spring;
class Triangle;

class MassPoint
{
public:

  MassPoint();
  ~MassPoint();

  double F_debug_medium[3];
  double F_debug_spring[3];

  short mIndex;



  double F_repulse[3];
  double F_springs_k[3];
  double dist_max;
  double store_dist_max;
  double store_old_dist_max;
  double F_max;
  double store_F_max;
  double store_old_F_max;

  double F_repulsive_old;
  double store_F_repuslive;
  double store_old_F_repulsive;
  double store_F[3];
  double store_old_F[3];



  double dist_point;
  double dist_triangle;

  double dist_move[3];
  double store_dist_move[3];
  double store_old_dist_move[3];
  double dist_v[3];
  
  //Data
  double position[3];         // position
  double v[3];                // velocity
  double F[3];                // force
  double norm_force;
  double a[3];                // acceleration
  double a_old[3];            // old acceleration
  double store_a_old[3];      // old acceleration
  double store_old_a_old[3];  // old acceleration
  double store_a_new[3];      // old acceleration
  double store_old_a_new[3];  // old acceleration
  double a_new[3];
  double point_mass;          // mass of this point
  double point_mass_invers;   // inverse of the mass

  double store_position[3];
  double store_old_position[3];
  double store_v[3];
  double store_old_v[3];
  double store_a[3];
  double store_old_a[3];

  std::vector <Spring*> mStoreNeighbourSprings;

  // norm of the force
  double norm_F[3];

  double log_move_max;
  double length_inside;
  double length_inside_old;

  // visualisation
  unsigned char r;
  unsigned char g;
  unsigned char b;
  double transparency;
  void SetColor(unsigned char r, unsigned char g, unsigned char b);

  // area of voronoi
  double mAreaVoronoi;

  // cell division
  unsigned int divide_number;//1 left, 2 right, 0 circle, 3 midd point

  //store neighbouring springs
  std::vector <Spring*> neighbourSprings;
  //store how many springs belong to surface (for voronoi calculation)
  int mSpringsOnSurface;


  // copy
  void copyfrom(MassPoint *mass_basic);

  // change length, used by reconstruction of real data
  void changeLength(double m_position[3], double l_c);
  // force profile of optical stretcher experiments
  double F_force_profile[3];

  bool self_interaction;

  unsigned int index;

  double radius;

  std::vector<Triangle*> mvpTriangles;   //stores the adjacent triangles
  std::vector<int> mvpPointsL;  //stores the 2 remaining points of the adjacent triangles
  std::vector<int> mvpPointsR;  //stores the 2 remaining points of the adjacent triangles
  std::vector<Spring*> mvpSprings;  //stores all springs of surface in correct order
};

#endif
