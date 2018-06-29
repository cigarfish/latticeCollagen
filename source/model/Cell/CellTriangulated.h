#ifndef CELLTRIANGULATED_H
#define CELLTRIANGULATED_H

#include "Cell.h"
#include "../Elements/ModelElementTriangulated.h"

#include <vector>

class MassPoint;
class Spring;
class Triangle;

#define pressure_last_number 100


class CellTriangulated :  public Cell, public ModelElementTriangulatedCell
{
public:
  CellTriangulated();

  //! Writes out an XML representation of this cell to the given QXmlStreamWriter
  void writeXML( QXmlStreamWriter * ) const {};

  double initialRadius;

  //Data
  double lengthAtWall;
  int count_force_point;
  bool can_divide;
  double volume_correction;
  double init_volume;

  bool store_can_divide;


  int pressure_low_since;
  bool pressure_high;
  double pressure_threshold;

  int count_stretch;

  bool able_to_move;

  double **A;

  //Elements (MassPoints, Springs, Triangles)
  std::vector <Spring*> spring;
  std::vector <Spring*> mStoreOldSpring;
  std::vector <Spring*> spring_divide;
  std::vector <Spring*> mStoreOldSpring_divide;
  std::vector <Spring*> spring_huell;
  std::vector <Spring*> spring_cytoscelett;
  std::vector <Spring*> spring_divide_cytoskeleton;
  std::vector <Spring*> mStoreOldSpring_divide_cytoskeleton;

  std::vector <Triangle*> triangle;
  std::vector <Triangle*> triangle_divide_ring;

  std::vector <MassPoint*> mass;
  std::vector <MassPoint*> mass_huell;



  //center and radii
  double radius_cell;

  //physical properties
  double nu_medium;//for enviroment
  double p;//inner stress
  double presure_inside;
  double mean_pressure;
  double k_hull;
  double nu_hull;
  double k_cytoskelett;
  double nu_cytoskelett;
  double eta_damper;
  double k2;

  double eta;
  double gamma;

  double point_mass;
  double point_mass_invers;

  unsigned int dim;
  unsigned int mass_number;

  int status;//0 normal, 1 divide circle//10 walls//11 gegen Wand
  int mStoreOldStatus;

  double time_relax;
  double relax_time;
  double meanSpringLength;

  //Optical stretcher experiments
  int stretch_long_start;
  int stretch_long_end;
  int stretch_short_start;
  int stretch_short_end;


  double v_reference;

  //cells as sinusoids
  int type;//0 cell, 1 sinusoid


  //create cells
  void triangulate( );
  void copyfrom( CellTriangulated *default_cell, double x, double y, double z, double r);
  void setMatrixA();

  //color
  void color_triangle_deform();


  //force and move
  void addInnerStress();
  void calcForce_force3d(double timeStep);
  void calcForce_froce3d_viscous(double timeStep);
  void calcForce_froce3d_spring(double timeStep);
  void projection_to_spring();
  void calcForce_force3d_second_damper();
  void force_add_to_points();
  void updateSecondDamper(double timeStep);
  void updateDamper_new(double timeStep);
  void calcForce_medium();
  void calcForce(double timeStep);
  void calcLength();
  void calcLengthWithNormDirection();
  void calcAcceleration(int number);
  double updateVelocity(int number, double delta_t);
  double updatePosition(int number, double delta_t);
  void addConstantForce(double x, double y, double z);
  void calcInteractionTetrahedral(CellTriangulated *cell_static);
  void calcInteractionTetrahedral_basic(CellTriangulated *cell_static);
  void clearForce();

  void setAndSolveSystem();

  void calcInteractionTetrahedral_2(CellTriangulated *cell_static);
  void calcInteractionTriangle(CellTriangulated *cell_static);
  void calcInteractionTriangle_basic(CellTriangulated *cell_static);
  bool calcInteraction_repulsive(CellTriangulated *cell_static);
  void checkTriangle();
  double calcVolume();

  //Interaction between cells and self interaction
  void calcInteraction(CellTriangulated *cell_static);
  void calcSelfInteraction();
  void calcSelfInteraction_all();
  void calcSelfInteraction_all2();
  bool calcSelfInteraction_reverseStep();
  void calcSelfInteraction_angle();

  bool check_overlap(CellTriangulated *cell_static);
  bool check_overlap_past(CellTriangulated *cell_static);
  //Remove last step; store before
  void removeLastStep();
  void storeLastStep();

  //Cell division
  void addDivideCircle(double v_x, double v_y, double v_z, int number, double nu);
  double divideCircleLength();
  void insertPoint(Triangle *triangle);
  void removeSpring( int i);
  double calcRadius();
  void projectToCell(CellTriangulated *originalCell, int side);

  void statusDepend( double timeStep, double time); 
  void updateCellCycle( double timeStep, double time);

  void createInternalStructureOfCell();//normaly calles from triangulate

  virtual ModelElement * Divide(CSModel *) {return NULL;};
};

#endif
