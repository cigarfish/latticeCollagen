///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Model3D.h                                                            //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012                                                                 //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef MODEL_3D_H
#define MODEL_3D_H

#include "../Cell/CellTriangulated.h"
#include "../BasicDatatypes/MassPoint.h"

#ifndef WIN32
  #include <zlib.h>
#endif

#include "CSModel.h"

#include "../BiologyLink.h"

#include <sstream>

class QXmlStreamWriter;
class QXmlStreamReader;

class QCSSimulationThread;

class BoundingBoxList;

class CSParameterContext;


class Model3D :	public CSModel
{
 // friend class NewCells;

protected:



  // Store parameter of default cell
  // Default radius of the cells at t=0 and right after cell division
  double defaultInitialCellRadius;

  // Default mass of mass point
  double default_PointMass;

/*
 * Default properties of springs
 */
  double default_k_huell;
  double default_nu_huell;
  double default_k_cytoskeleton;
  double default_nu_cytoskeleton;
  double default_damper;
  double default_k2;
  // Default pressure inside a cell
  double default_pressure_inside;
  // Default number of masspoint at surface 
  unsigned int default_Mass_Number;

  // Defined resistance of enviroment -- maybe its called also default
  double nu_medium;

  // Default triangulated cell
  CellTriangulated *defaultCellTriangulated;


/*
 * Model Parameters
 */
  //Volume of all Cells
  double volume_of_all_cells;
  //Debug for interaction
  int able_cell1;
  int able_cell2;
  //if move of single point > max_move --> reduce timestep
  double max_move;

  //time
  double defaultCellCycleTime;
  double defaultCellCycleSD;

  // State of the cell and choose experiment number
  // Define a special case of simulation
  unsigned int mode;

  bool conserve_volume;
  double pressure_threshold;
  bool force_profile;
  double force_push_two_cells;
  bool able_to_move;

  std::string output_prefix;
  std::string output_suffix;

  double mean_volume;

  double store_time;
  double store_time_old;
  double store_time_old_old;


  //! List of cells
  std::vector < CellTriangulated* > mStoreOldCells;
  std::vector < CellTriangulated* > mStoreNewCells;

  //! Alternative cell population:
  BoundingBoxList * cells2;

  vector <MassPoint*> watershed;

  // the timeStep will be calculated as std::pow( 2, log2timeStep )
  int mLog2timeStep;

  // simulation time limit
  int max_simulation_time;
  // Define frequency of output, every frequency. timeStep
  int frequency;


  CSParameterContext * mpParameters;
  void RegisterParameters();




/*
 * Boundary Box
 */
  //border
  double x_left,x_right,y_left,y_right,z_left,z_right;
  bool boundary_box;




/*
 * Povray
 */
  //povray camera
  double pov_cam_x;
  double pov_cam_y;
  double pov_cam_z;
  int print_sinusoids;
  //cut edge out
  bool mCutOut;
  int mCoutOutMax;
  bool mPovraySmooth;
  double timeBetweenOutputs;
  double currentTimeSinceLastOutput;
  bool enablePovrayOutput;
  unsigned int countpov;
  string outputPath;
  string inputPath;
  string watershedFitCells;
  bool compress_pov;
  double timeSinceLastOutput;





  // Movement of cells
  void calcForce();
  double updatePosition(int number, double delta_t);
  void calcAcceleration(int number);
  void updateVelocity(int number, double delta_t);
  void addConstantForce(double x, double y, double z);
  void clearForce();

  bool check_overlap();
  bool check_overlap_past();

  // Devide cells
  void cutCells(int number);

  // Do something dependend of mode
  void statusDepend();
  void updateCellCycle();

  // Calc interactions of other cells
  void interactions();
  void interactions2();

  //output
  void output();
  int mCounterOutput;

public:

  enum mMode {
      SingleCell = 102,
      StretchCell = 200,
      WatershedFit = 2001,
      PushTwoCells = 201,
      Pressure = 202,
      Cylinder = 301,
      Plane = 302,
      Wall = 203,
      Cube = 204
  };

  // simulation time limit
  double time_simulation;


  //! The individual model instantiation name (in case of multiple instances):
  std::string xmlName;
  static const std::string xmlType;

  //! List of cells
  std::vector < CellTriangulated* > cells;

  //! Included default cell parameters with biological units and methods for unit conversions, is used to get dimensionless default parameters for model
  BiologyLink *biolink;



  //! Adds and initializes a new cell
  void AddCell(double x, double y, double z, double r=0.5);




  double move(int number);
  void calc();
  void calcLength();
  void calcLengthWithNormDirection();
  // Change color after deformation
  void color_triangle_deform();

  void updateSecondDamper( double timeStep);
  void setMode(int i);

  void setAndSolveSystem();


  void calcForceOfDeformation();

  Model3D();
  ~Model3D(){};

  // Preliminary init
  void Reset();

  // Preliminary simulation step
  void Simulate();
  void SimulateInThread();
  void Simulate(int i);

  void removeLastStep();
  void storeLastStep();

  void setStretchPoints();
  void setObservationPoints();

  void setDefaultCell();

//  void setConserveVolume(bool box);
  int getDefaultNumberOfMassPoints();

  void writePovray( ofstream & pov);
  void writePovray_all();
  void writePovray_details( ofstream & pov);
  #if (!WIN32)
    void writePovray_details2_gzip( gzFile & povFile);
  #endif
  void writePovray_details_cut( ofstream & pov);
  void writeVTK( ofstream & pov );
  void writeVTKCell0( ofstream & pov );

  void calcVolume();

  CSParameterContext * GetParameters( std::string contextName ="" );
  void InitParameters( CSParameterContext * fromContext=0 );

  double GetSimulationProgress();

  void ClearAbleToMove();

  static void DefaultParameters( CSParameterContext * );

  // save to and load from XML:
  // save
  void writeXML( QXmlStreamWriter * ) const {};

  // load
  static Model3D * createFromXML( QXmlStreamReader * input,
                                  std::stringstream & errors,
                                  std::stringstream & warnings );

  void writeHDF5( H5::H5File * /*outputFile*/ ) const {};

  // method for reading the hdf5 data into an already allocated *Model
  // used after the createFromXML has created a certain Model and initialised
  // the parameters and, most importantly, its name.
  void readModelData( H5::H5File * /*inputFile*/,
                      std::stringstream & /*errors*/,
                      std::stringstream & /*warnings*/ ) {};

  // Do whatever else has to be done to be able to run a simulation.
  void SetupSimulation();
};

#endif
