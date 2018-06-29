///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Model3D.cpp                                                          //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012                                                                 //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "Model3D.h"

#define _USE_MATH_DEFINES

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <QXmlStreamReader>

#include "../BasicDatatypes/Spring.h"
#include "../BasicDatatypes/Triangle.h"
#include "../../gui/CSGLArena.h"
#include "../../gui/QCSSimulationThread.h"
#include "../../model/Cell/CellTriangulated.h"
#include "../../tools/math/mathematics.h"
#include "../../tools/model/BoundingBoxList.h"
#include "../../tools/parameters/CSParameterChoice.h"
#include "../../tools/parameters/CSParameterContextTemporary.h"

#include "../BiologyLink.h"


#if (!WIN32)
	#include <zlib.h>
#endif

const std::string Model3D::xmlType = "Model3D";

//initialize model
Model3D::Model3D() : CSModel (3) {

  mRandom.Init();

  mpParameters = NULL;
  RegisterParameters();

  //this->biolink = NULL;
  // Set initial pressure inside a cell
  default_pressure_inside = 1.0;//only at status 202

  // Create new cell population of triangulated cells

  // Set frequency of output
  frequency = 1000;

  //Starting time
  time = 0.;
  currentTimeSinceLastOutput = 0.;
  countpov = 0;
  cells2 =new BoundingBoxList(3,10000);
  compress_pov = false;
  mCutOut = false;
  mCoutOutMax = 20;
  mPovraySmooth = false;
  able_to_move = 1;

  this->conserve_volume = true;
  this->store_time_old_old = 0.;
  this->store_time_old = 0.;
  this->able_cell1 = 0;
  this->able_cell2 = 0;
  this->pov_cam_z = 0.;
  this->pov_cam_x = 0.;
  this->pov_cam_y = 0.;
  this->x_left = 0.;
  this->x_right = 0.;
  this->y_left = 0.;
  this->y_right = 0.;
  this->z_left = 0.;
  this->z_right = 0.;
  this->default_k_huell = 0.;
  this->max_simulation_time = 0.;
  this->default_Mass_Number = 0;
  this->mLog2timeStep = 0;
  this->mode = 0;
  this->timeBetweenOutputs = 0.;
  this->boundary_box = 0;
  this->pressure_threshold = 0.;
  this->default_PointMass = 0.;
  this->default_k2 = 0.;
  this->force_push_two_cells = 0.;
  this->volume_of_all_cells = 0.;
  this->default_nu_huell = 0.;
  this->mean_volume = 0.;
  this->nu_medium = 0.;
  this->default_nu_cytoskeleton = 0.;
  this->default_damper = 0.;
  this->enablePovrayOutput = 0.;
  this->force_profile = 0.;
  this->time_simulation = 0.;
  this->max_move = 0.;
  this->print_sinusoids = 0;
  this->store_time = 0.;
  this->default_k_cytoskeleton = 0.;
  this->defaultCellTriangulated = NULL;
  this->defaultInitialCellRadius = 0.;
  this->defaultCellCycleSD = 0.;
  this->defaultCellCycleTime = 0.;

  this->mCounterOutput = -1;
  this->timeSinceLastOutput = 0.;

}



//VTK output funktion
void Model3D::writeVTK( ofstream & vtkFile ){

  vtkFile.setf(std::ios_base::fixed);
  vtkFile.precision(5);

  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << "data" << std::endl;
  vtkFile << "ASCII" << std::endl;
  vtkFile << "DATASET POLYDATA" << std::endl;


  //count nodes over all cells
  int numberNodes = 0;
  for( unsigned int i = 0 ; i < this->cells.size() ; i++ )
    numberNodes += cells[i]->mass.size();


  vtkFile << "POINTS " << numberNodes << " float" << std::endl;

  //write all position out
  for( unsigned int i = 0 ; i < this->cells.size() ; i++ )
    for( unsigned int j = 0 ; j < this->cells[i]->mass.size() ; j++)
      vtkFile << this->cells[i]->mass[j]->position[0] << " " << this->cells[i]->mass[j]->position[1] << " " << this->cells[i]->mass[j]->position[2] << std::endl;

  //count triangles over all cells
  int numberTriangles = 0;
  for( unsigned int i = 0 ; i < this->cells.size() ; i++ )
    numberTriangles += cells[i]->triangle.size();

  vtkFile << "POLYGONS " << numberTriangles << " " << numberTriangles*4 << std::endl;

  //output all triangles
  int numberPointsBevor = 0;
  for( unsigned int i = 0 ; i < this->cells.size() ; i++ ){
    for( unsigned int j = 0 ; j < this->cells[i]->triangle.size() ; j++ ){
      vtkFile << "3 "
              << ( this->cells[i]->triangle[j]->points[0]->mIndex + numberPointsBevor ) << " "
              << ( this->cells[i]->triangle[j]->points[1]->mIndex + numberPointsBevor ) << " "
              << ( this->cells[i]->triangle[j]->points[2]->mIndex + numberPointsBevor ) << std::endl;
    }
    numberPointsBevor += this->cells[i]->mass.size();
  }



  //output additional information for points
  vtkFile << "POINT_DATA " << numberNodes << std::endl;
  vtkFile << "SCALARS fabs float" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;

  //write additinal information
  for( unsigned int i = 0 ; i < this->cells.size() ; i++ )
    for( unsigned int j = 0 ; j < this->cells[i]->mass.size() ; j++)
      vtkFile << norm(this->cells[i]->mass[j]->F) << std::endl;

  vtkFile << "SCALARS f_rep float" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;

  for( unsigned int i = 0 ; i < this->cells.size() ; i++ )
    for( unsigned int j = 0 ; j < this->cells[i]->mass.size() ; j++)
      vtkFile << (this->cells[i]->mass[j]->F_repulse[0]) << std::endl;

  vtkFile << "SCALARS Area float" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;

  for( unsigned int i = 0 ; i < this->cells.size() ; i++ )
    for( unsigned int j = 0 ; j < this->cells[i]->mass.size() ; j++)
      vtkFile << this->cells[i]->mass[j]->mAreaVoronoi << std::endl;

  vtkFile << "SCALARS countTriangles float" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;

  for( unsigned int i = 0 ; i < this->cells.size() ; i++ )
    for( unsigned int j = 0 ; j < this->cells[i]->mass.size() ; j++)
      vtkFile << (double)this->cells[i]->mass[j]->mvpTriangles.size() << std::endl;

  vtkFile << "VECTORS F float" << std::endl;
  for( unsigned int i = 0 ; i < this->cells.size() ; i++ )
    for( unsigned int j = 0 ; j < this->cells[i]->mass.size() ; j++)
        vtkFile << this->cells[i]->mass[j]->F[0] << " " << this->cells[i]->mass[j]->F[1] << " " << this->cells[i]->mass[j]->F[2] << " " << std::endl;


  vtkFile << "VECTORS Normal float" << std::endl;
  double ***m = (double***) malloc( sizeof(double**)*this->cells.size());
  for( int i = 0 ; i < this->cells[i]->mass.size() ; i++){
    m[i] = (double**) malloc( sizeof(double*)*this->cells[i]->mass.size());
    for( int j = 0 ; j < this->cells[i]->mass.size(); j++)
      m[i][j] = (double*) malloc( sizeof(double)*3);
  }

  for( unsigned int i = 0 ; i < this->cells.size() ; i++ ){
    for( unsigned int j = 0 ; j < this->cells[i]->mass.size() ; j++ ){
      m[i][j][0] = 0.;
      m[i][j][1] = 0.;
      m[i][j][2] = 0.;
    }
    for( unsigned int j = 0 ; j < this->cells[i]->triangle.size() ; j++){
      for( unsigned int l = 0 ; l < 3 ; l++ ){
        m[i][this->cells[i]->triangle[j]->points[l]->mIndex][0] =+ this->cells[i]->triangle[j]->mNormalvector[0]/3.;
        m[i][this->cells[i]->triangle[j]->points[l]->mIndex][1] =+ this->cells[i]->triangle[j]->mNormalvector[1]/3.;
        m[i][this->cells[i]->triangle[j]->points[l]->mIndex][2] =+ this->cells[i]->triangle[j]->mNormalvector[2]/3.;
      }
    }
  }
  for( unsigned int i = 0 ; i < this->cells.size() ; i++ )
    for( unsigned int j = 0 ; j < this->cells[i]->mass.size() ; j++)
        vtkFile << m[i][j][0] << " " << m[i][j][1] << " " << m[i][j][2] << " " << std::endl;

  for( int i = 0 ; i < this->cells.size() ; i++){
    for( int j = 0 ; j < this->cells[i]->mass.size() ; j++)
      free(m[i][j]);
    free(m[i]);
  }
  free(m);


}




//VTK output funktion
void Model3D::writeVTKCell0( ofstream & vtkFile ){

  vtkFile.setf(std::ios_base::fixed);
  vtkFile.precision(5);

  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << "data" << std::endl;
  vtkFile << "ASCII" << std::endl;
  vtkFile << "DATASET POLYDATA" << std::endl;


  //count nodes over all cells
  int numberNodes = this->cells[0]->mass.size();

  vtkFile << "POINTS " << numberNodes << " float" << std::endl;

  //write all position out
  for( unsigned int j = 0 ; j < this->cells[0]->mass.size() ; j++)
    vtkFile << this->cells[0]->mass[j]->position[0] << " " << this->cells[0]->mass[j]->position[1] << " " << this->cells[0]->mass[j]->position[2] << std::endl;

  //count triangles over all cells
  int numberTriangles = cells[0]->triangle.size();

  vtkFile << "POLYGONS " << numberTriangles << " " << numberTriangles*4 << std::endl;

  //output all triangles
  for( unsigned int j = 0 ; j < this->cells[0]->triangle.size() ; j++ ){
    vtkFile << "3 "
            << ( this->cells[0]->triangle[j]->points[0]->mIndex ) << " "
            << ( this->cells[0]->triangle[j]->points[1]->mIndex ) << " "
            << ( this->cells[0]->triangle[j]->points[2]->mIndex ) << std::endl;
  }



  //output additional information for points
  vtkFile << "POINT_DATA " << numberNodes << std::endl;
  /*
  vtkFile << "SCALARS fabs float" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;

  for( unsigned int j = 0 ; j < this->cells[0]->mass.size() ; j++)
    vtkFile << norm(this->cells[0]->mass[j]->F) << std::endl;

  vtkFile << "SCALARS f_rep float" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;

  for( unsigned int j = 0 ; j < this->cells[0]->mass.size() ; j++)
    vtkFile << (this->cells[0]->mass[j]->F_repulse[0]) << std::endl;
*/

  vtkFile << "SCALARS countTriangles float" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;

  for( unsigned int j = 0 ; j < this->cells[0]->mass.size() ; j++)
    vtkFile << (double)this->cells[0]->mass[j]->mvpTriangles.size() << std::endl;

  vtkFile << "VECTORS F float" << std::endl;
  for( unsigned int j = 0 ; j < this->cells[0]->mass.size() ; j++)
    vtkFile << this->cells[0]->mass[j]->F[0] << " " << this->cells[0]->mass[j]->F[1] << " " << this->cells[0]->mass[j]->F[2] << " " << std::endl;


  vtkFile << "VECTORS Normal float" << std::endl;

  double **m = (double**) malloc( sizeof(double*)*this->cells[0]->mass.size());
  for( int i = 0 ; i < this->cells[0]->mass.size() ; i++)
    m[i] = (double*) malloc( sizeof(double)*3);

  for( unsigned int j = 0 ; j < this->cells[0]->mass.size() ; j++ ){
    m[j][0] = 0.;
    m[j][1] = 0.;
    m[j][2] = 0.;
  }
  for( unsigned int j = 0 ; j < this->cells[0]->triangle.size() ; j++){
    for( unsigned int l = 0 ; l < 3 ; l++ ){
      m[this->cells[0]->triangle[j]->points[l]->mIndex][0] =+ this->cells[0]->triangle[j]->mNormalvector[0]/3.;
      m[this->cells[0]->triangle[j]->points[l]->mIndex][1] =+ this->cells[0]->triangle[j]->mNormalvector[1]/3.;
      m[this->cells[0]->triangle[j]->points[l]->mIndex][2] =+ this->cells[0]->triangle[j]->mNormalvector[2]/3.;
    }
  }
  for( unsigned int j = 0 ; j < this->cells[0]->mass.size() ; j++)
    vtkFile << m[j][0] << " " << m[j][1] << " " << m[j][2] << " " << std::endl;


  for( int i = 0 ; i < this->cells[0]->mass.size() ; i++)
    free(m[i]);
  free(m);



}



//povray output
void Model3D::writePovray_all(){

  std::stringstream pov;

  //add sinusoids if nessesary
  if( this->print_sinusoids == 2 )
    pov << "#include \"sinusoids.pov\" " << std::endl;


  //povray header
  pov << "global_settings { assumed_gamma 1}" << std::endl;
  pov << "light_source {  100*<-20, 20, -20> color rgb <1,1, 1> }\n";
  pov << "camera { location 50.0*<"
      << this->pov_cam_x << ","
      << this->pov_cam_z << ","
      << this->pov_cam_y << "> look_at <0, 0, 0>   right x*image_width/image_height angle 0.3}" << std::endl;
  pov << "background {color rgb <0.9,0.9,0.9>}" << std::endl;

  //all cells
  for( unsigned int c = 0 ; c < cells.size() ; c++){
    if( cells[c]->type == 0 ){
      if( !cells[c]->pressure_high ){
        for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){
          if( this->mCutOut )
            pov << "difference{" << std::endl;
          pov << "triangle{<"
              << cells[c]->triangle[i]->points[0]->position[0]
              << ","
              << cells[c]->triangle[i]->points[0]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[0]->position[1]
              << ">,<"
              << cells[c]->triangle[i]->points[1]->position[0]
              << ","
              << cells[c]->triangle[i]->points[1]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[1]->position[1]
              << ">,<"
              << cells[c]->triangle[i]->points[2]->position[0]
              << ","
              << cells[c]->triangle[i]->points[2]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[2]->position[1]
              << ">"
              << " texture { pigment { color rgb<1.0, 1.0, 0.0> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
            //<< " texture { pigment { color rgb<1.0, 1.0, 0.0> filter 0.5} finish {ambient 0.1 diffuse 0.9   phong 1}  }"//transparency
              << "}"<< std::endl;
          if( this->mCutOut ){
            pov << "}"<< std::endl;
            pov << "box{<0,0,0>,<"
                << this->mCoutOutMax << "," << this->mCoutOutMax << "," << -this->mCoutOutMax 
                << ">rotate <-25,15,-25>"
                << "texture {"
                << "pigment { color rgb<255./255., 115./255., 0> }"
                << "finish {ambient 0.1 diffuse 0.9   phong 1} }}}";
          }//if coutout
        }
      }
      else{
        for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){
          if( this->mCutOut )
            pov << "difference{" << std::endl;
          pov << "triangle{<" 
              << cells[c]->triangle[i]->points[0]->position[0]
              << ","
              << cells[c]->triangle[i]->points[0]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[0]->position[1]
              << ">,<"
              << cells[c]->triangle[i]->points[1]->position[0]
              << ","
              << cells[c]->triangle[i]->points[1]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[1]->position[1]
              << ">,<"
              << cells[c]->triangle[i]->points[2]->position[0]
              << ","
              << cells[c]->triangle[i]->points[2]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[2]->position[1]
              << ">"
              << " texture { pigment { color rgb<1.0, 0.5, 0.0> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
            //<< " texture { pigment { color rgb<1.0, 1.0, 0.0> filter 0.5} finish {ambient 0.1 diffuse 0.9   phong 1}  }"
              << "}"<< std::endl;
          if( this->mCutOut ){
            pov << "}"<< std::endl;
            pov << "box{<0,0,0>,<"
                << this->mCoutOutMax << "," << this->mCoutOutMax << "," << -this->mCoutOutMax 
                << ">rotate <-25,15,-25>"
                << "texture {"
                << "pigment { color rgb<255./255., 115./255., 0> }"
                << "finish {ambient 0.1 diffuse 0.9   phong 1} }}}";
          }//if coutout
        }//for
      }//else

      //all mass points
      for( unsigned int i = 0 ; i < cells[c]->mass_huell.size() ; i++){
        if( this->mCutOut )
          pov << "difference{" << std::endl;
        pov << "sphere{<" 
            << cells[c]->mass_huell[i]->position[0]
            << ","
            << cells[c]->mass_huell[i]->position[2]
            << "," 
            << cells[c]->mass_huell[i]->position[1]
            << ">,"
            << cells[c]->mass_huell[i]->radius
            << " texture { pigment { color rgb<"
            << cells[c]->mass_huell[i]->r/255. << ","
            << cells[c]->mass_huell[i]->g/255. << ","
            << cells[c]->mass_huell[i]->b/255.
            << "> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
            << "}"<< std::endl;
        if( this->mCutOut ){
          pov << "}"<< std::endl;
          pov << "box{<0,0,0>,<"
              << this->mCoutOutMax << "," << this->mCoutOutMax << "," << -this->mCoutOutMax 
              << ">rotate <-25,15,-25>"
              << "texture {"
              << "pigment { color rgb<255./255., 115./255., 0> }"
              << "finish {ambient 0.1 diffuse 0.9   phong 1} }}}";
        }//if coutout
      }

      //division ring
      for( unsigned int i = 0 ; i < cells[c]->spring_divide.size() ; i++){
        if( this->mCutOut )
          pov << "difference{" << std::endl;
        pov << "cylinder{<" 
            << cells[c]->spring_divide[i]->start->position[0]
            << ","
            << cells[c]->spring_divide[i]->start->position[2]
            << "," 
            << cells[c]->spring_divide[i]->start->position[1]
            << ">,<"
            << cells[c]->spring_divide[i]->end->position[0]
            << ","
            << cells[c]->spring_divide[i]->end->position[2]
            << "," 
            << cells[c]->spring_divide[i]->end->position[1]
            << ">"
            << cells[c]->spring_divide[i]->radius
            << " texture { pigment { color rgb<" 
            << cells[c]->spring_divide[i]->r/255. << ","
            << cells[c]->spring_divide[i]->g/255. << ","
            << cells[c]->spring_divide[i]->b/255. 
            << "> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
            << "}"<< std::endl;
        if( this->mCutOut ){
          pov << "}"<< std::endl;
          pov << "box{<0,0,0>,<"
              << this->mCoutOutMax << "," << this->mCoutOutMax << "," << -this->mCoutOutMax 
              << ">rotate <-25,15,-25>"
              << "texture {"
              << "pigment { color rgb<255./255., 115./255., 0> }"
              << "finish {ambient 0.1 diffuse 0.9   phong 1} }}}";
        }//if coutout
      }//for

      //springs at surface
      for( unsigned int i = 0 ; i < cells[c]->spring_huell.size() ; i++){
        if( (cells[c]->spring_huell[i]->start->position[0] != cells[c]->spring_huell[i]->end->position[0] ) || 
            (cells[c]->spring_huell[i]->start->position[1] != cells[c]->spring_huell[i]->end->position[1] ) || 
            (cells[c]->spring_huell[i]->start->position[2] != cells[c]->spring_huell[i]->end->position[2] ) ){
          if( this->mCutOut )
            pov << "difference{" << std::endl;
          pov << "cylinder{<" 
              << cells[c]->spring_huell[i]->start->position[0]
              << ","
              << cells[c]->spring_huell[i]->start->position[2]
              << "," 
              << cells[c]->spring_huell[i]->start->position[1]
              << ">,<"
              << cells[c]->spring_huell[i]->end->position[0]
              << ","
              << cells[c]->spring_huell[i]->end->position[2]
              << "," 
              << cells[c]->spring_huell[i]->end->position[1]
              << ">"
              << cells[c]->spring_huell[i]->radius
              << " texture { pigment { color rgb<" 
              << cells[c]->spring_huell[i]->r/255. << ","
              << cells[c]->spring_huell[i]->g/255. << ","
              << cells[c]->spring_huell[i]->b/255. 
              <<	"> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
              << "}"<< std::endl;
          if( this->mCutOut ){
            pov << "}"<< std::endl;
            pov << "box{<0,0,0>,<"
                << this->mCoutOutMax << "," << this->mCoutOutMax << "," << -this->mCoutOutMax 
               << ">rotate <-25,15,-25>"
               << "texture {"
               << "pigment { color rgb<255./255., 115./255., 0> }"
               << "finish {ambient 0.1 diffuse 0.9   phong 1} }}}";
          }//if coutout
        }//if
      }//for

      for( unsigned int i = 0 ; i < cells[c]->spring_cytoscelett.size() ; i++){
        if( this->mCutOut )
          pov << "difference{" << std::endl;
        pov << "cylinder{<" 
            << cells[c]->spring_cytoscelett[i]->start->position[0]
            << ","
            << cells[c]->spring_cytoscelett[i]->start->position[2]
            << "," 
            << cells[c]->spring_cytoscelett[i]->start->position[1]
            << ">,<"
            << cells[c]->spring_cytoscelett[i]->end->position[0]
            << ","
            << cells[c]->spring_cytoscelett[i]->end->position[2]
            << "," 
            << cells[c]->spring_cytoscelett[i]->end->position[1]
            << ">"
            << cells[c]->spring_cytoscelett[i]->radius
            << " texture { pigment { color rgb<" 
            << cells[c]->spring_cytoscelett[i]->r/255. << ","
            << cells[c]->spring_cytoscelett[i]->g/255. << ","
            << cells[c]->spring_cytoscelett[i]->b/255. 
            << "> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
            << "}"<< std::endl;
        if( this->mCutOut ){
          pov << "}"<< std::endl;
          pov << "box{<0,0,0>,<"
              << this->mCoutOutMax << "," << this->mCoutOutMax << "," << -this->mCoutOutMax 
              << ">rotate <-25,15,-25>"
              << "texture {"
              << "pigment { color rgb<255./255., 115./255., 0> }"
              << "finish {ambient 0.1 diffuse 0.9   phong 1} }}}";
        }//if coutout
      }

      //springs from force profile
      if( this->mode == StretchCell || mode == 2000 ){
        for( unsigned int i = 0 ; i < cells[c]->mass_huell.size() ; i++){
          double tmp[3] = { cells[c]->mass_huell[i]->position[0] + cells[c]->mass_huell[i]->F_force_profile[0],
                            cells[c]->mass_huell[i]->position[1] + cells[c]->mass_huell[i]->F_force_profile[1],
                            cells[c]->mass_huell[i]->position[2] + cells[c]->mass_huell[i]->F_force_profile[2]};
          if( this->mCutOut )
            pov << "difference{" << std::endl;
          pov << "cylinder{<" 
              << cells[c]->mass_huell[i]->position[0]
              << ","
              << cells[c]->mass_huell[i]->position[2]
              << "," 
              << cells[c]->mass_huell[i]->position[1]
              << ">,<"
              << tmp[0]
              << ","
              << tmp[2]
              << "," 
              << tmp[1]
              << ">"
              << cells[c]->spring[0]->radius
              << " texture { pigment { color rgb<"
              << 0 <<","
              << 1 <<","
              << 1 
              << "> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
              << "}"<< std::endl;
          if( this->mCutOut ){
            pov << "}"<< std::endl;
            pov << "box{<0,0,0>,<"
                << this->mCoutOutMax << "," << this->mCoutOutMax << "," << -this->mCoutOutMax 
                << ">rotate <-25,15,-25>"
                << "texture {"
                << "pigment { color rgb<255./255., 115./255., 0> }"
               << "finish {ambient 0.1 diffuse 0.9   phong 1} }}}";
          }//if coutout
        }//for
      }//if strech mode
    }//if help cell
  }//for all cells


  //add sinusoids if nessesary
  if( this->print_sinusoids == 2 )
    pov << "object{ sinusoids }" << std::endl;




  //create output file
  string filename;
  filename.append(outputPath+"/pov/");
  filename.append("complexCells"+output_prefix+"_");

  stringstream NumberString;
  NumberString.precision(10);
  NumberString << std::fixed << countpov/std::pow(10.,10) ;
  filename.append("1"+NumberString.str().substr(2, 100));
  countpov++;

  filename.append(".pov");


  if( this->compress_pov ){
    #if (!WIN32)
      filename.append(".gz");

      gzFile zipFile = gzopen(filename.c_str(), "wbt");
      gzwrite( zipFile, pov.str().c_str(), pov.str().size() );
      gzclose(zipFile);

    #else
      cout << "Error no Zlib is installed" << endl;
    #endif
  }//if compressed
  else{
    ofstream pov_file;
    pov_file.open(filename.c_str(),ios::out );
    pov_file << pov.rdbuf() << std::endl;
    pov_file.close();

  }//if else compressed

}


//povray and graphics and output
void Model3D::writePovray_details( ofstream & pov){

  //povray header
  pov << "global_settings { assumed_gamma 1}" << endl;
  pov << "light_source {  100*<-20, 20, -20> color rgb <1,1, 1> }\n";
  pov << "camera { location 50.0*<"
      << this->pov_cam_x << ","
      << this->pov_cam_z << ","
      << this->pov_cam_y << "> look_at <0, 0, 0>   right x*image_width/image_height angle 0.3}" << endl;
  pov << "background {color rgb <0.9,0.9,0.9>}" << endl;

  //all cells
  for( unsigned int c = 0 ; c < cells.size() ; c++){
    if( cells[c]->type == 0 ){
      if( !cells[c]->pressure_high ){
        for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){
          pov << "triangle{<"
              << cells[c]->triangle[i]->points[0]->position[0]
              << ","
              << cells[c]->triangle[i]->points[0]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[0]->position[1]
              << ">,<"
              << cells[c]->triangle[i]->points[1]->position[0]
              << ","
              << cells[c]->triangle[i]->points[1]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[1]->position[1]
              << ">,<"
              << cells[c]->triangle[i]->points[2]->position[0]
              << ","
              << cells[c]->triangle[i]->points[2]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[2]->position[1]
              << ">"
              << " texture { pigment { color rgb<1.0, 1.0, 0.0> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
            //<< " texture { pigment { color rgb<1.0, 1.0, 0.0> filter 0.5} finish {ambient 0.1 diffuse 0.9   phong 1}  }"//transparency
              << "}"<< endl;
        }
      }
      else{
        for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){
          pov << "triangle{<" 
              << cells[c]->triangle[i]->points[0]->position[0]
              << ","
              << cells[c]->triangle[i]->points[0]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[0]->position[1]
              << ">,<"
              << cells[c]->triangle[i]->points[1]->position[0]
              << ","
              << cells[c]->triangle[i]->points[1]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[1]->position[1]
              << ">,<"
              << cells[c]->triangle[i]->points[2]->position[0]
              << ","
              << cells[c]->triangle[i]->points[2]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[2]->position[1]
              << ">"
              << " texture { pigment { color rgb<1.0, 0.5, 0.0> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
            //<< " texture { pigment { color rgb<1.0, 1.0, 0.0> filter 0.5} finish {ambient 0.1 diffuse 0.9   phong 1}  }"
              << "}"<< endl;
        }//for
      }//else

      //all mass points
      for( unsigned int i = 0 ; i < cells[c]->mass_huell.size() ; i++){
        pov << "sphere{<" 
            << cells[c]->mass_huell[i]->position[0]
            << ","
            << cells[c]->mass_huell[i]->position[2]
            << "," 
            << cells[c]->mass_huell[i]->position[1]
            << ">,"
            << 0.007//cells[c]->mass_huell[i]->r
            << " texture { pigment { color rgb<"
            << cells[c]->mass_huell[i]->r/255. << ","
            << cells[c]->mass_huell[i]->g/255. << ","
            << cells[c]->mass_huell[i]->b/255.
            << "> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
            << "}"<< endl;
      }

      //division ring
      for( unsigned int i = 0 ; i < cells[c]->spring_divide.size() ; i++){
        pov << "cylinder{<" 
            << cells[c]->spring_divide[i]->start->position[0]
            << ","
            << cells[c]->spring_divide[i]->start->position[2]
            << "," 
            << cells[c]->spring_divide[i]->start->position[1]
            << ">,<"
            << cells[c]->spring_divide[i]->end->position[0]
            << ","
            << cells[c]->spring_divide[i]->end->position[2]
            << "," 
            << cells[c]->spring_divide[i]->end->position[1]
            << ">"
            << cells[c]->spring_divide[i]->radius
            << " texture { pigment { color rgb<" 
            << cells[c]->spring_divide[i]->r/255. << ","
            << cells[c]->spring_divide[i]->g/255. << ","
            << cells[c]->spring_divide[i]->b/255. 
            << "> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
            << "}"<< endl;
      }//for

      //springs at surface
      for( unsigned int i = 0 ; i < cells[c]->spring_huell.size() ; i++){
        if( (cells[c]->spring_huell[i]->start->position[0] != cells[c]->spring_huell[i]->end->position[0] ) || 
            (cells[c]->spring_huell[i]->start->position[1] != cells[c]->spring_huell[i]->end->position[1] ) || 
            (cells[c]->spring_huell[i]->start->position[2] != cells[c]->spring_huell[i]->end->position[2] ) ){
          pov << "cylinder{<" 
              << cells[c]->spring_huell[i]->start->position[0]
              << ","
              << cells[c]->spring_huell[i]->start->position[2]
              << "," 
              << cells[c]->spring_huell[i]->start->position[1]
              << ">,<"
              << cells[c]->spring_huell[i]->end->position[0]
              << ","
              << cells[c]->spring_huell[i]->end->position[2]
              << "," 
              << cells[c]->spring_huell[i]->end->position[1]
              << ">"
              << cells[c]->spring_huell[i]->radius
              << " texture { pigment { color rgb<" 
              << cells[c]->spring_huell[i]->r/255. << ","
              << cells[c]->spring_huell[i]->g/255. << ","
              << cells[c]->spring_huell[i]->b/255. 
              <<	"> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
              << "}"<< endl;
        }//if
      }//for

      for( unsigned int i = 0 ; i < cells[c]->spring_cytoscelett.size() ; i++){
        pov << "cylinder{<" 
            << cells[c]->spring_cytoscelett[i]->start->position[0]
            << ","
            << cells[c]->spring_cytoscelett[i]->start->position[2]
            << "," 
            << cells[c]->spring_cytoscelett[i]->start->position[1]
            << ">,<"
            << cells[c]->spring_cytoscelett[i]->end->position[0]
            << ","
            << cells[c]->spring_cytoscelett[i]->end->position[2]
            << "," 
            << cells[c]->spring_cytoscelett[i]->end->position[1]
            << ">"
            << cells[c]->spring_cytoscelett[i]->radius
            << " texture { pigment { color rgb<" 
            << cells[c]->spring_cytoscelett[i]->r/255. << ","
            << cells[c]->spring_cytoscelett[i]->g/255. << ","
            << cells[c]->spring_cytoscelett[i]->b/255. 
            << "> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
            << "}"<< endl;
      }

      //springs from force profile
      if( this->mode == StretchCell || mode == 2000 ){
        for( unsigned int i = 0 ; i < cells[c]->mass_huell.size() ; i++){
          double tmp[3] = { cells[c]->mass_huell[i]->position[0] + cells[c]->mass_huell[i]->F_force_profile[0],
                            cells[c]->mass_huell[i]->position[1] + cells[c]->mass_huell[i]->F_force_profile[1],
                            cells[c]->mass_huell[i]->position[2] + cells[c]->mass_huell[i]->F_force_profile[2]};

          pov << "cylinder{<" 
              << cells[c]->mass_huell[i]->position[0]
              << ","
              << cells[c]->mass_huell[i]->position[2]
              << "," 
              << cells[c]->mass_huell[i]->position[1]
              << ">,<"
              << tmp[0]
              << ","
              << tmp[2]
              << "," 
              << tmp[1]
              << ">"
              << cells[c]->spring[0]->radius
              << " texture { pigment { color rgb<"
              << 0 <<","
              << 1 <<","
              << 1 
              << "> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
              << "}"<< endl;
        }//for
      }//if strech mode
    }//if help cell
  }//for all cells
}


#if (!WIN32)
void Model3D::writePovray_details2_gzip( gzFile & povFile){

  std::stringstream pov;

  pov << "global_settings { assumed_gamma 1}" << endl;
  pov << "light_source {  100*<-20, 20, -20> color rgb <1,1, 1> }\n";
  pov << "camera { location 50.0*<"
      << this->pov_cam_x << ","
      << this->pov_cam_z << ","
      << this->pov_cam_y << "> look_at <0, 0, 0>   right x*image_width/image_height angle 0.3}" << endl;

  pov << "background {color rgb <0.9,0.9,0.9>}" << endl;

  //all triangles
  for( unsigned int c = 0 ; c < cells.size() ; c++){
    if( cells[c]->type == 0 ){
      if( !cells[c]->pressure_high ){
        for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){
          pov << "triangle{<" 
              << cells[c]->triangle[i]->points[0]->position[0]
              << ","
              << cells[c]->triangle[i]->points[0]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[0]->position[1]
              << ">,<"
              << cells[c]->triangle[i]->points[1]->position[0]
              << ","
              << cells[c]->triangle[i]->points[1]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[1]->position[1]
              << ">,<"
              << cells[c]->triangle[i]->points[2]->position[0]
              << ","
              << cells[c]->triangle[i]->points[2]->position[2]
              << "," 
              << cells[c]->triangle[i]->points[2]->position[1]
              << ">"
              << " texture { pigment { color rgb<1.0, 1.0, 0.0> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
            //<< " texture { pigment { color rgb<1.0, 1.0, 0.0> filter 0.5} finish {ambient 0.1 diffuse 0.9   phong 1}  }"
              << "}"<< endl;
        }//for
      }//if high
      else{
       for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){
         pov << "triangle{<" 
             << cells[c]->triangle[i]->points[0]->position[0]
             << ","
             << cells[c]->triangle[i]->points[0]->position[2]
             << "," 
             << cells[c]->triangle[i]->points[0]->position[1]
             << ">,<"
             << cells[c]->triangle[i]->points[1]->position[0]
             << ","
             << cells[c]->triangle[i]->points[1]->position[2]
             << "," 
             << cells[c]->triangle[i]->points[1]->position[1]
             << ">,<"
             << cells[c]->triangle[i]->points[2]->position[0]
             << ","
             << cells[c]->triangle[i]->points[2]->position[2]
             << "," 
             << cells[c]->triangle[i]->points[2]->position[1]
             << ">"
             << " texture { pigment { color rgb<1.0, 0.5, 0.0> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
           //<< " texture { pigment { color rgb<1.0, 1.0, 0.0> filter 0.5} finish {ambient 0.1 diffuse 0.9   phong 1}  }"
             << "}"<< endl;
        }//for
      }//if

      //all mass points
      for( unsigned int i = 0 ; i < cells[c]->mass_huell.size() ; i++){
        pov << "sphere{<" 
            << cells[c]->mass_huell[i]->position[0]
            << ","
            << cells[c]->mass_huell[i]->position[2]
            << "," 
            << cells[c]->mass_huell[i]->position[1]
            << ">,"
            << 0.007//cells[c]->mass_huell[i]->r
            << " texture { pigment { color rgb<"
            << cells[c]->mass_huell[i]->r/255. << ","
            << cells[c]->mass_huell[i]->g/255. << ","
            << cells[c]->mass_huell[i]->b/255.
            << "> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
            << "}"<< endl;
        }


      //division ring
      for( unsigned int i = 0 ; i < cells[c]->spring_divide.size() ; i++){
        pov << "cylinder{<" 
            << cells[c]->spring_divide[i]->start->position[0]
            << ","
            << cells[c]->spring_divide[i]->start->position[2]
            << "," 
            << cells[c]->spring_divide[i]->start->position[1]
            << ">,<"
            << cells[c]->spring_divide[i]->end->position[0]
            << ","
            << cells[c]->spring_divide[i]->end->position[2]
            << "," 
            << cells[c]->spring_divide[i]->end->position[1]
            << ">"
            << cells[c]->spring_divide[i]->radius
            << " texture { pigment { color rgb<" 
            << cells[c]->spring_divide[i]->r/255. << ","
            << cells[c]->spring_divide[i]->g/255. << ","
            << cells[c]->spring_divide[i]->b/255. 
            <<	"> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
            << "}"<< endl;
      }


      //springs at surface
      for( unsigned int i = 0 ; i < cells[c]->spring_huell.size() ; i++){
        if( (cells[c]->spring_huell[i]->start->position[0] != cells[c]->spring_huell[i]->end->position[0] ) || 
            (cells[c]->spring_huell[i]->start->position[1] != cells[c]->spring_huell[i]->end->position[1] ) || 
            (cells[c]->spring_huell[i]->start->position[2] != cells[c]->spring_huell[i]->end->position[2] ) ){
          pov << "cylinder{<" 
              << cells[c]->spring_huell[i]->start->position[0]
              << ","
              << cells[c]->spring_huell[i]->start->position[2]
              << "," 
              << cells[c]->spring_huell[i]->start->position[1]
              << ">,<"
              << cells[c]->spring_huell[i]->end->position[0]
              << ","
              << cells[c]->spring_huell[i]->end->position[2]
              << "," 
              << cells[c]->spring_huell[i]->end->position[1]
              << ">"
              << cells[c]->spring_huell[i]->radius
              << " texture { pigment { color rgb<" 
              << cells[c]->spring_huell[i]->r/255. << ","
              << cells[c]->spring_huell[i]->g/255. << ","
              << cells[c]->spring_huell[i]->b/255. 
              <<	"> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
              << "}"<< endl;
        }
      }


      for( unsigned int i = 0 ; i < cells[c]->spring_cytoscelett.size() ; i++){
        pov << "cylinder{<" 
            << cells[c]->spring_cytoscelett[i]->start->position[0]
            << ","
            << cells[c]->spring_cytoscelett[i]->start->position[2]
            << "," 
            << cells[c]->spring_cytoscelett[i]->start->position[1]
            << ">,<"
            << cells[c]->spring_cytoscelett[i]->end->position[0]
            << ","
            << cells[c]->spring_cytoscelett[i]->end->position[2]
            << "," 
            << cells[c]->spring_cytoscelett[i]->end->position[1]
            << ">"
            << cells[c]->spring_cytoscelett[i]->radius
            << " texture { pigment { color rgb<" 
            << cells[c]->spring_cytoscelett[i]->r/255. << ","
            << cells[c]->spring_cytoscelett[i]->g/255. << ","
            << cells[c]->spring_cytoscelett[i]->b/255. 
            << "> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
            << "}"<< endl;
      }

      //springs from force profile
      if( this->mode == StretchCell || mode == 2000 ){
        for( unsigned int i = 0 ; i < cells[c]->mass_huell.size() ; i++){
          double tmp[3] = { cells[c]->mass_huell[i]->position[0] + cells[c]->mass_huell[i]->F_force_profile[0],
                            cells[c]->mass_huell[i]->position[1] + cells[c]->mass_huell[i]->F_force_profile[1],
                            cells[c]->mass_huell[i]->position[2] + cells[c]->mass_huell[i]->F_force_profile[2]};
          pov << "cylinder{<" 
              << cells[c]->mass_huell[i]->position[0]
              << ","
              << cells[c]->mass_huell[i]->position[2]
              << "," 
              << cells[c]->mass_huell[i]->position[1]
              << ">,<"
              << tmp[0]
              << ","
              << tmp[2]
              << "," 
              << tmp[1]
              << ">"
              << cells[c]->spring[0]->radius
              << " texture { pigment { color rgb<"
              << 0 <<","
              << 1 <<","
              << 1 
              << "> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
              << "}"<< endl;
        }//for
      }//if
    }//if cell static
  }//all cells

  gzwrite( povFile, pov.str().c_str(), pov.str().size() );

}
#endif


void Model3D::writePovray_details_cut( ofstream & pov){

	double max = 20.;

	pov << "global_settings { assumed_gamma 1}" << endl;
	pov << "light_source {  100*<-20, 20, -20> color rgb <1,1, 1> }\n";
	pov << "camera { location 58.0*<"
		<< this->pov_cam_x << ","
		<< this->pov_cam_z << ","
		<< this->pov_cam_y << "> look_at <0, 0, 0>   right x*image_width/image_height angle 0.3}" << endl;

	pov << "background {color rgb <0.9,0.9,0.9>}" << endl;

	//all triangles
	
	for( unsigned int c = 0 ; c < cells.size() ; c++){
		if( cells[c]->type == 0 ){
		
		if( !cells[c]->pressure_high ){
		for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){
			pov << "difference{" << endl;
			pov << "triangle{<" 
				<< cells[c]->triangle[i]->points[0]->position[0]
				<< ","
				<< cells[c]->triangle[i]->points[0]->position[2]
				<< "," 
				<< cells[c]->triangle[i]->points[0]->position[1]
				<< ">,<"
				<< cells[c]->triangle[i]->points[1]->position[0]
				<< ","
				<< cells[c]->triangle[i]->points[1]->position[2]
				<< "," 
				<< cells[c]->triangle[i]->points[1]->position[1]
				<< ">,<"
				<< cells[c]->triangle[i]->points[2]->position[0]
				<< ","
				<< cells[c]->triangle[i]->points[2]->position[2]
				<< "," 
				<< cells[c]->triangle[i]->points[2]->position[1]
				<< ">"
				<< " texture { pigment { color rgb<1.0, 1.0, 0.0> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
				<< "}"<< endl;
				pov << "box{<0,0,0>,<"
					<< max << "," << max << "," << -max 
					<< ">rotate <-25,15,-25>"
					<< "texture {"
					<< "pigment { color rgb<255./255., 115./255., 0> }"
					<<	"finish {ambient 0.1 diffuse 0.9   phong 1} }}}";

		}
	
	}
	else{
		for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){
			pov << "difference{" << endl;
			pov << "triangle{<" 
				<< cells[c]->triangle[i]->points[0]->position[0]
				<< ","
				<< cells[c]->triangle[i]->points[0]->position[2]
				<< "," 
				<< cells[c]->triangle[i]->points[0]->position[1]
				<< ">,<"
				<< cells[c]->triangle[i]->points[1]->position[0]
				<< ","
				<< cells[c]->triangle[i]->points[1]->position[2]
				<< "," 
				<< cells[c]->triangle[i]->points[1]->position[1]
				<< ">,<"
				<< cells[c]->triangle[i]->points[2]->position[0]
				<< ","
				<< cells[c]->triangle[i]->points[2]->position[2]
				<< "," 
				<< cells[c]->triangle[i]->points[2]->position[1]
				<< ">"
				<< " texture { pigment { color rgb<1.0, 0.5, 0.0> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
				<< "}"<< endl;
				pov << "box{<0,0,0>,<"
					<< max << "," << max << "," << -max 
					<< ">rotate <-25,15,-25>"
					<< "texture {"
					<< "pigment { color rgb<255./255., 115./255., 0> }"
					<<	"finish {ambient 0.1 diffuse 0.9   phong 1} }}}";

		}
	
			}
		
		//all mass points
		for( unsigned int i = 0 ; i < cells[c]->mass_huell.size() ; i++){
			pov << "difference{" << endl;
			pov << "sphere{<" 
				<< cells[c]->mass_huell[i]->position[0]
				<< ","
				<< cells[c]->mass_huell[i]->position[2]
				<< "," 
				<< cells[c]->mass_huell[i]->position[1]
				<< ">,"
				<< 0.007//cells[c]->mass_huell[i]->r
				<< " texture { pigment { color rgb<"
					<< cells[c]->mass_huell[i]->r/255. << ","
					<< cells[c]->mass_huell[i]->g/255. << ","
					<< cells[c]->mass_huell[i]->b/255.
				<< "> } finish {ambient 0.1 diffuse 0.9   phong 1}  }"
				<< "}"<< endl;
				pov << "box{<0,0,0>,<"
					<< max << "," << max << "," << -max 
					<< ">rotate <-25,15,-25>"
					<< "texture {"
					<< "pigment { color rgb<"
					<< cells[c]->mass_huell[i]->r/255. << ","
					<< cells[c]->mass_huell[i]->g/255. << ","
					<< cells[c]->mass_huell[i]->b/255.
					<< ">}"
					<<	"finish {ambient 0.1 diffuse 0.9   phong 1} }}}";
		}



		//division ring
		for( unsigned int i = 0 ; i < cells[c]->spring.size() ; i++){
			pov << "difference{" << endl;
			pov << "cylinder{<" 
				<< cells[c]->spring[i]->start->position[0]
				<< ","
				<< cells[c]->spring[i]->start->position[2]
				<< "," 
				<< cells[c]->spring[i]->start->position[1]
				<< ">,<"
				<< cells[c]->spring[i]->end->position[0]
				<< ","
				<< cells[c]->spring[i]->end->position[2]
				<< "," 
				<< cells[c]->spring[i]->end->position[1]
				<< ">"
				//<< 0.003//cells[c]->mass_huell[i]->r
				<< cells[c]->spring[i]->radius
				<< " texture { pigment { color rgb<" 
					<< cells[c]->spring[i]->r/255. << ","
					<< cells[c]->spring[i]->g/255. << ","
					<< cells[c]->spring[i]->b/255. 
					<<	"> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
				<< "}"<< endl;
				pov << "box{<0,0,0>,<"
					<< max << "," << max << "," << -max 
					<< ">rotate <-25,15,-25>"
					<< "texture {"
					<< "pigment { color rgb<"
					<< cells[c]->spring[i]->r/255. << ","
					<< cells[c]->spring[i]->g/255. << ","
					<< cells[c]->spring[i]->b/255. 
					<< "> }"
					<<	"finish {ambient 0.1 diffuse 0.9   phong 1} }}}";
		}


		//springs from force profile
		if( this->mode == StretchCell || mode == WatershedFit ){
			for( unsigned int i = 0 ; i < cells[c]->mass_huell.size() ; i++){

			
				double tmp[3] = {	cells[c]->mass_huell[i]->position[0] + cells[c]->mass_huell[i]->F_force_profile[0],
									cells[c]->mass_huell[i]->position[1] + cells[c]->mass_huell[i]->F_force_profile[1],
									cells[c]->mass_huell[i]->position[2] + cells[c]->mass_huell[i]->F_force_profile[2]};

					pov << "cylinder{<" 
					<< cells[c]->mass_huell[i]->position[0]
					<< ","
					<< cells[c]->mass_huell[i]->position[2]
					<< "," 
					<< cells[c]->mass_huell[i]->position[1]
					<< ">,<"
					<< tmp[0]
					<< ","
					<< tmp[2]
					<< "," 
					<< tmp[1]
					<< ">"
					<< cells[c]->spring[0]->radius
					<< " texture { pigment { color rgb<"
						<< 0 <<","
						<< 1 <<","
						<< 1 
					<< "> }  finish {ambient 0.1 diffuse 0.9   phong 1} }"
					<< "}"<< endl;
			}




		}
	
	}
}


}


void Model3D::writePovray( ofstream &  pov){


	
		pov << "global_settings { assumed_gamma 1}" << std::endl;

//		pov << "light_source { <2, 2, 2>*500000 color rgb<1, 1, 1> shadowless}" << endl;
	//	pov << "light_source { <-2, -2, -2>*500000 color rgb<1, 1, 1> }" << endl;
//		pov << "light_source { <-2, 2, -2>*500000 color rgb<1, 1, 1> shadowless}" << endl;
	//	pov << "light_source { <2, -2, 2>*500000 color rgb<1, 1, 1> }" << endl;
//		pov << "light_source { <2, -2, -2>*500000 color rgb<1, 1, 1> shadowless}" << endl;
	//	pov << "light_source { <-2, -2, 2>*500000 color rgb<1, 1, 1> }" << endl;
	//	pov << "light_source { <-2, 2, -2>*500000 color rgb<1, 0, 1> }" << endl;
	//	pov << "light_source { <2, -2, -2>*500000 color rgb<1, 1, 0> }" << endl;

//		pov << "light_source { <0, 100000, 0> color rgb<1, 1, 1> shadowless}" << endl;

/*
		pov << "light_source {" << endl;
		pov << "<5, 100000, 5>, color rgb<1, 1, 1>" << endl;
		pov << "area_light"<<endl;
		pov << "<-10,0,0>, <0,0,-10>, 3, 3" << endl;
		pov << "adaptive 0  jitter   circular  orient" << endl;
		pov << "}"<< endl;
*/

//		   pov << "light_source { 1*x color rgb 0.7 area_light <8, 8, 0> <0, 8, 8>  4, 4 adaptive 0 jitter circular orient translate <40, 80, -40> }\n";
//        pov << "light_source { 1*x color rgb 0.7 area_light <8, 8, 0> <0, 8, 8>  4, 4 adaptive 0 jitter circular orient translate <-40, -80, 40> }\n";
        pov << "light_source {  100*<-20, 20, -20> color rgb <1,1, 1> }\n";

		//pov << "camera { location <0, 2, 0> look_at <0, 0, 0> }" << endl;
//pov << "camera { location <0, 0, -2> look_at <0, 0, 0> }" << endl;
pov << "camera { location 18.0*<0, 0, -25> look_at <0, 0, 0>   right x*image_width/image_height angle 0.3}" << std::endl;

     

		pov << "background {color rgb <0.9,0.9,0.9>}" << std::endl;


//mesh2

		for( unsigned int c = 0 ; c < cells.size() ; c++){

			if( cells[c]->type == 0 ){

			for( unsigned int i = 0 ; i < cells[c]->mass_huell.size() ; i++){
				cells[c]->mass_huell[i]->norm_F[0] = 0.;
				cells[c]->mass_huell[i]->norm_F[1] = 0.;
				cells[c]->mass_huell[i]->norm_F[2] = 0.;
			}	
	
//			pov << "zelle:\t" << c << endl;
//pov << "normalenvectoren:" << endl;

			for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){

				cells[c]->triangle[i]->points[0]->norm_F[0] += cells[c]->triangle[i]->mNormalvector[0]/3.;
				cells[c]->triangle[i]->points[0]->norm_F[1] += cells[c]->triangle[i]->mNormalvector[1]/3.;
				cells[c]->triangle[i]->points[0]->norm_F[2] += cells[c]->triangle[i]->mNormalvector[2]/3.;

/*
pov << i << ":\t" << cells[c]->triangle[i]->normalvector[0] ;
pov << "\t" << cells[c]->triangle[i]->normalvector[1];
pov << "\t" << cells[c]->triangle[i]->normalvector[2];
pov << "\t\t";
pov << "\t" << cells[c]->triangle[i]->points[0]->position[0];
pov << "\t" << cells[c]->triangle[i]->points[0]->position[1];
pov << "\t" << cells[c]->triangle[i]->points[0]->position[2];
pov << "\t" << cells[c]->triangle[i]->points[1]->position[0];
pov << "\t" << cells[c]->triangle[i]->points[1]->position[1];
pov << "\t" << cells[c]->triangle[i]->points[1]->position[2];
pov << "\t" << cells[c]->triangle[i]->points[2]->position[0];
pov << "\t" << cells[c]->triangle[i]->points[2]->position[1];
pov << "\t" << cells[c]->triangle[i]->points[2]->position[2];



pov << endl;
*/


			}









		pov << "mesh2 {" << std::endl;
		
		
		pov << "vertex_vectors {" << std::endl;
		pov << cells[c]->triangle.size()*3 << ","  << std::endl;
		for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){
			pov << "<" << cells[c]->triangle[i]->points[0]->position[0] << ","
					   << cells[c]->triangle[i]->points[0]->position[2] << ","
					   << cells[c]->triangle[i]->points[0]->position[1] << "> ";
			pov << "<" << cells[c]->triangle[i]->points[1]->position[0] << ","
					   << cells[c]->triangle[i]->points[1]->position[2] << ","
					   << cells[c]->triangle[i]->points[1]->position[1] << "> ";
			pov << "<" << cells[c]->triangle[i]->points[2]->position[0] << ","
					   << cells[c]->triangle[i]->points[2]->position[2] << ","
					   << cells[c]->triangle[i]->points[2]->position[1] << "> "
				<<std::endl;
		}
		pov << "}" << std::endl;


/*
		pov << "normal_vectors {" << endl;
		pov << cells[c]->triangle.size()*3 << ","  << endl;
		for( unsigned int i = 0 ; i < cells[0]->triangle.size() ; i++){
			pov << "<" << cells[c]->triangle[i]->points[0]->position[0] << ","
					   << cells[c]->triangle[i]->points[0]->position[2] << ","
					   << cells[c]->triangle[i]->points[0]->position[1] << "> ";
			pov << "<" << cells[c]->triangle[i]->points[1]->position[0] << ","
					   << cells[c]->triangle[i]->points[1]->position[2] << ","
					   << cells[c]->triangle[i]->points[1]->position[1] << "> ";
			pov << "<" << cells[c]->triangle[i]->points[2]->position[0] << ","
					   << cells[c]->triangle[i]->points[2]->position[2] << ","
					   << cells[c]->triangle[i]->points[2]->position[1] << "> "
				<<endl;
		}
		pov << "}" << endl;
*/
		pov << "normal_vectors {" << std::endl;
		pov << cells[c]->triangle.size()*3 << ","  << std::endl;
		for( unsigned int i = 0 ; i < cells[0]->triangle.size() ; i++){
			pov << "<" << cells[c]->triangle[i]->points[0]->norm_F[0] << ","
					   << cells[c]->triangle[i]->points[0]->norm_F[2] << ","
					   << cells[c]->triangle[i]->points[0]->norm_F[1] << "> ";
			pov << "<" << cells[c]->triangle[i]->points[1]->norm_F[0] << ","
					   << cells[c]->triangle[i]->points[1]->norm_F[2] << ","
					   << cells[c]->triangle[i]->points[1]->norm_F[1] << "> ";
			pov << "<" << cells[c]->triangle[i]->points[2]->norm_F[0] << ","
					   << cells[c]->triangle[i]->points[2]->norm_F[2] << ","
					   << cells[c]->triangle[i]->points[2]->norm_F[1] << "> "
				<<std::endl;
		}
		pov << "}" << std::endl;

		pov << "face_indices {" << std::endl;
		pov << cells[c]->triangle.size() << "," << std::endl;
		int f = 0;
		for( unsigned int i = 0 ; i < cells[c]->triangle.size() ; i++){
			pov << "<" << f << ","; f++;
			pov <<		  f << ","; f++;
			pov <<		  f << ">"; f++;
			pov << std::endl;
		}
		pov << "}" << std::endl;


		pov << "texture {" << std::endl;
		pov << "\t pigment { color rgb<255./255., 115./255., 0> }" << std::endl;
//		pov << "\t finish { ambient 0.2 diffuse 0.7 } }" << endl;
		pov << "\t finish {ambient 0.1 diffuse 0.9   phong 1} } no_shadow" << std::endl;
	
		
		pov << "}" << std::endl;

		}
		}

	

}
void Model3D::color_triangle_deform(){
	for( unsigned int i = 0 ; i < cells.size() ; i++ )
		cells[i]->color_triangle_deform();
}
void Model3D::output(){

  ofstream F;
  string outputFileName = outputPath + output_prefix + "_F_t.dat";
  F.open( outputFileName.c_str(), ios::out|ios::app );
  //F << time << "\t" << cells[0]->spring[101]->F[0] << "\t" <<   cells[0]->spring[101]->F[1] << "\t" <<  cells[0]->spring[101]->F[2] <<  endl;

  double a = -(cells[0]->spring[101]->l0_init - cells[0]->spring[101]->dlength - cells[0]->spring[101]->l0);
  double b =  cells[0]->spring[101]->l0_init - cells[0]->spring[101]->l0;
  //cells[0]->spring[101]->F[0];

	F	<< time << "\t" 
		<< cells[0]->spring[101]->dlength << "\t" 
		<< cells[0]->spring[101]->l0  << "\t" 
		<< cells[0]->spring[101]->l0_init << "\t" 
		<< a << "\t"
		<< b << "\t"
		<< cells[0]->spring[101]->F <<  std::endl;

	F.close();

	
	//volume

	/*
	ofstream output;
	output.open("../output/volume.txt",ios::out|ios::app );



	
	output	<< time 
			<< "\t"  
			<< cells[0]->volume//this->mean_volume
			<< "\t" 
			<< this->cells[0]->v_reference 
			<< "\t" 
			<< 4. / 3. * PI * cells[0]->radius_cell * cells[0]->radius_cell * cells[0]->radius_cell *cells[0]->volume_correction
			<< "\t" 
			<< this->cells[0]->radius_cell
			<< "\t"
			<< PI
			<< "\t"
			<< 4/3 
			<< endl;

	output.close();
	*/

	//presure
	ofstream V;
    outputFileName = outputPath + output_prefix + "_V_t.dat";
	V.open( outputFileName.c_str(), ios::out|ios::app );

  V << time << "\t" << volume_of_all_cells    << "\t" << this->cells[0]->volume << "\t" << this->cells[0]->radius_cell*this->cells[0]->radius_cell*this->cells[0]->radius_cell*4./3.*M_PI << "\t" <<this->cells[0]->volume/this->cells[0]->volume_correction;
  
  if( this->cells.size() > 1 )
  V <<  "\t" << this->cells[1]->volume ;

  V <<  std::endl;
	V.close();




	ofstream P;
    outputFileName = outputPath + output_prefix + "_P_t.dat";
	P.open( outputFileName.c_str(),  ios::out|ios::app );


	double mean = 0.;

	for( unsigned int i = 0 ; i < this->cells.size() ; i++ ){
		mean += cells[i]->presure_inside;	
	}
	mean /= cells.size();


	P << time << "\t" << mean    <<  std::endl;

	P.close();




	//radius
	
	


	ofstream R;
    outputFileName = outputPath + output_prefix +"_R.dat";
    R.open(outputFileName.c_str(), ios::out|ios::app);

	double center_of_mass[3] = {0.,0.,0.};
	for( unsigned int i = 0 ; i < cells.size() ; i++){
		center_of_mass[0] += cells[i]->mass[0]->position[0];
		center_of_mass[1] += cells[i]->mass[0]->position[1];
		center_of_mass[2] += cells[i]->mass[0]->position[2];
	}
	center_of_mass[0] /= cells.size();
	center_of_mass[1] /= cells.size();
	center_of_mass[2] /= cells.size();

	double r_gyr[3] = {0.,0.,0.};
	for( unsigned int i = 0 ; i < cells.size() ; i++){
		r_gyr[0] += (cells[i]->mass[0]->position[0] - center_of_mass[0])*(cells[i]->mass[0]->position[0] - center_of_mass[0]);
		r_gyr[1] += (cells[i]->mass[0]->position[1] - center_of_mass[1])*(cells[i]->mass[0]->position[1] - center_of_mass[1]);
		r_gyr[2] += (cells[i]->mass[0]->position[2] - center_of_mass[2])*(cells[i]->mass[0]->position[2] - center_of_mass[2]);
	}
	r_gyr[0] /= cells.size();
	r_gyr[1] /= cells.size();
	r_gyr[2] /= cells.size();

	r_gyr[0] = sqrt(r_gyr[0]);
	r_gyr[1] = sqrt(r_gyr[1]);
	r_gyr[2] = sqrt(r_gyr[2]);

	int left;
	int right;

	double max = 0.;

	for( unsigned int i = 0 ; i < (this->cells.size()-1) ; i++ ){
		for( unsigned int j = i+1 ; j < this->cells.size() ; j++){
			double tmp = dist(cells[i]->mass[0]->position,cells[j]->mass[0]->position);
			if ( tmp > max){
				left = i;
				right = j;
				max = tmp;
			}
		}
	}


	R << time << "\t" << this->cells.size() << "\t" << max+1    << "\t" << norm(r_gyr)<<  "\t" << log(timeStep) << "\t" <<std::endl;

	R.close();
	




	//special pressure
		if( mode == Pressure ){
	ofstream stretch;
    outputFileName = outputPath + output_prefix + "_pressure.dat";
	stretch.open( outputFileName.c_str(), ios::out|ios::app );

	stretch << time << "\t" << cells[0]->p << "\t" << cells[0]->radius_cell  << std::endl;

	stretch.close();

	}



		//optical stretcher
			if( mode == StretchCell ){
	ofstream stretch;
    outputFileName = outputPath + output_prefix + "_stretch.dat";
	stretch.open(outputFileName.c_str(),ios::out|ios::app );

	double stretch_long[3]= {	cells[0]->mass[ cells[0]->stretch_long_start]->position[0] - cells[0]->mass[ cells[0]->stretch_long_end]->position[0],
								cells[0]->mass[ cells[0]->stretch_long_start]->position[1] - cells[0]->mass[ cells[0]->stretch_long_end]->position[1],
								cells[0]->mass[ cells[0]->stretch_long_start]->position[2] - cells[0]->mass[ cells[0]->stretch_long_end]->position[2]};

	double stretch_short[3]= {	cells[0]->mass[ cells[0]->stretch_short_start]->position[0] - cells[0]->mass[ cells[0]->stretch_short_end]->position[0],
								cells[0]->mass[ cells[0]->stretch_short_start]->position[1] - cells[0]->mass[ cells[0]->stretch_short_end]->position[1],
								cells[0]->mass[ cells[0]->stretch_short_start]->position[2] - cells[0]->mass[ cells[0]->stretch_short_end]->position[2]};

	

	stretch << time << "\t" << norm(stretch_long) << "\t" << norm(stretch_short)  << std::endl;

	stretch.close();

	}




			//pushing two cells
			

	if( mode == PushTwoCells ){
	ofstream stretch;
    outputFileName = outputPath + output_prefix + "_push2cells.dat";
	stretch.open( outputFileName.c_str(), ios::out|ios::app );

	double push1[3]= {	cells[0]->mass[ cells[0]->stretch_long_start]->position[0] - cells[0]->mass[ cells[0]->stretch_long_end]->position[0],
						cells[0]->mass[ cells[0]->stretch_long_start]->position[1] - cells[0]->mass[ cells[0]->stretch_long_end]->position[1],
						cells[0]->mass[ cells[0]->stretch_long_start]->position[2] - cells[0]->mass[ cells[0]->stretch_long_end]->position[2]};

	double push2[3]= {	cells[1]->mass[ cells[1]->stretch_long_start]->position[0] - cells[1]->mass[ cells[1]->stretch_long_end]->position[0],
						cells[1]->mass[ cells[1]->stretch_long_start]->position[1] - cells[1]->mass[ cells[1]->stretch_long_end]->position[1],
						cells[1]->mass[ cells[1]->stretch_long_start]->position[2] - cells[1]->mass[ cells[1]->stretch_long_end]->position[2]};

	double push_g[3]= {	cells[0]->mass[ cells[0]->stretch_long_end]->position[0] - cells[1]->mass[ cells[1]->stretch_long_start]->position[0],
						cells[0]->mass[ cells[0]->stretch_long_end]->position[1] - cells[1]->mass[ cells[1]->stretch_long_start]->position[1],
						cells[0]->mass[ cells[0]->stretch_long_end]->position[2] - cells[1]->mass[ cells[1]->stretch_long_start]->position[2]};


	stretch << time << "\t" << norm(push1) << "\t" << norm(push2) << "\t" << norm(push_g) << "\t" << cells[0]->volume << "\t" << cells[1]->volume << std::endl;

	stretch.close();

	}



	ofstream timestream;
    outputFileName = outputPath + output_prefix + "_t.dat";
	timestream.open( outputFileName.c_str(), ios::out|ios::app );

	 std::time_t result = std::time(NULL);

	 timestream <<  time << "\t" << log(timeStep)/ log(2.) << "\t\t" << result << "\t" << this->cells[0]->status <<"\t" <<   std::asctime(std::localtime(&result)) ;//  <<  endl;

	timestream.close();

	//povray

	if( enablePovrayOutput ){

		currentTimeSinceLastOutput += timeStep;
		
		if( currentTimeSinceLastOutput >= timeBetweenOutputs){


      this->writePovray_all();

      /*
      

			string filename;
			filename.append(outputPath+"/pov/");//"../output/pov/");
			filename.append("complexCells"+output_prefix+"_");

			stringstream NumberString;
			NumberString.precision(10);
			NumberString << fixed << countpov/std::pow(10.,10) ;   
			filename.append("1"+NumberString.str().substr(2, 100));
			countpov++;

			filename.append(".pov");



            // todo:  bool parameter use_zlib: "Compress Files".
           
			if ( compress_pov )
            {
				#if (!WIN32)
                filename.append(".gz");

                gzFile zipFile = gzopen(filename.c_str(), "wbt");
				if( this->print_sinusoids == 2 ){
					    std::stringstream pov;
					pov << "#include \"sinusoids.pov\" " << endl;
					 gzwrite( zipFile, pov.str().c_str(), pov.str().size() );

				}
                writePovray_details2_gzip( zipFile );
				if( this->print_sinusoids == 2 ){
					std::stringstream pov;
					pov << "object{ sinusoids      }" << endl;
					 gzwrite( zipFile, pov.str().c_str(), pov.str().size() );
				}

                // gzwrite(zipFile, content.str().c_str(), content.str().size());

                gzclose(zipFile);
				#endif
            }
            else
            {
            
                ofstream pov;

                pov.open(filename.c_str(),ios::out );

               	if( print_sinusoids == 2 )
					pov << "#include \"sinusoids.pov\" " << endl;
    
			
                this->writePovray_details(pov);

							if( print_sinusoids == 2 )
					pov << "object{ sinusoids      }" << endl;


                pov.close();
            }

      */




			currentTimeSinceLastOutput =0.;




		}	
	}





}
/*
// set GUI parameter
void Model3D::setConserveVolume(bool box){
  conserve_volume = box;
}
*/
//set observation points
void Model3D::setObservationPoints(){

  double min_x = cells[0]->mass[0]->position[0];
  int min_x_pos;
  double max_x = cells[0]->mass[0]->position[0];
  int max_x_pos;

  for( unsigned int i = 0 ; i < cells[0]->mass.size() ; i++ ){
    if( cells[0]->mass[i]->position[0] < min_x){
      min_x_pos = i;
      min_x = cells[0]->mass[i]->position[0];
    }
    if( cells[0]->mass[i]->position[0] > max_x){
      max_x_pos = i;
      max_x = cells[0]->mass[i]->position[0];
    }
  }//for

  cells[0]->stretch_long_start = min_x_pos;
  cells[0]->stretch_long_end = max_x_pos;	

  cells[1]->stretch_long_start = min_x_pos;
  cells[1]->stretch_long_end = max_x_pos;	
}

//set modeof model
void Model3D::setMode(int i){
  mode = i;
}

//set initial values
void Model3D::setDefaultCell()
{
  CellTriangulated *c = new CellTriangulated();

    c->eta = 40.;
    c->gamma = 0.8;

    c->position.x = 0;
    c->position.y = 0;
    c->position.z = 0;
    c->radius_cell = defaultInitialCellRadius;
    c->k_hull = default_k_huell;
    c->k2 = default_k2;
    c->nu_hull = default_nu_huell;
    c->k_cytoskelett = default_k_cytoskeleton;
    c->nu_cytoskelett = default_nu_cytoskeleton;
    c->point_mass = default_PointMass;
    c->status = mode;
    c->dim = dimension;
    c->nu_medium = nu_medium;
    c->p = default_pressure_inside;
    c->mass_number = default_Mass_Number;
    c->initialRadius = defaultInitialCellRadius;
    c->relax_time = 1.;
    c->cycleTime = c->relax_time;
    c->eta_damper = default_damper;
    c->triangulate();

    c->pressure_threshold = pressure_threshold;
    c->type = 0;

    c->setMatrixA();

    c->mpRandom = &mRandom;

  defaultCellTriangulated = c;



}



// gives the default amount of mass points
int  Model3D::getDefaultNumberOfMassPoints(){
  return default_Mass_Number;
}

/*
 * Add a new complex cell - a default triangulated cell must be build before
 */
void Model3D::AddCell(double x, double y, double z, double r)
{
  CellTriangulated *c = new CellTriangulated();
  c->copyfrom(defaultCellTriangulated,x,y,z,r);
  
  // Set initial values
  c->initialRadius = defaultInitialCellRadius;
  c->dim = dimension;
  c->status = mode;

  // Add newly created cell to population
  cells.push_back(c);
  cells2->add(c);

  // Add cell also to area for visualisation
  mpArena->addObject( c->GLObject() );
}

//experiments
void Model3D::setStretchPoints(){


	int s = mRandom.GetRandomUniform01()*cells[0]->mass_huell.size();

	cells[0]->stretch_long_start = s;

	double max_value = 0;

	for( unsigned int i = 0 ; i < cells[0]->mass_huell.size() ; i++){

		double x_tmp[3] ={	cells[0]->mass_huell[i]->position[0] - cells[0]->mass_huell[s]->position[0],
							cells[0]->mass_huell[i]->position[1] - cells[0]->mass_huell[s]->position[1],
							cells[0]->mass_huell[i]->position[2] - cells[0]->mass_huell[s]->position[2]};

		double max_value_tmp = norm(x_tmp);

		if( max_value_tmp > max_value){
			max_value = max_value_tmp;
			cells[0]->stretch_long_end = i;
		}

	}//for


	double x[3] ={	cells[0]->mass_huell[s]->position[0] - cells[0]->mass[0]->position[0],
					cells[0]->mass_huell[s]->position[1] - cells[0]->mass[0]->position[1],
					cells[0]->mass_huell[s]->position[2] - cells[0]->mass[0]->position[2]};
		
	normvector(x);

/*
	cells[0]->mass_huell[s]->F_force_profile[0] = x[0];
 	cells[0]->mass_huell[s]->F_force_profile[1] = x[1];
	cells[0]->mass_huell[s]->F_force_profile[2] = x[2];

	cells[0]->mass_huell[cells[0]->strech_end]->F_force_profile[0] = -x[0];
	cells[0]->mass_huell[cells[0]->strech_end]->F_force_profile[1] = -x[1];
	cells[0]->mass_huell[cells[0]->strech_end]->F_force_profile[2] = -x[2];
*/

	double short_min = 1000;

	//calc force profile for each point
	for( unsigned i = 0 ; i < cells[0]->mass_huell.size() ; i++){

		double x_tmp[3] = {	cells[0]->mass_huell[i]->position[0] - cells[0]->mass[0]->position[0],
							cells[0]->mass_huell[i]->position[1] - cells[0]->mass[0]->position[1],
							cells[0]->mass_huell[i]->position[2] - cells[0]->mass[0]->position[2]};
		
		normvector(x_tmp);


		double cos_alpha = x_tmp[0]*x[0]+x_tmp[1]*x[1]+x_tmp[2]*x[2];
		double cos_alpha_q = cos_alpha*cos_alpha;


		double scalar = cos_alpha_q*1.5;



		cells[0]->mass_huell[i]->F_force_profile[0] = scalar * x_tmp[0]*0.5;//50;
		cells[0]->mass_huell[i]->F_force_profile[1] = scalar * x_tmp[1]*0.5;//50;
		cells[0]->mass_huell[i]->F_force_profile[2] = scalar * x_tmp[2]*0.5;//50;

		if( scalar < short_min ){

			cells[0]->stretch_short_start = i;
			short_min = scalar;
		}

	}
	
	
	
	max_value = 0;

	for( unsigned int i = 0 ; i < cells[0]->mass_huell.size() ; i++){

		double x_tmp[3] ={	cells[0]->mass_huell[i]->position[0] - cells[0]->mass_huell[cells[0]->stretch_short_start]->position[0],
							cells[0]->mass_huell[i]->position[1] - cells[0]->mass_huell[cells[0]->stretch_short_start]->position[1],
							cells[0]->mass_huell[i]->position[2] - cells[0]->mass_huell[cells[0]->stretch_short_start]->position[2]};

		double max_value_tmp = norm(x_tmp);

		if( max_value_tmp > max_value){
			max_value = max_value_tmp;
			cells[0]->stretch_short_end = i;
		}

	}//for





}
void Model3D::statusDepend(){
  for( unsigned int i = 0 ; i < cells.size() ; i++){

    this->cells[i]->statusDepend(timeStep , time);

    if( mode == WatershedFit ) {
      double force = 400;
      for( unsigned int j = 0 ; j < cells[i]->mass.size() ; j++ ){
        if( cells[i]->mass[j]->position[0] < this->x_left ){
          cells[i]->mass[j]->F[0] += force;
        }else if( cells[i]->mass[j]->position[0] > this->x_right ){
          cells[i]->mass[j]->F[0] -= force;
        }
        if( cells[i]->mass[j]->position[1] < this->y_left ){
          cells[i]->mass[j]->F[1] += force;
        }else if( cells[i]->mass[j]->position[1] > this->y_right ){
          cells[i]->mass[j]->F[1] -= force;
        }
        if( cells[i]->mass[j]->position[2] < this->z_left ){
          cells[i]->mass[j]->F[2] += force;
        }else if( cells[i]->mass[j]->position[2] > this->z_right ){
          cells[i]->mass[j]->F[2] -= force;
        }
      }
      cells[i]->mass[0]->position[1] = 0.;
      cells[i]->mass[0]->position[2] = 0.;
    }

/*	
	if( mode == 200 && !conserve_volume){
		double V = 0.;
		for( unsigned int l = 0 ; l < cells[i]->triangle.size() ; l++){
			cells[i]->triangle[l]->calcVolume(cells[i]->mass[0]->position);
			V += cells[i]->triangle[l]->Volume;
		}
		cells[i]->volume = V;

	}
	*/
  if( mode == Wall ){
    for( unsigned int i = 0 ; i < this->cells[0]->mass.size() ; i++){
      this->cells[0]->mass[i]->F[0] += this->force_push_two_cells;

      // this->mass_huell[i]->F[0] -= pow(0.0455/tmp_dist,12)-2*pow(0.0455/tmp_dist,6);
      this->cells[0]->mass[i]->F[0] -= pow(0.0455/(1 - this->cells[0]->mass[i]->position[0]),12);

      this->cells[0]->mass[i]->F_repulse[0] = pow(0.0455/(1 - this->cells[0]->mass[i]->position[0]),12);

      double a = 0;
    }
  }

  if( mode == PushTwoCells ){
    double dir[3];
    dir[0] = cells[0]->mass[0]->position[0] - cells[1]->mass[0]->position[0];
    dir[1] = cells[0]->mass[0]->position[1] - cells[1]->mass[0]->position[1];
    dir[2] = cells[0]->mass[0]->position[2] - cells[1]->mass[0]->position[2];

    normvector(dir);
    double f = this->force_push_two_cells/this->cells[0]->mass_huell.size();

    dir[0] *= f;
    dir[1] *= f;
    dir[2] *= f;

    for( unsigned int i = 0 ; i < this->cells[0]->mass_huell.size() ; i++){
/*
      cells[0]->mass_huell[i]->F[0] -= dir[0];
      cells[0]->mass_huell[i]->F[1] -= dir[1];
      cells[0]->mass_huell[i]->F[2] -= dir[2];

      cells[1]->mass_huell[i]->F[0] += dir[0];
      cells[1]->mass_huell[i]->F[1] += dir[1];
      cells[1]->mass_huell[i]->F[2] += dir[2];
*/

    	cells[0]->mass_huell[i]->F[0] -= f;
    	cells[1]->mass_huell[i]->F[0] += f;

    }


/*			cells[0]->mass[0]->F[0] -= dir[0]*f;
			cells[0]->mass[0]->F[1] -= dir[1]*f;
			cells[0]->mass[0]->F[2] -= dir[2]*f;
			cells[1]->mass[0]->F[0] += dir[0]*f;
			cells[1]->mass[0]->F[1] += dir[1]*f;
			cells[1]->mass[0]->F[2] += dir[2]*f;
*/


    cells[0]->mass[0]->position[1] = 0.;
    cells[0]->mass[0]->position[2] = 0.;
    cells[1]->mass[0]->position[1] = 0.;
    cells[1]->mass[0]->position[2] = 0.;

    }

  }

}
void Model3D::updateCellCycle(){
  for( unsigned int i = 0 ; i < cells.size() ; i++){
    cells[i]->updateCellCycle( timeStep , time);

    if( cells[i]->status > 100 ){

      if( cells[i]->can_divide )
        cutCells(i);
      if( cells[i]->status == 102 )
        cells[i]->time_relax+=timeStep;
    }//if regular cell cycle
  }//for all cells
}
//Model commands
void Model3D::Reset(){
std::cout << "here 100" << endl;
}
void Model3D::Simulate(int i){

  for( int j = 0 ; j<i ; j++ )
    this->Simulate();

}
void Model3D::Simulate(){


  int moveMode = 11;//2
  //1 Euler,2 Velocity Verlet(old)//3 VV (neu) // 4 Beeman // 11 first order - euler


//  this->storeLastStep();

  //clear everything
//  this->ClearAbleToMove();
  this->clearForce();

  //calc everythink
  this->calc();

  //do status depended thinks
  this->statusDepend();



  this->move(moveMode);//double move_max =
//  this->check_overlap();
/*
  //move_max = this->move(2);
  double a = defaultCellTriangulated->meanSpringLength*0.05;
  
//  while( (move_max > a) || (this->able_to_move == 0)){
  while( 0 == 1){
    this->ClearAbleToMove();

    F.open( outputFileName.c_str(), ios::out|ios::app );
    F <<"\t" ;//<< log(timeStep)/ log(2.) << endl; time <<
    F.close();

    this->output();

    //remove last step
    this->removeLastStep();

    //half of timestep
    this->timeStep *= 0.5;

    this->clearForce();
    this->statusDepend();

    //new move
    move_max = this->move(moveMode);

    bool co1 = this->check_overlap();
    bool co2 = this->check_overlap_past();

    bool sum = co1+co2;

  }//while
*/
  //after finaly one step move

  this->output();

  this->updateCellCycle();

  //set up output file
  this->timeBetweenOutputs = 0.01;
  this->enablePovrayOutput = true;
  if( enablePovrayOutput ){
    this->timeSinceLastOutput += this->timeStep;
    if( this->timeSinceLastOutput >= this->timeBetweenOutputs ){
      this->timeSinceLastOutput = 0.;
      this->mCounterOutput++;
      ofstream vtk;
      unsigned int countPrefix = this->mCounterOutput + 1000000;
      std::ostringstream outStream;
      outStream << countPrefix;
      string vtkOutputFileName = outputPath + output_prefix + "cells" + outStream.str() + ".vtk";
      vtk.open( vtkOutputFileName.c_str(), ios::out|ios::app );
      this->writeVTK( vtk );
      vtk.close();
    }
  }



  //for camera move
  //this->pov_cam_z += 5.*timeStep;

  this->time += this->timeStep;

  // if move small enough double timestep
//  if( move_max < defaultCellTriangulated->meanSpringLength*0.005 && (log(timeStep)/ log(2.) < -7.)   )
//   this->timeStep *= 2.;

}


//calc
void Model3D::calcForce(){
	for( unsigned int i = 0 ; i < cells.size() ; i++)
		cells[i]->calcForce(timeStep);
}
void Model3D::calcLength(){
  for( unsigned int i = 0 ; i < cells.size() ; i++)
    cells[i]->calcLength();
}
void Model3D::calcLengthWithNormDirection(){
  for( unsigned int i = 0 ; i < cells.size() ; i++)
    cells[i]->calcLengthWithNormDirection();
}
void Model3D::calcAcceleration(int number){
	for( unsigned int i = 0 ; i < cells.size() ; i++ )
			cells[i]->calcAcceleration(number);

}
void Model3D::calcVolume(){


  volume_of_all_cells = 0.;
  mean_volume = 0.;
  for( unsigned int j = 0 ; j < cells.size() ; j++){

    double V = cells[j]->calcVolume();//sums the volume of all pyramids
    volume_of_all_cells += V;
    mean_volume += V/cells.size();
		//*cells[j]->volume_correction;
/*
				double V = 0.;
		for( unsigned int i = 0 ; i < cells[0]->triangle.size() ;i++ ){
			cells[0]->triangle[i]->calcVolume(cells[0]->mass[0]->position);
			V += cells[0]->triangle[i]->Volume;
		}
		V_all += V/cells.size();
		cells[j]->volume = V;
		*/
  }
}
void Model3D::calc(){

  //Length of Springs
  this->calcLength();

  //set normalvectors
  for( unsigned int i = 0 ; i < this->cells.size() ; i++ ){
    for( unsigned int j = 0 ; j < this->cells[i]->triangle.size() ; j++ ){
      this->cells[i]->triangle[j]->setNormalVector();
    }
  }
  //Volume
  this->calcVolume();

  //angles in triangles
  for( unsigned int i = 0 ; i < this->cells.size() ; i++ ){
    for( unsigned int j = 0 ; j < this->cells[i]->triangle.size() ; j++ ){
      this->cells[i]->triangle[j]->mAngels[0] = angle(this->cells[i]->triangle[j]->springs[2]->direction,this->cells[i]->triangle[j]->springs[0]->direction);
      this->cells[i]->triangle[j]->mAngels[1] = angle(this->cells[i]->triangle[j]->springs[0]->direction,this->cells[i]->triangle[j]->springs[1]->direction);
      this->cells[i]->triangle[j]->mAngels[2] = angle(this->cells[i]->triangle[j]->springs[1]->direction,this->cells[i]->triangle[j]->springs[2]->direction);
    }
  }

  //calc voronoi area
  for( unsigned int i = 0 ; i < this->cells.size() ; i++ ){
    for( unsigned int j = 0 ; j < this->cells[i]->mass_huell.size() ; j++){

      this->cells[i]->mass_huell[j]->mAreaVoronoi = 0.;
      //sums over all remaining triangles
      for( unsigned int k = 0 ; k < (this->cells[i]->mass_huell[j]->mvpTriangles.size()-1) ; k++ ){
          double tmp = ( 1./( this->cells[i]->mass_huell[j]->mvpTriangles[k]->mAngels[ this->cells[i]->mass_huell[j]->mvpPointsL[k] ] )
        + 1./( this->cells[i]->mass_huell[j]->mvpTriangles[k+1]->mAngels[ this->cells[i]->mass_huell[j]->mvpPointsR[k+1] ] ))
        * this->cells[i]->mass_huell[j]->mvpSprings[ this->cells[i]->mass_huell[j]->mvpPointsR[k] ]->length
        * this->cells[i]->mass_huell[j]->mvpSprings[ this->cells[i]->mass_huell[j]->mvpPointsR[k] ]->length;
           this->cells[i]->mass_huell[j]->mAreaVoronoi += tmp;

      }
      this->cells[i]->mass_huell[j]->mAreaVoronoi /= 8.;
      double a = 0;
    }
  }

}

/*
 * store and remove timesteps by using dynamic step size
 */
void Model3D::removeLastStep(){
  for( unsigned int i = 0 ; i < cells.size() ; i++){
    cells[i]->removeLastStep();
    
 //   cells[i]->able_to_move = 1;
  }

  this->time = store_time;
//  this->store_time = this->store_time_old_old;

//  this->able_to_move = 1;

  /*
  // remove 1 step new cells
  for( unsigned int j = 0 ; j < this->mStoreNewCells.size() ; j++ ){
    cells2->remove(this->mStoreNewCells[j]);
    mpArena->removeObject( this->mStoreNewCells[j]->GLObject() );

    int tmp = 0;
    BoundingBoxList::iterator cellIterator = cells2->begin();
    while(dynamic_cast<CellTriangulated *>(* cellIterator) != this->mStoreNewCells[j] ){
      tmp++;
      ++cellIterator;
    }
    std::vector<CellTriangulated*>::iterator it = cells.begin()+tmp;
    cells.erase(it);

  }

  this->mStoreNewCells.clear();

  // add 1 step before removed cells
  for( unsigned int i = 0 ; i < this->mStoreOldCells.size() ; i++ ){
    cells.push_back(this->mStoreOldCells[i]);
    mpArena->addObject( this->mStoreOldCells[i]->GLObject() );
    cells2->add(this->mStoreOldCells[i]);
  }

  this->mStoreOldCells.clear();
  */

}
void Model3D::storeLastStep(){
  for( unsigned int i = 0 ; i < cells.size() ; i++)
    cells[i]->storeLastStep();

 // this->store_time_old_old = this->store_time_old;
 // this->store_time_old = this->store_time;
  this->store_time = time;

}


/*
 * add and calc forces
 */
void Model3D::calcForceOfDeformation(){
		this->calcLength();		//new length of springs
		this->calcForce();		//calculate force
		
		this->calcVolume();
	
		int i = 0;
			double V_diff = cells[i]->v_reference - cells[i]->volume;


			cells[i]->presure_inside = -V_diff;

			V_diff *= 500;
		
			for( unsigned int l = 0 ; l < cells[i]->triangle.size() ; l++){

				cells[i]->triangle[l]->points[0]->F[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
				cells[i]->triangle[l]->points[0]->F[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
				cells[i]->triangle[l]->points[0]->F[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;

				cells[i]->triangle[l]->points[1]->F[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
				cells[i]->triangle[l]->points[1]->F[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
				cells[i]->triangle[l]->points[1]->F[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;

				cells[i]->triangle[l]->points[2]->F[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
				cells[i]->triangle[l]->points[2]->F[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
				cells[i]->triangle[l]->points[2]->F[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;



			}

			double F_max = 0.;
			for( unsigned int i = 0 ;  i < cells[0]->mass_huell.size() ; i++){
				double F_tmp = norm( cells[0]->mass_huell[i]->F);
				if( F_tmp>F_max)
					F_max = F_tmp;
			}

			double center[3] = {cells[0]->mass[0]->position[0],cells[0]->mass[0]->position[1],cells[0]->mass[0]->position[2]};
			for( unsigned int i = 0 ; i < cells[0]->mass_huell.size() ; i++){
				
				double point[3]= {cells[0]->mass_huell[i]->position[0],cells[0]->mass_huell[i]->position[1],cells[0]->mass_huell[i]->position[2]};
				double dir[3] = {point[0]-center[0],point[1]-center[1],point[2]-center[2]};
				normvector(dir);
				
				double length = norm(cells[0]->mass_huell[i]->F);
				/*
				cells[0]->mass_huell[i]->F_force_profile[0] = dir[0] * cells[0]->mass_huell[i]->F[0];
				cells[0]->mass_huell[i]->F_force_profile[1] = dir[1] * cells[0]->mass_huell[i]->F[1];
				cells[0]->mass_huell[i]->F_force_profile[2] = dir[2] * cells[0]->mass_huell[i]->F[2];
				*/
				/*
				cells[0]->mass_huell[i]->F_force_profile[0] = dir[0] * length;
				cells[0]->mass_huell[i]->F_force_profile[1] = dir[1] * length;
				cells[0]->mass_huell[i]->F_force_profile[2] = dir[2] * length;
				*/
				cells[0]->mass_huell[i]->F_force_profile[0] = -cells[0]->mass_huell[i]->F[0]/F_max;
				cells[0]->mass_huell[i]->F_force_profile[1] = -cells[0]->mass_huell[i]->F[1]/F_max;
				cells[0]->mass_huell[i]->F_force_profile[2] = -cells[0]->mass_huell[i]->F[2]/F_max;

			}
}
void Model3D::addConstantForce(double x, double y, double z){
  for( unsigned int i = 0 ; i < cells.size() ; i++ )
    cells[i]->addConstantForce(x,y,z);
}
void Model3D::clearForce(){
  for( unsigned int i = 0 ; i < this->cells.size() ; i++)
    this->cells[i]->clearForce();
}


//additional damper
void Model3D::updateSecondDamper(double timeStep){
  for( unsigned int i = 0 ; i < cells.size() ; i++)
    cells[i]->updateSecondDamper(timeStep);
}



//update position and velocity
double Model3D::move(int number){

  double ret = 0.;

  this->ClearAbleToMove();

  if( number == 1){//euler

    this->calcLength();		//new length of springs
    this->calcForce();		//calculate force
    this->calcVolume();
    //this->statusDepend();	//status dependend action
    this->interactions2();	//interaction between cells
	//	this->updateDamper(timeStep);

    //move all points
    ret = this->updatePosition(number,timeStep);
    this->calcAcceleration(number);
    this->updateVelocity(number,timeStep);

  }else if( number == 2)	{//Velocity Verlet

    //new length of springs
    this->calcLength();
    //calculate force
    this->calcForce();

    //this->updateDamper(timeStep);
   this->calcVolume();
   //this->statusDepend();

    //after direct move
    this->interactions2();

    //this->updateDamper(timeStep);
   //if( this->able_to_move == 1){

    this->calcAcceleration(number);

    //move all points
    this->updateVelocity(number,timeStep);
    ret = this->updatePosition(number,timeStep);
    //}
  }else if( number == 3){
			// new (hopfully right) velocity verlert
			
			this->updateVelocity(1, timeStep*0.5);
			ret = this->updatePosition(1, timeStep);

			this->calcLength();
			this->calcForce();
			this->statusDepend();
			this->calcVolume();
			this->interactions2();

			this->calcAcceleration(number);

			this->updateVelocity(1, timeStep*0.5);

		}else if( number == 4){

			this->calcLength();		//new length of springs
			this->calcForce();		//calculate force
		
			this->calcVolume();
			this->statusDepend();	//status dependend action
			this->interactions2();	//interaction between cells

			this->calcAcceleration(number);

			ret = this->updatePosition(number, timeStep);
			this->updateVelocity(number, timeStep);//prediction

			this->calcForce();		//calculate force
			this->calcVolume();
			this->interactions2();	//interaction between cells
			
			//a_new
			this->calcAcceleration(number);


			this->updateVelocity(4,timeStep);//correction


  }else if( number == 11){//first order systems with euler

 //   this->calcLengthWithNormDirection();
    this->calcVolume();
    this->interactions2();	//interaction between cells
    this->setAndSolveSystem();
    ret = this->updatePosition(1, this->timeStep);

  }

    for( unsigned int i = 0 ; i < this->cells.size() ; i++){
      if( this->cells[i]->able_to_move == 0 )
        this->able_to_move = 0;
    }


  return ret;
}
double Model3D::updatePosition(int number, double delta_t){

  double max = 0.;

  for( unsigned int i = 0 ; i < cells.size() ; i++ ){
    double tmp_max = this->cells[i]->updatePosition(number, delta_t);
    if (tmp_max > max)
      max = tmp_max;
  }
  return max;

}
void Model3D::updateVelocity(int number, double delta_t){
  for( unsigned int i = 0 ; i < cells.size() ; i++ ){
    cells[i]->updateVelocity(number, delta_t);
  //  if( this->cells[i]->able_to_move == 0 )
   // 	this->able_to_move = 0;
  }
}

//first order dgl system
void Model3D::setAndSolveSystem(){

  for( unsigned int i = 0 ; i < this->cells.size() ; i++)
    this->cells[i]->setAndSolveSystem();

}


// check if 2 cells are overlap
bool Model3D::check_overlap(){


    BoundingBoxList::iterator cellIterator = cells2->begin();

    CellTriangulated * cell_1, * cell_2;

    while( (cellIterator != cells2->end()))
    {
        cell_1 = dynamic_cast<CellTriangulated *>( *cellIterator );

        CSListContainer< unsigned long > interactingCells = cell_1->mIntersectingList;

        unsigned long * otherIterator = interactingCells.begin();


      while ( (otherIterator != interactingCells.end())  )//&& (this->able_to_move == 1)
      {
          cell_2 = dynamic_cast<CellTriangulated *>( cells2->element( *otherIterator ) );

          bool bc1 = cell_1->check_overlap(cell_2);
          bool bc2 = cell_2->check_overlap(cell_1);

          //check overlap
          if( bc1 == 1 || bc2 == 2 ){
            this->able_to_move = 0;
            this->able_cell1 = bc1;
            this->able_cell2 = bc2;
            return 1;
          }

      ++otherIterator;
      }

      ++cellIterator;
  }

    return 0;


}
bool Model3D::check_overlap_past(){

    
    BoundingBoxList::iterator cellIterator = cells2->begin();

    CellTriangulated * cell_1, * cell_2;

    while( (cellIterator != cells2->end()))
    {
        cell_1 = dynamic_cast<CellTriangulated *>( *cellIterator );

        CSListContainer< ModelElement *> interactingCells = cells2->intersectingListByElement( cell_1 );

        ModelElement ** otherIterator = interactingCells.begin();


      while ( (otherIterator != interactingCells.end())  )//&& (this->able_to_move == 1)
      {
          cell_2 = dynamic_cast<CellTriangulated *>( *otherIterator );

          //check overlap
          if( cell_1->check_overlap_past(cell_2) == 1 || cell_2->check_overlap_past(cell_1)){
          //  this->able_to_move = 0;
            return 1;
          }

      ++otherIterator;
      }

      ++cellIterator;
  }

    return 0;


}
void Model3D::ClearAbleToMove(){

  this->able_to_move = 1;
  for( unsigned int i = 0 ; i < this->cells.size(); i++)
    this->cells[i]->able_to_move = 1;


}

//force


//interaction
void Model3D::interactions(){

	for( unsigned int i = 0 ; i < cells.size() ; i++){
		for( unsigned int j = 0 ; j < cells.size() ; j++){ 
			if( i != j )
				cells[i]->calcInteractionTetrahedral(cells[j]);


		}

	
		if( conserve_volume  ){
			/*
			double V = 0.;
			for( unsigned int l = 0 ; l < cells[i]->triangle.size() ; l++){
				cells[i]->triangle[l]->calcVolume(cells[i]->mass[0]->position);
				V += cells[i]->triangle[l]->Volume;
			}

			cells[i]->volume = V;//   *cells[i]->volume_correction;
			*/
			//double V_diff = cells[i]->v_reference - V;
			

			double V_diff = cells[i]->v_reference - cells[i]->volume;


			cells[i]->presure_inside = -V_diff;

			V_diff *= 500;
		
			for( unsigned int l = 0 ; l < cells[i]->triangle.size() ; l++){

				cells[i]->triangle[l]->points[0]->F[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
				cells[i]->triangle[l]->points[0]->F[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
				cells[i]->triangle[l]->points[0]->F[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;

				cells[i]->triangle[l]->points[1]->F[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
				cells[i]->triangle[l]->points[1]->F[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
				cells[i]->triangle[l]->points[1]->F[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;

				cells[i]->triangle[l]->points[2]->F[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
				cells[i]->triangle[l]->points[2]->F[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
				cells[i]->triangle[l]->points[2]->F[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;



			}
		

		
		}
	}
}
void Model3D::interactions2(){

  cells2->update();

  BoundingBoxList::iterator cellIterator = cells2->begin();

  this->able_to_move = 1;

  CellTriangulated * cell_1, * cell_2;

  while( (cellIterator != cells2->end()) ){
    cell_1 = dynamic_cast<CellTriangulated *>( *cellIterator );

    CSListContainer< unsigned long > interactingCells = cell_1->mIntersectingList;

    unsigned long * otherIterator = interactingCells.begin();

    while ( (otherIterator != interactingCells.end())  ){

      cell_2 = dynamic_cast<CellTriangulated *>( cells2->element( *otherIterator ) );

      if( cell_1->calcInteraction_repulsive(cell_2) == 0 || cell_2->calcInteraction_repulsive(cell_1) == 0)
        this->able_to_move = 0;

      ++otherIterator;
      }
    ++cellIterator;
  }

//TODO move to calc
conserve_volume=1;
  if( conserve_volume  ){
    for( unsigned int i = 0 ; i < cells.size() ; i++){

     double V_diff = cells[i]->v_reference - cells[i]->volume;

     //Debug Johannes
 //    fprintf(stderr, "%i\t: %2.3f \t %2.3f  %i \n", this->cells.size(), this->time , V_diff, this->mCounterOutput );

     cells[i]->presure_inside = V_diff;

     V_diff *= 400.0;

     for( unsigned int j = 0 ; j < this->cells[i]->mass.size() ; j++ ){
       this->cells[i]->mass[j]->F_repulse[0] = 0;
       this->cells[i]->mass[j]->F_repulse[1] = 0;
       this->cells[i]->mass[j]->F_repulse[2] = 0;
     }

     for( unsigned int l = 0 ; l < cells[i]->triangle.size() ; l++){

        cells[i]->triangle[l]->points[0]->F[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
        cells[i]->triangle[l]->points[0]->F[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
        cells[i]->triangle[l]->points[0]->F[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;

        cells[i]->triangle[l]->points[1]->F[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
        cells[i]->triangle[l]->points[1]->F[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
        cells[i]->triangle[l]->points[1]->F[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;

        cells[i]->triangle[l]->points[2]->F[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
        cells[i]->triangle[l]->points[2]->F[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
        cells[i]->triangle[l]->points[2]->F[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;

        cells[i]->triangle[l]->points[0]->F_repulse[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
        cells[i]->triangle[l]->points[0]->F_repulse[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
        cells[i]->triangle[l]->points[0]->F_repulse[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;

        cells[i]->triangle[l]->points[1]->F_repulse[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
        cells[i]->triangle[l]->points[1]->F_repulse[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
        cells[i]->triangle[l]->points[1]->F_repulse[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;

        cells[i]->triangle[l]->points[2]->F_repulse[0] += cells[i]->triangle[l]->mNormalvector[0]*V_diff;
        cells[i]->triangle[l]->points[2]->F_repulse[1] += cells[i]->triangle[l]->mNormalvector[1]*V_diff;
        cells[i]->triangle[l]->points[2]->F_repulse[2] += cells[i]->triangle[l]->mNormalvector[2]*V_diff;

      }//for all triangles
    }//for all cells


  }//if conserve volume


}


// cut 2 cells : remove (delete single cell) and creats two new cells
void Model3D::cutCells( int number ){

  CellTriangulated *original = cells[number];
  cells2->remove(cells[number]);

//  this->mStoreOldCells.push_back(cells[number]);

  mpArena->removeObject(cells[number]->GLObject() );
  std::vector<CellTriangulated*>::iterator it = cells.begin()+number;
  cells.erase(it);


  double midd_left[3] = {0.,0.,0.};
  double midd_right[3] = {0.,0.,0.};
  int midd_left_count = 0;
  int midd_right_count = 0;

  for( unsigned int i = 0 ; i < original->mass.size() ; i++ ){

    if( original->mass[i]->divide_number == 1 || original->mass[i]->divide_number == 4){
      midd_left[0] += original->mass[i]->position[0];
      midd_left[1] += original->mass[i]->position[1];
      midd_left[2] += original->mass[i]->position[2];
      midd_left_count++;
    }
    if( original->mass[i]->divide_number == 2 || original->mass[i]->divide_number == 4){
      midd_right[0] += original->mass[i]->position[0];
      midd_right[1] += original->mass[i]->position[1];
      midd_right[2] += original->mass[i]->position[2];
      midd_right_count++;
    }
  }//for

  midd_left[0] /= midd_left_count;
  midd_left[1] /= midd_left_count;
  midd_left[2] /= midd_left_count;
  midd_right[0] /= midd_right_count;
  midd_right[1] /= midd_right_count;
  midd_right[2] /= midd_right_count;

  //TODO calc max radius
  CellTriangulated *c1 = new CellTriangulated();
  c1->copyfrom(defaultCellTriangulated,midd_left[0],midd_left[1],midd_left[2],0.5);
  
  // Set initial radius
  c1->initialRadius = defaultInitialCellRadius;
  c1->status = 103;

  // Add newly created cell to population
  c1->projectToCell(original,1);
  cells.push_back(c1);
  mpArena->addObject( c1->GLObject() );

  CellTriangulated *c2 = new CellTriangulated();
  c2->copyfrom(defaultCellTriangulated,midd_right[0],midd_right[1],midd_right[2],0.5);
  // Set initial radius
  c2->initialRadius = defaultInitialCellRadius;
  c2->status = 103;

  // Add newly created cell to population
  c2->projectToCell(original,2);
  cells.push_back(c2);
  mpArena->addObject( c2->GLObject() );

  c1->mBoundingBox.xmax = midd_left[0];
  c1->mBoundingBox.ymax = midd_left[1];
  c1->mBoundingBox.zmax = midd_left[2];
  c1->mBoundingBox.xmin = midd_left[0];
  c1->mBoundingBox.ymin = midd_left[1];
  c1->mBoundingBox.zmin = midd_left[2];

  c2->mBoundingBox.xmax = midd_right[0];
  c2->mBoundingBox.ymax = midd_right[1];
  c2->mBoundingBox.zmax = midd_right[2];
  c2->mBoundingBox.xmin = midd_right[0];
  c2->mBoundingBox.ymin = midd_right[1];
  c2->mBoundingBox.zmin = midd_right[2];


  for( unsigned int i = 0 ; i < c1->mass.size() ; i++){
    if ( c1->mBoundingBox.xmin >  c1->mass[i]->position[0] )
      c1->mBoundingBox.xmin = c1->mass[i]->position[0];
    if ( c1->mBoundingBox.xmax <  c1->mass[i]->position[0] )
      c1->mBoundingBox.xmax = c1->mass[i]->position[0];
    if ( c1->mBoundingBox.ymin >  c1->mass[i]->position[1] )
      c1->mBoundingBox.ymin = c1->mass[i]->position[1];
    if ( c1->mBoundingBox.ymax <  c1->mass[i]->position[1] )
      c1->mBoundingBox.ymax = c1->mass[i]->position[1];
    if ( c1->mBoundingBox.zmin >  c1->mass[i]->position[2] )
      c1->mBoundingBox.zmin = c1->mass[i]->position[2];  
    if ( c1->mBoundingBox.zmax <  c1->mass[i]->position[2] )
      c1->mBoundingBox.zmax = c1->mass[i]->position[2];

    if ( c2->mBoundingBox.xmin >  c2->mass[i]->position[0] )
      c2->mBoundingBox.xmin = c2->mass[i]->position[0];
    if ( c2->mBoundingBox.xmax <  c2->mass[i]->position[0] )
      c2->mBoundingBox.xmax = c2->mass[i]->position[0];
    if ( c2->mBoundingBox.ymin >  c2->mass[i]->position[1] )
      c2->mBoundingBox.ymin = c2->mass[i]->position[1];
    if ( c2->mBoundingBox.ymax <  c2->mass[i]->position[1] )
      c2->mBoundingBox.ymax = c2->mass[i]->position[1];
    if ( c2->mBoundingBox.zmin >  c2->mass[i]->position[2] )
      c2->mBoundingBox.zmin = c2->mass[i]->position[2];
    if ( c2->mBoundingBox.zmax <  c2->mass[i]->position[2] )
      c2->mBoundingBox.zmax = c2->mass[i]->position[2];
  }


  	double tmp = 1./8.;
	c2->mBoundingBox.xmin -= tmp;
	c2->mBoundingBox.xmax += tmp;
	c2->mBoundingBox.ymin -= tmp;
	c2->mBoundingBox.ymax += tmp;
	c2->mBoundingBox.zmin -= tmp;
	c2->mBoundingBox.zmax += tmp;


	c1->mBoundingBox.xmin -= tmp;
	c1->mBoundingBox.xmax += tmp;
	c1->mBoundingBox.ymin -= tmp;
	c1->mBoundingBox.ymax += tmp;
	c1->mBoundingBox.zmin -= tmp;
	c1->mBoundingBox.zmax += tmp;
  /*
   for( unsigned int i = 0 ; i < c1->mass.size() ; i++){
     c1->mass[i]->position[0] -= 0.3;
     c2->mass[i]->position[0] += 0.3;
  }
  */

  
  double scale = 0.95;

  for( unsigned int i = 0 ;  i < c1->spring_cytoscelett.size() ; i++){

    double l1[3],l2[3];
    l1[0] = c1->spring_cytoscelett[i]->start->position[0]-c1->spring_cytoscelett[i]->end->position[0];
    l2[0] = c2->spring_cytoscelett[i]->start->position[0]-c2->spring_cytoscelett[i]->end->position[0];
    l1[1] = c1->spring_cytoscelett[i]->start->position[1]-c1->spring_cytoscelett[i]->end->position[1];
    l2[1] = c2->spring_cytoscelett[i]->start->position[1]-c2->spring_cytoscelett[i]->end->position[1];
    l1[2] = c1->spring_cytoscelett[i]->start->position[2]-c1->spring_cytoscelett[i]->end->position[2];
    l2[2] = c2->spring_cytoscelett[i]->start->position[2]-c2->spring_cytoscelett[i]->end->position[2];

    double l1_norm = norm(l1);
    double l2_norm = norm(l2);
    
    c1->spring_cytoscelett[i]->end->position[0] = c1->spring_cytoscelett[i]->start->position[0]+c1->spring_cytoscelett[i]->direction[0]*scale*l1_norm;
    c2->spring_cytoscelett[i]->end->position[0] = c2->spring_cytoscelett[i]->start->position[0]+c2->spring_cytoscelett[i]->direction[0]*scale*l2_norm;

    c1->spring_cytoscelett[i]->end->position[1] = c1->spring_cytoscelett[i]->start->position[1]+c1->spring_cytoscelett[i]->direction[1]*scale*l1_norm;
    c2->spring_cytoscelett[i]->end->position[1] = c2->spring_cytoscelett[i]->start->position[1]+c2->spring_cytoscelett[i]->direction[1]*scale*l2_norm;

    c1->spring_cytoscelett[i]->end->position[2] = c1->spring_cytoscelett[i]->start->position[2]+c1->spring_cytoscelett[i]->direction[2]*scale*l1_norm;
    c2->spring_cytoscelett[i]->end->position[2] = c2->spring_cytoscelett[i]->start->position[2]+c2->spring_cytoscelett[i]->direction[2]*scale*l2_norm;


    c1->spring_cytoscelett[i]->l0 = l1_norm*scale;
    c1->spring_cytoscelett[i]->l0_init = l1_norm*scale;

    c2->spring_cytoscelett[i]->l0 = l2_norm*scale;
    c2->spring_cytoscelett[i]->l0_init = l2_norm*scale;


  }
  c1->radius_cell *= scale;
  c2->radius_cell *= scale;
  c1->initialRadius *= scale;
  c2->initialRadius *= scale;
  for( unsigned int i = 0 ; i < c1->spring_huell.size() ; i++){
    
    double l1[3],l2[3];
    l1[0] = c1->spring_huell[i]->start->position[0]-c1->spring_huell[i]->end->position[0];
    l2[0] = c2->spring_huell[i]->start->position[0]-c2->spring_huell[i]->end->position[0];
    l1[1] = c1->spring_huell[i]->start->position[1]-c1->spring_huell[i]->end->position[1];
    l2[1] = c2->spring_huell[i]->start->position[1]-c2->spring_huell[i]->end->position[1];
    l1[2] = c1->spring_huell[i]->start->position[2]-c1->spring_huell[i]->end->position[2];
    l2[2] = c2->spring_huell[i]->start->position[2]-c2->spring_huell[i]->end->position[2];
    
    double l1_norm = norm(l1);
    double l2_norm = norm(l2);
    
    c1->spring_huell[i]->l0 = l1_norm;
    c1->spring_huell[i]->l0_init = l1_norm;

    c2->spring_huell[i]->l0 = l2_norm;
    c2->spring_huell[i]->l0_init = l2_norm;


  }
  

  cells2->add(c2);
  cells2->add(c1);
  /*
  double v1 = c1->calcVolume();
  double v2 = c2->calcVolume();

 c1->volume = v1;///this->defaultCellTriangulated->volume_correction;
  c2->volume = v2;///this->defaultCellTriangulated->volume_correction;
 c1->v_reference = v1;
 c2->v_reference = v2;
 */
  double r = 0.5;
 // c1->v_reference = v1;//4./3.*PI*(0.5*scale)*(0.5*scale)*(0.5*scale);
 // c2->v_reference = v2;//4./3.*PI*(0.5*scale)*(0.5*scale)*(0.5*scale);
  
  c1->volume_correction = 1.;
  c2->volume_correction = 1.;
  

 //   this->mStoreNewCells.push_back(c1);
 //  this->mStoreNewCells.push_back(c2);
  
  /*


  double shift = 0.04;

  for( unsigned int i = 0 ; i < c1->mass.size() ; i++ ){
    c1->mass[i]->position[0] -= shift;
    c2->mass[i]->position[0] += shift;

  }
  */

  for( unsigned int i = 0 ; i < c1->triangle.size() ; i++ ){
    c1->triangle[i]->setNormalVector();
    c2->triangle[i]->setNormalVector();
  }

}




/*
 * Simulation : Setup and simulion loop
 */
void Model3D::SetupSimulation(){

  if ( this->defaultCellTriangulated != NULL )
    return;

  mpSimulationThread->setUpdateInterval(1);

  timeStep = std::pow( 2., mLog2timeStep );

  std::cerr << "setting time step = " << timeStep << std::endl;

  if ( this->defaultCellTriangulated == NULL )
    this->setDefaultCell();

  if ( ! mpSimulationThread )
    mpSimulationThread = new QCSSimulationThread();

  switch ( mode ){
    case SingleCell:{
      setMode(102);
      AddCell( 0., 0., 0. );
      break;
    }
    case WatershedFit:{
      if( force_profile )

        setMode(2000);

        vector <vector <int> > cells_coord;//vector of vector of pixels
        vector <vector <int> > sinusoids_skeleton;
        vector <vector <int> > sinusoids;
        vector <int> sinusoids_skeleton_index;

        //input pixels of cells
        string line;
        string file =this->inputPath + this->watershedFitCells;
        ifstream myfile(file.c_str());

        if (myfile.is_open()){

          while ( myfile.good() ){
            getline (myfile,line);
            getline (myfile,line);

            vector <int> v_t;

            stringstream os(line);
            string temp;

            while (os >> temp) {//the stringstream makes temp a token
              int z = atoi(temp.c_str());
              v_t.push_back(z);
            }

            if( v_t.size() > 0)
              cells_coord.push_back(v_t);
            v_t.clear();

          }
        myfile.close();
      }

      else cout << "Unable to open file";

      //input pixels of cells
      file =this->inputPath + "sinusoids_skeleton.txt";
      ifstream myfile2(file.c_str());

      if (myfile2.is_open()){

        while ( myfile2.good() )
          {
            getline (myfile2,line);


            stringstream os2(line);
            string temp2;
            os2 >> temp2;
              int z;
              z = atoi(temp2.c_str());
			  sinusoids_skeleton_index.push_back(z);

            getline (myfile2,line);

            vector <int> v_t;

            stringstream os(line);
            string temp;

            while (os >> temp) {			  //the stringstream makes temp a token
              int z;
              z = atoi(temp.c_str());

              v_t.push_back(z);
            }

            if( v_t.size() > 0)
              sinusoids_skeleton.push_back(v_t);
            v_t.clear();

          }
        myfile2.close();
      }

    else cout << "Unable to open file";



	//input pixels of cells
	file =this->inputPath + "sinusoids.txt";
	ifstream myfile3(file.c_str());

    if (myfile3.is_open())
      {

        while ( myfile3.good() )
          {
            getline (myfile3,line);
            getline (myfile3,line);

            vector <int> v_t;

            stringstream os(line);
            string temp;

            while (os >> temp) {			  //the stringstream makes temp a token
              int z;
              z = atoi(temp.c_str());

              v_t.push_back(z);
            }

            if( v_t.size() > 0)
              sinusoids.push_back(v_t);
            v_t.clear();

          }
        myfile3.close();
      }

    else cout << "Unable to open file";












//end of reading pixels of cells


	/*
		double **sinusoids_skeleton_center = new double*[(int)sinusoids_skeleton.size()];
	 for( unsigned int i = 0 ; i < sinusoids_skeleton.size() ; i++)
      sinusoids_skeleton_center[i] = new double[3];
	 sinusoids_skeleton_center[i] = {0.,0.,0.};

	 	for( unsigned int j = 0 ; j < sinusoids_skeleton.size() ; j++ ){
			for( unsigned int j = 0 ; j <  sinusoids_skeleton[i].size() ; j=j+3 ){

			}
		}

		*/



	//translate to real x,y,z pixels
    int count_cells = (int) cells_coord.size();
    double **center_of_cells = new double*[count_cells];
    for( int i = 0 ; i < count_cells ; i++)
      center_of_cells[i] = new double[3];

	//calc mass of pixels
    double mass_center_of_all_cells[3] = {0.,0.,0.};

    for( unsigned int i = 0 ; i < cells_coord.size(); i++){

      double x = 0.;
      double y = 0.;
      double z = 0.;


      for( unsigned int j = 0 ; j <  cells_coord[i].size() ; j=j+3 ){

        x += cells_coord[i][j];
        y += cells_coord[i][j+1];
        z += cells_coord[i][j+2];

      }//over all voxels

	  x /= cells_coord[i].size()/3.;
	  y /= cells_coord[i].size()/3.;
	  z /= cells_coord[i].size()/3.;

      center_of_cells[i][0] = x;
      center_of_cells[i][1] = y;
      center_of_cells[i][2] = z;

      mass_center_of_all_cells[0] += x/count_cells;
      mass_center_of_all_cells[1] += y/count_cells;
      mass_center_of_all_cells[2] += z/count_cells;

    }//over all cells
	//end calc mass center of each cell


    double *norm_max_all = new double[count_cells];
    double max_of_all_norms = 0;

	//calc max length of each cell
    for( int i = 0 ; i < count_cells ; i++){

      double norm_max = 0.;
	  norm_max_all[i] = 0.;

	  //calc max from center of cell i to each point of this cell
      for( unsigned int j = 0 ; j <  cells_coord[i].size() ; j=j+3 ){

        double tmp[3] = {cells_coord[i][j]-center_of_cells[i][0], cells_coord[i][j+1]-center_of_cells[i][1], cells_coord[i][j+2]-center_of_cells[i][2]};
		double norm_tmp = norm(tmp);

        if( norm_tmp > norm_max )
          norm_max = norm_tmp;
      }//over all voxels


      norm_max_all[i] = norm_max;
      if( norm_max > max_of_all_norms)
        max_of_all_norms = norm_max;
		}






		 //output sinusoid pov file

				fstream c;

			string outputFileName = outputPath +  "pov/sinusoids.pov";
			c.open( outputFileName.c_str(), ios::out|ios::app );
			c << "			#declare  sinusoids =" << std::endl;
			c << "union{" << std::endl;
				for( unsigned int i = 0 ; i < sinusoids.size(); i++ ){

				c << "blob{  threshold 0.9 " << std::endl;
				for( unsigned int j = 0 ; j < sinusoids[i].size(); j=j+3 ){
				double x = sinusoids[i][j];
				double y = sinusoids[i][j+1];
				double z = sinusoids[i][j+2];
		c << "sphere { <" << x/max_of_all_norms-mass_center_of_all_cells[0]/max_of_all_norms << "," << z/max_of_all_norms-mass_center_of_all_cells[2]/max_of_all_norms << "," << y/max_of_all_norms-mass_center_of_all_cells[1]/max_of_all_norms << ">," << 1./max_of_all_norms << ",1}" << std::endl;

		}
		c <<" texture { pigment { color rgb<1,0,0>}	finish {ambient 0.1 diffuse 0.9   phong 1}} " ;
		c << "}" << std::endl;

		}
				c << "}" << std::endl;

	c.close();












	//add cells to model
	for(int i = 0; i < count_cells; i++){//1; i++){//
//for(int i = 0; i < 2; i++){//1; i++){//

      double x = center_of_cells[i][0];//-mass_center_of_all_cells[0]);///max_of_all_norms;
      double y = center_of_cells[i][1];//-mass_center_of_all_cells[1]);///max_of_all_norms;
      double z = center_of_cells[i][2];//-mass_center_of_all_cells[2]);///max_of_all_norms;
      //double r = pow(3./4./PI*cells_coord[i].size()/3.,1./3.);
	double r = norm_max_all[i];///max_of_all_norms;

	int tmp = (int) this->cells.size();
   	this->AddCell(x,y,z,r);
	this->cells[tmp]->type = 0;


		//from 1, 0 is nucleus
		for( unsigned int j = 1 ; j < this->cells[i]->mass.size() ; j++){

        for( unsigned int k = 0 ; k < cells_coord[i].size() ; k=k+3){




          double d[3];
          d[0] = this->cells[i]->mass[j]->position[0] - center_of_cells[i][0];
          d[1] = this->cells[i]->mass[j]->position[1] - center_of_cells[i][1];
          d[2] = this->cells[i]->mass[j]->position[2] - center_of_cells[i][2];
          normvector(d);


				double m[3];
				m[0] = cells_coord[i][k]   - center_of_cells[i][0];
				m[1] = cells_coord[i][k+1] - center_of_cells[i][1];
				m[2] = cells_coord[i][k+2] - center_of_cells[i][2];

				double t = dotmult(d,m);

				double tt = t*t;

				if( t >= 0 || tt <= 0.5 ){
					double m_  = tt-dotmult(m,m);




						if( m_*m_ <= 0.5 ){
							double a[3] = {static_cast<double>(cells_coord[i][k]),
                                           static_cast<double>(cells_coord[i][k+1]),
                                           static_cast<double>(cells_coord[i][k+2]) };
					this->cells[i]->mass[j]->position[0] = cells_coord[i][k];
					this->cells[i]->mass[j]->position[1] = cells_coord[i][k+1];
					this->cells[i]->mass[j]->position[2] = cells_coord[i][k+2];
			break;
						}
						else{
							if(this->cells[i]->mass[j]->position[0] < 0 ){
								double lambda = -center_of_cells[i][0]/d[0];
								this->cells[i]->mass[j]->position[0] = center_of_cells[i][0] + lambda*d[0];
								this->cells[i]->mass[j]->position[1] = center_of_cells[i][1] + lambda*d[1];
								this->cells[i]->mass[j]->position[2] = center_of_cells[i][2] + lambda*d[2];
							}//this->cells[i]->mass[j]->position[0] = 0;
							if(this->cells[i]->mass[j]->position[1] < 0 ){
								double lambda = -center_of_cells[i][1]/d[1];
								this->cells[i]->mass[j]->position[0] = center_of_cells[i][0] + lambda*d[0];
								this->cells[i]->mass[j]->position[1] = center_of_cells[i][1] + lambda*d[1];
								this->cells[i]->mass[j]->position[2] = center_of_cells[i][2] + lambda*d[2];
							}//this->cells[i]->mass[j]->position[1] = 0;
							if( this->cells[i]->mass[j]->position[2] < 0){
								double lambda = -center_of_cells[i][2]/d[2];
								this->cells[i]->mass[j]->position[0] = center_of_cells[i][0] + lambda*d[0];
								this->cells[i]->mass[j]->position[1] = center_of_cells[i][1] + lambda*d[1];
								this->cells[i]->mass[j]->position[2] = center_of_cells[i][2] + lambda*d[2];
							}//this->cells[i]->mass[j]->position[2] = 0;
							if(this->cells[i]->mass[j]->position[2] > 31){
								double lambda = (31-center_of_cells[i][2])/d[2];
								this->cells[i]->mass[j]->position[0] = center_of_cells[i][0] + lambda*d[0];
								this->cells[i]->mass[j]->position[1] = center_of_cells[i][1] + lambda*d[1];
								this->cells[i]->mass[j]->position[2] = center_of_cells[i][2] + lambda*d[2];
							}//this->cells[i]->mass[j]->position[2] = 31;


						}




					}






		}//k




		this->cells[i]->mass[j]->position[0] /= max_of_all_norms;
        this->cells[i]->mass[j]->position[0] -= mass_center_of_all_cells[0]/max_of_all_norms;
		this->cells[i]->mass[j]->position[1] /= max_of_all_norms;
 		this->cells[i]->mass[j]->position[1] -= mass_center_of_all_cells[1]/max_of_all_norms;
		this->cells[i]->mass[j]->position[2] /= max_of_all_norms;
		this->cells[i]->mass[j]->position[2] -= mass_center_of_all_cells[2]/max_of_all_norms;


			}

				this->cells[i]->mass[0]->position[0] /= max_of_all_norms;
				this->cells[i]->mass[0]->position[0] -= mass_center_of_all_cells[0]/max_of_all_norms;
				this->cells[i]->mass[0]->position[1] /= max_of_all_norms;
				this->cells[i]->mass[0]->position[1] -= mass_center_of_all_cells[1]/max_of_all_norms;
				this->cells[i]->mass[0]->position[2] /= max_of_all_norms;
				this->cells[i]->mass[0]->position[2] -= mass_center_of_all_cells[2]/max_of_all_norms;




		/*

			x = (sinusoids_skeleton[0][0]/max_of_all_norms - mass_center_of_all_cells[0]/max_of_all_norms);
			y = (sinusoids_skeleton[0][1]/max_of_all_norms - mass_center_of_all_cells[1]/max_of_all_norms);
			z = (sinusoids_skeleton[0][2]/max_of_all_norms - mass_center_of_all_cells[2]/max_of_all_norms);


		     this->AddCell(x,y,z,1./16.);

		*/





		for( unsigned int j = 0 ; j < cells_coord[i].size() ; j+=3 ){

			MassPoint *mp = new MassPoint();

			mp->position[0] = (cells_coord[i][j]/max_of_all_norms - mass_center_of_all_cells[0]/max_of_all_norms);
			mp->position[1] = (cells_coord[i][j+1]/max_of_all_norms - mass_center_of_all_cells[1]/max_of_all_norms);
			mp->position[2] = (cells_coord[i][j+2]/max_of_all_norms - mass_center_of_all_cells[2]/max_of_all_norms);

			mp->r = 0.6*255;
			mp->g = 0.6*255;
			mp->b = 0.6*255;

			mp->transparency = 0.1;

		//	this->cells[tmp]->mass.push_back(mp);

			//this->watershed.push_back(mp);//  ->cells[i]->mass.push_back(mp);
			}//over all voxels




			this->cells[tmp]->calcVolume();
			this->cells[tmp]->v_reference = this->cells[tmp]->volume;


			//V=4/3*pi*r^3
			//r=(V*3/4/pi)^(1/3)
			this->cells[tmp]->radius_cell = pow((this->cells[tmp]->volume*3./4./M_PI),(1./3.));

			double r_corr = this->cells[tmp]->radius_cell/r;
			for( unsigned int j = 0 ; j < this->cells[i]->spring.size() ; j++ ){

				this->cells[tmp]->spring[j]->l0 *= r_corr;
				this->cells[tmp]->spring[j]->l0_init = this->cells[tmp]->spring[j]->l0;
			}


			}//over cells

		for( unsigned int i = 0 ; i < sinusoids_skeleton.size() ; i++ ){
			for( unsigned int j = 0 ; j < sinusoids_skeleton[i].size() ; j=j+3 ){

			MassPoint *mp = new MassPoint();
			/*
			mp->position[0] = (sinusoids_skeleton[i][j]/max_of_all_norms - mass_center_of_all_cells[0]/max_of_all_norms);
			mp->position[1] = (sinusoids_skeleton[i][j+1]/max_of_all_norms - mass_center_of_all_cells[1]/max_of_all_norms);
			mp->position[2] = (sinusoids_skeleton[i][j+2]/max_of_all_norms - mass_center_of_all_cells[2]/max_of_all_norms);
		*/

			mp->position[0] = sinusoids_skeleton[i][j];
			mp->position[1] = sinusoids_skeleton[i][j+1];
			mp->position[2] = sinusoids_skeleton[i][j+2];

			mp->r = 0.0*255;
			mp->g = 0.0*255;
			mp->b = 0.0*255;

			mp->transparency = 1.0;

	//		this->cells[i]->mass.push_back(mp);

/*
		 int index;
			for( unsigned int k = 0 ; k < sinusoids_skeleton_index.size() ; k++){
				if( i == sinusoids_skeleton_index[k] ){
					index = k;
					break;
				}
			}
*/
		/*
				double dist[3];

				dist[0] = (sinusoids[i][0]/max_of_all_norms - mass_center_of_all_cells[0]/max_of_all_norms)-mp->position[0];
				dist[1] = (sinusoids[i][1]/max_of_all_norms - mass_center_of_all_cells[1]/max_of_all_norms)-mp->position[1];
				dist[2] = (sinusoids[i][2]/max_of_all_norms - mass_center_of_all_cells[2]/max_of_all_norms)-mp->position[2];

				double dist_norm_min = norm(dist);




			for( unsigned int l = 4 ; l < sinusoids[i].size() ; l=l+3){
				double dist[3];

				dist[0] = (sinusoids[i][l]/max_of_all_norms - mass_center_of_all_cells[0]/max_of_all_norms)-mp->position[0];
				dist[1] = (sinusoids[i][l+1]/max_of_all_norms - mass_center_of_all_cells[1]/max_of_all_norms)-mp->position[1];
				dist[2] = (sinusoids[i][l+2]/max_of_all_norms - mass_center_of_all_cells[2]/max_of_all_norms)-mp->position[2];

				double dist_norm = norm(dist);
				if( dist_norm < dist_norm_min)
					dist_norm_min = dist_norm;


			}
			*/


		int tmp = (int) this->cells.size();
		//this->AddCell(mp->position[0],mp->position[1],mp->position[2],dist_norm_min/max_of_all_norms);
	     this->AddCell(mp->position[0],mp->position[1],mp->position[2],1./4.*max_of_all_norms);
		 this->cells[tmp]->type = 1;




			for(unsigned int k = 1 ; k < this->cells[tmp]->mass.size() ; k++ ){

				        for( unsigned int l = 0 ; l < sinusoids[i].size() ; l=l+3){



							//hier weiter machen!!!!!!!!!!



			double s[3];
			/*
			s[0] = (sinusoids[i][l+0]/max_of_all_norms - mass_center_of_all_cells[0]/max_of_all_norms);
			s[1] = (sinusoids[i][l+1]/max_of_all_norms - mass_center_of_all_cells[1]/max_of_all_norms);
			s[2] = (sinusoids[i][l+2]/max_of_all_norms - mass_center_of_all_cells[2]/max_of_all_norms);
			*/
			s[0] = sinusoids[i][l+0];
			s[1] = sinusoids[i][l+1];
			s[2] = sinusoids[i][l+2];


          double d[3];
          d[0] = this->cells[tmp]->mass[k]->position[0] - this->cells[tmp]->mass[0]->position[0];
          d[1] = this->cells[tmp]->mass[k]->position[1] - this->cells[tmp]->mass[0]->position[1] ;
          d[2] = this->cells[tmp]->mass[k]->position[2] - this->cells[tmp]->mass[0]->position[2] ;
          double norm_d = normvector(d);


				double m[3];
				m[0] =  s[0] - this->cells[tmp]->mass[0]->position[0];
				m[1] =  s[1] - this->cells[tmp]->mass[0]->position[1];
				m[2] =  s[2] - this->cells[tmp]->mass[0]->position[2];
				double norm_m = norm(m);

				double t = dotmult(d,m);

				double tt = t*t;

				double diskriminante = tt-norm_m*norm_m+(1.5/max_of_all_norms)*(1.5/max_of_all_norms);


				if(t >= 0 || tt <= 0.5 ){
					double m_  = tt-dotmult(m,m);
	//			if( norm_d >= norm_m){
	//				if( diskriminante > 0){


						if( (m_*m_ <= 0.5) && (norm_d >norm_m) ){
//							double a[3] = {s[0], s[1], s[2]};
			/*
					this->cells[tmp]->mass[k]->position[0] = s[0];
					this->cells[tmp]->mass[k]->position[1] = s[1];
					this->cells[tmp]->mass[k]->position[2] = s[2];
			*/		/*
							this->cells[tmp]->mass[k]->r = 255;
							this->cells[tmp]->mass[k]->g = 255;
							this->cells[tmp]->mass[k]->b = 255;
						*/
			break;
						}


					}






		}//l
			}//all mass points
		for( unsigned int l = 0 ; l < this->cells[tmp]->mass.size() ; l++ ){




		this->cells[tmp]->mass[l]->position[0] /= max_of_all_norms;
        this->cells[tmp]->mass[l]->position[0] -= mass_center_of_all_cells[0]/max_of_all_norms;
		this->cells[tmp]->mass[l]->position[1] /= max_of_all_norms;
 		this->cells[tmp]->mass[l]->position[1] -= mass_center_of_all_cells[1]/max_of_all_norms;
		this->cells[tmp]->mass[l]->position[2] /= max_of_all_norms;
		this->cells[tmp]->mass[l]->position[2] -= mass_center_of_all_cells[2]/max_of_all_norms;

			}










			for( unsigned int j = 0 ; j < sinusoids[i].size() ; j=j+3 ){

			MassPoint *mp = new MassPoint();

			mp->position[0] = (sinusoids[i][j]/max_of_all_norms - mass_center_of_all_cells[0]/max_of_all_norms);
			mp->position[1] = (sinusoids[i][j+1]/max_of_all_norms - mass_center_of_all_cells[1]/max_of_all_norms);
			mp->position[2] = (sinusoids[i][j+2]/max_of_all_norms - mass_center_of_all_cells[2]/max_of_all_norms);

			mp->r = 0.0*255;
			mp->g = 0.0*255;
			mp->b = 0.0*255;

			mp->transparency = 0.8;

			//this->cells[0]->mass.push_back(mp);
//			this->cells[tmp]->mass.push_back(mp);
		 			}//over all voxels

			//this->watershed.push_back(mp);//  ->cells[i]->mass.push_back(mp);
			}//over all voxels


		}

		this->x_left = 0./ max_of_all_norms;
		this->x_left -= mass_center_of_all_cells[0]/max_of_all_norms;
		this->y_left = 0./ max_of_all_norms;
		this->y_left -= mass_center_of_all_cells[1]/max_of_all_norms;
		this->z_left = 0./ max_of_all_norms;
		this->z_left -= mass_center_of_all_cells[2]/max_of_all_norms;

		this->x_right = 200./ max_of_all_norms;
		this->x_right -= mass_center_of_all_cells[0]/max_of_all_norms;
		this->y_right = 200./ max_of_all_norms;
		this->y_right -= mass_center_of_all_cells[1]/max_of_all_norms;
		this->z_right = 31./ max_of_all_norms;
		this->z_right -= mass_center_of_all_cells[2]/max_of_all_norms;

			if( force_profile ){
				calcForceOfDeformation();

			}








			break;
		}
    case PushTwoCells:
		setMode(201);
        AddCell(0.6,0.0,0.0);
        AddCell(-0.6,0.0,0.0);
        setObservationPoints();
        break;
	case StretchCell:
		setMode(200);
        AddCell( 0., 0., 0. );
        setStretchPoints();
		break;
  case Wall:
    this->setMode(203);
    AddCell(0.0,0.0,0.0);
//    for( unsigned int i = 0 ; i < cells[0]->mass.size() ; i++)
//      cells[0]->mass[i]->v[0] = 10;
    break;
  default:
    setMode(102);
    AddCell( 0., 0., 0. );
    break;
  }


  enableSimulation = true;

  mpSimulationThread->setModel( this );


  mpSimulationThread->setMaxSteps( 1000000 );//TODO maschinen maximum


}
void Model3D::SimulateInThread(){

  if ( this->time <= this->time_simulation ){
    fprintf(stderr, "%2.3f \n", this->time );
    this->Simulate();
  }
  else
    enableSimulation = false;

}
double Model3D::GetSimulationProgress(){

  return (this->time / this->time_simulation);

}



/*
 * Parameter and create Model
 */
CSParameterContext * Model3D::GetParameters( std::string contextName ){

  if ( !mpParameters){
    RegisterParameters();
    return mpParameters;
  }

  if (!contextName.size()) return mpParameters;

  if ( contextName == name ) return mpParameters;

  return mpParameters->findContext(contextName);

}
void Model3D::RegisterParameters(){

  if ( mpParameters ) return;

  mpParameters = new CSParameterContext( this->name );

  // can we get the defaults from some predefined constants?
  std::vector<std::string> simualtion_mode;
   simualtion_mode.push_back("singleCell");
   simualtion_mode.push_back("stretchSingleCell");
   simualtion_mode.push_back("watershedFitCell");
   simualtion_mode.push_back("pushTwoCellsTogether");
   simualtion_mode.push_back("pressure");
   simualtion_mode.push_back("cylinder(1D)");
   simualtion_mode.push_back("twoPlates(2D)");
   simualtion_mode.push_back("against_a_wall");
   simualtion_mode.push_back("growth_in_a_cube");
  CSParameterChoice *mpSimulation_modeChoice = new CSParameterChoice(simualtion_mode, 0);
  mpParameters->addParameter("Simulation Mode", CSParameter::Choice, mpSimulation_modeChoice, "");

  mpParameters->addParameter("xml_name", CSParameter::String, new string("model3d_xml"), "");

  std::vector<std::string> dimension;
   dimension.push_back("1");
   dimension.push_back("2");
   dimension.push_back("3");
  CSParameterChoice *mpDimensionChoice = new CSParameterChoice(dimension, 2);
  mpParameters->addParameter("Dimension", CSParameter::Choice,mpDimensionChoice, "");

  // length scale
  mpParameters->addParameter("Cell Diameter", CSParameter::Double, new double(1.), "m");

  // cycle time
  mpParameters->addParameter("Cell Cycle Time", CSParameter::Double, &defaultCellCycleTime, "s");

  // std deviation of cycle times
  mpParameters->addParameter("Cycle Time SD", CSParameter::Double, &defaultCellCycleSD, "s");

  mpParameters->addParameter("Simulation Time", CSParameter::Double, &time_simulation, "" );
  mpParameters->addParameter("number_of_points", CSParameter::Int, &default_Mass_Number, "");
  mpParameters->addParameter("output_path", CSParameter::String, &outputPath, "");
  mpParameters->addParameter("input_path", CSParameter::String, new string("../../input/"), "");
  mpParameters->addParameter("watershedFitCells", CSParameter::String, new string("Cells_2.txt"), "");
  mpParameters->addParameter("k_spring_shell", CSParameter::Double, new double(80.), "");
  mpParameters->addParameter("eta_spring_shell", CSParameter::Double, new double(40.), "");
  mpParameters->addParameter("k_spring_cytoskeleton", CSParameter::Double, new double(40.), "");
  mpParameters->addParameter("eta_spring_cytoskeleton", CSParameter::Double, new double(9.), "");
  mpParameters->addParameter("k_second_damper", CSParameter::Double, new double(10.), "");
  mpParameters->addParameter("eta_second_damper", CSParameter::Double, new double(40.), "");
  mpParameters->addParameter("eta_medium", CSParameter::Double, new double(2.), "");
  mpParameters->addParameter("mass of mass point", CSParameter::Double, new double(1.), "");
  mpParameters->addParameter("max_steps", CSParameter::Int, new int(3000), "");
  mpParameters->addParameter("log_stepsize", CSParameter::Int, new int(-8), "");
  mpParameters->addParameter("initial_radius", CSParameter::Double, new double(0.5), "");
  mpParameters->addParameter("conserve_volume", CSParameter::Bool, new bool(true), "");
  mpParameters->addParameter("pressure_threshold", CSParameter::Double, new double(0.01), "");
  mpParameters->addParameter("force_push_two_cells", CSParameter::Double, new double(5.), "");


  CSParameterContext * outputContext = mpParameters->addContext( "Output" );
   outputContext->addParameter("output_prefix", CSParameter::String, new string("_s1000"), "");


  // generate the list of parameters that can be used as output file suffix:
  std::vector<CSParameter *> parms = mpParameters->getParameters();
  std::vector<CSParameter *>::iterator parmIt = parms.begin();
  std::vector<std::string> parmStrings;
  for ( parmIt = parms.begin(); parmIt != parms.end(); ++parmIt)
    parmStrings.push_back( (*parmIt)->name() );

  CSParameterChoice * suffixParameterChoice = new CSParameterChoice( parmStrings , 3);
  outputContext->addParameter("output_suffix (by parameter)", CSParameter::Choice, suffixParameterChoice, "" );




  CSParameterContext * output_povrayContext = outputContext->addContext( "Povray" );
   output_povrayContext->addParameter("output_to_pov", CSParameter::Bool, new bool(false), "");
   output_povrayContext->addParameter("time_beetwen_output", CSParameter::Double, new double(0.4), "");
   output_povrayContext->addParameter("compress_pov", CSParameter::Bool, new bool(true), "");
   output_povrayContext->addParameter("force_profile", CSParameter::Bool, new bool(false), "");
   output_povrayContext->addParameter("print_sinusoids", CSParameter::Int, new int(2), "");

  CSParameterContext * output_povray_camContext = output_povrayContext->addContext( "camera" );
   output_povray_camContext->addParameter("pov_cam_x", CSParameter::Double, new double(-25.0), "");
   output_povray_camContext->addParameter("pov_cam_y", CSParameter::Double, new double(-25.0), "");
   output_povray_camContext->addParameter("pov_cam_z", CSParameter::Double, new double(35.0), "");

  //set up in model3d.cpp -> createFromXML !!!!!!!!

  // fill the comboBoxXmlFileSuffixParameter
  // std::vector<CSParameter *> parms = mpParameterContext->getParameters();
  // for ( unsigned int i=0; i < parms.size(); ++i )
  //     comboBoxXmlFileSuffixParameter->addItem( QString(parms[i]->name().c_str()) );



}
void Model3D::DefaultParameters( CSParameterContext *emptyContext ){
  if ( !emptyContext ) return;


  // choose simulation mode
  std::vector<std::string> simualtion_mode;
   simualtion_mode.push_back("singleCell");
   simualtion_mode.push_back("stretchSingleCell");
   simualtion_mode.push_back("watershedFitCell");
   simualtion_mode.push_back("pushTwoCellsTogether");
   simualtion_mode.push_back("pressure");
   simualtion_mode.push_back("cylinder(1D)");
   simualtion_mode.push_back("twoPlates(2D)");
   simualtion_mode.push_back("against_a_wall");
   simualtion_mode.push_back("growth_in_a_cube");
  CSParameterChoice *mpSimulation_modeChoice = new CSParameterChoice(simualtion_mode, 0);
  emptyContext->setParameter("Simulation Mode", CSParameter::Choice, mpSimulation_modeChoice, "");

  //model name
  emptyContext->setParameter("xml_name", CSParameter::String, new string("model3d_xml"), "");

  std::vector<std::string> dimension;
   dimension.push_back("1");
   dimension.push_back("2");
   dimension.push_back("3");
   CSParameterChoice *mpDimensionChoice = new CSParameterChoice(dimension, 2);
  emptyContext->setParameter("Dimension", CSParameter::Choice,mpDimensionChoice, "");

  // length scale
  emptyContext->setParameter( "Cell Diameter",      CSParameter::Double, new double( 1. ), "");

  // cycle time
  emptyContext->setParameter( "Cell Cycle Time",    CSParameter::Double, new double( 1 ), "");

  // cycle time
  emptyContext->setParameter( "Cycle Time SD",    CSParameter::Double, new double( .1 ), "");

//  emptyContext->setParameter("dimension", CSParameter::Choice,mpDimensionChoice, "");
  emptyContext->setParameter("Simulation Time", CSParameter::Double, new double(10.), "" );
  emptyContext->setParameter("number_of_points", CSParameter::Int, new int(100), "");
  emptyContext->setParameter("output_path", CSParameter::String, new string("../../output/"), "");
  emptyContext->setParameter("input_path", CSParameter::String, new string("../../input/"), "");
  emptyContext->setParameter("watershedFitCells", CSParameter::String, new string("Cells_2.txt"), "");
  emptyContext->setParameter("k_spring_shell", CSParameter::Double, new double(80.), "");
  emptyContext->setParameter("eta_spring_shell", CSParameter::Double, new double(40.), "");
  emptyContext->setParameter("k_spring_cytoskeleton", CSParameter::Double, new double(40.), "");
  emptyContext->setParameter("eta_spring_cytoskeleton", CSParameter::Double, new double(9.), "");
  emptyContext->setParameter("k_second_damper", CSParameter::Double, new double(10.), "");
  emptyContext->setParameter("eta_second_damper", CSParameter::Double, new double(40.), "");
  emptyContext->setParameter("eta_medium", CSParameter::Double, new double(2.), "");
  emptyContext->setParameter("mass of mass point", CSParameter::Double, new double(1.), "");
  emptyContext->setParameter("log_stepsize", CSParameter::Int, new int(-8), "");
  emptyContext->setParameter("conserve_volume", CSParameter::Bool, new bool(true), "");
  emptyContext->setParameter("pressure_threshold", CSParameter::Double, new double(0.01), "");
  emptyContext->setParameter("force_push_two_cells", CSParameter::Double, new double(5.), "");

  CSParameterContext * outputContext = emptyContext->addContext( "Output" );
  outputContext->setParameter("output_prefix", CSParameter::String, new string("_s1000"), "");
  outputContext->addParameter("output_suffix", CSParameter::String, new string("byTime"), "");


  CSParameterContext * output_povrayContext = outputContext->addContext( "Povray" );
  output_povrayContext->setParameter("output_to_pov", CSParameter::Bool, new bool(false), "");
  output_povrayContext->setParameter("time_beetwen_output", CSParameter::Double, new double(0.4), "");
  output_povrayContext->setParameter("compress_pov", CSParameter::Bool, new bool(true), "");
  output_povrayContext->setParameter("force_profile", CSParameter::Bool, new bool(false), "");
  output_povrayContext->setParameter("print_sinusoids", CSParameter::Int, new int(2), "");

  CSParameterContext * output_povray_camContext = output_povrayContext->addContext( "camera" );
 // output_povray_camContext->addParameter("pov_cam_x", CSParameter::Double, new double(-25.0), "");
  output_povray_camContext->setParameter("pov_cam_x", CSParameter::Double, new double(-25.0), "");
  output_povray_camContext->setParameter("pov_cam_y", CSParameter::Double, new double(-25.0), "");
  output_povray_camContext->setParameter("pov_cam_z", CSParameter::Double, new double(35.0), "");

}

void Model3D::InitParameters( CSParameterContext *parms ){
  if (!parms){
    if ( !mpParameters )
      RegisterParameters();
    DefaultParameters( mpParameters );
    return;
  }

  std::vector<CSParameter *> parameters = parms->getParameters();

  std::vector<CSParameter *>::const_iterator parmsIt;

  for ( parmsIt=parameters.begin(); parmsIt!=parameters.end(); ++parmsIt ){
    if ( (*parmsIt)->name() == "Cell Diameter" ){
      this->defaultInitialCellRadius = (*parmsIt)->value()*0.5;
      double a = 0;
    }
    else if ( (*parmsIt)->name() == "Cell Cycle Time" )
      this->defaultCellCycleTime = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "Cycle Time SD" )
      this->defaultCellCycleSD = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "Simulation Mode" ){
      CSParameter * foundParm = mpParameters->findParameter( (*parmsIt)->name() );
      if ( foundParm )
        ((CSParameterChoice *) foundParm->dataPointer()) ->setCurrentIndex(((CSParameterChoice *)(*parmsIt)->dataPointer()) ->currentIndex());
      else
        std::cout << "Error in finding parameter \""
                  << (*parmsIt)->name()
                  << "\"\n";

      int mode = (*parmsIt)->value();
      switch( mode ){
        case 0:
          this->mode = 102;//SingleCell
          break;
        case 1:
          this->mode = 200;//StretchCell
          break;
        case 2:
          this->mode = 2001;//WatershedFit
          break;
        case 3:
          this->mode = 201;//PushTwoCells
          break;
        case 4:
          this->mode = 202;//Pressure
          break;
        case 5:
          this->mode = 301;//Cylinder
          break;
        case 6:
          this->mode = 302;//Plane
          break;
        case 7:
          this->mode = 203;//Wall
          break;
        case 8:
          this->mode = 204;//Cube
          break;
        default:
          this->mode = 102;//SingleCell
          break;
      }
    }
    else if ( (*parmsIt)->name() == "xml_name" )
      this->xmlName = (*parmsIt)->dataString();
    else if ( (*parmsIt)->name() == "Dimension" ){

      CSParameter * foundParm = mpParameters->findParameter( (*parmsIt)->name() );

      if ( foundParm )
        ((CSParameterChoice *) foundParm->dataPointer()) ->setCurrentIndex(((CSParameterChoice *)(*parmsIt)->dataPointer()) ->currentIndex());
      else
        std::cout << "Error in finding parameter \""
                  << (*parmsIt)->name()
                  << "\"\n";

      this->dimension = (*parmsIt)->value()+1;
    }
    else if ( (*parmsIt)->name() == "Simulation Time" )
      this->time_simulation = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "number_of_points" )
      this->default_Mass_Number = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "output_path" )
      this->outputPath = (*parmsIt)->dataString();
    else if ( (*parmsIt)->name() == "input_path" )
      this->inputPath = (*parmsIt)->dataString();
    else if ( (*parmsIt)->name() == "watershedFitCells" )
      this->watershedFitCells = (*parmsIt)->dataString();
    else if ( (*parmsIt)->name() == "k_spring_shell" )
      this->default_k_huell = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "eta_spring_shell" )
      this->default_nu_huell = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "k_spring_cytoskeleton" )
      this->default_k_cytoskeleton = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "eta_spring_cytoskeleton" )
      this->default_nu_cytoskeleton = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "k_second_damper" )
      this->default_k2 = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "eta_second_damper" )
      this->default_damper = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "eta_medium" )
      this->nu_medium = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "mass of mass point" )
      this->default_PointMass = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "log_stepsize" )
      this->mLog2timeStep = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "conserve_volume" )
      this->conserve_volume = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "pressure_threshold" )
      this->pressure_threshold = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "force_push_two_cells" )
      this->force_push_two_cells = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "output_prefix" )
      this->output_prefix = (*parmsIt)->dataString();
    else if ( (*parmsIt)->name() == "output_suffix" )
      this->output_suffix = (*parmsIt)->dataString();
    else if ( (*parmsIt)->name() == "output_to_pov" )
      this->enablePovrayOutput = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "time_beetwen_output" )
      this->timeBetweenOutputs = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "compress_pov" )
      this->compress_pov = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "force_profile" )
      this->force_profile = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "print_sinusoids" )
      this->print_sinusoids = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "pov_cam_x" )
      this->pov_cam_x = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "pov_cam_y" )
      this->pov_cam_y = (*parmsIt)->value();
    else if ( (*parmsIt)->name() == "pov_cam_z" )
      this->pov_cam_z = (*parmsIt)->value();
     else
       std::cerr << "ModelCellsSpherical::InitParameters:  Unknown parameter given\t"  << (*parmsIt)->name() << std::endl;
  }
}

Model3D * Model3D::createFromXML( QXmlStreamReader * xmlStream, std::stringstream & errors, std::stringstream & warnings ){
  Q_ASSERT( xmlStream->name() == "Model" && xmlStream->isStartElement() );
//TODO überarbeiten
  Model3D * model_3d = new Model3D();

  CSParameterContextTemporary * parms = NULL;

  while ( xmlStream->readNextStartElement() ){
    QStringRef elementName = xmlStream->name();

    if ( elementName == "parameters" )
    {
        parms = CSParameterContextTemporary::createFromXML( xmlStream );
    }
    else if ( elementName == "cells" )
    {
        // read in cell definitions
    }
    else
    {
        warnings << "\nIn Model3D::createFromXML():  Warning:  "
                 << "Unknown element type:\""
                 << elementName.toString().toStdString() << "\"\n";
        xmlStream->skipCurrentElement();
    }
  }

/*
    // t1m-debug
    std::stringstream parmDump;
    ((CSParameterContext *)parms)->dump( parmDump );

    std::cerr << parmDump.str();
    // ! t1m-debug
*/
  CSParameter * dummyParm = NULL;

  if ( parms ){
   // set the parameters from the context
   dummyParm = parms->findParameter("Simulation mode");
   if ( dummyParm ){
     CSParameterChoice *choice = (CSParameterChoice *)dummyParm->dataPointer();
     switch ( choice->currentIndex() ){
       case 0:
         model_3d->mode = SingleCell; break;
       case 1:
         model_3d->mode = StretchCell; break;
       case 2:
         model_3d->mode = WatershedFit; break;
       case 3:
         model_3d->mode = PushTwoCells; break;
       case 4:
         model_3d->mode = Pressure; break;
       case 5:
         model_3d->mode = Cylinder; break;
       case 6:
         model_3d->mode = Plane; break;
       case 7:
         model_3d->mode = Wall; break;
       case 8:
         model_3d->mode = Cube; break;
     }
   }

        dummyParm = parms->findParameter("dimension");
        if ( dummyParm )
        {
            CSParameterChoice *choice = (CSParameterChoice *)dummyParm->dataPointer();
            model_3d->dimension = choice->currentIndex() +1;
        }

        dummyParm = parms->findParameter("xml_name");
        if ( dummyParm )
        {
            model_3d->SetName( *(std::string *)dummyParm->dataPointer() );
            model_3d->xmlName = model_3d->name;
        }

        dummyParm = parms->findParameter("simulation time");
        if ( dummyParm )
            model_3d->time_simulation = *(double *)dummyParm->dataPointer();


        dummyParm = parms->findParameter("number_of_points");
        if ( dummyParm )
            model_3d->default_Mass_Number = *(int *)dummyParm->dataPointer();

        dummyParm = parms->findParameter("mass of mass point");
        if ( dummyParm )
            model_3d->default_PointMass = *(double *)dummyParm->dataPointer();

        dummyParm = parms->findParameter("initial_radius");
        if ( dummyParm )
            model_3d->defaultInitialCellRadius = *(double *)dummyParm->dataPointer();


        dummyParm = parms->findParameter("k_spring_shell");
        if ( dummyParm )
            model_3d->default_k_huell = *(double *)dummyParm->dataPointer();

        dummyParm = parms->findParameter("eta_spring_shell");
        if ( dummyParm )
            model_3d->default_nu_huell = *(double *)dummyParm->dataPointer();


        dummyParm = parms->findParameter("k_spring_cytoskeleton");
        if ( dummyParm )
            model_3d->default_k_cytoskeleton = *(double *)dummyParm->dataPointer();

        dummyParm = parms->findParameter("eta_spring_cytoskeleton");
        if ( dummyParm )
            model_3d->default_nu_cytoskeleton = *(double *)dummyParm->dataPointer();


        dummyParm = parms->findParameter("k_second_damper");
        if ( dummyParm )
            model_3d->default_k2 = *(double *)dummyParm->dataPointer();

        dummyParm = parms->findParameter("eta_second_damper");
        if ( dummyParm )
            model_3d->default_damper = *(double *)dummyParm->dataPointer();


        dummyParm = parms->findParameter("eta_medium");
        if ( dummyParm )
            model_3d->nu_medium = *(double *)dummyParm->dataPointer();


        // todo:  use as exponent (2^Dt)?  Or the actual step size?
        dummyParm = parms->findParameter("log_stepsize");
        if ( dummyParm )
            model_3d->mLog2timeStep = *(int *) dummyParm->dataPointer();


        dummyParm = parms->findParameter("conserve_volume");
        if ( dummyParm )
            model_3d->conserve_volume = *(bool *)dummyParm->dataPointer();

        std::string outputTrunk = "";
        dummyParm = parms->findParameter("output_prefix");
        if ( dummyParm )
			outputTrunk = *(std::string *)dummyParm->dataPointer();

		dummyParm = parms->findParameter("output_path");
        if ( dummyParm )
			model_3d->outputPath = *(std::string *)dummyParm->dataPointer();

		dummyParm = parms->findParameter("input_path");
        if ( dummyParm )
		model_3d->inputPath = *(std::string *)dummyParm->dataPointer();

		dummyParm = parms->findParameter("watershedFitCells");
        if ( dummyParm )
		model_3d->watershedFitCells = *(std::string *)dummyParm->dataPointer();

        dummyParm = parms->findParameter("output_suffix (by parameter)");
        if ( dummyParm )
        {
            CSParameterChoice * suffixChoice =
                (CSParameterChoice *) dummyParm->dataPointer();
            std::string suffixParmString = suffixChoice->currentString();
            CSParameter *suffixParm = parms->findParameter(suffixParmString);

            QString suffixParmQString = QString( suffixParmString.c_str() );
            suffixParmQString.replace( " ", "_" );

            outputTrunk += "-" + suffixParmQString.toStdString();

            if ( suffixParm )
                outputTrunk += "-" + suffixParm->dataString();
        }

        model_3d->output_prefix = outputTrunk;

        dummyParm = parms->findParameter("output_to_pov");
        if ( dummyParm )
            model_3d->enablePovrayOutput = *(bool *)dummyParm->dataPointer();

		dummyParm = parms->findParameter("time_beetwen_output");
        if ( dummyParm )
            model_3d->timeBetweenOutputs = *(double *)dummyParm->dataPointer();

		dummyParm = parms->findParameter("compress_pov");
        if ( dummyParm )
            model_3d->compress_pov = *(bool *)dummyParm->dataPointer();

		dummyParm = parms->findParameter("pressure_threshold");
        if ( dummyParm )
            model_3d->pressure_threshold = *(double *)dummyParm->dataPointer();

		dummyParm = parms->findParameter("force_profile");
        if ( dummyParm )
            model_3d->force_profile = *(bool *)dummyParm->dataPointer();


		dummyParm = parms->findParameter("print_sinusoids");
        if ( dummyParm )
            model_3d->print_sinusoids = *(int *)dummyParm->dataPointer();

		dummyParm = parms->findParameter("pov_cam_x");
        if ( dummyParm )
            model_3d->pov_cam_x = *(double *)dummyParm->dataPointer();

		dummyParm = parms->findParameter("pov_cam_y");
        if ( dummyParm )
            model_3d->pov_cam_y = *(double *)dummyParm->dataPointer();
		dummyParm = parms->findParameter("pov_cam_z");
        if ( dummyParm )
            model_3d->pov_cam_z = *(double *)dummyParm->dataPointer();

		dummyParm = parms->findParameter("force_push_two_cells");
		if ( dummyParm )
            model_3d->force_push_two_cells = *(double *)dummyParm->dataPointer();

        std::cerr << "Model3D::createFromXML():\n"
                  << "Simulation mode:\t" << model_3d->mode << std::endl
                  <<  "k_spring_shell\t" << model_3d->default_k_huell << std::endl
                  <<  "eta_spring_shell\t" << model_3d->default_nu_huell << std::endl
                  <<  "k_spring_cytoskeleton\t" << model_3d->default_k_cytoskeleton << std::endl
                  <<  "eta_spring_cytoskeleton\t" << model_3d->default_nu_cytoskeleton << std::endl
                  <<  "k_second_damper\t" << model_3d->default_k2 << std::endl
                  <<  "eta_second_damper\t" << model_3d->default_damper << std::endl
                  <<  "eta_medium\t" << model_3d->nu_medium << std::endl
                  <<  "mass of mass point\t" << model_3d->default_PointMass << std::endl
                  <<  "stepsize\t" << model_3d->mLog2timeStep << std::endl
                  <<  "initial_radius\t" << model_3d->defaultInitialCellRadius << std::endl
                  <<  "conserve_volume\t" << (model_3d->conserve_volume?"yes":"no") << std::endl
                  <<  "number_of_points\t" << model_3d->default_Mass_Number << std::endl
                  <<  "output_prefix\t" << model_3d->output_prefix << std::endl
                  <<  "xml_name\t" << model_3d->xmlName << std::endl
                  <<  "output_to_pov\t" << (model_3d->enablePovrayOutput?"yes":"no") << std::endl;

//        delete parms;
    }

    return model_3d;
}
