///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSGLVascularization.h                                                //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-06-05 19:10:37                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef CSGL_VASCULARIZATION_H
#define CSGL_VASCULARIZATION_H

#include "../CSGLObject.h"
#include "../../model/Model/ModelVascularization/VoronoiDiagram.hpp"

class CSVesselGraph;

class CSGLVascularization : public CSGLObject
{
public:
  CSGLVascularization( CSVesselGraph *CSvg );
  void draw();
  void setVoronoi( VoronoiDiagram *voronoi );


protected:
  double mRadius;
  CSVesselGraph *vg;
  VoronoiDiagram *voronoi;

public:

};


#endif // TRIANGULATED_CELL_H
