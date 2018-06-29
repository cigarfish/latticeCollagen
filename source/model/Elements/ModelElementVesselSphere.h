///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ModelElementVesselGraph.h                                            //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-09-07 20:12:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef MODEL_ELEMENT_VESSEL_SPHERE_H
#define MODEL_ELEMENT_VESSEL_SPHERE_H

#include "ModelElementSphere.h"
//#include "../model/BasicDatatypes/Graph.h"


#include <vector>
#include <sstream>

#include <Vector.h>

namespace H5 { class CompType; };

struct BoundingBox;

class CSGraphEdge;

class ModelElementVesselSphere : public ModelElementSphere
{
public:
  ModelElementVesselSphere(double x=0, double y=0, double z=0);

  static void HDF5DataFormat( H5::CompType & );
  static H5::CompType ParseHDF5DataFormat( H5::CompType & inputTypeDefinition,
                                           std::stringstream &,
                                           std::stringstream & warnings );

  unsigned int mIndex;

  double mYoungModulus;
  double mPoisonRatio;
  double mStiffness;
  double mPressure;
  double mConcentration;
  double mdConcentration;//concentration difference
  double mVolume;

  // added by Jieling
  int highlight; // 0: none; 1: ECM; 2: HSC

  std::vector<ModelElementVesselSphere*> mvpNeighbor;
  std::vector<CSGraphEdge*> mvpEdges;
  std::vector<CSGraphEdge*> mvpSegments;

  std::vector<ModelElementVesselSphere*> mvpBlackList;


  BoundingBox * boundingBox();
  BoundingBox * getBoundingBox();
  void setBoundingBox();

  enum VesselType {
      CentralVein = 1,
      PortalVein = 2,
      Sinusoid = 3,
      DeadEnd = 4
  };

  unsigned int mVesselType;


};

#endif //MODEL_ELEMENT_VESSEL_SPHERE_H
