///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Graph.h                                                              //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-09-05 11:21:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef GRAPH_H
#define GRAPH_H


#include <vector>
#include <../Elements/ModelElementVesselSphere.h>
#include "../../tools/model/BoundingBoxList.h"

class CSGraphEdge;
class CSGraphNode;

class CSGraphEdge
{
  public:
    ModelElementVesselSphere* mpStart;
    ModelElementVesselSphere* mpEnd;

    short int mIndex;

    CSGraphEdge();
    double length();
    double mLength0;
    double mDirection[3];

    double mRadius_start;
    double mRadius_end;

    double mVolume;

    double mViscosity;
    double mFlux;
    double mPermeability;
    double mInverseResistenceToFlow;

    void setViscosity();
    void setInverseResistanceToFlow();
    void calcFlux();
    double calcVolume();

};



class Graph
{
  public:
    Graph();
    void resize(double scale);
    void read( const char * filename, int i=1 );
    void readVTK(const char * filename);
    void readMXF(const char * filename);
    void setLength0();
    void setBoundingBoxList(BoundingBoxList *BoundingBoxList);
    void setBoundingBox();
    void shift(double x, int dim);
    double calcDimensions(double &lobule_radius, double &lobule_height);

    void initPressure();
    void calcPressure();
    void calcFlux();
    void calcVolume();
    void updateConcentration(double dt);

    // list of graph segments
    std::vector<CSGraphEdge *> mvSegments;
    // List of edges.
    // Differs from mvSegments as the edges between sampled are collected here.
    std::vector<CSGraphEdge *> mvEdge;
    std::vector<ModelElementVesselSphere*> mvNode;

    BoundingBoxList *mpBoundingBoxList;

    //bounding Box
    double mXmin;
    double mXmax;
    double mYmin;
    double mYmax;
    double mZmin;
    double mZmax;

    //physical properties
    double defaultYoungModulus;
    double defaultPoissonRatio;


};


#endif //GRAPH_H
