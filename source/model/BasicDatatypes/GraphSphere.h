///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphSphere.h                                                        //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-09-06 13:14:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef GRAPHSPHERE_H
#define GRAPHSPHERE_H

#include <Graph.h>

class GraphSphere: public Graph
{

  public:
    void SampleGraph(double distance);
    void ConnectNodes();
    void setStatic( unsigned int type );
    void setInitialOverlaps();
    void calcSpringForce(double k);
    void updateParametersForAllSpheres();
};


#endif //GRAPHSPHERE_H
