///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphSphere.cpp                                                      //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-09-06 13:16:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include <GraphSphere.h>
#include "../../tools/model/CSModelTools.h"
#include "../Elements/ModelElementSphere.h"
#include <cmath>
#include <cstdio>
#include <fstream>

void GraphSphere::SampleGraph(double distance)
{

  //list for new elements
  std::vector<CSGraphEdge*> new_edge;

  //resample all edges
  for( unsigned int i = 0 ; i < this->mvEdge.size() ; i++){

    double edgeLength = this->mvEdge[i]->length();

    //calc how many subdivisions
    int amount = std::ceil(edgeLength / distance);
    double sub_length = edgeLength / amount;
    double sub_radius = (this->mvEdge[i]->mpEnd->mRadius - this->mvEdge[i]->mpStart->mRadius) /amount;

    //choose type
    unsigned int type;
    if( this->mvEdge[i]->mpEnd->mVesselType == this->mvEdge[i]->mpStart->mVesselType )
      type = this->mvEdge[i]->mpEnd->mVesselType;
    else
      type = ModelElementVesselSphere::Sinusoid;

    ModelElementVesselSphere *nodeStartTmp = this->mvEdge[i]->mpStart;
    ModelElementVesselSphere *nodeEndTmp;

    //subdivide
    for( int j = 1 ; j < amount ; j++ ){

      CSGraphEdge *edge = new CSGraphEdge();
      ModelElementVesselSphere *node = new ModelElementVesselSphere(0,0,0);

      //set index
      node->mIndex = this->mvNode.size();

      node->setQuality(5,5);
//      this->mpBoundingBoxList->add(node);
      new_edge.push_back(edge);
      this->mvNode.push_back(node);
      node->mVesselType = type;

//      node->youngModulus = this->defaultYoungModulus;
//      node->poissonRatio = this->defaultPoissonRatio;

      //connect edges and nodes
      if( j == 1){//( amount - 1 ) ){
        edge->mpStart = nodeStartTmp;
        edge->mpEnd = node;
        this->mvEdge[i]->mpStart = node;
        nodeEndTmp = node;
      }else{
        edge->mpStart= nodeEndTmp;
        edge->mpEnd = node;
        this->mvEdge[i]->mpStart = node;
        nodeEndTmp = node;
      }

      //set radius
      node->mRadius = nodeStartTmp->mRadius + j * sub_radius;

      //set postion of additional node
      node->position.x = nodeStartTmp->position.x + this->mvEdge[i]->mDirection[0]*j*sub_length;
      node->position.y = nodeStartTmp->position.y + this->mvEdge[i]->mDirection[1]*j*sub_length;
      node->position.z = nodeStartTmp->position.z + this->mvEdge[i]->mDirection[2]*j*sub_length;
      
      node->setBoundingBox();
    }
  }

  for( unsigned int i = 0 ; i < new_edge.size() ; i++ ){
    //set index
    new_edge[i]->mIndex = this->mvEdge.size();
    this->mvEdge.push_back(new_edge[i]);

  }




}

void GraphSphere::ConnectNodes()
{
  for( unsigned int i = 0 ; i < this->mvEdge.size() ; i++ ){
    
    //add neighbours
    this->mvEdge[i]->mpStart->mvpNeighbor.push_back( this->mvEdge[i]->mpEnd );
    this->mvEdge[i]->mpEnd->mvpNeighbor.push_back( this->mvEdge[i]->mpStart );

    //add current edge
    this->mvEdge[i]->mpStart->mvpEdges.push_back(this->mvEdge[i]);
    this->mvEdge[i]->mpEnd->mvpEdges.push_back(this->mvEdge[i]);
  

  }
}

void GraphSphere::setStatic( unsigned int type ){

  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++ )
    if( this->mvNode[i]->mVesselType == type )
      this->mvNode[i]->mStatic = 1;

}

void GraphSphere::setInitialOverlaps(){

  fprintf(stderr, "set overlaps in vessel graph\n");

  for( unsigned int i = 0 ; i < (this->mvNode.size()-1) ; i++ ){
    for( unsigned int j = i+1 ; j < this->mvNode.size() ; j++){

      if( CSModelTools::OverlapSphericalElements( static_cast<ModelElementSphere *>(this->mvNode[i]) , static_cast<ModelElementSphere *>(this->mvNode[j]) )){
        this->mvNode[i]->mvpBlackList.push_back( this->mvNode[j] );
        this->mvNode[j]->mvpBlackList.push_back( this->mvNode[i] );
      }

    }
  }

  fprintf(stderr, "ready set overlaps in vessel graph\n");

}

void GraphSphere::calcSpringForce(double k){

  for( unsigned int i = 0 ; i < this->mvEdge.size() ; i++ ){

    double diffLengthDotK = (this->mvEdge[i]->length() - this->mvEdge[i]->mLength0)*k;

    double force[3];
    force[0] = this->mvEdge[i]->mDirection[0] * diffLengthDotK;
    force[1] = this->mvEdge[i]->mDirection[1] * diffLengthDotK;
    force[2] = this->mvEdge[i]->mDirection[2] * diffLengthDotK;

    this->mvEdge[i]->mpStart->directedForce.Add(force[0],force[1],force[2]);
    this->mvEdge[i]->mpEnd->directedForce.Add(-force[0],-force[1],-force[2]);

    this->mvEdge[i]->mpStart->accumulatedForceAbsolute += std::abs(diffLengthDotK);
    this->mvEdge[i]->mpEnd->accumulatedForceAbsolute += std::abs(diffLengthDotK);
  }

}

void GraphSphere::updateParametersForAllSpheres()
{
    for ( auto sphere: mvNode )
    {
        sphere->poissonRatio = defaultPoissonRatio;
        sphere->youngModulus = defaultYoungModulus;
    }
}
