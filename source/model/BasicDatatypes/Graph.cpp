///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Graph.cpp                                                            //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-09-05 11:21:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#include "Graph.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <string.h>
#include <stdio.h>

#ifdef _MSC_VER
#  define _USE_MATH_DEFINES
#  include <math.h>
#endif
#include <cmath>

#include <stdlib.h>

#include <locale.h>

#include <algorithm>
#include <iterator>
#include <set>

#include "../../tools/new_triangulation/LinearAlgebra.hpp"

Graph::Graph(){

 this->mpBoundingBoxList = NULL;



}

void Graph::setLength0(){
    for( unsigned int i = 0 ; i < this->mvEdge.size() ; i++ )
      this->mvEdge[i]->mLength0 = this->mvEdge[i]->length();
}

void Graph::readVTK(const char * filename){

    //read vtk graph
  std::fstream f;
  char cstring[256];
  cstring[0] = '0';
  f.open(filename, std::ios::in);

  //search postion
  bool search = 1;
  while( !f.eof() && search )
  {
    if (cstring[0] == 'P' && cstring[1] == 'O' && cstring[2] == 'I' && cstring[3] == 'N' && cstring[4] == 'T' && cstring[5] == 'S')
      search = 0;
    f.getline(cstring, sizeof(cstring));
  }

  //read position
  search = 1;
  while ( !f.eof() && search )
  {
    int index = 0;
    if ( cstring[0] == '\0' )
      search = 0;
    else{
      char *pch;
      pch = strtok(cstring," \t");

      ModelElementVesselSphere * ModelElementVesselSphere3;
      ModelElementVesselSphere * node = new ModelElementVesselSphere(0,0,0);
//      this->mpBoundingBoxList->add(node);
      ModelElementVesselSphere3 = node;
      ModelElementVesselSphere3->mIndex = mvNode.size()-1;
      mvNode.push_back(node);

      while( pch != NULL)
      {
        if( index == 3 )
        {
          index = 0;
          ModelElementVesselSphere * ModelElementVesselSphere2 = new ModelElementVesselSphere(0,0,0);
//          this->mpBoundingBoxList->add(node);
          mvNode.push_back(ModelElementVesselSphere2);
          ModelElementVesselSphere3 = ModelElementVesselSphere2;
          ModelElementVesselSphere3->mIndex = mvNode.size()-1;
        }
        int string_l = strlen(pch);
      char *tmp = new char[string_l];
      for( int i = 0 ; i < string_l ; i++){
    	  // if ( pch[i] == '.' )
    	  //     tmp[i] = ',';
    	  // else
    		  tmp[i] = pch[i];
      }
        ModelElementVesselSphere3->position.Set(atof(tmp),index);



        pch = strtok(NULL," \t");
        index++;
      }

      f.getline(cstring, sizeof(cstring));
    }
  }

  //index nodes
  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++){
    this->mvNode[i]->mIndex = i;
    this->mvNode[i]->setBoundingBox();
  }

  //search for connections
  search = 1;
  while( !f.eof() && search )
  {
    f.getline(cstring, sizeof(cstring));
    if (cstring[0] == 'E' && cstring[1] == 'D' && cstring[2] == 'G' && cstring[3] == 'E' && cstring[4] == 'S')
      search = 0;
  }

  //read connections
  search = 1;
  while( !f.eof() && search )
  {
    f.getline(cstring, sizeof(cstring));
    if (cstring[0] == 'E' && cstring[1] == 'D' && cstring[2] == 'G' && cstring[3] == 'E' && cstring[4] == '_')
      search = 0;
    else{
      char *pch;
      pch = strtok(cstring," \t");
      int start = atoi(pch);
      pch = strtok(NULL," \t");
      int end = atoi(pch);

      CSGraphEdge *edge = new CSGraphEdge();
       edge->mpStart = mvNode[start];
       edge->mpEnd = mvNode[end];
       mvSegments.push_back(edge);
       mvNode[start]->mvpSegments.push_back(edge);
       mvNode[end]->mvpSegments.push_back(edge);
    }
  }


  //index connections
  for( unsigned int i = 0 ; i < this->mvSegments.size() ; i++)
    this->mvSegments[i]->mIndex = i;


  f.getline(cstring, sizeof(cstring));
  f.getline(cstring, sizeof(cstring));


  //read radius of tubes
   search = 1;
   int index = 0;
   while ( !f.eof() && search )
   {
     f.getline(cstring, sizeof(cstring));
       if (cstring[0] == 'R' && cstring[1] == 'a' && cstring[2] == 'd' && cstring[3] == 'i' && cstring[4] == 'u'  && cstring[5] == 's')
         search = 0;
       else{
    	   char *pch;
    	   pch = strtok(cstring," \t");

    	   while( pch != NULL)
    	        {
    		     int string_l = strlen(pch);
    		      char *tmp = new char[string_l];
    		      for( int i = 0 ; i < string_l ; i++){
    		    	  // if ( pch[i] == '.' )
    		    	  //     tmp[i] = ',';
    		    	  // else
    		    		  tmp[i] = pch[i];
    		      }


    		        double t = atof(tmp);
                    this->mvSegments[index]->mRadius_start = t;

    		        index++;

    	        pch = strtok(NULL," \t");

    	        }
	     }

   }
   search = 1;
   while ( !f.eof() && search )
     {
       f.getline(cstring, sizeof(cstring));
         if (cstring[0] == 'V' && cstring[1] == 'e' && cstring[2] == 'r' && cstring[3] == 't' && cstring[4] == 'e'  && cstring[5] == 'x')
           search = 0;
     }
   //read radius of tubes
    search = 1;
    index = 0;
    while ( !f.eof() && search )
    {
      f.getline(cstring, sizeof(cstring));
        if (cstring[0] == 'V' && cstring[1] == 'e' && cstring[2] == 'r' && cstring[3] == 't' && cstring[4] == 'e'  && cstring[5] == 'x')
          search = 0;
        else{
     	   char *pch;
     	   pch = strtok(cstring," \t");

     	   while( pch != NULL)
     	        {
     		     int string_l = strlen(pch);
     		      char *tmp = new char[string_l];
     		      for( int i = 0 ; i < string_l ; i++){
     		    	  // if ( pch[i] == '.' )
     		    	  //     tmp[i] = ',';
     		    	  // else
     		    		  tmp[i] = pch[i];
     		      }


     		        double t = atof(tmp);

                    this->mvSegments[index]->mRadius_end = t;
     		       index++;

     	        pch = strtok(NULL," \t");

     	        }
 	     }

    }


    //read radius of tubes
     index = 0;
     search = 1;
     while ( !f.eof() && search )
     {
       f.getline(cstring, sizeof(cstring));

      	   char *pch;
      	   pch = strtok(cstring," \t");

      	   while( pch != NULL)
      	        {
      		     int string_l = strlen(pch);
      		      char *tmp = new char[string_l];
      		      for( int i = 0 ; i < string_l ; i++){
      		    	  // if ( pch[i] == '.' )
      		    	  //     tmp[i] = ',';
      		    	  // else
      		    		  tmp[i] = pch[i];
      		      }


      		        double t = atof(tmp);

      		       this->mvNode[index]->mRadius = t;
      		       index++;

      	        pch = strtok(NULL," \t");

      	        }


     }


  /*
  //read rest of file
  search = 1;
  while ( !f.eof() && search )
  {
    F << cstring << std::endl;
    f.getline(cstring, sizeof(cstring));
  }
  */

  f.close();

  std::copy( mvSegments.begin(), mvSegments.end(), std::back_inserter(mvEdge) );

  /*
  //debug output
  std::ofstream F;
  F.open("../../input/sin_graph1.txt_neu", std::ios::out|std::ios::app );

  F << "-------------" << std::endl;

  for( unsigned int i = 0 ; i < mvNode.size() ; i++){
    F << mvNode[i]->position.x << "\t" ;
    F << mvNode[i]->position.y << "\t" ;
    F << mvNode[i]->position.z << std::endl ;
  }

  F << "-------------" << std::endl;

  for( unsigned int i = 0 ; i < mvSegments.size() ; i++){
    F << mvSegments[i]->mpStart->mIndex << "\t" ;
    F << mvSegments[i]->mpEnd->mIndex << std::endl ;
  }
  */


}

void Graph::readMXF(const char * filename){

  setlocale (LC_ALL,"C");

  //read mxf graph
  std::fstream f;

  const int l = 512;
  char cstring[l];

  memset( cstring, 0, l );

  f.open(filename, std::ios::in);

  if ( !f.is_open() )
  {
    std::cerr << "Warning:  Couldn't open file \""
              << filename
              << "\" for reading." << std::endl;
    return;
  }

  if (!f.is_open())
    return;

  //search postion
  bool search = 1;
  while( !f.eof() && search )
  {
    if (cstring[0] == '<' &&
        cstring[1] == 'N' &&
        cstring[2] == 'o' &&
        cstring[3] == 'd' &&
        cstring[4] == 'e' &&
        cstring[5] == 'L' &&
        cstring[6] == 'i' &&
        cstring[7] == 's' &&
        cstring[8] == 't' &&
        cstring[9] == '>')
      search = 0;
    f.getline(cstring, sizeof(cstring));
  }

  //read position
  search = 1;
  while ( !f.eof() && search ){

    if (cstring[0] == '<' &&
        cstring[1] == '/' &&
        cstring[2] == 'N' &&
        cstring[3] == 'o' &&
        cstring[4] == 'd' &&
        cstring[5] == 'e' &&
        cstring[6] == 'L' &&
        cstring[7] == 'i' &&
        cstring[8] == 's' &&
        cstring[9] == 't' &&
        cstring[10] == '>')
          search = 0;
    else{
      char *pch;
      pch = strtok(cstring," \t");
      pch = strtok(NULL," ");
      pch = strtok(NULL,"\"");
      pch = strtok(NULL,"\"");


       ModelElementVesselSphere * node = new ModelElementVesselSphere(0,0,0);

       // // The following is handled in GraphSphere::updateParametersForAllSpheres()
       // node->youngModulus = this->defaultYoungModulus;
       // node->poissonRatio = this->defaultPoissonRatio;

       node->mIndex = mvNode.size();
       node->setQuality(5,5);
       mvNode.push_back(node);

//       int index = 0;

       //while( pch != NULL){
       for( int index = 0 ; index < 5 ; index++){
          if ( index >= 3)
          {
              if( index == 3)
                node->mRadius = atof(pch);
              else
                 node->mVesselType = (ModelElementVesselSphere::VesselType)atoi( pch );
          }
          else
            node->position.Set(atof(pch),index);
          pch = strtok(NULL,"\"");
          pch = strtok(NULL,"\"");
       //   index++;
       }
       node->setBoundingBox();

       // this->mpBoundingBoxList->add(node);

       f.getline(cstring, sizeof(cstring));

   // this->mpBoundingBoxList->update();

    }

  }


  //search connection
   search = 1;
   while( !f.eof() && search )
   {
     if (cstring[0] == '<' &&
         cstring[1] == 'S' &&
         cstring[2] == 'e' &&
         cstring[3] == 'g' &&
         cstring[4] == 'm' &&
         cstring[5] == 'e' &&
         cstring[6] == 'n' &&
         cstring[7] == 't' &&
         cstring[8] == 'L' &&
         cstring[9] == 'i' &&
         cstring[10] == 's' &&
         cstring[11] == 't' &&
         cstring[12] == '>')
       search = 0;
     f.getline(cstring, sizeof(cstring));
   }


  //read connection
  search = 1;
  while ( !f.eof() && search ){
    if (cstring[0] == '<' &&
        cstring[1] == '/' &&
        cstring[2] == 'S' &&
        cstring[3] == 'e' &&
        cstring[4] == 'g' &&
        cstring[5] == 'm' &&
        cstring[6] == 'e' &&
        cstring[7] == 'n' &&
        cstring[8] == 't' &&
        cstring[9] == 'L' &&
        cstring[10] == 'i' &&
        cstring[11] == 's' &&
        cstring[12] == 't' &&
        cstring[13] == '>')
      search = 0;
    else{
      char *pch;
      pch = strtok(cstring," \t");
      pch = strtok(NULL,"\"");
      pch = strtok(NULL,"\"");

      CSGraphEdge * edge = new CSGraphEdge();
       edge->mIndex = mvSegments.size();
       mvSegments.push_back(edge);

   //    while( pch != NULL){


          edge->mpStart = this->mvNode[atoi(pch)-1];
          	  int a = atoi(pch)-1;
          edge->mRadius_start = this->mvNode[a]->mRadius;
          pch = strtok(NULL,"\"");
          pch = strtok(NULL,"\"");

          edge->mpEnd = this->mvNode[atoi(pch)-1];
          	  int b = atoi(pch)-1;
          edge->mRadius_end = this->mvNode[b]->mRadius;

   //       pch = strtok(NULL,">");
   //       pch = strtok(NULL," ");

    //   }
       f.getline(cstring, sizeof(cstring));
    }

  }
  f.close();

  std::copy( mvSegments.begin(), mvSegments.end(), std::back_inserter( mvEdge ) );
}

void Graph::read( const char* filename , int i){

  switch( i ){
    case 2:
      readMXF(filename);
    break;
    default:
      readVTK(filename);
    break;
  }

}

void Graph::resize(double scale){

  //shift
  double center[3] = {0.,0.,0.};
  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++){
    center[0] += this->mvNode[i]->position.x;// / this->mvEdge.size();
    center[1] += this->mvNode[i]->position.y;// / this->mvEdge.size();
    center[2] += this->mvNode[i]->position.z;// / this->mvEdge.size();
  }
  center[0] /= this->mvNode.size();
  center[1] /= this->mvNode.size();
  center[2] /= this->mvNode.size();
  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++){

    this->mvNode[i]->position.x -= center[0];
    this->mvNode[i]->position.y -= center[1];
    this->mvNode[i]->position.z -= center[2];

    this->mvNode[i]->setBoundingBox();
  }


   /*
  //debug output
  std::ofstream F;
  F.open("../../input/sin_graph1.txt_neu_center", std::ios::out|std::ios::app );

  F << "-------------" << std::endl;


    F << center[0] << "\t" ;
    F << center[1] << "\t" ;
    F << center[2] << std::endl ;


  F << "-------------" << std::endl;
  */





  //scale and set GL-Object
  double tmp_scale = 1./scale;
  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++){
    this->mvNode[i]->position.Multiply(tmp_scale);
  }



}

void Graph::setBoundingBoxList(BoundingBoxList *BoundingBoxList){

	this->mpBoundingBoxList = BoundingBoxList;

    if (mpBoundingBoxList)
        for (auto node: mvNode)
            mpBoundingBoxList->add(node);
}

void Graph::setBoundingBox(){

  //get Bounding Box of Graph

  this->mXmin = 0.;
  this->mXmax = 0.;
  this->mYmin = 0.;
  this->mYmax = 0.;
  this->mZmin = 0.;
  this->mZmax = 0.;

  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++){

    if( this->mXmin > this->mvNode[i]->position.x )
      this->mXmin = this->mvNode[i]->position.x;
    if( this->mXmax < this->mvNode[i]->position.x )
      this->mXmax = this->mvNode[i]->position.x;
    if( this->mYmin > this->mvNode[i]->position.y )
      this->mYmin = this->mvNode[i]->position.y;
    if( this->mYmax < this->mvNode[i]->position.y )
      this->mYmax = this->mvNode[i]->position.y;
    if( this->mZmin > this->mvNode[i]->position.z )
      this->mZmin = this->mvNode[i]->position.z;
    if( this->mZmax < this->mvNode[i]->position.z )
      this->mZmax = this->mvNode[i]->position.z;

  }

}

void Graph::shift(double x, int dim){

  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++){
    this->mvNode[i]->position.Add(x,dim);
    this->mvNode[i]->setBoundingBox();
  }

}

double Graph::calcDimensions(double &lobule_radius, double &lobule_height){
/* over all blood vessel segments
  lobule_radius = 0.; double tmp_z_min = 0.; double tmp_z_max = 0.;
  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++){
    double calc_new_r = sqrt(this->mvNode[i]->position.x*this->mvNode[i]->position.x+this->mvNode[i]->position.y*this->mvNode[i]->position.y)+this->mvNode[i]->mRadius;
    if ( calc_new_r > lobule_radius)
    	lobule_radius = calc_new_r;
    double calc_z_min = this->mvNode[i]-> position.z - this->mvNode[i]->mRadius;
    if ( calc_z_min < tmp_z_min)
      tmp_z_min = calc_z_min;
    double calc_z_max = this->mvNode[i]-> position.z + this->mvNode[i]->mRadius;
    if ( calc_z_max > tmp_z_max)
      tmp_z_max = calc_z_max;
  }

  lobule_height = tmp_z_max-tmp_z_min;

  return (lobule_height*0.5 + tmp_z_min);
*/

  lobule_radius = 0.;
  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++){
    if( this->mvNode[i]->mVesselType != ModelElementVesselSphere::PortalVein )
      continue;
    double calc_new_r = sqrt(this->mvNode[i]->position.x*this->mvNode[i]->position.x+this->mvNode[i]->position.y*this->mvNode[i]->position.y);
    if ( calc_new_r > lobule_radius)
      lobule_radius = calc_new_r;
  }

  double tmp_z_min = 0.; double tmp_z_max = 0.;
  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++){
    if( this->mvNode[i]->mVesselType != ModelElementVesselSphere::CentralVein )
      continue;
    double calc_z_min = this->mvNode[i]-> position.z;
    if ( calc_z_min < tmp_z_min)
      tmp_z_min = calc_z_min;
    double calc_z_max = this->mvNode[i]-> position.z;
    if ( calc_z_max > tmp_z_max)
      tmp_z_max = calc_z_max;
  }
  lobule_height = tmp_z_max-tmp_z_min;

  return (lobule_height*0.5 + tmp_z_min);

}

void Graph::initPressure(){

  for( unsigned int i = 0 ; i < this->mvNode.size() ; i++){

    if( this->mvNode[i]->mVesselType == ModelElementVesselSphere::CentralVein ){
      this->mvNode[i]->mPressure = 0;
    }
    else if( this->mvNode[i]->mVesselType == ModelElementVesselSphere::PortalVein ){
      this->mvNode[i]->mPressure = 10;
    }
    else
      this->mvNode[i]->mPressure;
  }

}

void Graph::calcPressure(){

  //set viscosity
  for( unsigned int i = 0 ; i < this->mvEdge.size() ; i++ )
    this->mvEdge[i]->setViscosity();

  //set inverse of resistance to flow
    for( unsigned int i = 0 ; i < this->mvEdge.size() ; i++ )
      this->mvEdge[i]->setInverseResistanceToFlow();

  //matrix (sparse)G and vector b,x
  int size = (int)this->mvNode.size();
  SparseMatrix<double> *sG = new SparseMatrix<double>(size,size);
#ifdef _MSC_VER
  double *b = new double[size];
  double *x = new double[size];
#else
  double b[size];
  double x[size];
#endif
  for( int i = 0 ; i < size ; i++){
    b[i] = 0.;
    x[i] = this->mvNode[i]->mPressure;
  }

  // CONSTRUCT LINEAR SYSTEM
  for( int i = 0 ; i < size ; i++ ){

    sG->resetRow(i);

    if( this->mvNode[i]->mvpNeighbor.size() == 0 )
      fprintf( stderr, "ERROR: Node without neighbors!  %i \n", i);


    if(    this->mvNode[i]->mVesselType  == ModelElementVesselSphere::CentralVein
        || this->mvNode[i]->mVesselType  == ModelElementVesselSphere::PortalVein
        || this->mvNode[i]->mvpNeighbor.size() == 0 ){

      // BOUNDARY: ROOT
      sG->set(i,i, 1);
      if( this->mvNode[i]->mVesselType  == ModelElementVesselSphere::PortalVein )
        b[i] = 10;
      else
        b[i] = 0;


    }else{

      double G, sumG = 0;

      for( unsigned int v = 0 ; v < this->mvNode[i]->mvpNeighbor.size() ; v++ ){

        G = this->mvNode[i]->mvpEdges[v]->mInverseResistenceToFlow;

        sG->set(i,this->mvNode[i]->mvpNeighbor[v]->mIndex, -G);

        sumG += G;
      }

      sG->set(i,i,sumG);
      b[i] = 0.;
    }
  }

/*
//Debug Johanne
  std::ofstream F;
  F.open( "G.txt", std::ios::out|std::ios::app );
  for( int i = 0 ; i < size ; i++ ){
    for( int j = 0 ; j < size ; j++ ){
      F << sG->get(i,j) << "\t";
    }
    F << std::endl ;
  }
  std::ofstream B;
    B.open( "b.txt", std::ios::out|std::ios::app );
    for( int i = 0 ; i < size ; i++ ){
        B << b[i] << std::endl;
    }
//Debug Johannes - end
*/

  //solve System
  Solver<double> *S = new Solver<double>( size, Solver<double>::BiCGSTAB, 1e-15, 10000 );//BiCGSTAB;CG
  S->PreconditionJacobi( sG, b);
  S->solve( sG, b, x);

/*
//Debug Johannes
  std::ofstream X;
  X.open( "x.txt", std::ios::out|std::ios::app );
  for( int i = 0 ; i < size ; i++ )
    X << x[i] << std::endl;
//Debug Johannes end
*/

  //set pressure to nodes
  for( int i = 0 ; i < size ; i++ )
    if( this->mvNode[i]->mVesselType != ModelElementVesselSphere::PortalVein &&
        this->mvNode[i]->mVesselType != ModelElementVesselSphere::CentralVein ){
      this->mvNode[i]->mPressure = x[i];
      this->mvNode[i]->color.red = 1-x[i]/10.;
      this->mvNode[i]->color.blue = x[i]/10.;
//      fprintf( stderr , "%i : %i  %2.3f \n", i , this->mvNode[i]->mIndex, x[i] );
    }

  // FREE MEMORY
  delete S;
  delete sG;
#ifdef _MSC_VER
  delete[] b;
  delete[] x;
#endif
}

void Graph::calcFlux(){

  for( unsigned int i = 0 ; i < this->mvEdge.size() ; i++ )
    this->mvEdge[i]->calcFlux();

}

void Graph::calcVolume(){

  //calc Volume in edge and for each node
  for( unsigned int i = 0 ; i < this->mvEdge.size() ; i++ ){

    double v_half = this->mvEdge[i]->calcVolume() * 0.5;

    this->mvEdge[i]->mpStart->mVolume += v_half;
    this->mvEdge[i]->mpEnd->mVolume += v_half;

  }

}

void Graph::updateConcentration(double dt){

  //not tested!!!!

  int M = (int)this->mvNode.size();

  double border = 10;

  // Construct Matrix for vascular Marker Concentration
  for( int i = 0 ; i < M ; i++ ){
    if( this->mvNode[i]->mVesselType == ModelElementVesselSphere::CentralVein ||
        this->mvNode[i]->mVesselType == ModelElementVesselSphere::PortalVein ){

      if( this->mvNode[i]->mVesselType == ModelElementVesselSphere::PortalVein )
        this->mvNode[i]->mdConcentration = ( border - this->mvNode[i]->mConcentration );
      else{
        this->mvNode[i]->mdConcentration = this->mvNode[i]->mConcentration;
      }

    }else{

      this->mvNode[i]->mdConcentration = 0;
      float outFlow=0, inFlow=0, diffusion=0;

      float flux = 0.;

      // for all flow from/to neighbors
      for( unsigned int v = 0 ; v < this->mvNode[i]->mvpNeighbor.size() ; v++ ){

        float flow;

        ModelElementVesselSphere * origin;
        ModelElementVesselSphere * destination;

        if( this->mvNode[i]->mPressure > this->mvNode[i]->mvpNeighbor[v]->mPressure ){
        // flow to neighbor: f_{i+1/2}
          origin      = this->mvNode[i];
          destination = this->mvNode[i]->mvpNeighbor[v];
          flow = -fabs(this->mvNode[i]->mvpEdges[v]->mFlux);
        }else{
          // flow from neighbor: f_{i-1/2}
          origin      = this->mvNode[i]->mvpNeighbor[v];
          destination = this->mvNode[i];
          flow = fabs(this->mvNode[i]->mvpEdges[v]->mFlux);
        }

        // FLUX
        float courantNumber = flow * dt / this->mvNode[i]->mVolume;
        flux += flow * origin->mConcentration ;

      }

      this->mvNode[i]->mdConcentration = flux * dt / this->mvNode[i]->mVolume;

    }
  }


  for( int i=0 ; i < M ; i++ ){
    if( this->mvNode[i]->mVesselType == ModelElementVesselSphere::CentralVein )
      this->mvNode[i]->mConcentration = 0;
    else if( this->mvNode[i]->mVesselType == ModelElementVesselSphere::PortalVein )
      this->mvNode[i]->mConcentration = 10;
    else
      this->mvNode[i]->mConcentration += this->mvNode[i]->mdConcentration;
  }

}

CSGraphEdge::CSGraphEdge(){


}

void CSGraphEdge::calcFlux(){

  //calc pressure difference
  double deltaP = fabs( this->mpStart->mPressure - this->mpEnd->mPressure );

  //calc flux
  this->mFlux = deltaP * this->mInverseResistenceToFlow;

}

void CSGraphEdge::setInverseResistanceToFlow(){

  double radius = 0.5 * ( this->mpStart->mRadius + this->mpEnd->mRadius );

  this->mInverseResistenceToFlow = M_PI/8. * pow( radius, 4) / ( this->mViscosity * this->mLength0 );

}

void CSGraphEdge::setViscosity(){

  double radius = 0.5 * ( this->mpStart->mRadius + this->mpEnd->mRadius );

  float visc_45 = 6 * exp(-0.085*radius*2.) + 3.1 - 2.44 * exp(-0.06*pow(radius*2., 0.645));
  float visc_vivo = (1 + (visc_45-1)*pow(radius/(radius-0.55),2))*pow(radius/(radius-0.55),2);
  this->mViscosity = 4e-6 * visc_45; // kPa * s

}

double CSGraphEdge::calcVolume(){

  double radius = 0.5 * ( this->mpStart->mRadius + this->mpEnd->mRadius );
  this->mVolume = M_PI*radius*radius*this->mLength0;

  return this->mVolume;

}

double CSGraphEdge::length(){

  this->mDirection[0] = this->mpEnd->position.x - this->mpStart->position.x;
  this->mDirection[1] = this->mpEnd->position.y - this->mpStart->position.y;
  this->mDirection[2] = this->mpEnd->position.z - this->mpStart->position.z;

  double norm = sqrt(this->mDirection[0]*this->mDirection[0]+this->mDirection[1]*this->mDirection[1]+this->mDirection[2]*this->mDirection[2]);

  if ( norm == 0 )
  {
      this->mDirection[0] = 0;
      this->mDirection[1] = 0;
      this->mDirection[2] = 0;
  }
  else
  {
      this->mDirection[0] /= norm;
      this->mDirection[1] /= norm;
      this->mDirection[2] /= norm;
  }

  return norm;
}
