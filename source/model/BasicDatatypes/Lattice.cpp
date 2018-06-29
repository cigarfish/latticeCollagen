///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Lattice.cpp                                                          //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-10-11 16:22:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "Lattice.h"

#include <iostream>
#include <stdio.h>

#include <math.h>

#include "../../tools/model/CSModelTools.h"

#define cut(a) (a<0 ? floor(a) : ceil(a))

Lattice::Lattice(int xmin, int ymin, int zmin, int xl, int yl, int zl){

  // Allocate memory
  mpppvLattice = new latticeContent***[xl];
  for (int i = 0; i < xl; ++i) {
    mpppvLattice[i] = new latticeContent**[yl];

    for (int j = 0; j < yl; ++j){
      mpppvLattice[i][j] = new latticeContent*[zl];

    	for( int k = 0 ; k < zl ; k++)
    		mpppvLattice[i][j][k] = new latticeContent;
    }
  }

  int index = 0;
  for( int i = 0 ; i < xl ; i++)
	  for( int j = 0 ; j < yl ; j++)
		  for( int k = 0 ; k < zl ; k++)
		  {
			  this->mpppvLattice[i][j][k]->index = index;
			  index++;

		  }

  this->mXLenght = xl;
  this->mYLenght = yl;
  this->mZLenght = zl;

  this->mXminCut = xmin;
  this->mYminCut = ymin;
  this->mZminCut = zmin;



}

void
Lattice::Add( ModelElement *modelElement){

  int x = cut(modelElement->position.x)+mXShift;
  int y = cut(modelElement->position.y)+mYShift;
  int z = cut(modelElement->position.z)+mZShift;

//Debug
//  fprintf(stderr, "%i \t %i %i %i\n", modelElement->mType,x,y,z);


  this->mpppvLattice[x][y][z]->mpppvLattice.push_back(modelElement);
  modelElement->mpIndexLattice = this->mpppvLattice[x][y][z];


}

void
Lattice::Remove( ModelElement *modelElement){
//TODO



}


std::vector< ModelElement * >
Lattice::InteractWith(ModelElement *modelElement){

  std::vector<ModelElement *> returnList;

  int x = cut(modelElement->position.x)+mXShift;
  int y = cut(modelElement->position.y)+mYShift;
  int z = cut(modelElement->position.z)+mZShift;

  for( int i = x-1 ; i< x+1 ; i++){
	  for( int j = y-1; j<y+1; j++){
		  for( int k = z-1; k<z+1; k++){
			  for( unsigned int l = 0 ; l < this->mpppvLattice[i][j][k]->mpppvLattice.size() ; l++){

				  if( (CSModelTools::GetDistance3D(modelElement->position,this->mpppvLattice[i][j][k]->mpppvLattice[l]->position) - modelElement->currentRadius - this->mpppvLattice[i][j][k]->mpppvLattice[l]->currentRadius) < 0 )
					  if( modelElement->mIndex < this->mpppvLattice[i][j][k]->mpppvLattice[l]->mIndex )
     					  returnList.push_back(this->mpppvLattice[i][j][k]->mpppvLattice[l]);

			  }
		  }
	  }
  }


 return returnList;
}

bool
Lattice::check(ModelElement *modelElement){

	 int x = cut(modelElement->position.x)+mXShift;
	  int y = cut(modelElement->position.y)+mYShift;
	  int z = cut(modelElement->position.z)+mZShift;

	 if ( this->mpppvLattice[x][y][z] == modelElement->mpIndexLattice )
	  return true;

	 //remove
	 unsigned int i;
	 for( i=0 ; i <  this->mpppvLattice[x][y][z]->mpppvLattice.size() ; i++)
	   if( this->mpppvLattice[x][y][z]->mpppvLattice[i] == modelElement )
	     break;
	 this->mpppvLattice[x][y][z]->mpppvLattice.erase( this->mpppvLattice[x][y][z]->mpppvLattice.begin() + i);


	  //add
	  this->mpppvLattice[x][y][z]->mpppvLattice.push_back(modelElement);
	  modelElement->mpIndexLattice = this->mpppvLattice[x][y][z];

	  return false;

}

