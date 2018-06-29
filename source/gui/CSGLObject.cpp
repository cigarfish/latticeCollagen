///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSGLObject.cpp                                                       //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-06-04 23:49:23                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "CSGLObject.h"

#include <iostream>


/*!
  \brief Constructor
  \param position The initial position of the center of the object
  \param parent The parent object
*/
CSGLObject::CSGLObject( Vector3f * position )
  : mpPosition(position),
    mpColor(NULL)
{}


/*!
  \brief Destructor
  As all member variables are pointers delegated to this object, the destructor
  should not destroy anything the variables used in accessor (set/get) routines!
 */
CSGLObject::~CSGLObject()
{}


/*!
  \brief Drawing routine

  This function has to be reimplemented by every derived class.
  It should contain all GL routines to draw the shape of the object with the chosen color.
*/
void
CSGLObject::draw()
{}
