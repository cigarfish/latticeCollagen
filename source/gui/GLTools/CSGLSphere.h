///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSGLSphere.h                                                         //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-06-04 19:10:37                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef CS_GLTOOLS_CSGLSPHERE_H
#define CS_GLTOOLS_CSGLSPHERE_H

#include "../CSGLObject.h"

#include "../../model/BasicDatatypes/Vector.h"
#include "../../model/BasicDatatypes/Color.h"

/*!
  \brief A GL sphere
*/
class CSGLSphere : public CSGLObject
{
 public:
  CSGLSphere( Vector3f * where, ARGBColor *color, double * radius );
  
  void draw();
  void setQuality(int slice, int stacks);

 protected:
  double * mpRadius;
  int mSlices;
  int mStacks;
};


#endif //  CS_GLTOOLS_CSGLSPHERE_H
