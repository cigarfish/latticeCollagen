///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSGLTriangle.h                                                       //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2013-09-24 11:46:37                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef CS_GLTOOLS_CSGLTRIANGLE_H
#define CS_GLTOOLS_CSGLTRIANGLE_H

#include "../CSGLObject.h"

#include "../../model/BasicDatatypes/Vector.h"
#include "../../model/BasicDatatypes/Color.h"

/*!
  \brief A GL triangle
*/
class CSGLTriangle : public CSGLObject
{
 public:
  CSGLTriangle( Vector3f * where, ARGBColor *color, double (* points)[3][3] );

  void draw();

 protected:
  double (* mpPoints)[3][3];

};


#endif //  CS_GLTOOLS_CSGLTRIANGLE_H
