///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  BoundingBox.h                                                        //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-07-23 23:53:16                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef CS_BOUNDING_BOX_H
#define CS_BOUNDING_BOX_H

struct BoundingBox
{
  double xmin, xmax;
  double ymin, ymax;
  double zmin, zmax;
};

#endif // CS_BOUNDING_BOX_H
