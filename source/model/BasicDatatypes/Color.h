///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Color.h                                                              //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-06-04 22:15:00                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef COLOR_H
#define COLOR_H


struct RGBColor
{
  unsigned char red;
  unsigned char green;
  unsigned char blue;
};


struct ARGBColor
{
  float alpha;
  float red;
  float green;
  float blue;
};


#endif //COLOR_H
