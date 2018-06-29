////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  isnan.h                                                       //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2014-11-28 12:41:43                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_ISNAN_H
#define CS_ISNAN_H

#include <cmath>

#if (WIN32)
#  define isnan(a) _isnan(a)
#  define isinf(a) !_finite(a)
#else
#  undef isnan
#  undef isinf
#  define isnan(a) (a!=a)
#  define isinf(a) std::isinf(a)
#endif


#endif // CS_ISNAN_H
