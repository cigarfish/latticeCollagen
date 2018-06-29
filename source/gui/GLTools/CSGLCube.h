////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSGLCube.h                                                    //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2016-10-14 14:24:49                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris;                                        //
//      Hoehme Lab, Universitaet Leipzig.                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_GLTOOLS_CSGLCUBE_H
#define CS_GLTOOLS_CSGLCUBE_H

#include "../CSGLObject.h"

class CSGLCube : public CSGLObject
{
public:
    CSGLCube( Vector3f *position, ARGBColor *color, double * sideLength );

    void draw();

protected:
    double * mpSideLength;
};


#endif // CS_GLTOOLS_CSGLCUBE_H
