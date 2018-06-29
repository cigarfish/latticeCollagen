////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSGLCube.cpp                                                  //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2016-10-14 14:29:36                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris;                                        //
//      Hoehme Lab, Universitaet Leipzig.                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "CSGLCube.h"

#ifdef __BUILD_MAC__
# include <OpenGL/OpenGL.h>
# include <GLUT/glut.h>
#else
# include <GL/glut.h>
#endif


CSGLCube::CSGLCube(Vector3f* position, ARGBColor* color, double* sideLength)
    : CSGLObject( position ),
      mpSideLength( sideLength )
{
    mpColor = color;
}

void
CSGLCube::draw()
{
    // if ( !mVisible )
    //     return;

    glPushMatrix();

    glColor4d( mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha );
    glTranslated( mpPosition->x+*mpSideLength/2, mpPosition->y+*mpSideLength/2, mpPosition->z+*mpSideLength/2 );
    glutSolidCube( *mpSideLength );

    glPopMatrix();
}
