///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSSphere.cpp                                                         //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-06-04 22:00:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "CSGLSphere.h"

#include "../../model/BasicDatatypes/Vector.h"
#include "../../model/BasicDatatypes/Color.h"

#ifdef __BUILD_MAC__
# include <OpenGL/OpenGL.h>
# include <GLUT/glut.h>
#else
# include <GL/glut.h>
#endif


/*!
  \brief Constructor
  \param where Initial position of the center of the sphere
  \param radius Radius of the sphere
  \param parent The parent object
*/
CSGLSphere::CSGLSphere( Vector3f * where, ARGBColor * color, double * radius )
  : CSGLObject(where)
{
  this->mpPosition = where;
  this->mpRadius   = radius;
  this->mpColor    = color;

  this->mSlices = 20;
  this->mStacks = 20;
}


/*!
  \brief Drawing routine

  Reimplemented from GLObject.  Draws a GLUT sphere at mPosition with
  radius *mpRadius and color *mpColor.
*/
void
CSGLSphere::draw()
{
  glPushMatrix();

  glColor4d( mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha );
  // Let's do this with glut, not the Geometry/Patch mechanism.
  glTranslated(mpPosition->x, mpPosition->y, mpPosition->z);
  glutSolidSphere(*mpRadius, this->mSlices, this->mStacks);
  //glutSwapBuffers();

  glPopMatrix();
}


/*!
  \brief Drawing routine

  Reimplemented from GLObject.  Set quality values of sphere //default 20 20
*/
void
CSGLSphere::setQuality(int slice, int stacks)
{
  this->mSlices  = slice;
  this->mStacks = stacks;
}

