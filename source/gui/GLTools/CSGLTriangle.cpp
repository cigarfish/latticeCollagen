///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSGLTriangle.cpp                                                     //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2013-09-24 11:51:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "CSGLTriangle.h"

#include "../../model/BasicDatatypes/Vector.h"
#include "../../model/BasicDatatypes/Color.h"

#include "../../tools/math/mathematics.h"

#ifdef __BUILD_MAC__
# include <OpenGL/OpenGL.h>
# include <GLUT/glut.h>
# define glutSolidCylinder(radio, altura, slices, stacks) gluCylinder(gluNewQuadric(), radio, radio, altura, slices, stacks)
#else
# include <GL/freeglut.h>
#endif

/*!
  \brief Constructor
  \param where Initial position of the center of the sphere
  \param radius Radius of the sphere
  \param parent The parent object
*/
CSGLTriangle::CSGLTriangle( Vector3f * where , ARGBColor * color, double (* points)[3][3] )
  : CSGLObject(where)
{
  this->mpPosition = where;
  this->mpColor    = color;
  this->mpPoints    = points;

}


/*!
  \brief Drawing routine

  Reimplemented from GLObject.  Draws a GLUT sphere at mPosition with
  radius *mpRadius and color *mpColor.
*/
void
CSGLTriangle::draw()
{


  glPushMatrix();
   glBegin(GL_TRIANGLES);  //tells OpenGL that we're going to start drawing triangles
   glColor4d( mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha );
   glVertex3f((*mpPoints)[0][0],(*mpPoints)[0][1],(*mpPoints)[0][2]);  //specifies the first vertex of our triangle
   glColor4d( mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha );
   glVertex3f((*mpPoints)[1][0],(*mpPoints)[1][1],(*mpPoints)[1][2]);  //specifies the first vertex of our triangle
   glColor4d( mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha );
   glVertex3f((*mpPoints)[2][0],(*mpPoints)[2][1],(*mpPoints)[2][2]);  //specifies the first vertex of our triangle
   glEnd();                //tells OpenGL that we've finished drawing
  glPopMatrix();

  glPushMatrix();
     glBegin(GL_TRIANGLES);  //tells OpenGL that we're going to start drawing triangles
     glColor4d( mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha );
     glVertex3f((*mpPoints)[0][0],(*mpPoints)[0][1],(*mpPoints)[0][2]);  //specifies the first vertex of our triangle
     glColor4d( mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha );
     glVertex3f((*mpPoints)[2][0],(*mpPoints)[2][1],(*mpPoints)[2][2]);  //specifies the first vertex of our triangle
     glColor4d( mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha );
     glVertex3f((*mpPoints)[1][0],(*mpPoints)[1][1],(*mpPoints)[1][2]);  //specifies the first vertex of our triangle
     glEnd();                //tells OpenGL that we've finished drawing
    glPopMatrix();


/*
  double normalVector[3];
  double m[3];
  double m2[3];
  crossProduct(mPoints[0],mPoints[1],mPoints[2],normalVector);

  mean(mPoints[0],mPoints[1],mPoints[2], m);
  m2[0] = m[0] + normalVector[0];
  m2[1] = m[1] + normalVector[1];
  m2[2] = m[2] + normalVector[2];

  glPushMatrix();
   glColor4d( mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha );
   glTranslated(m2[0],m2[1], m2[2]);
   glutSolidSphere(0.1, 10, 10);
  glPopMatrix();

*/

/*
  glPushMatrix();
   glColor4d( mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha );

   double tmp_rot[3];
   double tmp_r = rotateCylinder( tmp_rot , m ,m2, 1 );
//   double tmp_r = rotateCylinder(tmp_rot,cell->spring[i]);
   glTranslatef( m[0],
                 m[1],
                 m[2] );

   glRotatef(tmp_rot[0],tmp_rot[1], tmp_rot[2], 0.0);

    glutSolidCylinder(0.1,tmp_r, 5, 5);//radius 0.005

  glPopMatrix();

*/

}
