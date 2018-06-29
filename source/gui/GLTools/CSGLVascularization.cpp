///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSGLVascularization.cpp                                              //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-06-05 19:10:37                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include <QVector3D>

#include "CSGLVascularization.h"

#include "../../model/Model/ModelVascularization/CSVesselGraph.h"


#include "../tools/math/mathematics.h"

#ifdef __BUILD_MAC__
# include <OpenGL/OpenGL.h>
# include <GLUT/glut.h>
# define glutSolidCylinder(radio, altura, slices, stacks) gluCylinder(gluNewQuadric(), radio, radio, altura, slices, stacks)
#else
# include <GL/freeglut.h>
#endif


CSGLVascularization::CSGLVascularization( CSVesselGraph *CSvg )
  : CSGLObject(NULL)
{
  vg = CSvg;
  mRadius = 0;
  voronoi = NULL;
}

void
CSGLVascularization::setVoronoi( VoronoiDiagram *voronoi )
{
  this->voronoi = voronoi;
}



/*!
  \brief Drawing routine

  Reimplemented from GLObject.  Draws a GLUT sphere at mPosition with radius mRadius and color mColor.
*/
void
CSGLVascularization::draw()
{
  double scale = 10.;

/*
  for( int i = 0 ; i < vg->countVesselNodes ; i++){

    glPushMatrix();
    glColor4d( 1.,1.,1., 0.5 );

    // Let's do this with glut, not the Geometry/Patch mechanism.
    glTranslated( vg->vesselNodes[i]->position[0]/scale, vg->vesselNodes[i]->position[1]/scale, vg->vesselNodes[i]->position[2]/scale);
    glutSolidSphere(1./scale, 5, 5);
    //glutSwapBuffers();

    glPopMatrix();

  }
*/

  for( int i = 0 ; i < vg->countVesselSegments; i++){
    glPushMatrix();

    glColor4d( 1.,
               1.-((vg->vesselSegments[i]->vesselNodes[0]->marker+vg->vesselSegments[i]->vesselNodes[1]->marker)/14.),
               1.-((vg->vesselSegments[i]->vesselNodes[0]->marker+vg->vesselSegments[i]->vesselNodes[1]->marker)/14.),
               1.);     //sets the current colour to red


    double tmp_rot[3];

    double tmp_r = rotateCylinder(tmp_rot,vg->vesselSegments[i]->vesselNodes[0]->position,vg->vesselSegments[i]->vesselNodes[1]->position,scale);

    glTranslatef( vg->vesselSegments[i]->vesselNodes[0]->position[0]/scale,
                  vg->vesselSegments[i]->vesselNodes[0]->position[1]/scale,
                  vg->vesselSegments[i]->vesselNodes[0]->position[2]/scale );

    glRotatef(tmp_rot[0],tmp_rot[1], tmp_rot[2], 0.0);

    glutSolidCylinder(vg->vesselSegments[i]->radius/scale/voronoi->LATTICE_CONSTANT,tmp_r, 5, 5);//radius 0.005

    glPopMatrix();

  }

  scale = 10./voronoi->LATTICE_CONSTANT;




  for( int i = 0 ; i <  this->voronoi->getCountVertices(); i++){

    glPushMatrix();

    glColor4d( 1-this->voronoi->get(i)->conc/10.  ,1.,1., 0.2 );

    POSITION_T *t;
    t = this->voronoi->get(i)->pos();

    // Let's do this with glut, not the Geometry/Patch mechanism.
    glTranslated( ((double)t[0])/scale/voronoi->LATTICE_CONSTANT, (double)t[1]/voronoi->LATTICE_CONSTANT/scale, (double)t[2]/voronoi->LATTICE_CONSTANT/scale );
    glutSolidSphere(voronoi->LATTICE_CONSTANT/40, 4, 4);
    //glutSwapBuffers();

    glPopMatrix();
  }



}
