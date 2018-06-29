/*
 *  TriangulatedCell.cpp
 *  deformable Cells
 *
 *  Created by Johannes Neitsch on 05/06/2012.
 *  Copyright 2012 The Drasdo Group, IZBI Leipzig. All rights reserved.
 *
 */

#include <QVector3D>

#include "TriangulatedCell.h"

#include "../../model/Cell/CellTriangulated.h"
#include "../tools/math/mathematics.h"

#include "../BasicDatatypes/MassPoint.h"
#include "../BasicDatatypes/Spring.h"
#include "../BasicDatatypes/Triangle.h"

#ifdef __BUILD_MAC__
# include <OpenGL/OpenGL.h>
# include <GLUT/glut.h>
# define glutSolidCylinder(radio, altura, slices, stacks) gluCylinder(gluNewQuadric(), radio, radio, altura, slices, stacks)
#else
# include <GL/freeglut.h>
#endif


TriangulatedCell::TriangulatedCell( CellTriangulated *cell_ )
  : CSGLObject(&cell_->position)
{
	cell = cell_;
}



/*!
  \brief Drawing routine

  Reimplemented from GLObject.  Draws a GLUT sphere at mPosition with radius mRadius and color mColor.
*/
void
TriangulatedCell::draw(){

  for( unsigned int i = 0 ; i < cell->spring.size(); i++){
    glPushMatrix();

    glColor4d( cell->spring[i]->r/255.,
               cell->spring[i]->g/255.,
               cell->spring[i]->b/255.,
               cell->spring[i]->transparency);     //sets the current colour to red


    double tmp_rot[3];
    double tmp_r = rotateCylinder(tmp_rot,cell->spring[i]);

    glTranslatef( cell->spring[i]->start->position[0],
                  cell->spring[i]->start->position[1],
                  cell->spring[i]->start->position[2] );

    glRotatef(tmp_rot[0],tmp_rot[1], tmp_rot[2], 0.0);

    glutSolidCylinder(cell->spring[i]->radius,tmp_r, 5, 5);//radius 0.005

    glPopMatrix();

  }
  
  if( cell->status == 200 || cell->status == 2000){//200 optical stretcher, 2000 force of deformation after watershed in tissue
    for( unsigned int i = 0 ; i < cell->mass_huell.size() ;i++){

      glPushMatrix();

      glColor4d( 0,1,1,1);     //sets the current colour to red

      double tmp_rot[3];
      double tmp[3] = { cell->mass_huell[i]->position[0]+cell->mass_huell[i]->F_force_profile[0],
                        cell->mass_huell[i]->position[1]+cell->mass_huell[i]->F_force_profile[1],
                        cell->mass_huell[i]->position[2]+cell->mass_huell[i]->F_force_profile[2]};

      double tmp_r = rotateCylinder(tmp_rot,cell->mass_huell[i]->position,tmp);
//        double tmp_r = rotateCylinder(tmp_rot,cell->spring[i]);

      glTranslatef( cell->mass_huell[i]->position[0],
                    cell->mass_huell[i]->position[1],
                    cell->mass_huell[i]->position[2] );

      glRotatef(tmp_rot[0],tmp_rot[1], tmp_rot[2], 0.0);

      glutSolidCylinder(0.005,tmp_r, 5, 5);//radius 0.005

      glPopMatrix();

    }
  }


  for( unsigned int i = 0 ; i < cell->mass.size() ;i++){

    glPushMatrix();
/*if( i < 49 ){
    glColor4d( cell->mass[i]->r/255.,
               cell->mass[i]->g/255.,
               cell->mass[i]->b/255.,
               cell->mass[i]->transparency);     //sets the current colour to red
//}
/*
else{
    glColor4d( 1,
               cell->mass[i]->g/255.,
               cell->mass[i]->b/255.,
               cell->mass[i]->transparency);     //sets the current colour to red
}
*/
    glTranslated( cell->mass[i]->position[0],
                  cell->mass[i]->position[1],
                  cell->mass[i]->position[2]);
  //Debug-Johannes
//if( i != 49)
    glutSolidSphere(this->cell->mass[i]->radius, 5, 5);
//else
//	glutSolidSphere(0.1, 5, 5);
    glPopMatrix();

  }

	if( !cell->pressure_high){ 
	for( unsigned int i = 0 ; i <cell->triangle.size() ; i++){//
	//	if( cell->triangle[i]->interaction != 0){
        glPushMatrix();

		glBegin(GL_TRIANGLES);  //tells OpenGL that we're going to start drawing triangles
		
		glColor4d(cell->triangle[i]->r/255.,cell->triangle[i]->g/255.,cell->triangle[i]->b/255., cell->triangle[i]->transparency);     //sets the current colour to red
		
		glVertex3f(cell->triangle[i]->points[0]->position[0],cell->triangle[i]->points[0]->position[1],cell->triangle[i]->points[0]->position[2]);  //specifies the first vertex of our triangle
  
		glColor4d(cell->triangle[i]->r/255.,cell->triangle[i]->g/255.,cell->triangle[i]->b/255., cell->triangle[i]->transparency);     //sets the current colour to red
		
		glVertex3f(cell->triangle[i]->points[1]->position[0],cell->triangle[i]->points[1]->position[1],cell->triangle[i]->points[1]->position[2]);  //specifies the first vertex of our triangle
  
		glColor4d(cell->triangle[i]->r/255.,cell->triangle[i]->g/255.,cell->triangle[i]->b/255., cell->triangle[i]->transparency);     //sets the current colour to red
	
		glVertex3f(cell->triangle[i]->points[2]->position[0],cell->triangle[i]->points[2]->position[1],cell->triangle[i]->points[2]->position[2]);  //specifies the first vertex of our triangle
        glEnd();                //tells OpenGL that we've finished drawing

		glPopMatrix();
//	}
  










  }
	}else{

		for( unsigned int i = 0 ; i < cell->triangle.size(); i++){
	      glPushMatrix();
	//	  if( cell->triangle[i]->interaction != 0){
		glBegin(GL_TRIANGLES);  //tells OpenGL that we're going to start drawing triangles
		
		glColor4d(1.,0.,0., cell->triangle[i]->transparency);     //sets the current colour to red
	
		glVertex3f(cell->triangle[i]->points[0]->position[0],cell->triangle[i]->points[0]->position[1],cell->triangle[i]->points[0]->position[2]);  //specifies the first vertex of our triangle
  
			

		glColor4d(1.,0.,0., cell->triangle[i]->transparency);     //sets the current colour to red
	

		glVertex3f(cell->triangle[i]->points[1]->position[0],cell->triangle[i]->points[1]->position[1],cell->triangle[i]->points[1]->position[2]);  //specifies the first vertex of our triangle
  
			

		glColor4d(1.,0.,0., cell->triangle[i]->transparency);     //sets the current colour to red
	

		glVertex3f(cell->triangle[i]->points[2]->position[0],cell->triangle[i]->points[2]->position[1],cell->triangle[i]->points[2]->position[2]);  //specifies the first vertex of our triangle
        glEnd();                //tells OpenGL that we've finished drawing

		glPopMatrix();
	//	  }
	}
	}


  //draw normalvector of triangles
  
  	for( unsigned int i = 0 ; i < cell->triangle.size(); i++){
        glPushMatrix();

       glColor4d( 0,
                   1,
                   1,
                   1);     //sets the current colour to red

       double m[3];

       mean(this->cell->triangle[i]->points[0]->position,
            this->cell->triangle[i]->points[1]->position,
            this->cell->triangle[i]->points[2]->position,
            m);


       double z[3];
       z[0] = m[0] + this->cell->triangle[i]->mNormalvector[0];
       z[1] = m[1] + this->cell->triangle[i]->mNormalvector[1];
       z[2] = m[2] + this->cell->triangle[i]->mNormalvector[2];

        double tmp_rot[3];
        double tmp_r = rotateCylinder(tmp_rot,m,z,10);

        glTranslatef( m[0],
                      m[1],
                      m[2] );

        glRotatef(tmp_rot[0],tmp_rot[1], tmp_rot[2], 0.0);
		
        glutSolidCylinder(0.004,tmp_r, 5, 5);//radius 0.005

        glPopMatrix();

	}
  



}
