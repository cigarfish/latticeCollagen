#include "CSGLBar.h"

#include <cmath>

#include "../../model/BasicDatatypes/Vector.h"
#include "../../model/BasicDatatypes/Color.h"


#ifdef __BUILD_MAC__
# include <OpenGL/OpenGL.h>
# include <GLUT/glut.h>
# define glutSolidCylinder(radio, altura, slices, stacks) gluCylinder(gluNewQuadric(), radio, radio, altura, slices, stacks)
#else
# include <GL/freeglut.h>
# include <GL/glut.h>
#endif


CSGLBar::CSGLBar(Vector3f * endpoint1, Vector3f * endpoint2, ARGBColor * color)
	: CSGLObject(endpoint1)
{
	this->mpPosition = endpoint1;
	this->mpColor = color;

	this->p1 = endpoint1;
	this->p2 = endpoint2;

	this->mSlices = 20;
	this->mStacks = 20;

	//this->mpRadius = 0.005; // default radius of the cylinder is 0.005
}

void
CSGLBar::draw()
{
	glPushMatrix();

	glColor4d(mpColor->red, mpColor->green, mpColor->blue, mpColor->alpha);
	
	double tmp_rot[3];
	
	double vx = this->p2->x - this->p1->x; // x2 - x1
	double vy = this->p2->y - this->p1->y; // y2 - y1
	double vz = this->p2->z - this->p1->z; // z2 - z1

	//std::cout << "	-> GLBar p1: " << p1->x << ", " << p1->y << ", " << p1->z
	//		  << " p2: " << p2->x << ", " << p2->y << ", " << p2->z << std::endl;
	
	if (vz == 0)
		vz = .0001;

	double tmp_r = sqrt(vx*vx + vy*vy + vz*vz);
	double ax = 57.2957795*acos(vz / tmp_r);
	//std::cout << "	-> GLBar tmp_r: " << tmp_r << ", ax: " << ax << std::endl;
	double rx = -vy / tmp_r;
	double ry =  vx / tmp_r;
	if (vz < 0.0)
	{
		ax = -ax;
		rx = -rx;
		ry = -ry;
	}
	//double rx = -vy*vz;
	//double ry = vx*vz;

	tmp_rot[0] = ax;
	tmp_rot[1] = rx;
	tmp_rot[2] = ry;
	//std::cout << "	-> GLBar tmp_rot: " << tmp_rot[0]  << ", " << tmp_rot[1] << ", " << tmp_rot[2] << std::endl;

	glTranslatef(this->p1->x,
				 this->p1->y,
				 this->p1->z);

	glRotatef(tmp_rot[0], tmp_rot[1], tmp_rot[2], 0.0);

	glutSolidCylinder(0.05, tmp_r, 5, 5);//radius 0.005

	glPopMatrix();
}