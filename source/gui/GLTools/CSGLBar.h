#pragma once

// added by Jieling, 30.05.2017
// draw bars connecting two points

#ifndef CS_GLTOOLS_CSGLBAR_H
#define CS_GLTOOLS_CSGLBAR_H

#include "../CSGLObject.h"

#include "../../model/BasicDatatypes/Vector.h"
#include "../../model/BasicDatatypes/Color.h"

class CSGLBar : public CSGLObject
{
public:
	CSGLBar(Vector3f * endpoint1, Vector3f * endpoint2, ARGBColor * color);

	void draw();
	void changeP1(Vector3f *P) { p1 = P; }
	void changeP2(Vector3f *P) { p2 = P; }
	//void setRaduis(double r) { mpRadius = r; }
	//void setQuality(int slice, int stacks);

protected:
	Vector3f *p1;
	Vector3f *p2;
	double mpRadius;
	int mSlices;
	int mStacks;
};

#endif
