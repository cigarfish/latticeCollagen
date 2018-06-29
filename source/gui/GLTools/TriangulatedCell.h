/*
 *  TriangulatedCell.h
 *  deformable Cells
 *
 *  Created by Johannes Neitsch on 05/06/2012.
 *  Copyright 2012 The Drasdo Group, IZBI Leipzig. All rights reserved.
 *
 */

#ifndef TRIANGULATED_CELL_H
#define TRIANGULATED_CELL_H

#include "../CSGLObject.h"

class CellTriangulated;


class TriangulatedCell : public CSGLObject
{
public:
	TriangulatedCell( CellTriangulated *cell );
	void draw();


protected:
	double mRadius;


public:
	CellTriangulated *cell;

};


#endif // TRIANGULATED_CELL_H
