////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSVoxel.h                                                     //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2016-10-14 11:52:45                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//      Hoehme Lab, Universitaet Leipzig.                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef CS_VOXEL_H
#define CS_VOXEL_H

#include "ModelElement.h"

#include "../../gui/GLTools/CSGLCube.h"


class CSVoxel : public ModelElement
{
public:
    CSVoxel() : CSVoxel(0.,0.,0.) {};

    CSVoxel( double x, double y, double z, double sideLength=1 )
        : ModelElement( x, y, z, ModelElement::TypeVoxel ),
          contains(ModelElement::TypeUnspecified),
          visited(false),
          mSideLength(sideLength)
        {
            mpGLObject = new CSGLCube( &position, &color, &mSideLength);
            mBoundingBox.xmin = position.x;
            mBoundingBox.xmax = position.x + mSideLength;

            mBoundingBox.ymin = position.y;
            mBoundingBox.ymax = position.y + mSideLength;

            mBoundingBox.zmin = position.z;
            mBoundingBox.zmax = position.z + mSideLength;
        };

    virtual ~CSVoxel() {};

    double size() const { return mSideLength; };
    void setSize( const double & length ) { mSideLength = length; };

    virtual BoundingBox * boundingBox()
        {
            return &mBoundingBox;
        }

    void Reset() { contains = ModelElement::TypeUnspecified; visited = false; };

    ModelElement::Type contains;

    bool visited;
protected:
    double mSideLength;
};

#endif // CS_VOXEL_H
