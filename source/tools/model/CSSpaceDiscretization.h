////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:                                                                //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2017-02-10 15:29:24                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                        INRIA, Paris-Rocquencourt;                          //
//        and Hoehme Lab, Universitaet Leipzig.                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef CS_SPACE_DISCRETIZATION_H
#define CS_SPACE_DISCRETIZATION_H

#include "../../model/Elements/ModelElement.h"
#include "../../tools/model/BoundingBoxList.h"

#include <array>

class CSModel;
class CSVoxel;


class CSSpaceDiscretization
{
public:
    CSSpaceDiscretization( CSModel*, double voxelSize = 1. );
    CSSpaceDiscretization( const CSSpaceDiscretization & );
    CSSpaceDiscretization( CSSpaceDiscretization && );

    CSSpaceDiscretization & operator=( const CSSpaceDiscretization & rhs );

    ~CSSpaceDiscretization();

    CSVoxel * getVoxel( size_t x, size_t y, size_t z );
    CSVoxel * getVoxel( double posx, double posy, double posz );

    void setVoxelSize( const double & edgeLength ) { mVoxelSize=edgeLength; };
    const double & voxelSize() const { return mVoxelSize; };

    void update( bool updateLimits=false );

    const std::array<size_t,3> & dimensions() const { return mDimensions; };

    void destroy(CSVoxel * voxel);


protected:
    size_t indexOf( size_t x, size_t y, size_t z ) const
        {
            return z*mDimensions[0]*mDimensions[1] + y*mDimensions[0] + x;
        };

    std::array<size_t,3> mDimensions;
    BoundingBoxList *mpBBListCopy;
    Vector3f mMinCorner;
    double mVoxelSize;
    CSVoxel **mpVoxels;
    CSModel *mpModel;
};


#endif // CS_SPACE_DISCRETIZATION_H
