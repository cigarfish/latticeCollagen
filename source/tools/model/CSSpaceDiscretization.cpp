////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSSpaceDiscretization.c                                       //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2017-02-10 16:56:24                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                        INRIA, Paris-Rocquencourt;                          //
//        and Hoehme Lab, Universitaet Leipzig.                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "CSSpaceDiscretization.h"

#include "../../model/Model/CSModel.h"
#include "../../model/Elements/CSVoxel.h"

#include "BoundingBoxList.h"

#include <cstring>


CSSpaceDiscretization::CSSpaceDiscretization( CSModel* model, double voxelSize )
    : mpBBListCopy( model->mpBBList ),
      mVoxelSize(voxelSize),
      mpVoxels(nullptr),
      mpModel(model)
{
    update(true);
}


CSSpaceDiscretization::CSSpaceDiscretization( const CSSpaceDiscretization & other )
    : mDimensions( other.mDimensions ),
      mpBBListCopy( other.mpBBListCopy ),
      mMinCorner( other.mMinCorner ),
      mVoxelSize(other.mVoxelSize),
      mpModel(other.mpModel)
{// ToDo
    // copy the voxels in mpVoxels
}

CSSpaceDiscretization::CSSpaceDiscretization( CSSpaceDiscretization && other )
    : mpBBListCopy( other.mpBBListCopy ),
      mVoxelSize(other.mVoxelSize),
      mpVoxels(std::move(other.mpVoxels)),
      mpModel(other.mpModel)
{// ToDo
    other.mpBBListCopy = nullptr;
    other.mpVoxels     = nullptr;
}

CSSpaceDiscretization &
CSSpaceDiscretization::operator=( const CSSpaceDiscretization & rhs )
{
    mDimensions  = rhs.mDimensions;
    mpBBListCopy = rhs.mpBBListCopy;
    mMinCorner   = rhs.mMinCorner;
    mVoxelSize   = rhs.mVoxelSize;
    mpModel      = rhs.mpModel;
    // ToDo: copy mpVoxels
    return *this;
}


CSSpaceDiscretization::~CSSpaceDiscretization()
{
    if ( mpVoxels )
        for ( size_t i=0; i<mDimensions[0]*mDimensions[1]*mDimensions[2]; ++i )
            if (mpVoxels[i]) delete mpVoxels[i];

    delete[] mpVoxels;
}


CSVoxel *
CSSpaceDiscretization::getVoxel( size_t x, size_t y, size_t z )
{
    if ( x >= mDimensions[0] || y >= mDimensions[1] || z >= mDimensions[2] )
        return NULL;

    size_t i = indexOf( x, y, z );

    if ( !mpVoxels[i] )
    {
        // calculate position of voxel
        double posx = mMinCorner[0] + x*mVoxelSize;
        double posy = mMinCorner[1] + y*mVoxelSize;
        double posz = mMinCorner[2] + z*mVoxelSize;
        mpVoxels[i] = new CSVoxel(posx, posy, posz, mVoxelSize);
    }

    return mpVoxels[i];
}


CSVoxel *
CSSpaceDiscretization::getVoxel( double posx, double posy, double posz )
{
    size_t xext = static_cast<size_t>(std::floor( ( posx - mMinCorner[0] ) / mVoxelSize ));
    size_t yext = static_cast<size_t>(std::floor( ( posy - mMinCorner[1] ) / mVoxelSize ));
    size_t zext = static_cast<size_t>(std::floor( ( posz - mMinCorner[2] ) / mVoxelSize ));

    if ( xext >= mDimensions[0] || yext >= mDimensions[1] || zext >= mDimensions[2] )
        return NULL;

    size_t i = indexOf( xext, yext, zext );

    if ( !mpVoxels[i] )
    {
        // align position to voxelization
        double posx_corr = xext * mVoxelSize + mMinCorner[0];
        double posy_corr = yext * mVoxelSize + mMinCorner[1];
        double posz_corr = zext * mVoxelSize + mMinCorner[2];


        mpVoxels[i] = new CSVoxel( posx_corr, posy_corr, posz_corr, mVoxelSize );
    }

    return mpVoxels[i];
}


void
CSSpaceDiscretization::update( bool updateLimits )
{
    if (mpModel)
        if ( mpModel->mpBBList && updateLimits )
        {
            if ( mpVoxels )
                for ( size_t i=0; i<mDimensions[0]*mDimensions[1]*mDimensions[2]; ++i )
                    if (mpVoxels[i]) delete mpVoxels[i];

            mpBBListCopy = new BoundingBoxList( * mpModel->mpBBList );

            // assuming, mpBBList already is updated... else we should do:
            mpBBListCopy->update();

            // get world dimensions
            BoundingBox worldLimits( mpBBListCopy->absoluteLimits() );

            double extent = std::floor( std::fabs( worldLimits.xmin/mVoxelSize ) );
            if ( worldLimits.xmin <= 0. )
                extent *= -1;
            else
                extent += 1.;
            worldLimits.xmin = extent * mVoxelSize;
            long leftLimit = (long)extent;

            extent = std::ceil(fabs( worldLimits.xmax/mVoxelSize ) );
            if ( worldLimits.xmax <= 0. )
            {
                extent -= 1;
                extent *= -1;
            }
            worldLimits.xmax = extent * mVoxelSize;
            mDimensions[0] = (int)std::fabs(extent - leftLimit);


            extent = std::floor( std::fabs( worldLimits.ymin/mVoxelSize ) );
            if ( worldLimits.ymin <= 0. )
                extent *= -1;
            else
                extent += 1.;
            worldLimits.ymin = extent * mVoxelSize;
            leftLimit = (long)extent;

            extent = std::ceil(fabs( worldLimits.ymax/mVoxelSize ) );
            if ( worldLimits.ymax <= 0. )
            {
                extent -= 1;
                extent *= -1;
            }
            worldLimits.ymax = extent * mVoxelSize;
            mDimensions[1] = (int)std::fabs(extent - leftLimit);


            extent = std::floor( std::fabs( worldLimits.zmin/mVoxelSize ) );
            if ( worldLimits.zmin <= 0. )
                extent *= -1;
            else
                extent += 1.;
            worldLimits.zmin = extent * mVoxelSize;
            leftLimit = (long)extent;

            extent = std::ceil(fabs( worldLimits.zmax/mVoxelSize ) );
            if ( worldLimits.zmax <= 0. )
            {
                extent -= 1;
                extent *= -1;
            }
            worldLimits.zmax = extent * mVoxelSize;
            mDimensions[2] = (int)std::fabs(extent - leftLimit);

            size_t total = mDimensions[0]*mDimensions[1]*mDimensions[2];

            mpVoxels = new CSVoxel*[total] {};
//            mpVoxels = (CSVoxel **) malloc(total*sizeof(CSVoxel*));
//            memset( mpVoxels, 0, total*sizeof(CSVoxel*));

            mMinCorner.x = worldLimits.xmin;
            mMinCorner.y = worldLimits.ymin;
            mMinCorner.z = worldLimits.zmin;
        }
}


void
CSSpaceDiscretization::destroy(CSVoxel * voxel)
{
    size_t xext = static_cast<size_t>(std::floor( ( voxel->position[0] - mMinCorner[0] ) / mVoxelSize ));
    size_t yext = static_cast<size_t>(std::floor( ( voxel->position[0] - mMinCorner[1] ) / mVoxelSize ));
    size_t zext = static_cast<size_t>(std::floor( ( voxel->position[0] - mMinCorner[2] ) / mVoxelSize ));

    // do not delete if voxel is not in mpVoxels!
    if ( xext > mDimensions[0] || yext > mDimensions[1] || zext > mDimensions[2] )
        return;

    size_t index = indexOf( xext, yext, zext );

    if ( mpVoxels[index] )
    {
        if ( voxel != mpVoxels[index] )
            return;

        delete voxel;
        mpVoxels[index] = nullptr;
    }
}
