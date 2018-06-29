////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  LobuleLesionVoxelization.cpp                                  //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2017-04-12 15:50:07                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                        INRIA, Paris-Rocquencourt;                          //
//        and Hoehme Lab, Universitaet Leipzig.                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "LobuleLesionVoxelization.h"

#include "../../Interactions/CSInteractionHertz.h"

#include <queue>


void
LobuleLesionVoxelization::exec()
{
    resetVoxels();
    mvRegion.clear();

    ModelCellsSpherical * model = static_cast<ModelCellsSpherical *>(mpModel);
    std::cout <<  "# starting area discretization" << std::endl;
    //const double voxelSize = .4;

    delete mpBBListCopy;
    mpBBListCopy = new BoundingBoxList( *mpModel->mpBBList );

    // the queue for the region growing algorithm
    std::queue<CSVoxel *> q;

    CSInteractionHertz overlapTest;

    std::vector<CSVoxel *> voxels;

    const double radius = 1.5 + model->biolink->scaleBiologyToInternal( 1e-6 * model->mKillZoneRadius, BiologyLink::ScaleLength );
    for ( double zext = -model->mlobule_height/2 + mVoxelSize;
          zext < model->mlobule_height/2 - mVoxelSize;
          zext += mVoxelSize )
        for ( double yext = -radius; yext < radius; yext += mVoxelSize )
            for ( double xext = -radius; xext < radius; xext += mVoxelSize )
            {
                if ( radius*radius > xext*xext+yext*yext )
                {
                    if ( CSVoxel * voxel = getVoxel(xext,yext,zext) )
                        mpBBListCopy->add( voxel );
                }
            }

    mpBBListCopy->update();

    for ( size_t i=0; i < mDimensions[0]*mDimensions[1]*mDimensions[2]; ++i )
    {
        if ( CSVoxel *voxel = mpVoxels[i] )
            for ( auto contact: voxel->mIntersectingList )
                overlapTest( voxel, mpBBListCopy->element(contact) );
    }

    // put all pixels on the z-axis into the work queue for the space filling
    for ( double zext = -model->mlobule_height/2;
          zext < model->mlobule_height/2;
          zext += mVoxelSize )
    {
        CSVoxel * voxel = getVoxel(0., 0., zext);
        q.push( voxel );
    }

    const std::vector<std::array<double, 3>> searchDirections = {{{-1.5*mVoxelSize, 0., 0.}},
                                                                 {{ 1.5*mVoxelSize, 0., 0.}},
                                                                 {{0.,  1.5*mVoxelSize, 0.}},
                                                                 {{0., -1.5*mVoxelSize, 0.}},
                                                                 {{0., 0., -1.5*mVoxelSize}},
                                                                 {{0., 0.,  1.5*mVoxelSize}}};
    double minz=0., maxz=0.;
    while ( !q.empty() )
    {
        CSVoxel * voxel = q.front();
        q.pop();
        Vector3f pos = voxel->position;
        Vector3f dir;

        if ( pos.z < minz )
            minz = pos.z;
        else if ( pos.z > maxz )
            maxz = pos.z;

        for ( auto direction: searchDirections )
        {
            dir.x = direction.data()[0];
			dir.y = direction.data()[1];
			dir.z = direction.data()[2];
            Vector3f next = pos + dir;

            if ( next.x*next.x+next.y*next.y >= radius*radius
                 || next.z < -model->mlobule_height/2 + 2*mVoxelSize
                 || next.z >  model->mlobule_height/2 - 2*mVoxelSize)
                 continue;

            CSVoxel * nextvoxel = getVoxel( next.x, next.y, next.z );

            if (!nextvoxel)
                continue;

            if ( !nextvoxel->visited /*hasAttribute("visited")*/ )
            {
                nextvoxel->visited = true; //setAttribute("visited");

                if ( nextvoxel->contains != ModelElement::TypeCellSpherical &&
                     nextvoxel->contains != ModelElement::TypeCellSphericalPolar )
                {
                    q.push(nextvoxel);
                    mvRegion.push_back(nextvoxel);
                    // put the voxels into the mpArena
                    model->mpArena->addObject(nextvoxel->GLObject());
                }
            }
        }
    }

    int numstacks = 1 + (maxz - minz)/mVoxelSize;
    const double areaUnitInMeters = model->biolink->getLengthInMicrometers(mVoxelSize)*model->biolink->getLengthInMicrometers(mVoxelSize)* 1e-12;
    // std::cout << "# voxels in region: " << mvRegion.size()
    //           << "   number of levels: " << numstacks
    //           << "   voxel area (edge length^2): "
    //           << areaUnitInMeters
    //           << " m" << std::endl;
    // std::cout << "# average lesion area: " << areaUnitInMeters * mvRegion.size() / numstacks << " m^2" <<std::endl;
//    mObservables.lesionArea = areaUnitInMeters * mvRegion.size() / numstacks;
    mLesionArea = areaUnitInMeters * mvRegion.size() / numstacks;
//    std::cout <<  "# finished area discretization" << std::endl;
}


void
LobuleLesionVoxelization::resetVoxels()
{
    for ( size_t i = 0; i < mDimensions[0]*mDimensions[1]*mDimensions[2]; ++i )
        if ( mpVoxels[i] )
            mpVoxels[i]->visited = false;
}
