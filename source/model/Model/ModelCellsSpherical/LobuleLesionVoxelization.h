////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  LobuleLesionVoxelization.h                                    //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2017-04-05 15:58:24                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_LOBULE_LESION_VOXELIZATION_H
#define CS_LOBULE_LESION_VOXELIZATION_H

#include "../../../tools/model/CSSpaceDiscretization.h"

#include "ModelCellsSpherical.h"


class LobuleLesionVoxelization : public CSSpaceDiscretization
{
public:
    LobuleLesionVoxelization( ModelCellsSpherical * model, double voxelSize=1. )
        : CSSpaceDiscretization( static_cast<CSModel *>(model), voxelSize ),
          mLesionArea(0.)
    {};

    ~LobuleLesionVoxelization() { /* ToDo */ };

    LobuleLesionVoxelization( const LobuleLesionVoxelization & other )
        : CSSpaceDiscretization(other)
    { /* ToDo */ };

    LobuleLesionVoxelization( const LobuleLesionVoxelization && other )
        : CSSpaceDiscretization(std::move(other))
    { /* ToDo */ };

    LobuleLesionVoxelization & operator=( const CSSpaceDiscretization & ) { /* ToDo */ };

    void exec();

    double lesionArea() const { return mLesionArea; };

protected:

    void resetVoxels();
    std::vector<CSVoxel *> mvRegion;
    double mLesionArea;
};

#endif // CS_LOBULE_LESION_VOXELIZATION_H
