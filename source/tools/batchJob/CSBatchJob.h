////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSBatchJob.h                                                  //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-01-23 15:20:34                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef CS_BATCH_JOB_H
#define CS_BATCH_JOB_H

#include <string>
#include <vector>

class CSModel;

#if defined( CS_BUILD_IMAGEPROCESSING )
class JobManager;
#endif


class CSBatchJob
{
 public:
    enum Mode { ModeSimulationThread=0, ModeJobManager };
    CSBatchJob( const std::string & inputFileName );

    int run();

    bool isValid() {return mValidFlag;};

 protected:
    std::vector<CSModel *> mModels;
#if defined( CS_BUILD_IMAGEPROCESSING )
    JobManager * mpJobManager;
#endif
    bool mValidFlag;
};

#endif // CS_BATCH_JOB_H
