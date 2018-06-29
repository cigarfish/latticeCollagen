////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  ImageProcessingWorker.h                                       //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2014-05-30 18:22:50                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_IMAGEPROCESSING_WORKER_H
#define CS_IMAGEPROCESSING_WORKER_H

#include <QObject>

#include "ImageProcessing.h"

class CSParameterContext;


class ImageProcessingWorker : public QObject
{
    Q_OBJECT

    friend class ImageProcessing;

public:
    ImageProcessingWorker( CSParameterContext * rootContext, ImageProcessing *parent )
        : QObject(), mpCaller(parent), mpParameterContextRoot(rootContext), mpPipelineContextCopy(NULL) {}

    void setPipeline( ImageProcessing::Pipeline pipeline ) { mPipelineID = pipeline; }

public slots:
    void execute();

signals:
    void finished();
    void aborted(QString reason);
    void displayResult(void * pipeline, int dm, bool displayGraphAlreadyDisplayed);

private:

    ImageProcessing             * mpCaller;
    CSParameterContext          * mpParameterContextRoot;
    CSParameterContextTemporary * mpPipelineContextCopy;
    ImageProcessing::Pipeline     mPipelineID;
};

#endif // CS_IMAGEPROCESSING_WORKER_H
