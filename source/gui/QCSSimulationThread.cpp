///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSSimulationThread.cpp                                              //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-06-26 17:34:39                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "QCSSimulationThread.h"

#include "../model/Model/CSModel.h"
#include "../../Core.h"

#include <limits>


QCSSimulationThread::QCSSimulationThread()
    : QThread(),
      mStepCounter(0),
      mMaxSteps(std::numeric_limits<unsigned long int>::max()),
      mUpdate(1),
      mStopped(false),
      mPaused(false),
      mContinued(false)
{}


QCSSimulationThread::~QCSSimulationThread()
{}


void
QCSSimulationThread::run()
{
  mStopped = false;

  if ( !mContinued )
      {
          mContinued = true;
          mStepCounter=0;
      }

  while (mpModel->enableSimulation)  
  {
      if ( mStepCounter > mMaxSteps )
          mStopped = true;

      if (mPaused)
          {
              emit paused();
              mWaitCondition.wait(&mWaitMutex);
          }

      if ( mStopped )
          return;

      mpModel->SimulateInThread();

      if ( ! (mStepCounter % mUpdate) )
          emit updateStep(mStepCounter);

      ++mStepCounter;
    }

  emit finishedNormally();
}
