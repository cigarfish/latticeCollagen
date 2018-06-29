///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSSimulationThread.h                                                //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-06-26 15:59:06                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef Q_CS_SIMULATION_THREAD_H
#define Q_CS_SIMULATION_THREAD_H

#include <QThread>
#include <QMutex>
#include <QWaitCondition>

class CSModel;


class QCSSimulationThread : public QThread
{
   Q_OBJECT

 public:
   QCSSimulationThread();
   ~QCSSimulationThread();

   void setModel( CSModel * model )
   { mpModel = model; };

   void setMaxSteps( int steps )
   { mMaxSteps = steps; };

   void setUpdateInterval( int interval )
   { mUpdate = interval; }

   bool isPaused() const
   { return mPaused; };

   int currentStep() const
   { return mStepCounter; };

 protected:
   void run();


 public slots:
   void stop() { mStopped = true; resume(); };
   void pause(bool pause=true)
   {
     if ( mPaused != pause )
       {
         if ( ! mPaused )
             {
                 mPaused = true;
                 mWaitMutex.lock();
             }
         else
           resume();
       }
   };

   void resume()
   {
       mPaused = false;
       // should not hurt, if not locked
       mWaitMutex.unlock();
       mWaitCondition.wakeAll();
   };

   void reset()
   {
       // stop the thread
       mStopped = true;
       resume();

       mContinued = false;
       mStepCounter = 0;
   };

 signals:
   void updateStep(int);
   void paused();
   void finishedNormally();

 protected:
   QMutex   mWaitMutex;
   QWaitCondition mWaitCondition;

   CSModel  * mpModel;
   unsigned long mStepCounter;
   unsigned long    mMaxSteps;
   unsigned long      mUpdate;

   bool   mStopped;
   bool    mPaused;
   bool mContinued;
};

#endif //Q_CS_SIMULATION_THREAD_H
