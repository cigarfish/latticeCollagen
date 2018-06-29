///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  complexCells.h                                                       //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2013-10-24                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef COMPLEXCELLS_H

#include "ui_tabComplexCells.h"

#include <string>

class CSParameterContext;
class CSModel;
class Model3D;
class QDebugStream;
class QCSParameterModel;
class QCSSimulationThread;
class QCSGLDisplay;


class complexCells : public QWidget, private Ui::complexCells
{
  Q_OBJECT

public:

  complexCells(QWidget *parent=0);
  ~complexCells();

  void initFromModel( CSModel * );


private slots:

  // Parameterization
  void ResetParametersToDefaults();
  void PushParametersToModel();
  void ParametersChanged();

  // Simulation
  void startSimulationClicked();
  void pauseResumeSimulation();
  void resetSimulation();
  void buttonAbortSimulationButtonClicked();

  void threadUpdate();
  void doAfterThreadIsFinished();

  void createNewDisplay();

private:

  void ConnectModelToGUI();

  void GUIStatePristine();
  void GUIStateInitialized();
  void GUIStateRunning();
  void GUIStatePaused();
  void GUIStateFinished();

signals:

  void abortThread();
  void pauseResumeThread(bool);

private:

  std::string mName;

  Model3D *mpModel;

  bool mIsReset;
  bool mStateParametersChanged;

  std::vector<QCSGLDisplay *> mvDisplays;

  CSParameterContext *mpParameterContext;
  QCSParameterModel * mpParameterModel;

  bool mThreadPaused;

};


#endif //COMPLEXCELLS_H




