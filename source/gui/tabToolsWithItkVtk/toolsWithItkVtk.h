///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Vascularization.h                                                    //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes@neitsch.de>                               //
//    Created:  2012-12-20                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef TOOLSWITHITKVTK_H
#define TOOLSWITHITKVTK_H
#include "ui_toolsWithItkVtk.h"
 
#include <vector>

class toolsWithItkVtk;
class QCSSimulationThread;
class QCSGLDisplay;




class toolsWithItkVtk : public QWidget, private Ui::toolsWithItkVtk
{
    Q_OBJECT
public:
  toolsWithItkVtk(QWidget *parent=0);
  ~toolsWithItkVtk();
private slots:
  void load();
private:
  QCSSimulationThread *mpSimulationThread;
  bool mReset;
  std::vector<QCSGLDisplay *> mDisplays;
  
};
#endif
