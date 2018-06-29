///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSVTKDisplayRender.h                                                //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes@neitsch.de>                               //
//    Created:  2013-08-28 15:23:24                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include <QWidget>
#include <QTextEdit>


#if defined( CS_BUILD_IMAGEPROCESSING )

  #include "QVTKWidget.h"
  #include "QCSVTKDisplay.h"
  
  #include <vtkRenderWindow.h>

#endif

#ifndef Q_CS_VTK_DISPLAY_RENDER_H
#define Q_CS_VTK_DISPLAY_RENDER_H

class QCSVTKDisplayRender : public QCSVTKDisplay
{
   Q_OBJECT
public:
  QCSVTKDisplayRender(QWidget * controlWidget=0, QWidget * parent=0);
 // ~QCSVTKDisplayRender();
private slots:
  void toPovray();

public:
  vtkRenderWindow * renWin;

};

#endif //Q_CS_VTK_DISPLAY_RENDER_H
