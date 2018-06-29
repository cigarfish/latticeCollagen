///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSVTKDisplay.h                                                      //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-03-24 23:46:24                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef Q_CS_VTK_DISPLAY_H
#define Q_CS_VTK_DISPLAY_H

#include <QWidget>
#include <QTextEdit>
#include "QVTKWidget.h"

class QMenuBar;
class QHBoxLayout;
class QSplitter;
class vtkRenderWindow;
class vtkRenderWindowInteractor;


class QCSVTKDisplay : public QWidget
{
  Q_OBJECT

 public:
  QCSVTKDisplay(QWidget * controlWidget=0, QWidget * parent=0);
  ~QCSVTKDisplay();

  void addControlWidget( QWidget * cWidget );
  QWidget * getControlWidget() const { return mpControlWidget; };

  QVTKWidget * getVTKWidget() const { return mpQVTKWidget; };

  void showConsole() const { mpConsole->show(); };
  void hideConsole() const { mpConsole->hide(); };

  vtkRenderWindow * getRenderWindow() const { return mpQVTKWidget->GetRenderWindow(); };
  void setRenderWindow( vtkRenderWindow * renwin ) { mpQVTKWidget->SetRenderWindow(renwin); };

  vtkRenderWindowInteractor * getInteractor() const { return mpQVTKWidget->GetInteractor(); };

 private slots:
  void testTriggered();

 protected:
  QMenuBar *mpMenuBar;
  QSplitter * mpSplitter;
  QHBoxLayout * mpLayout;
  QWidget * mpControlWidget;
  QTextEdit * mpConsole;
  QVTKWidget * mpQVTKWidget;
};

#endif //Q_CS_VTK_DISPLAY_H
