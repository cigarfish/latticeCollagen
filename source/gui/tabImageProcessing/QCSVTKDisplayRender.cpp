///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSVTKDisplayRender.cpp                                              //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes@neitsch.de>                               //
//    Created:  2013-08-28 15:23:24                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#if defined( CS_BUILD_IMAGEPROCESSING )

  #include "QCSVTKDisplayRender.h"
  #include "QCSVTKDisplay.h"
  #include <QShortcut>

  //#include "vtkPOVExporter.h"
 #include "../tabToolsWithItkVtk/vtkPOVExporterWireFrame.h"
#endif


QCSVTKDisplayRender::QCSVTKDisplayRender( QWidget * controlWidget, QWidget * parent )
  : QCSVTKDisplay(controlWidget,parent)
{
  QShortcut * shortcut = new QShortcut(QKeySequence(QObject::tr("p")), this );
  shortcut = new QShortcut(QKeySequence(Qt::SHIFT + Qt::Key_P), this );
  connect( shortcut, SIGNAL(activated()), this, SLOT(toPovray()));

 
}


void QCSVTKDisplayRender::toPovray(){
  
  vtkPOVExporterWireFrame *povexp = vtkPOVExporterWireFrame::New();
  
  povexp->SetRenderWindow(this->renWin);
  povexp->SetFileName("TestPOVExporter.pov");
  cout << "Writing file TestPOVExporter.pov..." << endl;
  
  povexp->Write();
  cout << "Done writing file TestPOVExporter.pov..." << endl;
  
  povexp->Delete();
  
}




