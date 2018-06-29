////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  TiSimDemoMainWindow.cpp                                       //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2015-02-22 22:47:26                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include <QtGui>
#ifndef __BUILD_MAC__
# include <qapplication.h>
#endif

#include "../QCSGLDisplay.h"

#include "TiSimDemoMainWindow.h"

#include "TiSimDemoCentralWidget.h"


TiSimDemoMainWindow::TiSimDemoMainWindow(QWidget *parent, Qt::WindowFlags flags)
    : QMainWindow(parent, flags),
      mpActiveDisplay(NULL)
{
    resize(1000, 600);

    setCentralWidget(new TiSimDemoCentralWidget(this));

    QAction * actionQuit = new QAction("&Quit",this);
    actionQuit->setShortcuts(QKeySequence::Quit);
    connect(actionQuit, SIGNAL(triggered()), qApp, SLOT(quit()));
}


TiSimDemoMainWindow::~TiSimDemoMainWindow()
{}


/*!
  \brief Set the active display
  \param display new active display
*/
void
TiSimDemoMainWindow::setDisplay(QCSGLDisplay *display)
{
    if (!display)
        return;

    // first, allow only one display
    if (mpActiveDisplay != NULL)
        delete mpActiveDisplay;

    mpActiveDisplay = display;
}


/*!
  \brief Get the active display
 */
QCSGLDisplay *
TiSimDemoMainWindow::getDisplay()
{
    return mpActiveDisplay;
}


void
TiSimDemoMainWindow::closeEvent( QCloseEvent * event )
{
    QMessageBox message;

#ifndef CS_TI_QUANT_ONLY
    message.setText("You are about to quit TiSim Demo.");
#else
    message.setText("You are about to quit TiQuant.");
#endif
    message.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel );

    int close = message.exec();

    if ( close == QMessageBox::Cancel )
    {
        event->ignore();
        return;
    }

    qApp->quit();
}
