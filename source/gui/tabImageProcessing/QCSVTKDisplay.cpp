///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSVTKDisplay.cpp                                                    //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-03-24 23:55:49                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "QCSVTKDisplay.h"

#include <QtGui>

#include "QVTKWidget.h"


QCSVTKDisplay::QCSVTKDisplay( QWidget * controlWidget, QWidget * parent )
  : QWidget(parent),
    mpControlWidget(controlWidget)
{
  mpLayout = new QHBoxLayout(this);

  setLayout(mpLayout);

  mpSplitter = new QSplitter( Qt::Vertical, this );

  mpLayout->addWidget( mpSplitter );

  mpQVTKWidget = new QVTKWidget(this);
  mpQVTKWidget->resize( 800, 600 );
  mpSplitter->setSizePolicy( QSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding) );

  mpSplitter->insertWidget( 0, mpQVTKWidget );

  if ( mpControlWidget )
    mpLayout->addWidget(mpControlWidget);

  mpConsole = new QTextEdit(this);
  mpSplitter->insertWidget( 1, mpConsole ); // TODO: span over both columns, if there's a controlWidget
  mpSplitter->setCollapsible( 1, true );
  mpConsole->hide();

  mpMenuBar = new QMenuBar(this);

  // QMenu * viewMenu = mpMenuBar->addMenu("View");
  // QAction * viewTest = new QAction( tr("&Test"), this );
  // connect( viewTest, SIGNAL(triggered()), this, SLOT(testTriggered()) );
  // viewMenu->addAction(viewTest);

  connect( mpConsole, SIGNAL(textChanged()), mpConsole, SLOT(show()) );

  this->setFocusProxy(mpQVTKWidget);
}


QCSVTKDisplay::~QCSVTKDisplay()
{
    delete mpMenuBar;
    delete mpConsole;
    delete mpQVTKWidget;
    delete mpLayout;
}


void
QCSVTKDisplay::testTriggered()
{
//  mpConsole->append("Test menu item triggered");
}


void
QCSVTKDisplay::addControlWidget( QWidget * widget )
{
  if (mpControlWidget)
    delete mpControlWidget;

  mpControlWidget = widget;
  mpControlWidget->setSizePolicy( QSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred) );
  mpLayout->addWidget(widget);
  mpQVTKWidget->setMinimumSize(QSize(mpControlWidget->height(),mpControlWidget->height()));
  this->setFocusProxy(mpQVTKWidget);
}
