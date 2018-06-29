////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  TiSimDemoCentralWidget.cpp                                    //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2015-02-22 22:56:49                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "TiSimDemoCentralWidget.h"

#include "../QCS3DTabWidget.h"
#include "monolayerDemo.h"


TiSimDemoCentralWidget::TiSimDemoCentralWidget(QWidget *parent)
    : QTabWidget(parent)
{
    QSizePolicy spolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    setSizePolicy(spolicy);

    QWidget * pageView3D = new QCS3DTabWidget(this);

    QWidget * cellSphericalTab = new monolayerDemo(this);

    addTab( pageView3D, "Visualization" );
    addTab( cellSphericalTab, "Spherical Cell Model" );
}


TiSimDemoCentralWidget::~TiSimDemoCentralWidget()
{}
