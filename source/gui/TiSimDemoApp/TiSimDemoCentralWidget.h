////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  TiSimDemoCentralWidget.h                                      //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2015-02-22 22:54:49                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef TISIM_DEMO_CENTRALWIDGET_H
#define TISIM_DEMO_CENTRALWIDGET_H

#include <QTabWidget>


class TiSimDemoCentralWidget : public QTabWidget
{
public:
    TiSimDemoCentralWidget(QWidget *parent);

    ~TiSimDemoCentralWidget();
};


#endif // TISIM_DEMO_CENTRALWIDGET_H
