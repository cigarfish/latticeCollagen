////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  TiSimDemoMainWindow.h                                         //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2015-02-22 22:38:28                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef TISIM_DEMO_MAINWINDOW_H
#define TISIM_DEMO_MAINWINDOW_H

#include <QMainWindow>

class QCSGLDisplay;


class TiSimDemoMainWindow : public QMainWindow
{
    Q_OBJECT

    friend class TiSimDemoCentralWidget;

public:
    TiSimDemoMainWindow(QWidget * parent =0, Qt::WindowFlags flags=0);

    ~TiSimDemoMainWindow();

    QCSGLDisplay * getDisplay();

protected:
    void setDisplay( QCSGLDisplay * );
    void closeEvent( QCloseEvent * );

    QCSGLDisplay * mpActiveDisplay;
    QList<QCSGLDisplay *> mDisplays;

    QString mDataDirName;
    bool mSaved;
};

#endif // TISIM_DEMOSMAINWINDOW_H
