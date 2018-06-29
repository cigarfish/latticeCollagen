////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  main.cpp                                                      //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2015-02-22 22:30:12                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#define _CS_MAIN_

#include <qapplication.h>

#ifdef __BUILD_MAC__
# include <GLUT/glut.h>
#else
# include <GL/glut.h>
#endif

#include "TiSimDemoMainWindow.h"

#include "../../Core.h"


int main(int argc, char **argv)
{
    core = new Core();
    core->init();

    // Init and Start QT
    QApplication mainApp(argc, argv);

#if ! defined( __BUILD_MAC__ )
    glutInit(&argc,argv);
#endif

    QMainWindow * mainDisplay = new TiSimDemoMainWindow();

    mainDisplay->show();

    return mainApp.exec();
}
