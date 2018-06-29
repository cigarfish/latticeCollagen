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

#ifndef VASCULARIZATION_H
#define VASCULARIZATION_H
#include "ui_Vascularization.h"
 
#include <vector>

class Vascularization;
class QCSSimulationThread;
class QCSGLDisplay;

class Model_Vascularization;
class Cube;

class Vascularization : public QWidget, private Ui::Vascularization
{
    Q_OBJECT
public:
    Vascularization(QWidget *parent=0);
    ~Vascularization();
private slots:
	void init();
	void createNewDisplay();
	void load();
	void play();
private:
	Model_Vascularization *vascularization;
    QCSSimulationThread *mpSimulationThread;
	bool mReset;
    std::vector<QCSGLDisplay *> mDisplays;

};
#endif
