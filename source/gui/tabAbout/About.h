///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  About.h                                                              //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-06-10                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ABOUT_H
#define ABOUT_H
#include "ui_About.h"
 
#include <vector>



class About : public QWidget, private Ui::About
{
    Q_OBJECT
public:
  About(QWidget *parent=0);
  ~About();
};
#endif
