///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSGLArena.h                                                          //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-06-04 23:17:12                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef CS_GL_ARENA_H
#define CS_GL_ARENA_H

#include <QMutex>

#include "CSGLObject.h"
#include "QCSGLDisplay.h"
#include <vector>
#include <deque>

#include "../model/BasicDatatypes/Color.h"

class QCSGLDisplay;

/*!
  \brief The container for a GLDisplay's objects

  This class is used to collect and draw all GLDisplay's GLObjects.
*/
class CSGLArena : public CSGLObject
{
 public:
  CSGLArena();
  ~CSGLArena();

  void draw();
  
  void addObject(CSGLObject *newObject);
  void removeObject(CSGLObject *obsoleteObject);

  void clear();

  void registerDisplay( QCSGLDisplay * disp )
  { 
    // todo: search if display is registered already;
    mDisplays.push_back(disp);
  };

 private:
  std::deque<CSGLObject *> mObjectList;
  std::vector<QCSGLDisplay *> mDisplays;
  QMutex mMutex;
};

#endif // CS_GL_ARENA_H
