///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSGLArena.cpp                                                        //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-06-04 23:21:08                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "CSGLArena.h"


#include "CSGLArena.h"
#include "CSGLObject.h"

#include <algorithm>

/*!
  \brief Constructor
  \param parent The parent object
*/
CSGLArena::CSGLArena()
  : CSGLObject(new Vector3f(0.,0.,0.))
{
  mObjectList.clear();
  mDisplays.clear();
}


CSGLArena::~CSGLArena()
{}


/**
 * @brief Draw all GLObjects in the Arena
 *
 * This function will go through the list of objects and call the draw() routine
 *
 */
void
CSGLArena::draw()
{
  mMutex.lock();

  std::deque<CSGLObject *>::const_iterator listPtr = mObjectList.begin();
  std::deque<CSGLObject *>::const_iterator listEnd = mObjectList.end();
  
  while ( listPtr != listEnd )
  {
    if ( (*listPtr)->isTransparent() )
      glDepthMask(GL_FALSE);

    (*listPtr)->draw();
    ++listPtr;

    glDepthMask(GL_TRUE);
  }

  mMutex.unlock();
}


/*!
  \brief Add a GLObject to the arena
  \param newObject The GLObject to add

  Adds the newObject to the vector of objects.
*/ 
void
CSGLArena::addObject(CSGLObject *newObject)
{
  mMutex.lock();

  if ( newObject->isTransparent() )
    mObjectList.push_back(newObject);
  else
    mObjectList.push_front(newObject);

  mMutex.unlock();
}


/*!
  \brief Remove a GLObject to the arena
  \param obsoleteObject The GLObject to be removed

  Erases the obsoleteObject from the vector of objects.
*/ 
void
CSGLArena::removeObject(CSGLObject *obsoleteObject)
{
  mMutex.lock();

  std::deque<CSGLObject *>::iterator oIt;

  oIt = std::find( mObjectList.begin(), mObjectList.end(), obsoleteObject );

  // did we find it?
  if ( oIt != mObjectList.end() )
    mObjectList.erase( oIt );

  mMutex.unlock();
}


/*!
  \brief Remove all GLObjects from arena
*/ 
void
CSGLArena::clear()
{
    mMutex.lock();
	mObjectList.clear();
    mMutex.unlock();
}

