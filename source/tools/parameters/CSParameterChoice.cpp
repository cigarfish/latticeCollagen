///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSParameterChoice.cpp                                                //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-05-24 15:59:49                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "CSParameterChoice.h"

#include <iostream>


CSParameterChoice::CSParameterChoice( const char* list[], int numElements, int index )
    : mCurrentIndex(index),
      mMaster(false)
{
  for ( unsigned int i=0; i<numElements; ++i )
    mList.push_back(list[i]);
}


CSParameterChoice::CSParameterChoice( const std::vector<std::string> & list, int index )
    : mCurrentIndex(index),
      mMaster(false)
{
  for ( unsigned int i=0; i<list.size(); ++i )
    mList.push_back( list[i] );
}


//! Read in an array of strings which is limited with an empty string.
CSParameterChoice::CSParameterChoice( const std::string list[], int index )
    : mCurrentIndex(index),
      mMaster(false)
{
  for ( unsigned int i=0; !list[i].empty(); ++i )
    mList.push_back( list[i] );
}


void CSParameterChoice::setChoices( const std::vector<std::string> & list, int index, bool overwrite )
{
  mCurrentIndex = index;

  if ( overwrite )
    mList.clear();

  for ( unsigned int i=0; i<list.size(); ++i )
    mList.push_back( list[i] );
}
