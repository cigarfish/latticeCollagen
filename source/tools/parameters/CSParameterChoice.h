///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSParameterChoice.h                                                  //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-05-24 15:46:44                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef CS_PARAMETER_CHOICE_H
#define CS_PARAMETER_CHOICE_H

#include <string>
#include <vector>


class CSParameterChoice
{
 public:
  CSParameterChoice( const char * list[], int numElements, int index=0 );
  CSParameterChoice( const std::vector<std::string> & list, int index=0 );
  CSParameterChoice( const std::string list[], int index=0 );

  ~CSParameterChoice() { mList.clear(); };

  void setCurrentIndex( int index ) { mCurrentIndex = index; };
  void setChoices( const std::vector<std::string> & list, int index=0, bool overwrite=true );
  unsigned int currentIndex() const { return mCurrentIndex; };

  std::string currentString() const { return mList[mCurrentIndex]; };
  std::vector<std::string> choices() { return mList; };

  // flag/methods for influencing activation of sub-contexts
  // mMaster == true means that the activation of one of the first mList.size()
  // contexts of the parent context is influenced by mCurrentIndex.
  // cf QCSParameterModel.cpp
  void setControlSubContexts( bool master=true ) { mMaster = master; };
  bool controlsSubContexts() const {return mMaster;};

 private:
  std::vector<std::string> mList;
  int mCurrentIndex;
  bool mMaster;
};

#endif // CS_PARAMETER_CHOICE_H
