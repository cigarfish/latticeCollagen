///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSParameterTreeItem.h                                                //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-03-19 12:40:12                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef CS_PARAMETER_TREE_ITEM_H
#define CS_PARAMETER_TREE_ITEM_H

#include <string>


class CSParameterTreeItem
{
 public:
  enum NodeType { Context, Parameter }
;
  CSParameterTreeItem( std::string name, NodeType type, CSParameterTreeItem *parent=NULL )
    : mName(name),
      mType(type),
      mpParent(parent),
      mVisibility(true)
  {};

  std::string name() const { return mName; };
  NodeType type() const { return mType; };
  CSParameterTreeItem * parent() const { return mpParent; };

  void setName( std::string name ) { mName = name; };
  void setParent( CSParameterTreeItem * parent ) { mpParent = parent; };

  bool isVisible() const { return mVisibility; };
  void setVisible( bool visible = true ) { mVisibility = visible; };

 protected:
  std::string mName;
  NodeType mType;
  CSParameterTreeItem * mpParent;
  bool mVisibility;
};


#endif //CS_PARAMETER_TREE_ITEM_H
