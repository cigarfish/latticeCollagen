///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSParameterModel.cpp                                                 //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-03-16 12:33:15                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include <cassert>

#include "QCSParameterModel.h"
#include "../tools/parameters/CSParameterContext.h"
#include "../tools/parameters/CSParameter.h"

#include <iostream>

QCSParameterModel::QCSParameterModel(CSParameterContext *root, QObject * parent)
  : QAbstractItemModel(parent),
    mpRoot(root)
{
    connect( this, SIGNAL( choiceChanged(const QModelIndex &) ), this, SLOT( changeActivatedSubContexts(const QModelIndex &) ) );
}


int
QCSParameterModel::rowCount( const QModelIndex &parent ) const
{
  if ( parent.column() > 0 )
    return 0;

  CSParameterContext * parentContext;

  if ( !parent.isValid() )
    parentContext = mpRoot;
  else
    parentContext = static_cast<CSParameterContext *>(parent.internalPointer());

  if (parentContext->type() == CSParameterTreeItem::Parameter)
    return 0;

  int numRows = parentContext->getChildren().size() + parentContext->getParameters().size();
  return numRows;
}


int
QCSParameterModel::columnCount( const QModelIndex & /*parent*/ ) const
{
  return 3;
}


Qt::ItemFlags
QCSParameterModel::flags( const QModelIndex &index ) const
{
  if ( !index.isValid() )
    return 0;

  Qt::ItemFlags itemFlags = Qt::ItemIsSelectable;

  CSParameterTreeItem *object = (CSParameterTreeItem *) index.internalPointer();

  if ( object->isVisible() )
      itemFlags |= Qt::ItemIsEnabled;

  if ( object->type() == CSParameterTreeItem::Context )
      return itemFlags;

  if ( index.column() == 1 )
    {
      itemFlags |= Qt::ItemIsEditable;

      // the internalPointer must be a CSParameter, as this returned already for a
      // CSParameterContext
      if ( ((CSParameter *) object)->dataType() == CSParameter::Bool )
        itemFlags |= Qt::ItemIsUserCheckable;
    }

  return itemFlags;
}


QVariant
QCSParameterModel::headerData( int section, Qt::Orientation orientation, int role ) const
{
  if ( role == Qt::DisplayRole )
    if ( orientation == Qt::Horizontal )
      {
        switch (section) {
        case 0:
          return QString("Parameter");
        case 1:
          return QString("Value");
        case 2:
          return QString("Unit");
        }
      }

  return QVariant();
}


QVariant
QCSParameterModel::data( const QModelIndex &index, int role ) const
{
  if ( !index.isValid() )
    return QVariant();

  if ( index.row() > rowCount(parent(index)) )
    return QVariant();

  if ( (role != Qt::DisplayRole) && (role != Qt::CheckStateRole) )
    return QVariant();

  int column = index.column();

  CSParameterTreeItem *indexItem = (CSParameterTreeItem *)index.internalPointer();

  if ( role == Qt::DisplayRole )
    {
      if (indexItem->type() == CSParameterTreeItem::Parameter)
        {
          CSParameter * dummyParameter = static_cast<CSParameter *>( indexItem );

          if (!dummyParameter) // should not happen -- raise exception?!
            return QVariant();

          if ( column == 0 )
            return QVariant( dummyParameter->name().c_str() );

          // Booleans get handled via the Qt::CheckStateRole and delegate
          if ( dummyParameter->dataType() == CSParameter::Bool )
            {
              // bool yesno = * (bool *) dummyParameter->dataPointer();
              // return QVariant( yesno ? Qt::Checked : Qt::Unchecked );
              return QVariant();
            }

          if ( column == 1 )
            return QVariant( dummyParameter->dataString().c_str() );
          if ( column == 2 )
            return QVariant( dummyParameter->unit().c_str() );
        }
      else
        {
          CSParameterContext * dummyContext = static_cast<CSParameterContext *>( indexItem );

          if (!dummyContext) // should not happen -- raise exception?!
            return QVariant();

          if ( column == 0 )
            return QVariant( dummyContext->name().c_str() );
        }
    }

  if ( role == Qt::CheckStateRole && column == 1 )
    if ( indexItem->type() == CSParameterTreeItem::Parameter )
      {
        CSParameter * parm = static_cast<CSParameter *>(indexItem);
        if ( parm->dataType() == CSParameter::Bool )
          {
            bool yesno = * (bool *) parm->dataPointer();
            return QVariant( yesno ? Qt::Checked : Qt::Unchecked );
          }
      }

  return QVariant();
}


bool
QCSParameterModel::setData( const QModelIndex &index, QVariant &value, int role )
{
  if ( !index.isValid() )
    return false;

  if ( index.column() != 1 )
    return false;

  CSParameterTreeItem * ptr =
    static_cast<CSParameterTreeItem *>(index.internalPointer());

  if ( ptr->type() != CSParameterTreeItem::Parameter )
    return false;

  CSParameter * parm = static_cast<CSParameter *>(ptr);

  void * data = parm->dataPointer();

  bool changed = false;
  bool okFlag = false;

  switch ( parm->dataType() )
    {
    case  CSParameter::Bool:
        if ( role == Qt::CheckStateRole )
        {
            parm->setData( value.toBool() );
            okFlag = true;
            changed = true;
        }
        break;

    case CSParameter::Int:
    {
        int val = value.toInt( &okFlag );
        if ( val != * static_cast<int *>(data) )
        {
            * static_cast<int *>(data) = val;
            changed = true;
        }
    }
    break;

    case CSParameter::Long:
    {
        long val  = (long) value.toInt( &okFlag );
        if ( val != * static_cast<long *>(data) )
        {
            * static_cast<long *>(data) = val;
            changed = true;
        }
    }
    break;

    case CSParameter::Float:
    {
        float val = value.toFloat( &okFlag );
        if ( val != * static_cast<float *>(data) )
        {
            * static_cast<float *>(data) = val;
            changed = true;
        }
    }
    break;

    case CSParameter::Double:
    {
        double val = value.toDouble( &okFlag );
        if ( val != * static_cast<double *>(data) )
        {
            *static_cast<double *>(data) = val;
            changed = true;
        }
    }
    break;

    case CSParameter::String:
    case CSParameter::DirName:
    case CSParameter::FileName:
    {
        std::string val = value.toString().toStdString();
        if ( val != * static_cast<std::string *>(data) )
        {
            static_cast<std::string *>(data)->clear();
            static_cast<std::string *>(data)->append( val );
            changed = true;
        }
        okFlag = true;
    }
    break;

    case CSParameter::RangeInt:
        break;

    case CSParameter::RangeLong:
        break;

    case CSParameter::RangeFloat:
        break;

    case CSParameter::RangeDouble:
        break;

    case CSParameter::Choice:
    {
        unsigned int val = value.toUInt(&okFlag);
        if ( val != static_cast<CSParameterChoice *>(data)->currentIndex() )
        {
            static_cast<CSParameterChoice *>(data)->setCurrentIndex( val );
            if ( static_cast<CSParameterChoice *>(data)->controlsSubContexts() )
                emit choiceChanged(index);
            changed = true;
        }
    }
    break;
    }

  if ( changed )
      emit QAbstractItemModel::dataChanged( index, index );

  return okFlag;
}



QModelIndex
QCSParameterModel::index( int row, int col, const QModelIndex &parent ) const
{
  if ( row < 0 || col < 0 )
    return QModelIndex();

  // if ( !hasIndex(row, col, parent) )
  //   return QModelIndex();

  CSParameterContext * parentContext;

  if ( !parent.isValid() )
    parentContext = mpRoot;
  else
    parentContext = static_cast<CSParameterContext *>(parent.internalPointer());

  int numParms = parentContext->getParameters().size();

  if ( row < numParms )
    {
      CSParameter * parm = parentContext->getParameters()[row];
      QModelIndex index =  createIndex( row, col, parm );
      return index;
    }
  else
    {
      int numChildren = parentContext->getChildren().size();

      if ( row > numChildren + numParms )
        return QModelIndex();

      CSParameterContext * context = parentContext->getChildren()[row-numParms];

      return createIndex( row, col, context );
    }   
}


QModelIndex
QCSParameterModel::parent( const QModelIndex & index ) const
{
  if ( !index.isValid() )
    return QModelIndex();

  CSParameterTreeItem * thisItem = static_cast<CSParameterTreeItem *>( index.internalPointer() );
  
  // if ( thisItem == mpRoot )
  //   return QModelIndex();

  CSParameterContext * parentContext = static_cast<CSParameterContext *>( thisItem->parent() );

  if ( parentContext == mpRoot )
    return QModelIndex();

  assert(parentContext);

  CSParameterContext * grandParent = parentContext->parent();

  std::vector<CSParameterContext *> uncles = grandParent->getChildren();
  std::vector<CSParameterContext *>::iterator it = std::find(uncles.begin(),uncles.end(),parentContext);
  int offset = grandParent->getParameters().size();
  int row = std::distance( uncles.begin(), it );
  int parentsRow = offset + row;

  return createIndex( parentsRow, 0, parentContext );
}


void
QCSParameterModel::changeActivatedSubContexts(const QModelIndex &index)
{
    CSParameter * item = static_cast<CSParameter *>( index.internalPointer() );

    // Assuming the index is at the CSParameter::Choice entry that triggered this slot:
    CSParameterChoice * choice = static_cast<CSParameterChoice *>( item->dataPointer() );

    unsigned int currentIndex = choice->currentIndex();
    unsigned int numChoices = choice->choices().size();

    // the subContexts:
    std::vector<CSParameterContext *> children = item->parent()->getChildren();
    unsigned int numChildren = children.size();

    if ( numChildren < numChoices )
    {
        std::cout << "Warning:  QCSParameterModel::changeActivatedSubContexts():\n"
                  << "          You have fewer subcontexts than choices for \n"
                  << "           Parameter " << item->name() << std::endl
                  << "           Context " << item->parent()->name() << std::endl;
    }

    bool visibility;
    CSParameterContext *child;
    unsigned int i = 0;
    while ( i < numChildren )
    {
        child = children[i];
        visibility = false;
        if ( i == currentIndex )
            visibility =true;

        child->setVisible( visibility );
        ++i;
    }
}
