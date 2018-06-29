///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSParameterModel.h                                                  //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-03-14 12:44:26                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef QCS_PARAMETER_MODEL_H
#define QCS_PARAMETER_MODEL_H

#include <QAbstractItemModel>

class CSParameterContext;


class QCSParameterModel : public QAbstractItemModel
{
  Q_OBJECT

 public:
  QCSParameterModel( CSParameterContext *root, QObject * parent=0 );

  int rowCount( const QModelIndex &parent = QModelIndex() ) const;
  int columnCount( const QModelIndex &parent = QModelIndex() ) const;
  Qt::ItemFlags flags( const QModelIndex &parent = QModelIndex() ) const;

  QVariant headerData( int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;
  QVariant data( const QModelIndex &index, int role = Qt::DisplayRole ) const;
  bool setData( const QModelIndex &index, QVariant &value, int role );

  QModelIndex index( int row, int col, const QModelIndex &parent ) const;
  QModelIndex parent( const QModelIndex & index ) const;

  CSParameterContext* rootContext() const { return mpRoot; };

 signals:
  void choiceChanged( const QModelIndex & );
  void dataChanged(int);


 private slots:
  void changeActivatedSubContexts( const QModelIndex & );


 private:
  CSParameterContext * mpRoot;
};

#endif // QCS_PARAMETER_MODEL_H
