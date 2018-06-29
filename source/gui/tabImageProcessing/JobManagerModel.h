////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  JobManagerModel.h                                             //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2012-11-18 00:19:11                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef JOBMANAGER_MODEL_H
#define JOBMANAGER_MODEL_H

#include <QAbstractItemModel>


class CSParameterContext;
class JobManager;


class JobManagerModel : public QAbstractItemModel
{
  Q_OBJECT

  typedef CSParameterContext Job;

 public:
  JobManagerModel( JobManager *root, QObject * parent=0 );

  int rowCount( const QModelIndex &parent = QModelIndex() ) const;
  int columnCount( const QModelIndex &parent = QModelIndex() ) const;
  Qt::ItemFlags flags( const QModelIndex &parent = QModelIndex() ) const;

  QVariant headerData( int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;
  QVariant data( const QModelIndex &index, int role = Qt::DisplayRole ) const;
    //bool setData( const QModelIndex &index, const QVariant &value, int role );

  QModelIndex index( int row, int col, const QModelIndex &parent ) const;
  QModelIndex parent( const QModelIndex & index ) const;

  bool removeRows( int row, int count, const QModelIndex & parent = QModelIndex() );

public slots:
  void dataAboutToChange() { beginResetModel(); };
  void dataChanged() { endResetModel(); };


 private:
  JobManager * mpJobManager;
  Job * mpRoot;

};

#endif // JOBMANAGER_MODEL_H
