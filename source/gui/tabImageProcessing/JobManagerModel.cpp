////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  JobManagerModel.cpp                                           //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2012-11-18 00:25:36                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "JobManagerModel.h"
#include "../../tools/JobManager.h"
#include "../../tools/parameters/CSParameterContext.h"
// t1m-debug
#include <iostream>
// !t1m-debug

JobManagerModel::JobManagerModel(JobManager *manager, QObject * parent)
  : QAbstractItemModel(parent),
    mpJobManager( manager )
{
    connect( manager, SIGNAL( aboutToChange() ), this, SLOT( dataAboutToChange() ) );
    connect( manager, SIGNAL( changed() ), this, SLOT( dataChanged() ) );
}


int
JobManagerModel::rowCount( const QModelIndex &parent ) const
{
    return mpJobManager->Queue()->getChildren().size();
}


int
JobManagerModel::columnCount( const QModelIndex & /*parent*/ ) const
{
  return 3;
}


Qt::ItemFlags
JobManagerModel::flags( const QModelIndex &index ) const
{
  if ( !index.isValid() )
    return 0;

  return  Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}


QVariant
JobManagerModel::headerData( int section, Qt::Orientation orientation, int role ) const
{
  if ( role == Qt::DisplayRole )
    if ( orientation == Qt::Horizontal )
      {
        switch (section) {
        case 0:
          return QString("ID");
        case 1:
          return QString("Type");
        case 2:
          return QString("Parameters");
        }
      }

  return QVariant();
}


QVariant
JobManagerModel::data( const QModelIndex &index, int role ) const
{
  if ( !index.isValid() )
    return QVariant();

  if ( index.row() > rowCount(parent(index)) )
    return QVariant();

  if ( (role != Qt::DisplayRole) && (role != Qt::CheckStateRole) )
    return QVariant();

  int column = index.column();

  Job *indexJob = (Job *)index.internalPointer();

  if ( role == Qt::DisplayRole )
    {
        switch ( column )
        {
        case 0:
            return QVariant( QString( indexJob->name().c_str() ) );
        case 1: 
        {
            CSParameter * typeParm = indexJob->findParameter("type");
            if ( ! typeParm )
                return QVariant();

            int * type = (int *) typeParm->dataPointer();
            return QVariant( JobManager::JobTypeString[*type].c_str() );
        }
        case 2:
        {
            QStringList parmsList;
            const std::vector<CSParameter *> & parms = indexJob->getParameters();
            std::vector<CSParameter *>::const_iterator parmIt = parms.begin();
            while ( parmIt != parms.end() )
            {
                QString name = QString( (*parmIt)->name().c_str() );
                if ( name == "type" )
                    break;
                QString dataString = QString( (*parmIt)->dataString().c_str() );
                parmsList << (name + "=" + dataString);
                ++parmIt;
            }

            //parmsList << "Test" << "Test" << "Test";
            const std::vector<CSParameterContext *> & ctxts = indexJob->getChildren();
            std::vector<CSParameterContext *>::const_iterator ctxtIt = ctxts.begin();
            while ( ctxtIt != ctxts.end() )
            {
                //parmsList << QString( (*ctxtIt)->name().c_str() );
                const std::vector<CSParameter *> parms = (*ctxtIt)->getParameters();
                parmIt = parms.begin();
                while ( parmIt != parms.end() )
                {
                    QString string = QString( (*parmIt)->name().c_str() )
                        + "=" + QString( (*parmIt)->dataString().c_str() );
                    parmsList << string;
                    ++parmIt;
                }

                ++ctxtIt;
            }
            // t1m-debug
            //QString strng = parmsList.join("\n");
            //std::cout << strng.toStdString() << std::endl;
            // !t1m-debug

            return parmsList.join("\n");
        }
        default:
            return QVariant();
        }

    }

  return QVariant();
}


QModelIndex
JobManagerModel::index( int row, int col, const QModelIndex &parent ) const
{
  if ( row < 0 || col < 0 || col > 3 )
    return QModelIndex();

  if ( !hasIndex(row, col, parent) )
     return QModelIndex();

  Job * job = mpJobManager->Queue()->getChildren()[row];

  return createIndex( row, col, job );
}


QModelIndex
JobManagerModel::parent( const QModelIndex & index ) const
{
    return QModelIndex();
}


bool
JobManagerModel::removeRows( int row, int count, const QModelIndex & parent )
{
    std::vector<Job *> jobs = mpJobManager->Queue()->getChildren();
    std::vector<Job *> deleteJobs;
    for ( int line=row, num=0; num < count; ++line, ++num )
    {
        deleteJobs.push_back(jobs[line]);
    }
    std::vector<Job *>::iterator delIt;
    for (delIt=deleteJobs.begin(); delIt != deleteJobs.end(); ++delIt)
        mpJobManager->DeleteJob( (CSParameterContextTemporary *)*delIt );
    return true;
}
