///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSParameter.h                                                        //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-03-14 17:07:16                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef CS_PARAMETER_H
#define CS_PARAMETER_H

#include <string>
#include <map>

#include "CSParameterTreeItem.h"
#include "CSParameterChoice.h"

class CSParameterContext;
class QXmlStreamWriter;
class QXmlStreamReader;


class CSParameter : public CSParameterTreeItem
{
 public:
  enum DataType { Bool,
                  Int,
                  Long,
                  Float,
                  Double,
                  String,
                  DirName,
                  FileName,
                  RangeInt,
                  RangeLong,
                  RangeFloat,
                  RangeDouble,
                  Choice };

  static const char * DataTypeString[];

  CSParameter( std::string name,
               CSParameter::DataType parmType,
               void * data = NULL,
               const std::string & unit = "",
               CSParameterContext * parent = NULL );

  //! Copy-constructor
  CSParameter( const CSParameter * otherParm );

  DataType dataType() const {return mDataType;};


  void setData(CSParameter::DataType parmType, void * data, bool force=false) 
  {
    if ( force || !mpData )
      {
        mDataType = parmType;
        mpData    =     data;
      }
  };


  void setData(bool data)
  {
    if ( mDataType == CSParameter::Bool )
      {
        if ( mpData )
          * (bool *) mpData = data;
        else
          mpData = (void *)(new bool(data));
      }
  };


  void setData(int data)
  {
    if ( mDataType == CSParameter::Choice )
      {
        if ( mpData )
          static_cast<CSParameterChoice *>(mpData)->setCurrentIndex( data );
      }
    else if ( mDataType == CSParameter::Int )
        setData( (long int) data );
  };


  void setData(long int data)
  {
    if ( mDataType == CSParameter::Long || mDataType == CSParameter::Int )
      {
        if ( mpData )
          * (long int*) mpData = data;
        else
          mpData = (void *)(new long(data));
      }
  };


  void setData( float data )
  {
    if ( mDataType == CSParameter::Float )
      {
        if ( mpData )
          * (float *) mpData = data;
        else
          mpData = (void *)(new float(data));
      }
  };


  void setData( double data )
  {
    if ( mDataType == CSParameter::Double )
      {
        if ( mpData )
          * (double *) mpData = data;
        else
          mpData = (void *)(new double(data));
      }
  };


  void setData( std::string data )
  {
    if ( mDataType == CSParameter::String ||
         mDataType == CSParameter::DirName ||
         mDataType == CSParameter::FileName )
      {
        if ( mpData )
          {
            ((std::string *) mpData)->clear();
            ((std::string *) mpData)->append( data );
          }
        else
          mpData = (void *)(new std::string(data));
      }
  };


  void setData( int start, int end )
  {
    if ( mDataType == CSParameter::RangeInt )
      {
        if ( mpData )
          {
            ((std::pair<int,int> *)mpData)->first = start;
            ((std::pair<int,int> *)mpData)->second = end;
          }
        else
          mpData = (void *)(new std::pair<int, int>(start,end));
      }
  };


  void setData( long start, long end )
  {
    if ( mDataType == CSParameter::RangeLong )
      {
        if ( mpData )
          {
            ((std::pair<long,long> *)mpData)->first = start;
            ((std::pair<long,long> *)mpData)->second = end;
          }
        else
          mpData = (void *)(new std::pair<long, long>(start,end));
      }
  };


  void setData( float start, float end )
  {
    if ( mDataType == CSParameter::RangeFloat )
      {
        if ( mpData )
          {
            ((std::pair<float, float> *)mpData)->first = start;
            ((std::pair<float, float> *)mpData)->second = end;
          }
        else
          mpData = (void *)(new std::pair<float, float>(start, end));
      }
  };


  void setData( double start, double end )
  {
    if ( mDataType == CSParameter::RangeDouble )
      {
        if ( mpData )
          {
            ((std::pair<double, double> *)mpData)->first = start;
            ((std::pair<double, double> *)mpData)->second = end;
          }
        else
          mpData = (void *)(new std::pair<double, double>(start, end));
      }
  };


  void setData( const char * list[], int numElements, int index )
  {
    if ( mDataType == CSParameter::Choice )
      {
        if ( ! mpData )
          mpData = (void *)( new CSParameterChoice( list, numElements, index ) );
      }
  };


  void setData( const std::vector<std::string> list, int index )
  {
    if ( mDataType == CSParameter::Choice )
      if ( !mpData )
        mpData = (void *)( new CSParameterChoice( list, index ) );
  };


  CSParameter & operator=( CSParameter & otherParm )
  {
    void *data = otherParm.dataPointer();
    CSParameter::DataType type = otherParm.dataType();

    switch (type)
      {
      case Bool:
          setData( * (bool *) data );
          break;
      case Int:
      case Long:
          setData( * (long *) data );
          break;
      case Float:
      case Double:
          setData( * (double *) data );
          break;
      case String:
      case DirName:
      case FileName:
          setData( * (std::string *) data );
          break;
      case RangeInt:
      {
          std::pair<int,int> * dummy = (std::pair<int,int> *) data;
          setData( dummy->first, dummy->second );
          break;
      }
      case RangeLong:
      {
          std::pair<long, long> * dummy = (std::pair<long,long> *) data;
          setData( dummy->first, dummy->second );
      }
      break;
      case RangeFloat:
      {
          std::pair<float, float> * dummy = (std::pair<float,float> *) data;
          setData( dummy->first, dummy->second );
      }
      break;
      case RangeDouble:
      {
          std::pair<double, double> * dummy = (std::pair<double,double> *) data;
          setData( dummy->first, dummy->second );
      }
      break;
      case Choice:
      {
          CSParameterChoice * dummy = (CSParameterChoice *) data;
          setData( dummy->choices(), dummy->currentIndex() );
      }
      break;

      }

    return *this;
  }


  double value()
  {
      switch ( mDataType )
      {
      case CSParameter::Bool:
          return (double) ( *(bool *) mpData);
      case CSParameter::Double:
          return (double) ( *(double *) mpData);
      case CSParameter::Float:
          return (double) ( *(float *) mpData);
      case CSParameter::Int:
          return (double) ( *(int *) mpData);
      case CSParameter::Long:
          return (double) ( *(long *) mpData);
      case Choice:
          return (double) ( (CSParameterChoice *) mpData)->currentIndex() ;
      case String:
      case DirName:
      case FileName:
      case RangeInt:
      case RangeLong:
      case RangeFloat:
      case RangeDouble:
          return 0.;
      }

      return 0.;
  };


  void setUnit( const std::string unit ) { mUnit = unit; };

  std::pair<CSParameter::DataType, void *> data() const
  { return std::pair<CSParameter::DataType, void *>(mDataType, mpData); };

  void * dataPointer() const { return mpData; };

  std::string dataString() const;
  std::string unit() const;

  void setParent(CSParameterContext * parent) { CSParameterTreeItem::setParent( (CSParameterTreeItem *)parent ); };
  CSParameterContext * parent() const { return (CSParameterContext *) mpParent; };

  // XML output and input methods
  void writeXML( QXmlStreamWriter * ) const;
  static CSParameter * createFromXML( QXmlStreamReader * );

  void setAttribute( const std::string &attr, const std::string &val ) { mAttributes[attr]=val; };
  const std::string getAttribute( const std::string &attr ) const
  {
      std::map<std::string,std::string>::const_iterator found
          = mAttributes.find(attr);
      if ( found == mAttributes.end() )
          return std::string();
      return found->second;
  };
 protected:
  void *      mpData;
  DataType mDataType;
  std::string  mUnit;
  std::map<std::string,std::string> mAttributes;
};

#endif // CS_PARAMETER_H
