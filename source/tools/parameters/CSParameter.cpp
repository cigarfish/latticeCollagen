///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSParameter.cpp                                                      //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-03-14 17:37:56                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include <sstream>

#include "CSParameter.h"
#include "CSParameterContext.h"

#include <QtCore>

// t1m-debug
#include <iostream>
// ! t1m-debug

const char * CSParameter::DataTypeString[] = {
    "Bool",
    "Int",
    "Long",
    "Float",
    "Double",
    "String",
    "DirName",
    "FileName",
    "RangeInt",
    "RangeLong",
    "RangeFloat",
    "RangeDouble",
    "Choice"
};




CSParameter::CSParameter( std::string name,
                          CSParameter::DataType type,
                          void * data,
                          const std::string & unit,
                          CSParameterContext * parent )
  : CSParameterTreeItem( name, CSParameterTreeItem::Parameter, parent ),
    mpData(data),
    mDataType(type),
    mUnit(unit)
{}


// CSParameter::CSParameter()
// {}


CSParameter::CSParameter( const CSParameter * otherParm )
    : CSParameterTreeItem( otherParm->name(), CSParameterTreeItem::Parameter ),
      mpData( otherParm->dataPointer() ),
      mDataType( otherParm->dataType() ),
      mUnit( otherParm->unit() )
{}


std::string
CSParameter::dataString() const
{
  std::stringstream output;

  switch (mDataType)
    {
    case CSParameter::Bool:
      output << ((*(bool *)mpData)?"yes":"no");
      break;

    case CSParameter::Int:
      output << *((int *) mpData);
      break;

    case CSParameter::Long:
      output << *((long *) mpData);
      break;

    case CSParameter::Float:
      output << *((float *) mpData);
      break;

    case CSParameter::Double:
      output << *((double *) mpData);
      break;

    case CSParameter::String:
    case CSParameter::DirName:
    case CSParameter::FileName:
      output << *((std::string *) mpData);
      break;

    case CSParameter::RangeInt:
      output << ((std::pair<int,int> *) mpData)->first
             << "-"
             << ((std::pair<int,int> *) mpData)->second;
      break;

    case CSParameter::RangeLong:
      output << ((std::pair<long,long> *) mpData)->first
             << "-"
             << ((std::pair<long,long> *) mpData)->second;
      break;

    case CSParameter::RangeFloat:
      output << ((std::pair<float,float> *) mpData)->first
             << "-"
             << ((std::pair<float,float> *) mpData)->second;
      break;
    case CSParameter::RangeDouble:
      output << static_cast< std::pair<double,double> *>( mpData )->first
             << "-"
             << static_cast< std::pair<double,double> *>( mpData )->second;
      break;

    case CSParameter::Choice:
      output << static_cast<CSParameterChoice *>( mpData )->currentString();
      break;
    }

  return output.str();
}


std::string
CSParameter::unit() const
{
  if ( mDataType != String   &&
       mDataType != DirName  &&
       mDataType != FileName &&
       mDataType != Choice )
    {
      return mUnit;
    }
  else
    return std::string();
}


void
CSParameter::writeXML( QXmlStreamWriter * xml ) const
{
    xml->writeStartElement( "Parameter" );
    xml->writeAttribute( "name", this->name().c_str() );

    CSParameter::DataType dataType = this->dataType();

    xml->writeAttribute( "type", CSParameter::DataTypeString[dataType] );

    std::string valueString;

    switch (dataType)
    {
    case CSParameter::Bool:
    case CSParameter::Int:
    case CSParameter::Long:
    case CSParameter::Float:
    case CSParameter::Double:
    case CSParameter::String:
    case CSParameter::DirName:
    case CSParameter::FileName:
        valueString =  this->dataString();
        break;
    case CSParameter::RangeInt:
    case CSParameter::RangeLong:
    case CSParameter::RangeFloat:
    case CSParameter::RangeDouble:
        valueString = this->dataString();
        break;
    case CSParameter::Choice:
        //concatenate the options with a colon delimiter
        std::stringstream concatChoices;
        CSParameterChoice * ptr =
            (CSParameterChoice *)this->dataPointer();
        std::vector<std::string> choices = ptr->choices();

        for( unsigned int i=0; i < choices.size(); ++i )
        {
            // if ( i == ptr->currentIndex() )
            //     concatChoices << "*";
            concatChoices << choices[i];
            // if ( i == ptr->currentIndex() )
            //     concatChoices << "*";
            concatChoices << ":";
        }
        concatChoices << ptr->currentIndex();
        valueString = concatChoices.str();
    }

    xml->writeAttribute( "value", valueString.c_str() );

    xml->writeAttribute( "unit", this->unit().c_str() );
    xml->writeEndElement();
}


CSParameter *
CSParameter::createFromXML( QXmlStreamReader * xml )
{
    CSParameter::DataType parmType;
    void * parmPointer;
    std::string parmName;
    std::string parmUnit;

    parmName = xml->attributes().value("name").toString().toStdString();
    parmUnit = xml->attributes().value("unit").toString().toStdString();

    std::string parmTypeString = xml->attributes().value("type").toString().toStdString();
    QString parmValueString = xml->attributes().value("value").toString();

    if ( parmTypeString == "Bool" )
    {
        parmType = CSParameter::Bool;
        std::string yesno = parmValueString.toStdString();
        parmPointer = (void *) new bool( (yesno=="yes") ? true : false );
    }
    else if ( parmTypeString == "Int" )
    {
        parmType = CSParameter::Int;
        parmPointer = (void *) new long(parmValueString.toLong());
    }
    else if ( parmTypeString == "Long" )
    {
        parmType = CSParameter::Long;
        parmPointer = (void *) new long(parmValueString.toLong());

    }
    else if ( parmTypeString == "Float" )
    {
        parmType = CSParameter::Float;
        parmPointer = (void *) new double(parmValueString.toDouble());
    }
    else if ( parmTypeString == "Double" )
    {
        parmType = CSParameter::Double;
        parmPointer = (void *) new double(parmValueString.toDouble());
    }
    else if ( parmTypeString == "String" )
    {
        parmType = CSParameter::String;
        parmPointer = (void *) new std::string(parmValueString.toStdString());
    }
    else if ( parmTypeString == "DirName" )
    {
        parmType = CSParameter::DirName;
        parmPointer = (void *) new std::string(parmValueString.toStdString());
    }
    else if ( parmTypeString == "FileName" )
    {
        parmType = CSParameter::FileName;
        parmPointer = (void *) new std::string(parmValueString.toStdString());
    }
    // Range* not yet supported
    else if ( parmTypeString == "Choice" )
    {
        parmType = CSParameter::Choice;
        QStringList choices = parmValueString.split(QChar(':'));
        int currentIndex = choices.last().toInt();
        std::vector<std::string> choicesString;
        // size() should be greater than 1, maybe assert()?
        for ( int i=0; i<choices.size()-1; ++i )
            choicesString.push_back( choices.at(i).toStdString() );
        parmPointer = (void *) new CSParameterChoice( choicesString, currentIndex );
    }
    else
    {
        std::cout << "unknown parameter type\n";
        return NULL;
    }

    // while ( xml->readNextStartElement() )
    // {
    //     if ( xml->name() == "parameter" )
    //         readParameter(xml);
    //     else
    //         xml->skipCurrentElement();
    // }

    // std::cout << "parameter: "
    //           << parmName
    //           << " = "
    //           << parmValueString.toStdString() << " " << parmUnit
    //           << std::endl;

    return new CSParameter( parmName, parmType, parmPointer, parmUnit );
}
