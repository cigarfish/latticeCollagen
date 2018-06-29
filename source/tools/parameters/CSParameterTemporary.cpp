////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSParameterTemporary.cpp                                      //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2012-11-17 21:10:44                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "CSParameterTemporary.h"


CSParameterTemporary::CSParameterTemporary( CSParameter * otherParm )
    : CSParameter( otherParm->name(), otherParm->dataType(), NULL, otherParm->unit() )
{
    // copy the data into a local pointer
    void * otherData = otherParm->dataPointer();

    switch ( otherParm->dataType() )
    {
    case CSParameter::Bool:
        setData( *(bool *) otherData );
        break;
    case CSParameter::Int:
    case CSParameter::Long:
        setData( *(long *) otherData );
        break;
    case CSParameter::Float:
    case CSParameter::Double:
        setData( *(double *) otherData );
        break;
    case CSParameter::String:
    case CSParameter::DirName:
    case CSParameter::FileName:
        setData( *(std::string *) otherData );
        break;
    case CSParameter::RangeInt:
    case CSParameter::RangeLong:
    {
        std::pair<long, long> * longRange =
            (std::pair<long, long> *) otherData;
        setData( longRange->first, longRange->second );
    }
    break;
    case CSParameter::RangeFloat:
    case CSParameter::RangeDouble:
    {
        std::pair<double, double> * doubleRange =
            (std::pair<double, double> *) otherData;
        setData( doubleRange->first, doubleRange->second );
    }
    break;
    case CSParameter::Choice:
        CSParameterChoice *choice = static_cast<CSParameterChoice *>( otherData );
        setData( choice->choices(), choice->currentIndex() );
        break;
    }
}


CSParameterTemporary::~CSParameterTemporary()
{
    switch ( mDataType )
    {
    case CSParameter::Bool:
        delete (bool *) mpData;
        break;
    case CSParameter::Int:
    case CSParameter::Long:
        delete (long *) mpData;
        break;
    case CSParameter::Float:
    case CSParameter::Double:
        delete (double *) mpData;
        break;
    case CSParameter::String:
    case CSParameter::DirName:
    case CSParameter::FileName:
        delete (std::string *) mpData;
        break;
    case CSParameter::RangeInt:
    case CSParameter::RangeLong:
        delete (std::pair<long, long> *) mpData;
        break;
    case CSParameter::RangeFloat:
    case CSParameter::RangeDouble:
        delete (std::pair<double, double> *) mpData;
        break;
    case CSParameter::Choice:
        delete (CSParameterChoice *) mpData;
        break;
    }
}
