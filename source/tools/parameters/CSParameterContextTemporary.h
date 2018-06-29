////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSParameterContextTemporary.h                                 //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2012-11-17 20:28:45                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_PARAMETER_CONTEXT_TEMPORARY_H
#define CS_PARAMETER_CONTEXT_TEMPORARY_H

#include "CSParameterContext.h"


class CSParameterContextTemporary : public CSParameterContext
{
public:
    CSParameterContextTemporary( CSParameterContext * );
    CSParameterContextTemporary( const std::string & );

    ~CSParameterContextTemporary();

    static CSParameterContextTemporary * createFromXML( QXmlStreamReader * xml );
};

#endif //  CS_PARAMETER_CONTEXT_TEMPORARY_H
