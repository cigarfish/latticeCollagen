////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSParameterTemporary.h                                        //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2012-11-17 21:07:21                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////



#ifndef CS_PARAMETER_TEMPORARY_H
#define CS_PARAMETER_TEMPORARY_H

#include "CSParameter.h"


class CSParameterTemporary : public CSParameter
{
public:
    CSParameterTemporary( CSParameter * );

    ~CSParameterTemporary();
};

#endif // CS_PARAMETER_TEMPORARY_H
