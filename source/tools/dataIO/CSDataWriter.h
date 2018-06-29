////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSDataWriter.h                                                //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-06-20 17:06:18                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CSDATAWRITER_H
#define CSDATAWRITER_H

#include <string>

class CSModel;


class CSDataWriter
{
public:
    CSDataWriter( const std::string & bundleName );

    void exec( CSModel * model =0 );

private:
    std::string mDirectoryName;
};

#endif // CSDATAWRITER_H
