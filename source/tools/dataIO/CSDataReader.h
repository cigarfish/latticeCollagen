////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSDataReader.h                                                //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-06-20 17:28:50                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CSDATAREADER_H
#define CSDATAREADER_H

#include <string>
#include <vector>
#include <sstream>

class CSModel;


class CSDataReader
{
public:
    CSDataReader( const std::string & bundleName );

    std::vector<CSModel *> exec( std::stringstream & errors,
                               std::stringstream & warnings );

private:
    std::string mDirectoryName;
};

#endif // CSDATAREADER_H
