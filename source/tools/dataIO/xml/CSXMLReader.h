////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSXMLReader.h                                                 //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2012-11-10 19:30:52                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////



#ifndef CS_XML_READER_H
#define CS_XML_READER_H

#include <string>
#include <vector>
#include <sstream>

class CSModel;
class QXmlStreamReader;


class CSXMLReader
{
public:
    CSXMLReader( const std::string & xmlFileName );

    std::vector<CSModel *> exec( std::stringstream & errors,
                               std::stringstream & warnings );

private:
    std::vector<CSModel *> readCellSysXML( QXmlStreamReader *,
                                         std::stringstream & errors,
                                         std::stringstream & warnings );

    std::string mInputFileName;
};

#endif // CS_XML_READER_H
