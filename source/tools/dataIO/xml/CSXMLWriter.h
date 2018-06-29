///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSXMLWriter.h                                                       //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-08-09 16:37:01                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef Q_CS_XML_WRITER_H
#define Q_CS_XML_WRITER_H

class QIODevice;

#include <string>

class CSModel;
class QXmlStreamWriter;


class CSXMLWriter
{
public:
    CSXMLWriter( const std::string &ouptutFileName );

    void exec();
    void exec( CSModel *model );

private:
    void writeXMLHeader( QXmlStreamWriter * writer );
    std::string mOutputFileName;
};

#endif // Q_CS_XML_WRITER_H
