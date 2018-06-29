///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSXMLWriter.cpp                                                      //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-08-09 17:01:19                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include <QtCore>

#include "CSXMLWriter.h"
#include "../../Core.h"
#include "../../model/Model/CSModel.h"
#include "../parameters/CSParameterContext.h"

#include <iostream>


CSXMLWriter::CSXMLWriter( const std::string &fileName )
  : mOutputFileName(fileName)
{}


void
CSXMLWriter::exec()
{
    QFile outputHandle( QString(mOutputFileName.c_str()) );

    if ( outputHandle.open( QFile::WriteOnly | QFile::Text ) )
    {
        QXmlStreamWriter writer( &outputHandle );

        writeXMLHeader( &writer );

        std::map<std::string, CSModel *>::const_iterator modelIt;

        for ( modelIt = core->models.begin(); modelIt != core->models.end(); ++modelIt )
            modelIt->second->writeXML(&writer);

        writer.writeEndDocument();
    }
    else
        std::cerr << "CSXMLWriter:  Could not open file " << mOutputFileName << std::endl;
}


void
CSXMLWriter::exec( CSModel *model )
{
    QFile outputHandle( QString(mOutputFileName.c_str()) );

    if ( outputHandle.open( QFile::WriteOnly | QFile::Text ) )
    {
        QXmlStreamWriter writer( &outputHandle );

        writeXMLHeader( &writer );

        model->writeXML( &writer );

        writer.writeEndDocument();
    }
    else
        std::cerr << "CSXMLWriter:  Could not open file " << mOutputFileName << std::endl;   
}


void
CSXMLWriter::writeXMLHeader( QXmlStreamWriter * writer )
{
    writer->setAutoFormatting(true);

    writer->writeStartDocument();

    writer->writeDTD("<!DOCTYPE CellSysXML>");
    writer->writeStartElement("CellSysXML");
    writer->writeAttribute("version", "0.01");
}
