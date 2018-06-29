////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSXMLReader.cpp                                               //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2012-11-10 21:34:16                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CSXMLReader.h"

#include <QtCore>
#include <iostream>

#include "../../model/Model/CSModel.h"


CSXMLReader::CSXMLReader( const std::string & input )
    : mInputFileName( input )
{}


std::vector<CSModel *>
CSXMLReader::exec( std::stringstream & errors,
                   std::stringstream & warnings )
{
    QFile inputHandle( QString(mInputFileName.c_str()) );

    std::vector<CSModel *> models;

    if ( inputHandle.open( QFile::ReadOnly | QFile::Text ) )
    {
        QXmlStreamReader reader( &inputHandle );

        if ( reader.readNextStartElement() )
        {
            if ( reader.name() == "CellSysXML" )
            {
                // std::cout << "got a CellSysXML version "
                //           << reader.attributes().value("version").toString().toStdString()
                //           << std::endl;
                models = readCellSysXML( &reader, errors, warnings );
            }
            else
            {
                errors << "The file \"" << mInputFileName
                       << "\" is not a CellSysXML file.";
                reader.raiseError(QObject::tr("The file is not a CellSysXML file."));
            }
        }
        else
        {
            errors << "Error:  " << reader.errorString().toStdString() << std::endl;
        }
    }
    else
        errors << "CSXMLReader:  Error:  Could not open file " << mInputFileName << std::endl;

    return models;
}


std::vector<CSModel *>
CSXMLReader::readCellSysXML( QXmlStreamReader * xml,
                             std::stringstream & errors,
                             std::stringstream & warnings )
{
    Q_ASSERT( xml->isStartElement() && xml->name() == "CellSysXML" );

    std::vector<CSModel *> models;

    while ( !xml->atEnd()  )
    {
        CSModel * model = NULL;

        if ( xml->readNextStartElement() )
        {   
            //std::cout << xml->name().toString().toStdString() << std::endl;
            if ( xml->name() == "Model" )
            {
                model = CSModel::ModelFromXML( xml, errors, warnings );
            }
            else
            {
                // std::cout << "got element:  " << xml->name().toString().toStdString() << std::endl;
                xml->skipCurrentElement();
            }
        }

        if (model)
            models.push_back(model);
    }

    return models;
}

