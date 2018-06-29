////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSBatchJob.cpp                                                //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-01-23 18:59:25                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CSBatchJob.h"

#include "../../model/Model/CSModel.h"
#include "../../model/Model/ModelCellsTriangulated/Model3D.h"
#include "../../model/Model/ModelCellsSpherical/ModelCellsSpherical.h"

#if defined( CS_BUILD_IMAGEPROCESSING )
#  include "../JobManager.h"
#endif

#include "../dataIO/CSDataReader.h"
#include "../dataIO/xml/CSXMLReader.h"

#include <string>
#include <iostream>

#include <QFile>
#include <QFileInfo>
#include <QXmlStreamReader>

    
CSBatchJob::CSBatchJob( const std::string & inputPath )
    :
#if defined( CS_BUILD_IMAGEPROCESSING )
      mpJobManager(NULL),
#endif
      mValidFlag( true )
{
    std::stringstream errors;   // stream for  errors  that may occur in the reading process
    std::stringstream warnings; // stream for warnings that may occur in the reading process

    // see if file exists
    QString qPathString( inputPath.c_str() ); 
    QFile input( qPathString );


    if ( !input.exists() )
    {
        std::cerr << "\nIn CSBatchJob:  Error:  input file/bundle does not exist.\n";
        mValidFlag = false;
        return;
    }

    QFileInfo finfo( input );
    if ( !finfo.isDir() )
    {
        if ( !input.open(QFile::ReadOnly) )
        {
            std::cerr << "\nIn CSBatchJob:  Error:  input file not readable, check permissions\n";
            mValidFlag = false;
            return;
        }

        QXmlStreamReader xmlInput( &input );

        // is it a valid xml file
        // read header:
        if ( xmlInput.readNextStartElement() )
        {
#if defined( CS_BUILD_IMAGEPROCESSING )
            if ( xmlInput.name() == "ParameterContext" && xmlInput.attributes().value("name") == "JobManager" )
            {
                mpJobManager = new JobManager();
                mpJobManager->InitWithJobQueueFile( inputPath );
                input.close();
                return;
            }
#endif

            if ( xmlInput.name() != "CellSysXML" )
            {
                std::cerr << "\nIn CSBatchJob:  Error: input file is not a CellSysXML file.\n";
                mValidFlag = false;
                input.close();
                return;
            }
        }
        else
        {
            std::cerr << "\nIn CSBatchJob:  Error:  input file not in XML format.\n";
            std::cerr << xmlInput.errorString().toStdString() << std::endl;
            mValidFlag = false;
            input.close();
            return;
        }
        // Read single (legacy) CellSysXML file
        CSXMLReader xmlReader( inputPath );
        mModels = xmlReader.exec( errors, warnings );
    }
    else // if we read a data bundle
    {
        CSDataReader bundleReader( inputPath );
        mModels = bundleReader.exec( errors, warnings );
    }

    if ( errors.str().size() )
    {
        std::cerr << errors.str();
        mValidFlag = false;
    }
    else
        std::cerr << warnings.str();

    std::cerr << std::endl;
}


int
CSBatchJob::run()
{
    if ( !mValidFlag )
        return -1;

#if defined( CS_BUILD_IMAGEPROCESSING )
    if ( mpJobManager )
    {
        mpJobManager->Start();
        return 0;
    }
#endif

    if ( !mModels.size() )
        return 3;

    int returnCode = 0;
    std::vector<CSModel *>::iterator modelIt;
    for ( modelIt = mModels.begin(); modelIt != mModels.end(); ++modelIt )
    {
        (*modelIt)->SetupSimulation();

        if ( (*modelIt)->Run(true) )
            std::cerr << "CSBatchJob::run():  An error occured in the simulation.\n";
        // else if ( (modelIt)->state() )
        // {
        //     std::cerr << "CSBatchJob::run():  Error occurred when simulating model \""
        //               << (*modelIt)->name << "\"" << std::endl;
        //     returnCode = 4;
        // }
    }

    return returnCode;
}
