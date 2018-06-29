////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSDataReader.cpp                                              //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-06-21 18:30:22                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "CSDataReader.h"

#include "xml/CSXMLReader.h"
#include "hdf5/CSHDF5Reader.h"

// for checking the files' existence:
#include <fstream>


CSDataReader::CSDataReader( const std::string & bundleName )
    : mDirectoryName( bundleName )
{}


std::vector<CSModel *>
CSDataReader::exec( std::stringstream & errors, std::stringstream & warnings )
{
    std::vector<CSModel *> models;

    std::string xmlFile = mDirectoryName + "/CSModelDescription.xml";

    {
        std::ifstream file(xmlFile.c_str());
        bool exists = file.is_open();

        if ( !exists )
        {
            errors << "Error in CSDataReader:  Directory \""
                   << mDirectoryName << "\" is not a CellSys data bundle." << std::endl;
            return models;
        }

        file.close();
    }

    std::string hdf5File = mDirectoryName + "/CSModelData.hdf5";

    {
        std::ifstream file(hdf5File.c_str());
        bool exists = file.is_open();

        if ( !exists )
        {
            errors << "Error in CSDataReader:  Directory \""
                   << mDirectoryName << "\" does not a contain HDF5 data." << std::endl;
            return models;
        }

        file.close();
    }

    CSXMLReader xmlReader( xmlFile );

    models = xmlReader.exec( errors, warnings );

    CSHDF5Reader hdf5Reader( hdf5File );

    std::vector<CSModel *>::iterator modelIter;

    for ( modelIter=models.begin(); modelIter!=models.end(); ++modelIter )
        hdf5Reader.readModelData( *modelIter, errors, warnings );

    return models;
}
