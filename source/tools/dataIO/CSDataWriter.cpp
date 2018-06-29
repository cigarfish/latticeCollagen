////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSDataWriter.cpp                                              //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-06-20 17:11:46                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////



#include "CSDataWriter.h"
#include "xml/CSXMLWriter.h"
#include "hdf5/CSHDF5Writer.h"

#include "../../model/Model/CSModel.h"

#include <string>


CSDataWriter::CSDataWriter( const std::string & bundleName )
    : mDirectoryName(bundleName)
{}


void
CSDataWriter::exec( CSModel * model )
{
    std::string hdf5File = mDirectoryName + "/CSModelData.hdf5";
    CSHDF5Writer hdf5Writer( hdf5File );

    if (model)
        hdf5Writer.exec(model);
    else
        hdf5Writer.exec();

    std::string xmlFile = mDirectoryName + "/CSModelDescription.xml";
    CSXMLWriter xmlWriter( xmlFile );

    if (model)
        xmlWriter.exec(model);
    else
        xmlWriter.exec();
}
