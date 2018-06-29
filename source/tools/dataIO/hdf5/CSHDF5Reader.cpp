////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSHDF5Reader.cpp                                              //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-06-20 15:28:01                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "CSHDF5Reader.h"

#include "../../../model/Model/CSModel.h"

#include <H5Cpp.h>


CSHDF5Reader::CSHDF5Reader( const std::string & inputFileName )
    : mInputFileName( inputFileName )
{
    mpHDF5File = new H5::H5File( mInputFileName, H5F_ACC_RDONLY );
}


CSHDF5Reader::~CSHDF5Reader()
{
    mpHDF5File->close();
    delete mpHDF5File;
}


void
CSHDF5Reader::readModelData( CSModel * emptyModel,
                             std::stringstream & errors,
                             std::stringstream & warnings )
{
    if ( !emptyModel )
        return;
//    try {
    emptyModel->readModelData( mpHDF5File, errors, warnings );
    emptyModel->Reset(false);
//    }
//    catch {
//    }
}
