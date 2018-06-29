////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSHDF5Writer.cpp                                              //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-05-27 15:10:28                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "CSHDF5Writer.h"

#include "../../Core.h"
#include "../../model/Model/CSModel.h"
#include "../../model/Cell/CellSpherical.h"

#include <map>

#include <H5Cpp.h>


CSHDF5Writer::CSHDF5Writer( const std::string & fileName )
    : mOutputFileName( fileName )
{
    mpHDF5File = new H5::H5File( mOutputFileName, H5F_ACC_TRUNC );
}


CSHDF5Writer::~CSHDF5Writer()
{
    mpHDF5File->close();
    delete mpHDF5File;
}


void
CSHDF5Writer::exec( CSModel *model )
{
    if (model)
        model->writeHDF5( mpHDF5File );
    else
    {
        std::map<std::string, CSModel *> models = core->models;

        std::map<std::string, CSModel *>::const_iterator modelIterator;

        for ( modelIterator = models.begin(); modelIterator != models.end(); ++modelIterator )
            modelIterator->second->writeHDF5( mpHDF5File );
    }
}
