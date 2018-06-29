////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSHDF5Reader.h                                                //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-06-20 15:24:34                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_HDF5_READER_H
#define CS_HDF5_READER_H

#include <string>
#include <sstream>

class CSModel;
namespace H5 {
    class H5File;
}


class CSHDF5Reader
{
public:
    CSHDF5Reader( const std::string & fileName );
    ~CSHDF5Reader();

    void readModelData( CSModel * modelToBePopulated,
                        std::stringstream & errors,
                        std::stringstream & warnings );

private:
    std::string mInputFileName;
    H5::H5File * mpHDF5File;
};

#endif // CS_HDF5_READER_H
