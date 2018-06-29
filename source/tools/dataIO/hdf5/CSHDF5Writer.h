////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSHDF5Writer.h                                                //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-05-27 15:02:49                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_HDF5_WRITER_H
#define CS_HDF5_WRITER_H

#include <string>

class CSModel;
namespace H5 {
    class H5File;
}


class CSHDF5Writer
{
public:
    CSHDF5Writer( const std::string & outputFileName );
    ~CSHDF5Writer();

    void exec( CSModel *model =0 );

private:
    std::string mOutputFileName;
    H5::H5File * mpHDF5File;
};

#endif // CS_HDF5_WRITER_H
