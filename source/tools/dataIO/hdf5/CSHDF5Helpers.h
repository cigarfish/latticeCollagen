////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSHDF5Helpers.h                                               //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-06-03 14:44:29                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_HDF5_HELPERS_H
#define CS_HDF5_HELPERS_H


#define H5_DOUBLE_ARRAY( _name, _length )                           \
    const hsize_t _name ## _Len = _length;                         \
    H5::ArrayType _name ( H5::PredType::NATIVE_DOUBLE, 1, &_name ## _Len ); {}

#define H5_LONG_ARRAY( _name, _length )                             \
    const hsize_t _name ## _Len = _length;                         \
    H5::ArrayType _name ( H5::PredType::NATIVE_LONG, 1, &_name ## _Len ); {}

#define H5_BOOL_ARRAY( _name, _length )                             \
    const hsize_t _name ## _Len = _length;                         \
    H5::ArrayType _name ( H5::PredType::NATIVE_B8, 1, &_name ## _Len ); {}



#endif // CS_HDF5_HELPERS_H
