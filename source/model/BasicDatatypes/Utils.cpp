///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Utils.cpp                                                            //
//                                                                                   //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                       //
//    Created:  2016-04-15                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "Utils.h"

BasicDatatypes::Tensor outer(const Vector3f& lhs, const Vector3f& rhs)
{
    BasicDatatypes::Tensor ret;
    ret(0, 0) = lhs.x * rhs.x;
    ret(1, 1) = lhs.y * rhs.y;
    ret(2, 2) = lhs.z * rhs.z;

    ret(0, 1) = lhs.x * rhs.y;
    ret(0, 2) = lhs.x * rhs.z;

    ret(1, 0) = lhs.y * rhs.x;
    ret(1, 2) = lhs.y * rhs.z;

    ret(2, 0) = lhs.z * rhs.x;
    ret(2, 1) = lhs.z * rhs.y;

    return ret;
}

BasicDatatypes::SymTensor outer(const Vector3f& vec)
{
    BasicDatatypes::SymTensor ret;
    ret(0, 0) = vec.x * vec.x;
    ret(1, 1) = vec.y * vec.y;
    ret(2, 2) = vec.z * vec.z;

    ret(0, 1) = vec.x * vec.y;
    ret(0, 2) = vec.x * vec.z;
    ret(1, 2) = vec.y * vec.z;

    return ret;
}

