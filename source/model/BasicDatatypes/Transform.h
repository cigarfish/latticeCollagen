///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Transform.h                                                          //
//                                                                                   //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                       //
//    Created:  2016-04-15                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <functional>
#include <cmath>
#include <numeric>
#include <utility>
#include <initializer_list>
#include <stdexcept>
#include "tools/concepts/concepts.h"
#include "tools/utils/macros.h"
#include "vector_detail.h"
#include "Matrix.h"

namespace BasicDatatypes {

using BasicDatatypes::Matrix;

template <class T>
class Transform
{
public:
    static constexpr std::size_t N = 3;

    Transform(const Matrix<T, N, N>& _mat)
        : mat(_mat)
    {}

    Transform(Matrix<T, N, N>&& _mat)
        : mat(_mat)
    {}

    Transform(const Matrix<T, N, N>& _mat, const Matrix<T, N, N>& _inv)
        : mat(_mat), inv(_inv)
    {}

    Transform(Matrix<T, N, N>&& _mat, Matrix<T, N, N>&& _inv)
        : mat(_mat), inv(_inv)
    {}

    VectorClass<T, 3> transform(const VectorClass<T, 3>& other) const;
    VectorClass<T, 3> operator()(const VectorClass<T, 3>& other) const;

    Transform inverse() const
    {
        return inv;
    }

    // COMMONLY USED TRANSFORMS
    static Transform identity();

    static Transform rotate_x(T angle);
    static Transform rotate_y(T angle);
    static Transform rotate_z(T angle);
    static Transform rotate(VectorClass<T, 3> axis, T angle);

    static Transform scale(T f);
    static Transform scale(T fx, T fy, T fz);

private:
    Matrix<T, N, N> mat, inv;
};

template <class T>
VectorClass<T, 3> Transform<T>::transform(const VectorClass<T, 3>& other) const
{
    return mat * other;
}

template <class T>
VectorClass<T, 3> Transform<T>::operator()(const VectorClass<T, 3>& other) const
{
    return transform(other);
}

template <class T>
Transform<T> Transform<T>::identity()
{
    auto matrix = Matrix<T, 3, 3>::identity();
    return Transform<T>(matrix, matrix);
}

template <class T>
Transform<T> Transform<T>::rotate_x(T angle)
{
    return Transform<T>(Matrix<T, 3, 3>::rotate_x(angle),
                        Matrix<T, 3, 3>::rotate_x(-angle));
}

template <class T>
Transform<T> Transform<T>::rotate_y(T angle)
{
    return Transform<T>(Matrix<T, 3, 3>::rotate_y(angle),
                        Matrix<T, 3, 3>::rotate_y(-angle));
}

template <class T>
Transform<T> Transform<T>::rotate_z(T angle)
{
    return Transform<T>(Matrix<T, 3, 3>::rotate_z(angle),
                        Matrix<T, 3, 3>::rotate_z(-angle));
}

template <class T>
Transform<T> Transform<T>::rotate(VectorClass<T, 3> axis, T angle)
{
    return Transform<T>(Matrix<T, 3, 3>::rotate(axis, angle),
                        Matrix<T, 3, 3>::rotate(axis, -angle));
}

template <class T>
Transform<T> Transform<T>::scale(T f)
{
    return Transform<T>(Matrix<T, 3, 3>::scale(f),
                        Matrix<T, 3, 3>::scale(1. / f));
}

template <class T>
Transform<T> Transform<T>::scale(T fx, T fy, T fz)
{
    return Transform<T>(Matrix<T, 3, 3>::rotate(fx, fy, fz),
                        Matrix<T, 3, 3>::rotate(1. / fx, 1. / fy, 1. / fz));
}

} //end namespace
