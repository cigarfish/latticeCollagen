///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Matrix.h                                                             //
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

using namespace std::placeholders;

namespace BasicDatatypes {

template <typename Derived, typename T, std::size_t N>
    void swap(BaseVector<Derived, T, N>& lhs, BaseVector<Derived, T, N>& rhs);

template <typename T, std::size_t ORDER>
struct Matrix_init
{
    using type = std::initializer_list<typename Matrix_init<T, ORDER-1>::type >;
};

template <typename T>
struct Matrix_init<T, 1>
{
    using type = std::initializer_list<T>;
};

template<typename T>
struct Matrix_init<T, 0>;

template<typename T, std::size_t ORDER>
using Matrix_initializer = typename Matrix_init<T, ORDER>::type;

template <class T, class DATA>
inline void insert_flat(const T* first, const T* last, DATA& data)
{
    data.insert(data.cback(), first, last);
}

template <class T, class DATA>
inline void insert_flat(const std::initializer_list<T>* first, const std::initializer_list<T>* last,
                        DATA& data)
{
    for (;first!=last;++first)
        insert_flat(first->begin(), first->end(), data);
}

template <class T, class DATA>
inline void insert_flat(const std::initializer_list<T>& list, DATA& data)
{
    insert_flat(list.begin(), list.end(), data);
}

//
// TODO
//
// allow access to most common used members by name
// ie. add something like
//
// union {
//     T data[9];
//     T xx, xy, xz, ...;
//     };
//
template <typename Derived, typename T, std::size_t N, std::size_t M>
class MatrixBase
{
public:
    static constexpr std::size_t ORDER = 2;

    /*constexpr*/ MatrixBase()
        : MatrixBase(T(0))
    {}

    MatrixBase(T (&data)[N]) // avoid pointer decay
    {
        std::copy(std::begin(data), std::end(data), begin());
    }

    explicit MatrixBase(const MatrixBase& other)
    {
        std::copy(other.begin(), other.end(), begin());
    }

    explicit MatrixBase(MatrixBase&& other)
    {
        //swap(other, *this);
    }

    explicit MatrixBase(const T value)
    {
        set(value);
    }

    //MatrixBase(std::initializer_list<T> values)
    //    : MatrixBase(T(0))
    //{
    //    ASSERT(values.size() <= N * M, "Initializer list with size " +
    //           std::to_string(values.size()) + " is larger than the container" +
    //           " size " + std::to_string(N) + ".");
    //    std::copy(values.begin(), values.end(), begin());
    //}

    MatrixBase(Matrix_initializer<T, ORDER> values)
        : MatrixBase(T(0))
    {
        // HACK
        current_end = begin();
        ASSERT(values.size() <= N * M, "Initializer list with size " +
               std::to_string(values.size()) + " is larger than the container" +
               " size " + std::to_string(N) + ".");
        insert_flat(values, *this);
    }

    MatrixBase& operator=(const MatrixBase& other)
    {
        std::copy(other.begin(), other.end(), begin());
        return *this;
    }

    MatrixBase& operator=(MatrixBase&& matrix)
    {
        swap(matrix, *this);
        return *this;
    }

    MatrixBase& operator=(const T value)
    {
        set(value);
        return *this;
    }

    //using iterator = typename std::vector<T>::iterator;
    //using const_iterator = typename std::vector<T>::const_iterator;

    MatrixBase& operator+=(const MatrixBase& rhs)
    {
        return operator_base(rhs, std::plus<T>());
    }

    MatrixBase& operator-=(const MatrixBase& rhs)
    {
        return operator_base(rhs, std::minus<T>());
    }

    MatrixBase& operator*=(const T& scalar)
    {
        std::transform(begin(), end(), begin(),
                       std::bind1st(std::multiplies<T>(), scalar));
        return *this;
    }

    MatrixBase& operator/=(const T& scalar)
    {
        std::transform(begin(), end(), begin(),
                       std::bind2nd(std::divides<T>(), scalar));
        return *this;
    }

    ~MatrixBase() {}

    T* begin() { return static_cast<Derived*>(this)->begin(); }
    const T* begin() const { return static_cast<const Derived*>(this)->begin(); }

    T* end() { return static_cast<Derived*>(this)->end(); }
    const T* end() const { return static_cast<const Derived*>(this)->end(); }

    const T* cbegin() { return begin(); }
    const T* cend() const { return end(); }

    T* data() { return begin(); }
    const T* data() const { return begin(); }

    constexpr std::size_t size() const
    {
        return static_cast<const Derived*>(this)->size();
    }

    std::size_t memory_size() const
    {
        return static_cast<const Derived*>(this)->memory_size();
    }

    T& get(const std::size_t i, const std::size_t j)
    {
        return operator[](index(i, j));
    }

    const T& get(const std::size_t i, const std::size_t j) const
    {
        return operator[](index(i, j));
    }

    T& operator()(const std::size_t i, const std::size_t j)
    {
        return get(i, j);
    }

    const T& operator()(const std::size_t i, const std::size_t j) const
    {
        return get(i, j);
    }

    // Access operator
    T& operator[](const std::size_t i) { return *(begin() + i); }
    const T& operator[](const std::size_t i) const { return *(begin() + i); }

    T& at(const std::size_t i)
    {
        if (!(i < size()))
            throw std::out_of_range("Index " + std::to_string(i) + " larger than "
                                    + std::to_string(size()));

        return operator[](i);
    }

    const T& at(const std::size_t i) const
    {
        if (!(i < size()))
            throw std::out_of_range("Index " + std::to_string(i) + " larger than "
                                    + std::to_string(size()));

        return operator[](i);
    }

    T& front() { return operator[](0); }
    const T& front() const { return operator[](0); }

    T& back()
    {
        return operator[](size() - 1);
    }

    const T& back() const
    {
        return operator[](size() - 1);
    }

    T* cback()
    {
        return current_end;
    }

    const T* cback() const
    {
        return current_end;
    }

    void swap(MatrixBase& other)
    {
        std::swap_ranges(begin(), end(), other.begin());
    }

    void set(const T value)
    {
        //std::cout << "begin:" << begin() << " end:" << end() << " val:"
        //    << value << std::endl;
        //for (std::size_t i = 0; i < size(); i++)
        //{
        //    (this)->at(i) = value;
        //}
        // WHY DOESNT THIS WORK!!
        std::fill(begin(), end(), value);
    }

    constexpr void Reset()
    {
        set(T(0));
    }

    // HACK
    T* insert(T * pos, const T* begin, const T* end)
    {
        auto this_begin = pos;
        for(;begin!=end; ++begin)
        {
            *this_begin = *begin;
            this_begin++;
        }

        current_end = this_begin;
        return pos;
    }

    void print()
    {
        std::cout << "[";
        for (std::size_t i = 0; i < N; ++i)
        {
            if (i != 0) std::cout << " ";
            std::cout << "[";
            for (std::size_t j = 0; j < M; ++j)
            {
                std::cout << get(i, j);
                if (j < M-1) std::cout << " ";
            }
            std::cout << "]";
            if (i < N - 1) std::cout << ", " << std::endl;
        }
        std::cout << "]" << std::endl;
    }

private:

    const std::size_t index(std::size_t i, std::size_t j) const
    {
        return static_cast<const Derived*>(this)->index(i, j);
    }

    template <class Function>
    MatrixBase& operator_base(const MatrixBase& rhs, Function f)
    {
        std::transform(rhs.begin(), rhs.end(), begin(), begin(), f);
        return *this;
    }

    // HACK
    T * current_end;
};

// swap
template <typename Derived, typename T, std::size_t N, std::size_t M>
void swap(MatrixBase<Derived, T, N, M>& lhs, MatrixBase<Derived, T, N, M>& rhs)
{
    lhs.swap(rhs);
}

template <typename T, std::size_t N, std::size_t M>
class Matrix : public MatrixBase<Matrix<T, N, M>, T, N, M>
{
public:
    using MatrixBase<Matrix, T, N, M>::MatrixBase;
    using MatrixBase<Matrix<T, 3, 3>, T, 3, 3>::ORDER;

    // constructors
    Matrix()
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>()
    {}

    Matrix(T (&data)[2])
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>(data)
    {}

    Matrix(const Matrix& vector)
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>(vector)
    {}

    Matrix(Matrix&& vector)
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>(vector)
    {}

    explicit Matrix(const T value)
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>(value)
    {}

    Matrix(Matrix_initializer<T, ORDER> values)
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>(values)
    {}

    Matrix& operator=(const Matrix& vector)
    {
        MatrixBase<Matrix<T, 3, 3>, T, 3, 3>::operator = (vector);
        return *this;
    }

    Matrix& operator=(Matrix&& vector)
    {
        MatrixBase<Matrix<T, 3, 3>, T, 3, 3>::operator = (vector);
        return *this;
    }

    Matrix& operator=(const T value)
    {
        MatrixBase<Matrix<T, 3, 3>, T, 3, 3>::operator = (value);
        return *this;
    }

    //Matrix(T _x, T _y)
    //    : x(_x), y(_y)
    //{}

    ~Matrix() {}

    T data[N * M];

    std::size_t index(const std::size_t i, const std::size_t j) const
    {
        return (i * N + j);
    }

    T* begin() { return std::begin(data); }
    const T* begin() const { return std::begin(data); }

    T* end() { return std::end(data); }
    const T* end() const { return std::end(data); }

    constexpr std::size_t size() const { return N * M; }
    std::size_t memory_size() const { return sizeof(*this); }

    // COMMONLY USED STUFF
    static Matrix zero();
    static Matrix identity();
};

template <typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> Matrix<T, N, M>::zero()
{
    Matrix<T, N, M> m;
    return m;
}

template <typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> Matrix<T, N, M>::identity()
{
    Matrix<T, N, M> m;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < M; ++j)
            if (i == j)
                m.get(i, j) = 1;

    return m;
}

// specialize this
template <typename T>
class Matrix<T, 3, 3> : public MatrixBase<Matrix<T, 3, 3>, T, 3, 3>
{
public:
    using MatrixBase<Matrix<T, 3, 3>, T, 3, 3>::ORDER;

    // constructors
    Matrix()
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>()
    {}

    Matrix(T (&data)[2])
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>(data)
    {}

    Matrix(const Matrix& vector)
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>(vector)
    {}

    Matrix(Matrix&& vector)
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>(vector)
    {}

    explicit Matrix(const T value)
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>(value)
    {}

    Matrix(Matrix_initializer<T, ORDER> values)
        : MatrixBase<Matrix<T, 3, 3>, T, 3, 3>(values)
    {}

    Matrix& operator=(const Matrix& vector)
    {
        MatrixBase<Matrix<T, 3, 3>, T, 3, 3>::operator = (vector);
        return *this;
    }

    Matrix& operator=(Matrix&& vector)
    {
        MatrixBase<Matrix<T, 3, 3>, T, 3, 3>::operator = (vector);
        return *this;
    }

    Matrix& operator=(const T value)
    {
        MatrixBase<Matrix<T, 3, 3>, T, 3, 3>::operator = (value);
        return *this;
    }

    ~Matrix() {}

    // DATA
    union {
        T data[9];
        struct {T xx, xy, xz, yx, yy, yz, zx, zy, zz; };
        // Add vectors
    };

    std::size_t index(const std::size_t i, const std::size_t j) const
    {
        return (i * 3 + j);
    }

    T* begin() { return std::begin(data); }
    const T* begin() const { return std::begin(data); }

    T* end() { return std::end(data); }
    const T* end() const { return std::end(data); }

    constexpr std::size_t size() const { return 9; }
    std::size_t memory_size() const { return sizeof(*this); }

    // COMMONLY USED MATRICES
    static Matrix zero();
    static Matrix identity();

    static Matrix rotate_x(T angle);
    static Matrix rotate_y(T angle);
    static Matrix rotate_z(T angle);
    static Matrix rotate(VectorClass<T, 3> axis, T angle);

    static Matrix scale(T f);
    static Matrix scale(T fx, T fy, T fz);
};

template <typename T>
Matrix<T, 3, 3> Matrix<T, 3, 3>::zero()
{
    return {{0,0,0}, {0,0,0}, {0,0,0}};
}

template <typename T>
Matrix<T, 3, 3> Matrix<T, 3, 3>::identity()
{
    return {{1,0,0}, {0,1,0}, {0,0,1}};
}

template <class T>
Matrix<T, 3, 3> Matrix<T, 3, 3>::rotate_x(T angle)
{
    T c = std::cos(angle);
    T s = std::sin(angle);

    return {{1.0, 0.0, 0.0},
            {0.0, c  ,  -s},
            {0.0, s  ,   c}};
}

template <class T>
Matrix<T, 3, 3> Matrix<T, 3, 3>::rotate_y(T angle)
{
    T c = std::cos(angle);
    T s = std::sin(angle);

    return {{c  , 0.0,   s},
            {0.0, 1.0, 0.0},
            {-s , 0.0,   c}};
}

template <class T>
Matrix<T, 3, 3> Matrix<T, 3, 3>::rotate_z(T angle)
{
    T c = std::cos(angle);
    T s = std::sin(angle);

    return {{  c,  -s, 0.0},
            {  s,   c, 0.0},
            {0.0, 0.0, 1.0}};
}

template <class T>
Matrix<T, 3, 3> Matrix<T, 3, 3>::rotate(VectorClass<T, 3> axis, T angle)
{
    T c = std::cos(angle);
    T s = std::sin(angle);
    T cc = 1.0 - c;

    T t1 = axis.x * axis.y * cc;
    T t2 = axis.x * axis.z * cc;
    T t3 = axis.y * axis.z * cc;

    T u1 = axis.x * s;
    T u2 = axis.y * s;
    T u3 = axis.z * s;

    return {{axis.x * axis.x * cc + c, t1 - u3, t2 + u2},
            {t1 + u3, axis.y * axis.y * cc + c, t3 - u1},
            {t2 - u2, t3 + u1, axis.z * axis.z * cc + c}};
}

template <class T>
Matrix<T, 3, 3> Matrix<T, 3, 3>::scale(T f)
{
    return {{  f, 0.0, 0.0},
            {0.0,   f, 0.0},
            {0.0, 0.0,   f}};
}

template <class T>
Matrix<T, 3, 3> Matrix<T, 3, 3>::scale(T fx, T fy, T fz)
{
    return {{ fx, 0.0, 0.0},
            {0.0,  fy, 0.0},
            {0.0, 0.0,  fz}};
}

template <typename T, std::size_t N>
class SymmetricMatrix : public MatrixBase<SymmetricMatrix<T, N>, T, N, N>
{
public:
    using MatrixBase<SymmetricMatrix, T, N, N>::MatrixBase;

    SymmetricMatrix()
        : MatrixBase<SymmetricMatrix, T, N, N>()
    {}

    std::size_t index(const std::size_t i, const std::size_t j) const
    {
        return (i <= j) ? offset(i, j) : offset(j, i);
    }

    // DATA
    union {
        T data[6];
        struct {T xx, xy, xz, yx, yy, zz; };
        // Add vectors
    };

private:
    std::size_t offset(const std::size_t i, const std::size_t j) const
    {
        return (i * N - (i - 1) * i / 2 + j - i);
    }

};

template <typename Derived, typename T, std::size_t N>
T trace(MatrixBase<Derived, T, N, N> matrix)
{
    T trace {0};
    for (std::size_t i = 0; i < N; ++i)
        trace += matrix.get(i, i);
    return trace;
}

template <typename Derived, typename T, std::size_t N, std::size_t M>
MatrixBase<Derived, T, N, M>
operator*(const MatrixBase<Derived, T, N, M>& m, const T factor)
{
    MatrixBase<Derived, T, N, M> ret {m};
    ret *= factor;
    return ret;
}

template <typename Derived, typename T, std::size_t N, std::size_t M>
MatrixBase<Derived, T, N, M>
operator*(const T factor, const MatrixBase<Derived, T, N, M>& m)
{
    MatrixBase<Derived, T, N, M> ret {m};
    ret *= factor;
    return ret;
}

template <typename Derived, typename T, std::size_t N, std::size_t M>
MatrixBase<Derived, T, N, M>
operator/(const MatrixBase<Derived, T, N, M>& m, const T factor)
{
    MatrixBase<Derived, T, N, M> ret {m};
    ret /= factor;
    return ret;
}

template <typename Derived, typename T, std::size_t N, std::size_t M>
MatrixBase<Derived, T, N, M>
operator/(const T factor, const MatrixBase<Derived, T, N, M>& m)
{
    MatrixBase<Derived, T, N, M> ret {m};
    ret /= factor;
    return ret;
}

template <typename Derived, typename T, std::size_t N, std::size_t M>
VectorClass<T, N>
operator*(const MatrixBase<Derived, T, N, M>& m, const VectorClass<T, M>& vector)
{
    VectorClass<T, N> r;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < M; ++j)
            r[i] += m(i, j) * vector[j];
    return r;
}

template <typename Derived, typename T, std::size_t N, std::size_t M>
VectorClass<T, M>
operator*(const VectorClass<T, N>& vector, const MatrixBase<Derived, T, N, M>& m)
{
    VectorClass<T, N> r;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < N; ++j)
            r[i] += m(j, i) * vector[j];
    return r;
}

using Tensor = Matrix<double, 3, 3>;
using SymTensor = SymmetricMatrix<double, 3>;
using SymmetricTensor = SymmetricMatrix<double, 3>;

} // end namespace BasicDatatypes
