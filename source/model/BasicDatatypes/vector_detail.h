#ifndef VECTOR_DETAIL_H
#define VECTOR_DETAIL_H

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <locale>
#include <functional>
#include <cmath>
#include <numeric>
#include <utility>
#include <initializer_list>
#include <stdexcept>
#include <bitset>

#include "tools/dataIO/CSVReader.h"
#include "tools/concepts/concepts.h"
#include "tools/utils/macros.h"

using namespace std::placeholders;
using std::bitset;

template <std::size_t N>
class BoolVector : public bitset<N>
{
public:
    using bitset<N>::all;
    using bitset<N>::flip;

    BoolVector()
        : bitset<N>()
    {}

    template< class CharT, class Traits, class Alloc>
    BoolVector(const std::basic_string<CharT, Traits, Alloc>& str,
               typename std::basic_string<CharT, Traits, Alloc>::size_type pos = 0,
               typename std::basic_string<CharT, Traits, Alloc>::size_type n =
               std::basic_string<CharT, Traits, Alloc>::npos,
               CharT zero = CharT('0'), CharT one = CharT('1'))
        : bitset<N>(str, n)
    {}

    template< class CharT >
    BoolVector(const CharT * str,
               typename std::basic_string<CharT>::size_type n =
               std::basic_string<CharT>::npos,
               CharT zero = CharT('0'), CharT one = CharT('1'))
        : bitset<N>(str, n)
    {}

    explicit operator bool() const { return all(); }
};

template <std::size_t N>
bool All(const BoolVector<N>& vector)
{
    return vector.all();
}

template <std::size_t N>
bool None(const BoolVector<N>& vector)
{
    return vector.none();
}

template <std::size_t N>
bool Any(const BoolVector<N>& vector)
{
    return vector.any();
}

template <std::size_t N>
BoolVector<N> operator~(const BoolVector<N>& vector)
{
    BoolVector<N> ret {vector};
    ret.flip();
    return ret;
}

template <std::size_t N>
BoolVector<N> operator!(const BoolVector<N>& vector)
{
    BoolVector<N> ret {vector};
    ret.flip();
    return ret;
}

template <typename Derived, typename T, std::size_t N>
struct BaseVector
{
    using reverse_iterator       = typename std::reverse_iterator<T*>;
    using const_reverse_iterator = typename std::reverse_iterator<const T*>;

    constexpr BaseVector()
        : BaseVector(T(0))
    {}

    BaseVector(T (&data)[N]) // avoid decay to pointer
    {
        std::copy(std::begin(data), std::end(data), begin());
    }

    explicit BaseVector(const BaseVector& vector)
    {
        std::copy(vector.cbegin(), vector.cend(), begin());
    }

    explicit BaseVector(BaseVector&& vector)
    {
        swap(vector, *this);
    }

    explicit BaseVector(const T value)
    {
        set(value);
    }

    template <class InputIt>
    BaseVector(InputIt first, InputIt last, std::size_t pos = 0)
        : BaseVector(T(0))
    {
        auto d_first = begin() + pos;
        while (first != last && d_first != end())
            *d_first++ = *first++;

        // This checks whether more data was provided using the InputIterator if
        // so report an error!
        ASSERT(std::distance(first, last) == 0, "Input iterator with size " +
               std::to_string(std::distance(first, last) + N) +
               " is larger than the container" + " size " +
               std::to_string(N) + ".");
    }

    BaseVector(std::initializer_list<T> ilist)
        : BaseVector(T(0)) // this make sure that if values.size() < N we have initialized values
    {
        ASSERT(ilist.size() <= N, "Initializer list with size " +
               std::to_string(ilist.size()) + " is larger than the container" +
               " size " + std::to_string(N) + ".");
        std::copy(ilist.begin(), ilist.end(), begin());
    }

    template <typename U,
              typename = Enable_if<All(Convertible<U, T>())> >
    void assign(const U& value)
    {
        std::fill(begin(), end(), value);
    }

    // Note this function will only take the first N values
    // as InputIt may not be a RandomAccessIterator we cant use std::distance
    template <class InputIt>
    void assign(InputIt first, InputIt last)
    {
        auto d_first = begin();
        while (first != last && d_first != end())
            *d_first++ = *first++;

        // This checks whether more data was provided using the InputIterator if
        // so report an error!
        ASSERT(std::distance(first, last) == 0, "Input iterator with size " +
               std::to_string(std::distance(first, last) + N) +
               " is larger than the container" + " size " +
               std::to_string(N) + ".");
    }

    void assign(std::initializer_list<T> ilist)
    {
        ASSERT(ilist.size() <= N, "Initializer list with size " +
               std::to_string(ilist.size()) + " is larger than the container" +
               " size " + std::to_string(N) + ".");
        std::copy(ilist.begin(), ilist.end(), begin());
    }

    BaseVector& operator=(const BaseVector& vector)
    {
        std::copy(vector.cbegin(), vector.cend(), begin());
        return *this;
    }

    BaseVector& operator=(BaseVector&& vector)
    {
        swap(vector, *this);
        return *this;
    }

    BaseVector& operator=(const T value)
    {
        set(value);
        return *this;
    }

    T* begin() { return static_cast<Derived*>(this)->begin(); }
    const T* begin()  const { return static_cast<Derived*>(this)->begin(); }
    const T* cbegin() const { return static_cast<const Derived*>(this)->cbegin(); }

    reverse_iterator rbegin() { return static_cast<Derived*>(this)->rbegin(); }
    const_reverse_iterator crbegin() const { return static_cast<const Derived*>(this)->crbegin(); }

    T* end() { return static_cast<Derived*>(this)->end(); }
    const T* end()  const { return static_cast<Derived*>(this)->end(); }
    const T* cend() const { return static_cast<const Derived*>(this)->cend(); }

    reverse_iterator rend() { return static_cast<Derived*>(this)->rend(); }
    const_reverse_iterator crend() const { return static_cast<const Derived*>(this)->crend(); }

    T* data() { return begin(); }
    const T* data() const { return cbegin(); }

    constexpr std::size_t size() const
    {
        return static_cast<const Derived*>(this)->size();
    }

    std::size_t memory_size() const
    {
        return static_cast<const Derived*>(this)->memory_size();
    }

    // Access operator
    T& operator[](const std::size_t i) { return *(begin() + i); }
    const T& operator[](const std::size_t i) const { return *(cbegin() + i); }

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

    BaseVector& operator+=(const BaseVector& rhs)
    {
        return operator_base(rhs, std::plus<T>());
    }

    BaseVector& operator-=(const BaseVector& rhs)
    {
        return operator_base(rhs, std::minus<T>());
    }

    BaseVector& operator*=(const T& scalar)
    {
        return operator_scalar_base(scalar, std::multiplies<T>());
    }

    BaseVector& operator/=(const T& scalar)
    {
        return operator_scalar_base(scalar, std::divides<T>());
    }

    BaseVector& operator*=(const BaseVector& rhs)
    {
        return operator_base(rhs, std::multiplies<T>());
    }

    BaseVector& operator/=(const BaseVector& rhs)
    {
        return operator_base(rhs, std::divides<T>());
    }

    void swap(BaseVector& other)
    {
        std::swap_ranges(begin(), end(), other.begin());
    }

    void set(const T value)
    {
        std::fill(begin(), end(), value);
    }
	
    constexpr void Reset()
    {
        set(T(0));
    }

	auto Norm2() const -> decltype(std::sqrt(std::inner_product(cbegin(), cend(), cbegin(), T(0))))
	{
		return std::sqrt(std::inner_product(cbegin(), cend(), cbegin(), T(0)));
	}

    auto Norm() const -> decltype(Norm2())
    {
        return Norm2();
    }

    T Normalize()
    {
        auto norm {Norm2()};
        if (!std::isnormal(norm))
            throw std::overflow_error("Divide by zero exception!");

        *this /= norm;
        return norm;
    }

private:

    template <class Function>
    BaseVector& operator_base(const BaseVector& rhs, Function f)
    {
        std::transform(begin(), end(), rhs.cbegin(), begin(), f);
        return *this;
    }

    template <class Function>
    BaseVector& operator_scalar_base(const T scalar, Function f)
    {
        std::transform(begin(), end(), begin(), std::bind2nd(f, scalar));
        return *this;
    }

};

template <typename T, std::size_t N>
struct VectorClass : public BaseVector<VectorClass<T, N>, T, N>
{
    using reverse_iterator       = typename std::reverse_iterator<T*>;
    using const_reverse_iterator = typename std::reverse_iterator<const T*>;

    // constructors
    constexpr VectorClass()
        : BaseVector<VectorClass<T, N>, T, N>()
    {}

    constexpr explicit VectorClass(T (&data)[N])
        : BaseVector<VectorClass<T, N>, T, N>(data)
    {}

    explicit VectorClass(const VectorClass& vector)
        : BaseVector<VectorClass<T, N>, T, N>(vector)
    {}

    explicit VectorClass(VectorClass&& vector)
        : BaseVector<VectorClass<T, N>, T, N>(vector)
    {}

    explicit VectorClass(const T value)
        : BaseVector<VectorClass<T, N>, T, N>(value)
    {}

    VectorClass(std::initializer_list<T> values)
        : BaseVector<VectorClass<T, N>, T, N>(values)
    {}

    template <class InputIt>
    VectorClass(InputIt first, InputIt last, std::size_t pos = 0)
        : BaseVector<VectorClass<T, N>, T, N>(first, last, pos)
    {}

    VectorClass& operator=(const VectorClass& vector)
    {
        BaseVector<VectorClass<T, N>, T, N>::operator = (vector);
        return *this;
    }

    VectorClass& operator=(VectorClass&& vector)
    {
        BaseVector<VectorClass<T, N>, T, N>::operator = (vector);
        return *this;
    }

    VectorClass& operator=(const T value)
    {
        BaseVector<VectorClass<T, N>, T, N>::operator = (value);
        return *this;
    }

    ~VectorClass() {}

    // Iterator access
    T* begin() { return std::begin(data); }
    const T* begin()  const { return std::begin(data); }
    const T* cbegin() const { return std::cbegin(data); }

    reverse_iterator rbegin() { return std::rbegin(data); }
    const_reverse_iterator crbegin() const { return std::crbegin(data); }

    T* end() { return std::end(data); }
    const T* end()  const { return std::end(data); }
    const T* cend() const { return std::cend(data); }

    reverse_iterator rend() { return std::rend(data); }
    const_reverse_iterator crend() const { return std::crend(data); }

    constexpr std::size_t size() const { return N; }
    std::size_t memory_size() const { return sizeof(*this); }

    T data[N];
};

template<typename T>
struct VectorClass<T, 2> : public BaseVector<VectorClass<T, 2>, T, 2>
{
    using reverse_iterator       = typename std::reverse_iterator<T*>;
    using const_reverse_iterator = typename std::reverse_iterator<const T*>;

    // constructors
    VectorClass()
        : BaseVector<VectorClass<T, 2>, T, 2>()
    {}

    VectorClass(T (&data)[2])
        : BaseVector<VectorClass<T, 2>, T, 2>(data)
    {}

    VectorClass(const VectorClass& vector)
        : BaseVector<VectorClass<T, 2>, T, 2>(vector)
    {}

    VectorClass(VectorClass&& vector)
        : BaseVector<VectorClass<T, 2>, T, 2>(vector)
    {}

    explicit VectorClass(const T value)
        : BaseVector<VectorClass<T, 2>, T, 2>(value)
    {}

    VectorClass(std::initializer_list<T> values)
        : BaseVector<VectorClass<T, 2>, T, 2>(values)
    {}

    VectorClass& operator=(const VectorClass& vector)
    {
        BaseVector<VectorClass<T, 2>, T, 2>::operator = (vector);
        return *this;
    }

    VectorClass& operator=(VectorClass&& vector)
    {
        BaseVector<VectorClass<T, 2>, T, 2>::operator = (vector);
        return *this;
    }

    VectorClass& operator=(const T value)
    {
        BaseVector<VectorClass<T, 2>, T, 2>::operator = (value);
        return *this;
    }

    VectorClass(T _x, T _y)
        : x(_x), y(_y)
    {}

    ~VectorClass() {}

    union {
        T data[2];
        struct { T x, y; };
    };

    T* begin() { return std::begin(data); }
    const T* begin()  const { return std::begin(data); }
    const T* cbegin() const { return std::cbegin(data); }

    reverse_iterator rbegin() { return std::rbegin(data); }
    const_reverse_iterator crbegin() const { return std::crbegin(data); }

    T* end() { return std::end(data); }
    const T* end()  const { return std::end(data); }
    const T* cend() const { return std::cend(data); }

    reverse_iterator rend() { return std::rend(data); }
    const_reverse_iterator crend() const { return std::crend(data); }

    constexpr std::size_t size() const { return 2; }
    std::size_t memory_size() const { return sizeof(*this); }
};

template<typename T>
struct VectorClass<T, 3> : public BaseVector<VectorClass<T, 3>, T, 3>
{
    using reverse_iterator       = typename std::reverse_iterator<T*>;
    using const_reverse_iterator = typename std::reverse_iterator<const T*>;

    // constructors
    //constexpr VectorClass()
	VectorClass()
        : BaseVector<VectorClass<T, 3>, T, 3>()
    {}

    explicit VectorClass(T (&data)[3])
        : BaseVector<VectorClass<T, 3>, T, 3>(data)
    {}

    VectorClass(const VectorClass& vector)
        : BaseVector<VectorClass<T, 3>, T, 3>(vector)
    {}

    VectorClass(VectorClass&& vector)
        : BaseVector<VectorClass<T, 3>, T, 3>(vector)
    {}

    explicit VectorClass(const T value)
        : BaseVector<VectorClass<T, 3>, T, 3>(value)
    {}

    VectorClass(std::initializer_list<T> values)
        : BaseVector<VectorClass<T, 3>, T, 3>(values)
    {}

    template <class InputIt>
    VectorClass(InputIt first, InputIt last)
        : BaseVector<VectorClass<T, 3>, T, 3>(first, last)
    {}

    VectorClass& operator=(const VectorClass& vector)
    {
        BaseVector<VectorClass<T, 3>, T, 3>::operator = (vector);
        return *this;
    }

    VectorClass& operator=(VectorClass&& vector)
    {
        BaseVector<VectorClass<T, 3>, T, 3>::operator = (vector);
        return *this;
    }

    VectorClass& operator=(const T value)
    {
        BaseVector<VectorClass<T, 3>, T, 3>::operator = (value);
        return *this;
    }

    VectorClass(T _x, T _y, T _z)
        : x(_x), y(_y), z(_z)
    {}

    void DEPRECATED(Multiply(const T value))
    {
        *this *= value;
    }

    void DEPRECATED(Set(T x, T y, T z))
    {
        *this = VectorClass(x, y, z);
    }

    void DEPRECATED(Set(const VectorClass& vector))
    {
        *this = VectorClass(vector);
    }

    void DEPRECATED(Add(T x, T y, T z))
    {
        *this += VectorClass(x, y, z);
    }

    void DEPRECATED(Add(const VectorClass& vector))
    {
        *this += vector;
    }

    void DEPRECATED(Add(T d_ , int index))
    {
        switch( index )
        {
            case 0:
                x += d_;
            break;
            case 1:
                y += d_;
            break;
            case 2:
                z += d_;
            break;
        }
    }

    void DEPRECATED(Set(T x_, int i))
    {
        switch( i)
        {
            case 0:
                x = x_;
            break;
            case 1:
                y = x_;
            break;
            case 2:
                z = x_;
             break;
        }
    }

    ~VectorClass() {}

    union {
        T data[3];
        struct { T x, y, z; };
        struct { T r, g, b; };
        struct { T red, green, blue; };
        VectorClass<T, 2> xy;
    };

    T* begin() { return std::begin(data); }
    const T* begin()  const { return std::begin(data); }
    const T* cbegin() const { return std::cbegin(data); }

    reverse_iterator rbegin() { return std::rbegin(data); }
    const_reverse_iterator crbegin() const { return std::crbegin(data); }

    T* end() { return std::end(data); }
    const T* end()  const { return std::end(data); }
    const T* cend() const { return std::cend(data); }

    reverse_iterator rend() { return std::rend(data); }
    const_reverse_iterator crend() const { return std::crend(data); }

    constexpr std::size_t size() const { return 3; }
    std::size_t memory_size() const { return sizeof(*this); }
};

template<typename T>
struct VectorClass<T, 4> : public BaseVector<VectorClass<T, 4>, T, 4>
{
    using reverse_iterator       = typename std::reverse_iterator<T*>;
    using const_reverse_iterator = typename std::reverse_iterator<const T*>;

    // constructors
    constexpr VectorClass()
        : BaseVector<VectorClass<T, 4>, T, 4>()
    {}

    explicit VectorClass(T (&data)[4])
        : BaseVector<VectorClass<T, 4>, T, 4>(data)
    {}

    VectorClass(const VectorClass& vector)
        : BaseVector<VectorClass<T, 4>, T, 4>(vector)
    {}

    VectorClass(VectorClass&& vector)
        : BaseVector<VectorClass<T, 4>, T, 4>(vector)
    {}

    explicit VectorClass(const T value)
        : BaseVector<VectorClass<T, 4>, T, 4>(value)
    {}

    VectorClass(std::initializer_list<T> values)
        : BaseVector<VectorClass<T, 4>, T, 4>(values)
    {}

    template <class InputIt>
    VectorClass(InputIt first, InputIt last, std::size_t pos = 0)
        : BaseVector<VectorClass<T, 4>, T, 4>(first, last, pos)
    {}

    VectorClass& operator=(const VectorClass& vector)
    {
        BaseVector<VectorClass<T, 4>, T, 4>::operator = (vector);
        return *this;
    }

    VectorClass& operator=(VectorClass&& vector)
    {
        BaseVector<VectorClass<T, 4>, T, 4>::operator = (vector);
        return *this;
    }

    VectorClass& operator=(const T value)
    {
        BaseVector<VectorClass<T, 4>, T, 4>::operator = (value);
        return *this;
    }

    VectorClass(T _x, T _y, T _z, T _w)
        : x(_x), y(_y), z(_z), w(_w)
    {}

    ~VectorClass() {}

    union {
        T data[4];
        struct { T x, y, z, w; };
        struct { T r, g, b, a; };
        struct { T red, green, blue, alpha; };
        VectorClass<T, 2> xy;
        VectorClass<T, 3> xyz;
        VectorClass<T, 3> rgb;
    };

    T* begin() { return std::begin(data); }
    const T* begin()  const { return std::begin(data); }
    const T* cbegin() const { return std::cbegin(data); }

    reverse_iterator rbegin() { return std::rbegin(data); }
    const_reverse_iterator crbegin() const { return std::crbegin(data); }

    T* end() { return std::end(data); }
    const T* end()  const { return std::end(data); }
    const T* cend() const { return std::cend(data); }

    reverse_iterator rend() { return std::rend(data); }
    const_reverse_iterator crend() const { return std::crend(data); }

    constexpr std::size_t size() const { return 4; }
    std::size_t memory_size() const { return sizeof(*this); }
};

// swap
template <typename Derived, typename T, std::size_t N>
void swap(BaseVector<Derived, T, N>& lhs, BaseVector<Derived, T, N>& rhs)
{
    lhs.swap(rhs);
}

// Numerical methods of vectors
template <typename T, std::size_t N>
auto Norm2(const VectorClass<T, N>& vector)
{
    return std::sqrt(std::inner_product(vector.cbegin(), vector.cend(),
                                        vector.cbegin(), T(0)));
}

template <typename T, std::size_t N>
auto Norm(const VectorClass<T, N>& vector)
{
    return Norm2(vector);
}

template <typename T, std::size_t N>
auto Norm2Squared(const VectorClass<T, N>& vector)
{
    return std::inner_product(vector.cbegin(), vector.cend(),
                              vector.cbegin(), T(0));
}

template <typename T, std::size_t N>
auto Normalize(const VectorClass<T, N>& vector)
{
    VectorClass<T, N> ret {vector};

    auto norm {Norm(vector)};
    if (!std::isnormal(norm))
        throw std::overflow_error("Divide by zero exception!");

    ret /= norm;
    return ret;
}

template <typename T, std::size_t N>
auto UnitVector(const VectorClass<T, N>& vector)
{
    return Normalize(vector);
}

// If condition true select lhs, otherwise rhs
template <typename T, std::size_t N>
VectorClass<T, N> select(const BoolVector<N>& condition,
                         const VectorClass<T, N>& lhs,
                         const VectorClass<T, N>& rhs)
{
    VectorClass<T, N> ret;
    for (std::size_t i = 0; i < N; ++i)
        ret[i] = (condition[i]) ? lhs[i] : rhs[i];
    return ret;
}

template <typename T,
          typename = Enable_if<Floating_Point<T>()> >
VectorClass<T, 2> UnitVector(const T theta)
{
    return {cos(theta), sin(theta)};
}

template <typename T,
          typename = Enable_if<Floating_Point<T>()> >
VectorClass<T, 3> UnitVector(const T theta, const T phi)
{
    double t = clip(theta, 0., 2. * M_PI);
    double p = clip(phi,   0., M_PI);
    return {cos(t) * sin(p), sin(t) * sin(p), cos(p)};
}

template <typename T, std::size_t N>
auto Dot(const VectorClass<T, N>& lhs, const VectorClass<T, N>& rhs)
{
    return std::inner_product(lhs.cbegin(), lhs.cend(), rhs.cbegin(), T(0));
}

template <typename T, typename Transform, std::size_t N>
VectorClass<T, N> apply(const Transform& transform,
                        const VectorClass<T, N>& vector)
{
    return transform(vector);
}

template <typename T>
T Cross(const VectorClass<T, 2>& lhs, const VectorClass<T, 2>& rhs)
{
    return lhs.x * rhs.y - lhs.y * rhs.x;
}

template <typename T>
VectorClass<T, 3> Cross(const VectorClass<T, 3>& lhs,
                        const VectorClass<T, 3>& rhs)
{
    return {lhs.y * rhs.z - lhs.z * rhs.y,
            lhs.z * rhs.x - lhs.x * rhs.z,
            lhs.x * rhs.y - lhs.y * rhs.x};
}

// mathematical operators
template <typename T, std::size_t N>
VectorClass<T, N> operator+(const VectorClass<T, N>& lhs,
                            const VectorClass<T, N>& rhs)
{
    VectorClass<T, N> ret {lhs};
    ret += rhs;
    return ret;
}

template <typename T, std::size_t N>
VectorClass<T, N> operator-(const VectorClass<T, N>& lhs,
                            const VectorClass<T, N>& rhs)
{
    VectorClass<T, N> ret {lhs};
    ret -= rhs;
    return ret;
}

template <typename T, std::size_t N>
VectorClass<T, N> operator-(const VectorClass<T, N>& vector)
{
    VectorClass<T, N> ret;
    std::transform(vector.cbegin(), vector.cend(), ret.begin(), std::negate<T>());
    return ret;
}

template <typename T, std::size_t N>
VectorClass<T, N> operator*(const VectorClass<T, N>& lhs,
                            const VectorClass<T, N>& rhs)
{
    VectorClass<T, N> ret {lhs};
    ret *= rhs;
    return ret;
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
VectorClass<T, N> operator*(const VectorClass<T, N>& lhs, const U& rhs)
{
    VectorClass<T, N> ret {lhs};
    ret *= rhs;
    return ret;
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
VectorClass<T, N> operator*(const U& lhs, const VectorClass<T, N>& rhs)
{
    return rhs * lhs;
}

template <typename T, std::size_t N>
VectorClass<T, N> operator/(const VectorClass<T, N>& lhs,
                            const VectorClass<T, N>& rhs)
{
    VectorClass<T, N> ret {lhs};
    ret /= rhs;
    return ret;
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
VectorClass<T, N> operator/(const VectorClass<T, N>& lhs, const U& rhs)
{
    VectorClass<T, N> ret {lhs};
    ret /= rhs;
    return ret;
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
VectorClass<T, N> operator/(const U& lhs, const VectorClass<T, N>& rhs)
{
    VectorClass<T, N> ret (lhs);
    ret /= rhs;
    return ret;
}

// other math functions
template <typename T, std::size_t N>
VectorClass<T, N> abs(const VectorClass<T, N>& vector)
{
    VectorClass<T, N> ret;
    std::transform(vector.cbegin(), vector.cend(), ret.begin(),
                   static_cast<T (*)(T)>(&std::abs));
    return ret;
}

template <typename T, std::size_t N>
auto minComponent(const VectorClass<T, N>& vector)
{
    auto res = std::min_element(vector.cbegin(), vector.cend());
    return std::distance(std::cbegin(vector), res);
}

/*template <typename T, std::size_t N>
auto min(const VectorClass<T, N>& vector)
{
    return *std::min_element(vector.cbegin(), vector.cend());
}*/

/*template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
auto min(const VectorClass<T, N>& vector, const U value)
{
    VectorClass<T, N> ret;
    // min has too many overloads the compiler struggles
    T const& (*min) (T const&, T const&) = std::min<T>;
    std::transform(vector.cbegin(), vector.cend(), ret.begin(),
                   std::bind(min, _1, T(value)));
    return ret;
}*/

template <typename T, std::size_t N>
auto maxComponent(const VectorClass<T, N>& vector)
{
    auto res = std::max_element(vector.cbegin(), vector.cend());
    return std::distance(std::cbegin(vector), res);
}

/*template <typename T, std::size_t N>
auto max(const VectorClass<T, N>& vector)
{
    return *std::max_element(vector.cbegin(), vector.cend());
}*/

/*template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
auto max(const VectorClass<T, N>& vector, const U value)
{
    VectorClass<T, N> ret;
    // min has too many overloads the compiler struggles
    T const& (*max) (T const&, T const&) = std::max<T>;
    std::transform(vector.cbegin(), vector.cend(), ret.begin(),
                   std::bind(max, _1, T(value)));
    return ret;
}*/

// Computes (x * y) + z
template <typename T, std::size_t N,
          typename = Enable_if<Floating_Point<T>()> >
auto fma(const VectorClass<T, N>& x, const VectorClass<T, N>& y,
         const VectorClass<T, N>& z)
{
    // can we make this a AVX?
    VectorClass<T, N> ret;
    for (std::size_t i = 0; i < N; i++)
        ret[i] = std::fma(x[i], y[i], z[i]);
    return ret;
}

template <typename T, std::size_t N,
          typename = Enable_if<Floating_Point<T>()> >
bool isfinite(const VectorClass<T, N>& vector)
{
    bool (*isfinite) (const T) = std::isfinite;
    return std::all_of(vector.cbegin(), vector.cend(), std::bind(isfinite, _1));
}

/*template <typename T, std::size_t N,
          typename = Enable_if<Floating_Point<T>()> >
bool isnan(const VectorClass<T, N>& vector)
{
    bool (*isnan) (const T) = std::isnan;
    return std::all_of(vector.cbegin(), vector.cend(), std::bind(isnan, _1));
}*/

/*template <typename T, std::size_t N,
          typename = Enable_if<Floating_Point<T>()> >
bool isinf(const VectorClass<T, N>& vector)
{
    bool (*isinf) (const T) = std::isinf;
    return std::all_of(vector.cbegin(), vector.cend(), std::bind(isinf, _1));
}*/

/*template <typename T, std::size_t N,
          typename = Enable_if<Floating_Point<T>()> >
bool isnormal(const VectorClass<T, N>& vector)
{
    bool (*isnormal) (const T) = std::isnormal;
    return std::all_of(vector.cbegin(), vector.cend(), std::bind(isnormal, _1));
}*/

template <typename T>
T clip(const T& n, const T& lower_bound, const T& upper_bound)
{
    return std::max(lower_bound, std::min(n, upper_bound));
}

template <typename T, std::size_t N>
VectorClass<T, N> lerp(const VectorClass<T, N>& start,
                       const VectorClass<T, N> end, double t)
{
    return start + (end - start) * clip(t, 0., 1.);
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
auto clamp(const VectorClass<T, N>& vector,
           const U lower_bound, const U upper_bound)
{
    VectorClass<T, N> ret;
    std::transform(vector.cbegin(), vector.cend(), ret.begin(),
                   std::bind(clip<T>, _1, T(lower_bound), T(upper_bound)));
    return ret;
}

template <typename T, std::size_t N>
auto saturate(const VectorClass<T, N>& vector)
{
    return clamp(vector, T(0), T(1));
}

// mathematical comparisions
template <typename T, std::size_t N>
BoolVector<N> operator==(const VectorClass<T, N>& vec1,
                         const VectorClass<T, N>& vec2)
{
    BoolVector<N> ret;
    for (std::size_t i = 0; i < N; ++i)
        ret[i] = (vec1[i] == vec2[i]);

    return ret;
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
bool operator==(const VectorClass<T, N>& vector, const U value)
{
    return all_of(vector.cbegin(), vector.cend(),
                  std::bind1st(std::equal_to<T>(), value));
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
bool operator==(const U value, const VectorClass<T, N>& vector)
{
    return (vector == value);
}

template <typename T, std::size_t N>
BoolVector<N> operator!=(const VectorClass<T, N>& vec1,
                         const VectorClass<T, N>& vec2)
{
    return !(vec1 == vec2);
}

// should return false if a single value does not equal value
template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
bool operator!=(const VectorClass<T, N>& vector, const U value)
{
    return any_of(vector.cbegin(), vector.cend(),
                  std::bind1st(std::not_equal_to<T>(), value));
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
bool operator!=(const U value, const VectorClass<T, N>& vector)
{
    return (vector != value);
}

template <typename T, std::size_t N>
BoolVector<N> operator<(const VectorClass<T, N>& vec1,
                        const VectorClass<T, N>& vec2)
{
    BoolVector<N> ret;
    for (std::size_t i = 0; i < N; ++i)
        ret[i] = (vec1[i] < vec2[i]);

    return ret;
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
bool operator<(const VectorClass<T, N>& vector, const U value)
{
    return all_of(vector.cbegin(), vector.cend(),
                  std::bind2nd(std::less<T>(), T(value)));
}

template <typename T, std::size_t N>
BoolVector<N> operator<=(const VectorClass<T, N>& vec1,
                         const VectorClass<T, N>& vec2)
{
    return !(vec1 > vec2);
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
bool operator<=(const VectorClass<T, N>& vector, const U value)
{
    return !(vector > T(value));
}

template <typename T, std::size_t N>
BoolVector<N> operator>(const VectorClass<T, N>& vec1,
                        const VectorClass<T, N>& vec2)
{
    return vec2 < vec1;
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
bool operator>(const VectorClass<T, N>& vector, const U value)
{
    return all_of(vector.cbegin(), vector.cend(),
                  std::bind2nd(std::greater<T>(), T(value)));
}

template <typename T, std::size_t N>
BoolVector<N> operator>=(const VectorClass<T, N>& vec1,
                         const VectorClass<T, N>& vec2)
{
    return !(vec1 < vec2);
}

template <typename T, typename U, std::size_t N,
          typename = Enable_if<All(Convertible<U, T>())> >
bool operator>=(const VectorClass<T, N>& vector, const U value)
{
    return !(vector < T(value));
}

// computes the determinant of matrix constructed of vectors a, b, c as rows
template <typename T>
auto determinant(const VectorClass<T, 3>& a, const VectorClass<T, 3>& b,
                 const VectorClass<T, 3>& c)
{
    return a.x * (b.y * c.z - c.y * b.z) -
           a.y * (b.x * c.z - c.x * b.z) +
           a.z * (b.x * c.y - b.y * c.x);
}

template <typename T>
auto Angle(const VectorClass<T, 2>& a, const VectorClass<T, 2>& b)
{
    auto dot = Dot(a, b);
    auto crs = Cross(a, b);
    return atan2(crs, dot);
}

template <typename T>
auto Angle(const VectorClass<T, 3>& a, const VectorClass<T, 3>& b)
{
    auto dot = Dot(a, b);
    auto crs = Norm(Cross(a, b));
    return atan2(crs, dot);
}

template <typename T>
auto Angle2(const VectorClass<T, 3>& a, const VectorClass<T, 3>& b)
{
    auto n1 = Norm(a);
    auto n2 = Norm(b);
    auto vec1 = n2 * a - n1 * b;
    auto vec2 = n2 * a + n1 * b;
    return 2. * atan2(Norm(vec1), Norm(vec2));
}

// input and output
template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& stream,
                         const BaseVector<VectorClass<T, N>, T, N>& vec)
{
    std::ostringstream ss;

    ss << "(";
    for (auto p = vec.cbegin(); p!=vec.cend() - 1; ++p)
        ss << *p << ", ";
    ss << *(vec.cend() - 1) << ")";

    stream << ss.str();

    return stream;
}

template<typename T, std::size_t N>
std::istream& operator>>(std::istream& stream,
                         BaseVector<VectorClass<T, N>, T, N>& vec)
{
    stream.imbue(std::locale(std::locale(), new csv_reader()));
    vec.assign(std::istream_iterator<T>(stream), std::istream_iterator<T>());
    return stream;
}

#endif
