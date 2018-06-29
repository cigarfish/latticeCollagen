
#ifndef VECTOR_H_
#define VECTOR_H_

//#include <cstdlib>

// Import Basic Vector library
#include "vector_detail.h"
#include "Matrix.h"
#include "Transform.h"

// Aligned Data types required for the solvers
#include "AlignedVector.h"

template <typename T, typename = Enable_if<Floating_Point<T>()> >
auto Rotate(const VectorClass<T, 3>& vector, VectorClass<T, 3> axis, const T angle)
{
	// TODO in this case we dont have to assemble the inverse
	auto rotate = BasicDatatypes::Transform<T>::rotate(axis, angle);
	return apply(rotate, vector);
}

template <typename T, typename = Enable_if<Floating_Point<T>()> >
auto Rotate_x(const VectorClass<T, 3>& vector, const T angle)
{
	auto rotate = BasicDatatypes::Transform<T>::rotate_x(angle);
	return apply(rotate, vector);
}

template <typename T, typename = Enable_if<Floating_Point<T>()> >
auto Rotate_y(const VectorClass<T, 3>& vector, const T angle)
{
	auto rotate = BasicDatatypes::Transform<T>::rotate_y(angle);
	return apply(rotate, vector);
}

template <typename T, typename = Enable_if<Floating_Point<T>()> >
auto Rotate_z(const VectorClass<T, 3>& vector, const T angle)
{
	auto rotate = BasicDatatypes::Transform<T>::rotate_z(angle);
	return apply(rotate, vector);
}

// useful shortcuts
using Vector2 = VectorClass<double, 2>;
using Vector3 = VectorClass<double, 3>;
using Vector4 = VectorClass<double, 4>;

using Vector2d = VectorClass<double, 2>;
using Vector3d = VectorClass<double, 3>;
using Vector4d = VectorClass<double, 4>;

// uncomment when old vector class can be removed
// this should be float
using Vector2f = VectorClass<double, 2>;
using Vector3f = VectorClass<double, 3>;
using Vector4f = VectorClass<double, 4>;

using Vector2i = VectorClass<int, 2>;
using Vector3i = VectorClass<int, 3>;
using Vector4i = VectorClass<int, 4>;

using Vector2l = VectorClass<long, 2>;
using Vector3l = VectorClass<long, 3>;
using Vector4l = VectorClass<long, 4>;

using Vector2u = VectorClass<unsigned int, 2>;
using Vector3u = VectorClass<unsigned int, 3>;
using Vector4u = VectorClass<unsigned int, 4>;

using Vector2ul = VectorClass<unsigned long, 2>;
using Vector3ul = VectorClass<unsigned long, 3>;
using Vector4ul = VectorClass<unsigned long, 4>;

using Points = std::pair<Vector3d, Vector3d>;

/*
struct Vector2f 
{
public:

  double data[2];

  double &x, &y;

  // add an assert for index<2 in here!
  double & operator[](const size_t index) { return data[index]; };

  Vector2f();
  Vector2f(double x, double y);
  Vector2f(const Vector2f&);

  void Set(const Vector2f &);
  void Set(double x, double y);

  void Add(double dx, double dy);

  //! Scale the vector to size 1.
  //! \return Norm of the vector before the normalization.
  double Normalize();
  double Norm() const;

  void Multiply(double factor);
};

struct Vector3f 
{
public:

  double data[3];

  double &x,&y,&z;

  // add an assert for index<3 in here!
  double & operator[](const size_t index) { return data[index]; };

  Vector3f();
  Vector3f(double x, double y, double z);
  Vector3f(const Vector3f&);

  void Set(const Vector3f &);
  void Set(double x, double y, double z);
  void Set(double pos[3]);
  void Set(double x, int i);
  inline Vector3f & operator=( const Vector3f & otherVector );

  void Add(double dx, double dy, double dz);
  void Add(Vector3f v);
  void Add(double d_ , int index);

  //! Another version of the Add(Vector3f) method;
  inline Vector3f & operator+=(const Vector3f &otherVector);

  inline Vector3f   operator +(const Vector3f &otherVector) const;

  //! Return a reversed version of this Vector3f
  inline Vector3f   operator -() const;

  inline Vector3f   operator -(const Vector3f &otherVector) const;

  inline Vector3f & operator-=(const Vector3f &otherVector);

  //! Scale the vector to size 1.
  //! \return Norm of the vector before the normalization.
  double Normalize();
  double Norm() const;
  void Multiply(double factor);

  //! Another version of the Multiply() method.
  inline Vector3f & operator*=(double factor);

  inline Vector3f & operator/=(double factor);

  inline Vector3f   operator *(double factor) const;

  //! Rotate the vector around an axis at a certain angle
  //  leaving the vector norm constant.
  //
  //    \param rotationAxis A normalized Vector3f.
  //    \param angle The angle in rads to rotate the vector (right-hand rule).
  //
  void Rotate( Vector3f & rotationAxis, double angle );

};

inline
Vector3f & Vector3f::operator=( const Vector3f & otherVector )
{
  this->Set(otherVector);
  return *this;
}

inline
Vector3f & Vector3f::operator+=(const Vector3f &otherVector)
{
  this->x += otherVector.x;
  this->y += otherVector.y;
  this->z += otherVector.z;
  return *this;
}

inline
Vector3f Vector3f::operator +(const Vector3f &otherVector) const
{
  Vector3f tmp(*this);
  return tmp+=otherVector;
}

inline
Vector3f Vector3f::operator-() const
{
  return Vector3f( -this->x, -this->y, -this->z );
}

inline
Vector3f Vector3f::operator-( const Vector3f & otherVector ) const
{
  return (-otherVector) + *this;
}

inline
Vector3f & Vector3f::operator-=( const Vector3f & otherVector )
{
  this->x -= otherVector.x;
  this->y -= otherVector.y;
  this->z -= otherVector.z;
  return *this;
}

inline
Vector3f & Vector3f::operator*=(double factor)
{
  this->x *= factor;
  this->y *= factor;
  this->z *= factor;
  return *this;
}

inline
Vector3f & Vector3f::operator/=(double factor)
{
  this->x /= factor;
  this->y /= factor;
  this->z /= factor;
  return *this;
}

inline
Vector3f Vector3f::operator *(double factor) const
{
  Vector3f tmp(*this);
  return tmp*=factor;
}

inline
Vector3f operator*( double factor, const Vector3f & vector )
{
  return vector*factor;
}
*/
#endif 

/* VECTOR_H_ */
