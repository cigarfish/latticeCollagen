/**
 * @file	Vector.cpp
 * @date	09.01.2012
 * @author	Johannes Neitsch (Position.cpp)
 * @brief	3D Vector (double precision)
 *
 * extended description
 *
 */

/** short funtion discription
 * @param	paramter_name	parameter_description
 * @return	return_description
 *
 * extended description
 *
 */

#include "Vector.h" 
#include "math.h"

#include <cmath>

std::ostream& operator<<(std::ostream& stream, const Points& points)
{
	stream << "(" << points.first << ", " << points.second << ")";
	return stream;
}


/*
// Default constructor
Vector3f::Vector3f()
  : x(*data),
    y(*(data+1)),
    z(*(data+2))
{
  x = 0.;
  y = 0.;
  z = 0.;
}

// Constructor with set
Vector3f::Vector3f(double x, double y, double z)
  : x(*data),
    y(*(data+1)),
    z(*(data+2))
{
  Set(x,y,z);
}

// copy constructor
Vector3f::Vector3f(const Vector3f &otherVector)
  : x(*data),
    y(*(data+1)),
    z(*(data+2))
{
  this->Set(otherVector);
}


// Set
void Vector3f::Set(const Vector3f &otherVector)
{
  x = otherVector.x;
  y = otherVector.y;
  z = otherVector.z;
}

void Vector3f::Set(double x_, double y_, double z_) 
{
  x = x_;
  y = y_;
  z = z_;
}

void Vector3f::Set(double pos[3]) 
{
  x = pos[0];
  y = pos[1];
  z = pos[2];
}

void Vector3f::Set(double x_, int i) 
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

double Vector3f::Normalize()
{

  double norm = sqrt(x * x + y * y + z * z);
  if (norm != 0){
    x /= norm;
    y /= norm;
    z /= norm;
  }
  else
  {
	// Error: To do
  }

  return norm;
}

double Vector3f::Norm() const
{
  return sqrt(x * x + y * y + z * z);
}

void Vector3f::Multiply(double factor)
{
  x *= factor;
  y *= factor;
  z *= factor;
}

void Vector3f::Add(const double dx, const double dy, const double dz)
{
  x += dx;
  y += dy;
  z += dz;
}

void Vector3f::Add(Vector3f v)
{
  x += v.x;
  y += v.y;
  z += v.z;
}

void Vector3f::Add(double d_ , int index)
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


void Vector3f::Rotate( Vector3f & rotationAxis, double angle )
{
  double sinAngle = std::sin(angle);
  double cosAngle = std::cos(angle);

  // todo normalize rotationAxis; and don't do it in-place!

  double newX = this->x * cosAngle + sinAngle * ( rotationAxis.y * this->z - rotationAxis.z * this->y );
  double newY = this->y * cosAngle + sinAngle * ( rotationAxis.z * this->x - rotationAxis.x * this->z);
  double newZ = this->z * cosAngle + sinAngle * ( rotationAxis.x * this->y - rotationAxis.y * this->x );

  this->Set( newX, newY, newZ );
}


//
// Vector2f
//
// Default constructor
Vector2f::Vector2f()
  : x( *data ),
    y( *(data+1) )
{
  x = 0.;
  y = 0.;
}


// Constructor with set
Vector2f::Vector2f(double x, double y)
  : x( *data ),
    y( *(data+1) )
{
  Set(x,y);
}


// Copy constructor
Vector2f::Vector2f(const Vector2f &otherVector)
  : x( *data ),
    y( *(data+1) )
{
  Set(otherVector);
}


// Set
void
Vector2f::Set(const Vector2f &otherVector)
{
  x = otherVector.x;
  y = otherVector.y;
}


void Vector2f::Set(double x_, double y_) 
{
  x = x_;
  y = y_;
}

double Vector2f::Normalize()
{
  double norm = sqrt(x * x + y * y);
  if (norm != 0){
    x /= norm;
    y /= norm;

  }
  else
  {
	// Error: To do
  }
  return norm;
}

double Vector2f::Norm() const
{
  return sqrt(x * x + y * y);
}

void Vector2f::Multiply(double factor)
{
  x *= factor;
  y *= factor;
}

void Vector2f::Add(const double dx, const double dy)
{
  x += dx;
  y += dy;

}
*/
