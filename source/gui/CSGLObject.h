///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSGLObject.h                                                         //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-06-04 22:19:34                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef CS_GL_OBJECT_H
#define CS_GL_OBJECT_H

#include "../model/BasicDatatypes/Vector.h"
#include "../model/BasicDatatypes/Color.h"


/*!
  \brief Base class for 3D Objects

  The base class for 3D objects to be drawn on the GLDisplay.
  Defines the necessary member functions for GL drawing in Qt and color setting:\n
  void draw() \n
  void setColor(QtColor color);
*/
class CSGLObject //: public QObject
{
 public:
  CSGLObject( Vector3f * position );
  ~CSGLObject();

  /*!
    \brief The draw routine

    Reimlement to draw the actual shape with GL routines.
  */
  virtual void draw() =0;

  /*!
    \brief Set the object's color

    Only sets the member variable mColor;
  */
  void setColor(ARGBColor * color) {mpColor = color;};
  ARGBColor * getColor() const { return mpColor; };

  /*!
    \brief Query transparency setting

    This function queries mpColor's alpha channel - if it exists.  If mpColor
    is still NULL, the default of non-transparency is assumed.
  */
  bool isTransparent() const
  {
      if (mpColor)
          return (mpColor->alpha!=1.);

      return false;
  };

 protected:
  ARGBColor   * mpColor;
  Vector3f * mpPosition;
};

#endif // CS_GL_OBJECT_H
