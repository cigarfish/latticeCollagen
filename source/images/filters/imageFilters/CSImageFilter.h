///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSImageFilter.h                                                      //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2011-10-21 16:42:54                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef CS_IMAGEFILTER_H
#define CS_IMAGEFILTER_H

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "vtkImageData.h"


template<class ImageTypePointer>
class CSImageFilter
{
 public:
  CSImageFilter()
    : mpInputImage(NULL),mpOutputImage(NULL) {};

  virtual ~CSImageFilter() {};

  void SetInput(ImageTypePointer indata)
  { mpInputImage = indata; }

  ImageTypePointer GetInput()
  { return mpInputImage; }

  ImageTypePointer GetOutput()
  {
    if (!mpOutputImage)
      filter();

    return mpOutputImage ? mpOutputImage : mpInputImage;
  }

 protected:
  /*!
    This is the function in which to define the filter.
  */
  virtual void filter() {};

  ImageTypePointer mpInputImage;
  ImageTypePointer mpOutputImage;
};

typedef unsigned char                         itkPixelComponentType;
typedef itk::RGBPixel<itkPixelComponentType>  itkRGBPixelType;
typedef itk::Image< itkRGBPixelType, 2 >      itkImageType;

typedef CSImageFilter<vtkImageData *> CSVTKImageFilter;
typedef CSImageFilter<itkImageType::Pointer> CSITKImageFilter;

#endif // CSIMAGEFILTER_H
