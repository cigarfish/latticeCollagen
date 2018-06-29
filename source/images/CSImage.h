///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSImage.h                                                            //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2011-10-20 14:44:34                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef CS_IMAGE_H
#define CS_IMAGE_H


#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkImageToVTKImageFilter.h"

#include "filters/imageFilters/CSImageFilter.h"

class vtkImageData;

typedef itk::Image< itk::RGBPixel<unsigned char>, 2 >   itkImageType;
typedef itk::ImageToVTKImageFilter<itkImageType>      itkToVTKFilter;


class CSImage
{
public:

    CSImage();
    ~CSImage();

    void loadFromFile( const char * fileName );
    void loadFromFile( const std::string fileName);

    void addITKFilter( CSITKImageFilter * itkFilter );
    void addVTKFilter( CSVTKImageFilter * vtkFilter );

    vtkImageData * getImageData()
    {
        return mpVTKImage;
    };

    void setImageData( vtkImageData * newImage );
    void setImageData( itkImageType * newImage );

    void update();

private:
    itkImageType::Pointer mpITKImage;
    vtkImageData * mpVTKImage;
    itkToVTKFilter::Pointer mpConnector;

    std::vector<CSITKImageFilter *> mITKFilters;
    std::vector<CSVTKImageFilter *> mVTKFilters;
};


#endif /* CS_IMAGE_H */
