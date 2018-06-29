///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  exampleFilter.h                                                      //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2011-10-27 18:39:09                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef EXAMPLE_ITK_FILTER_H
#define EXAMPLE_ITK_FILTER_H

#include "CSImageFilter.h"

#include "itkRGBPixel.h"

#include "itkImageAdaptor.h"
#include "itkNthElementImageAdaptor.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBinaryThresholdImageFilter.h"
//#include "itkImageToVTKImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

class ExampleITKFilter : public CSITKImageFilter
{
public:
    ExampleITKFilter();
    ~ExampleITKFilter();

    void Update();

private:
    void filter();

    //----------DEFINE-ALL-THE-DEPLOYED-TYPES-----------------------------------------------------------------------------------
    typedef unsigned char                       PixelComponentType;
    typedef itk::RGBPixel<PixelComponentType>   RGBPixelType;

    typedef itk::Image<RGBPixelType, 2>         RGBImageType;
    typedef itk::Image<PixelComponentType, 2>       ScalarImageType;

    typedef itk::NthElementImageAdaptor<RGBImageType, PixelComponentType> RGBChannelImageAdaptorType;

    typedef itk::ImageFileReader<RGBImageType>      RGBReaderType;
    typedef itk::ImageFileReader<ScalarImageType>   ScalarReaderType;
    typedef itk::ImageFileWriter<RGBImageType>      RGBWriterType;
    typedef itk::ImageFileWriter<ScalarImageType>   ScalarWriterType;

    typedef itk::BinaryThresholdImageFilter<ScalarImageType, ScalarImageType>               ThresholdFilterType;
    //  typedef itk::ImageToVTKImageFilter<RGBImageType>                                            ITKVTKRGBConnectorType;
    typedef itk::RescaleIntensityImageFilter<RGBChannelImageAdaptorType, ScalarImageType>   RescalIntensityFilterType;
    //-------------------------------------------------------------------------------------------------------------------------

    //----------DEFINE-ALL-THE-DEPLOYED-READER--WRITER--FILTER----------------------------------------------------------------
    RGBReaderType::Pointer mp_rgbReader;
    RGBChannelImageAdaptorType::Pointer mp_redChannel;
    RGBChannelImageAdaptorType::Pointer mp_greenChannel;
    RGBChannelImageAdaptorType::Pointer mp_blueChannel;
    RescalIntensityFilterType::Pointer mp_rescaleFilterRed;
    RescalIntensityFilterType::Pointer mp_rescaleFilterGreen;
    RescalIntensityFilterType::Pointer mp_rescaleFilterBlue;
    ThresholdFilterType::Pointer mp_thresFilterRed;
    ThresholdFilterType::Pointer mp_thresFilterGreen;
    ThresholdFilterType::Pointer mp_thresFilterBlue;
    ScalarWriterType::Pointer mp_scalarWriter;
};

#endif //EXAMPLE_ITK_FILTER.H
