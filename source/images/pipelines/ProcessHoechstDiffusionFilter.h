///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ProcessHoechstDiffusionFilter.h                                      //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-12-06                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef PROCESSHOECHSTDIFFUSIONFILTER_H_
#define PROCESSHOECHSTDIFFUSIONFILTER_H_

#include "BasePipeline.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBluePixelAccessor.h"
#include "itkImageAdaptor.h"
#include "itkImageToImageFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkRedPixelAccessor.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"


class ProcessHoechstDiffusionFilter : public BasePipeline<3>
{
    typedef itk::RGBPixel<CScalarPixelType>     CRGBPixelType;
    typedef itk::Image<unsigned char, 2>        CScalarImageType;
    typedef itk::Image<short, 2>                SScalarImageType;
    typedef itk::Image<float, 2>                FScalarImageType;
    typedef itk::Image<CRGBPixelType, 2>        CRGBImageType;


    class BlueChannelPixelAccessor
    {
    public:
        typedef CRGBPixelType       InternalType;
        typedef CScalarPixelType    ExternalType;

        static CScalarPixelType Get( const CRGBPixelType &input )
        {
            return static_cast<CScalarPixelType>( input.GetBlue() );
        }
    };

    class RedChannelPixelAccessor
    {
    public:
        typedef CRGBPixelType       InternalType;
        typedef CScalarPixelType    ExternalType;

        static CScalarPixelType Get( const CRGBPixelType &input )
        {
            return static_cast<CScalarPixelType>( input.GetRed() );
        }
    };

    typedef itk::ImageAdaptor<CRGBImageType, BlueChannelPixelAccessor>   BlueAdaptorType;
    typedef itk::ImageAdaptor<CRGBImageType, RedChannelPixelAccessor>    RedAdaptorType;

    typedef itk::BinaryBallStructuringElement<CScalarImageType::PixelType, 2>     StructuringElementType;

    typedef itk::ImageFileReader<CRGBImageType>         CRGBReaderType;
    typedef itk::ImageFileWriter<CScalarImageType>      CScalarWriterType;
    typedef itk::ImageFileWriter<FScalarImageType>      FScalarWriterType;
    typedef itk::ImageFileWriter<CRGBImageType>         RGBVoWriterType;

    typedef itk::BinaryThresholdImageFilter<CScalarImageType, CScalarImageType>                                     ThresholdFilterType;
    typedef itk::BinaryMorphologicalClosingImageFilter<CScalarImageType, CScalarImageType, StructuringElementType>  ClosingImageFilterType;
    typedef itk::ImageFileReader<CRGBImageType>                                                                     CRGBVoReaderType;
    typedef itk::RescaleIntensityImageFilter<BlueAdaptorType, CScalarImageType>                                     BlueRescaleIntensityImageFilterType;
    typedef itk::RescaleIntensityImageFilter<RedAdaptorType, CScalarImageType>                                      RedRescaleIntensityImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarImageType>                                                 ImageToShapeLabelMapFilterType;
    typedef itk::ShapeOpeningLabelMapFilter<ImageToShapeLabelMapFilterType::OutputImageType>                        ShapeOpeningLabelMapFilterType;
    typedef itk::LabelMapToLabelImageFilter<ImageToShapeLabelMapFilterType::OutputImageType, CScalarImageType>      LabelMapToLabelImageFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarImageType, CScalarImageType, CRGBImageType>                         LabelOverlayImageFilterType;
    typedef itk::SignedMaurerDistanceMapImageFilter<CScalarImageType, FScalarImageType>                             DistanceMapImageFilterType;
    typedef itk::RescaleIntensityImageFilter<FScalarImageType, SScalarImageType>                                    RescaleFloatToShortImageFilterType;

public:
    ProcessHoechstDiffusionFilter();
    virtual ~ProcessHoechstDiffusionFilter();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);

    void BuildFilename();
    void WriteDistBlueValueTableFile();


    std::string mDatasetID;

    std::string mFullFilename;
    std::string mPath;
    std::string mFilename;
    std::string mFileExtension;

    unsigned char mThreshold;
    StructuringElementType mClosingStructuringElement;
    unsigned int mMinimalSpheroidSize;

    double mTimeStep;
    double mVoxSpacing[2];

    int mFrame;
    std::vector<float> mDists;
    std::vector<CScalarPixelType> mBlueValues;
};

#endif /* PROCESSHOECHSTDIFFUSIONFILTER_H_ */
