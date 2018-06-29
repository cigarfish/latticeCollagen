///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  BasePipeline.h                                                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-15                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef BASEPIPELINE_H_
#define BASEPIPELINE_H_

#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>

#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIO.h>
#endif

#include <QFileInfo>
#include <QString>

#include "../../tools/parameters/CSParameterChoice.h"

class CSParameterContext;

template< unsigned int VImageDimension > class BasePipeline
{
public:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

protected:
    typedef unsigned char                                   CScalarPixelType;
    typedef float                                           FScalarPixelType;
    typedef long                                            LScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType>                 RGBCPixelType;

    typedef itk::Image<CScalarPixelType, VImageDimension>   CScalarVoImageType;
    typedef itk::Image<FScalarPixelType, VImageDimension>   FScalarVoImageType;
    typedef itk::Image<LScalarPixelType, VImageDimension>   LScalarVoImageType;
    typedef itk::Image<RGBCPixelType, VImageDimension>      RGBCVoImageType;

    typedef typename CScalarVoImageType::Pointer            CScalarImagePointerType;
    typedef typename FScalarVoImageType::Pointer            FScalarImagePointerType;
    typedef typename LScalarVoImageType::Pointer            LScalarImagePointerType;
    typedef typename RGBCVoImageType::Pointer               RGBCImagePointerType;


    typedef typename CScalarVoImageType::RegionType         CScalarRegionType;
    typedef typename CScalarVoImageType::SpacingType        CScalarSpacingType;
    typedef typename CScalarVoImageType::SizeType           CScalarSizeType;
    typedef typename CScalarVoImageType::IndexType          CScalarIndexType;


    typedef itk::ImageFileReader<CScalarVoImageType>        ScalarVoReaderType;
    typedef itk::ImageFileReader<RGBCVoImageType>           RGBVoReaderType;
    typedef itk::ImageFileWriter<CScalarVoImageType>        ScalarVoWriterType;

    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                       ImageToShapeLabelMapFilterType;
    typedef itk::BinaryThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>         ThresholdFilterType1;
    typedef itk::BinaryThresholdImageFilter<LScalarVoImageType, CScalarVoImageType>         ThresholdFilterType2;
    typedef typename ImageToShapeLabelMapFilterType::OutputImageType                        ShapeLabelMapType;
    typedef itk::LabelMapToLabelImageFilter<ShapeLabelMapType, LScalarVoImageType>          LabelMapToLabelImageFilterType;
    typedef itk::MinimumMaximumImageCalculator<CScalarVoImageType>                          CMinMaxCalculatorType;
    typedef itk::SignedMaurerDistanceMapImageFilter<CScalarVoImageType, FScalarVoImageType> SignedMaurerDistanceMapImageFilterType;

public:
    BasePipeline() { m_paramContext = NULL; };
    virtual ~BasePipeline() {};

    void SetParameterContext(CSParameterContext *paramContext)
    {
        m_paramContext = paramContext;
    };

    static void ReadImage(std::string filename, CScalarImagePointerType &image, CScalarSpacingType spacing);
    static void ReadImage(std::string filename, RGBCImagePointerType &image, CScalarSpacingType spacing);
    static void ReadImageLayer(std::string filename, CScalarImagePointerType &image, CScalarSpacingType spacing, CScalarPixelType foreground);
    static void CreateImage(CScalarImagePointerType &image, CScalarRegionType region, CScalarPixelType constant, CScalarSpacingType spacing);
    static void CreateImage(FScalarImagePointerType &image, CScalarRegionType region, FScalarPixelType constant, CScalarSpacingType spacing);
    static void BuildDistanceMap(CScalarImagePointerType &inImage, FScalarImagePointerType &distMap, CScalarSpacingType spacing, bool distanceToSegmentation);
    static void BuildNonZeroDistanceMap(CScalarImagePointerType &inImage, FScalarImagePointerType &distMap, FScalarPixelType constant, CScalarSpacingType spacing, bool distanceToSegmentation);
    static void BuildLabelImage(CScalarImagePointerType &inImage, LScalarImagePointerType &labelImage, CScalarSpacingType spacing, bool fullyConneceted);
    static void LabelToBinImage(CScalarImagePointerType inLabelImage, unsigned char labelToExtract, CScalarImagePointerType &outBinImage);
    static void LabelToBinImage(LScalarImagePointerType inLabelImage, unsigned long labelToExtract, CScalarImagePointerType &outBinImage);

    static unsigned int GetNumberOfDimensions(std::string file);
    static long double ComputeVolumeOfRegions(std::string filename, CScalarSpacingType spacing, long double voxelVolume, CScalarPixelType foreground = 255);
    static long double ComputeVolumeOfRegions(CScalarImagePointerType &image, long double voxelVolume, CScalarPixelType foreground = 255);

    virtual void Update() = 0;

protected:
    virtual void ParseParameterContext() = 0;
    virtual void WriteLogFile(std::string timeStamp) = 0;

    CSParameterContext *m_paramContext;
};

#include "BasePipeline.tpp"

#endif /* BASEPIPELINE_H_ */
