///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentBileNetworkOnDPPIV60x.h                                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SEGMENTBILENETWORKONDPPIV60X_H_
#define SEGMENTBILENETWORKONDPPIV60X_H_

#include "BasePipeline.h"

#include "itkRGBPixel.h"

#include "itkLabelMap.h"

#include "itkBinaryBallStructuringElement.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include "../filters/imageFilters/itkAdaptiveOtsuThresholdingImageFilter.h"
#include "../filters/imageFilters/itkBinaryThinningImageFilter3D.h"
#include "../filters/imageFilters/MultiChannelBinaryThresholdImageFilter.h"


class SegmentBileNetworkOnDPPIV60x : public BasePipeline<3>
{
    typedef itk::RGBPixel<CScalarPixelType>             RGBPixelType;

    typedef itk::Image<RGBPixelType,3>                  RGBVoImageType;

    typedef itk::ImageFileWriter<RGBVoImageType>        RGBVoWriterType;

    typedef itk::BinaryBallStructuringElement<CScalarVoImageType::PixelType, 3>     StructuringElementType;

    typedef itk::AdaptiveOtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                               AdaptiveOtsuThresholdImageFilterType;
    typedef itk::BinaryDilateImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>                BinaryDilateFilterType;
    typedef itk::BinaryMorphologicalOpeningImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>  OpeningImageFilterType;
    typedef itk::GrayscaleDilateImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>             GreyscaleDilateImageFilterType;
    typedef itk::GrayscaleErodeImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>              GreyscaleErodeImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                                                   ImageToShapeLabelMapFilterType;
    typedef itk::BinaryThinningImageFilter3D<CScalarVoImageType,CScalarVoImageType>                                     Thinning3DImageFilterType;
    typedef itk::BinaryThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                     ThresholdFilterType;
    typedef itk::InvertIntensityImageFilter<CScalarVoImageType>                                                         InvertIntensityImageFilterType;
    typedef itk::MedianImageFilter<CScalarVoImageType, CScalarVoImageType >                                             MedianImageFilterType;
    typedef itk::LabelMapToLabelImageFilter<ImageToShapeLabelMapFilterType::OutputImageType, CScalarVoImageType>        LabelMapToLabelImageFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarVoImageType, CScalarVoImageType, RGBVoImageType>                        LabelOverlayImageFilterType;
    typedef itk::OrImageFilter<CScalarVoImageType, CScalarVoImageType, CScalarVoImageType>                              OrImageFilterType;
    typedef itk::OtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                       OtsuThresholdImageFilterType;
    typedef itk::ShapeOpeningLabelMapFilter<ImageToShapeLabelMapFilterType::OutputImageType>                            ShapeOpeningLabelMapFilterType;
    typedef itk::SubtractImageFilter<CScalarVoImageType, CScalarVoImageType>                                            SubtractImageFilterType;
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter<CScalarVoImageType>                                        HoleFillingImageFilterType;

public:
    SegmentBileNetworkOnDPPIV60x();
    virtual ~SegmentBileNetworkOnDPPIV60x();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteDataSetSummary();

    QString m_fullFilenameChannelDPPIV;
    QFileInfo m_infoFullFilenameChannelDPPIV;
    std::string m_pathChannelDPPIV;
    std::string m_filenameChannelDPPIV;
    std::string m_fileExtensionChannelDPPIV;

    QString m_fullFilenameChannelLastStep;
    QFileInfo m_infoFullFilenameChannelLastStep;
    std::string m_pathChannelLastStep;
    std::string m_filenameChannelLastStep;
    std::string m_fileExtensionChannelLastStep;

    QString m_fullFilenameNecroticRegion;
    QFileInfo m_infoFullFilenameNecroticRegion;
    std::string m_pathNecroticRegion;
    std::string m_filenameNecroticRegion;
    std::string m_fileExtensionNecroticRegion;

    bool m_hasNR;

    QString m_fullFilenameSegCV;
    QFileInfo m_infoFullFilenameSegCV;
    QString m_fullFilenameSegPV;
    QFileInfo m_infoFullFilenameSegPV;

    bool m_withCVMask;
    bool m_withPVMask;

    CScalarVoImageType::SpacingType m_spacing;

    unsigned int m_entryPoint;

    int m_medianRadius;
    int m_greyscaleOpeningRadius;

    int m_thresholdingMode;                //0 - adaptive otsu, 1 - normal otsu, 2 - manual
    CScalarVoImageType::SizeType m_adapOtsuRadius;
    int m_adapOtsuSamplePoints;
    CScalarPixelType m_otsuThreshold;
    CScalarPixelType m_lowerThreshold;
    CScalarPixelType m_upperThreshold;

    CScalarVoImageType::SizeType m_holeFilling1NeighborhoodRadius;
    unsigned int m_holeFilling1MajThreshold;

    CScalarVoImageType::SizeType m_inverseHoleFillingNeighborhoodRadius;
    unsigned int m_inverseHoleFillingMajThreshold;

    CScalarVoImageType::SizeType m_openingNeighborhoodRadius;

    CScalarVoImageType::SizeType m_maskVeinRadius;

    unsigned int m_minimalBileSize;

    float m_overlayOpacity;

    bool m_saveEverything;                   //1 - save everything, 0 - save only essentials
    std::string m_logFilenameSave;
    std::string m_filenameSave;
    std::string m_finalSaveSuffixes[3];
};

#endif /* SEGMENTBILENETWORKONDPPIV60X_H_ */
