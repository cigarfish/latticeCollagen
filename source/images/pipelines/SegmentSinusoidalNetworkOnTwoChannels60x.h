///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentSinusoidalNetworkOnTwoChannels60x.h                           //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SEGMENTSINUSOIDALNETWORKONTWOCHANNELS60X_H_
#define SEGMENTSINUSOIDALNETWORKONTWOCHANNELS60X_H_

#include "BasePipeline.h"

#include "itkRGBPixel.h"

#include "itkLabelMap.h"

#include "itkBinaryBallStructuringElement.h"

#include "itkAndImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#if (ITK_VERSION_MAJOR < 4 || ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR < 4 || ITK_VERSION_MAJOR > 4)
#include "itkCompose2DVectorImageFilter.h"
#else
#include "itkComposeImageFilter.h"
#endif
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkNearestNeighborExtrapolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include "../filters/imageFilters/itkAdaptiveOtsuThresholdingImageFilter.h"
#include "../filters/imageFilters/itkBinaryThinningImageFilter3D.h"
#include "../filters/imageFilters/MultiChannelBinaryThresholdImageFilter.h"
#include "../filters/imageFilters/IterativeCavityFillingImageFilter.h"
#include "../filters/imageFilters/CavityFillingImageFilter.h"


class SegmentSinusoidalNetworkOnTwoChannels60x : public BasePipeline<3>
{
    typedef itk::Vector<CScalarPixelType, 2>            CVector2DPixelType;
    typedef itk::RGBPixel<CScalarPixelType>             RGBPixelType;
    typedef double                                      CoordRepType;

    typedef itk::Image<CVector2DPixelType, 3>           Vector2DVoImageType;
    typedef itk::Image<RGBPixelType,3>                  RGBVoImageType;

    typedef itk::BinaryBallStructuringElement<CScalarVoImageType::PixelType, 3>     StructuringElementType;

    typedef itk::ImageFileWriter<RGBVoImageType>        RGBVoWriterType;

    typedef itk::AdaptiveOtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                               AdaptiveOtsuThresholdImageFilterType;
    typedef itk::AndImageFilter<CScalarVoImageType, CScalarVoImageType, CScalarVoImageType>                             AndImageFilterType;
    typedef itk::BinaryDilateImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>                DilateImageFilterType;
    typedef itk::BinaryErodeImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>                 ErodeImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                                                   ImageToShapeLabelMapFilterType;
    typedef itk::BinaryMorphologicalClosingImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>  ClosingImageFilterType;
    typedef itk::BinaryMorphologicalOpeningImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>  OpeningImageFilterType;
    typedef itk::BinaryThinningImageFilter3D<CScalarVoImageType,CScalarVoImageType>                                     Thinning3DImageFilterType;
    typedef itk::BinaryThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                     ThresholdFilterType;
#if (ITK_VERSION_MAJOR < 4 || ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR < 4 || ITK_VERSION_MAJOR > 4)
    typedef itk::Compose2DVectorImageFilter<CScalarVoImageType, Vector2DVoImageType>                                    ComposeVectorImageFilterType;
#else
    typedef itk::ComposeImageFilter<CScalarVoImageType, Vector2DVoImageType>                                            ComposeVectorImageFilterType;
#endif
    typedef itk::GrayscaleDilateImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>             GrayscaleDilateImageFilterType;
    typedef itk::GrayscaleErodeImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>              GrayscaleErodeImageFilterType;
    typedef itk::InvertIntensityImageFilter<CScalarVoImageType>                                                         InvertIntensityImageFilterType;
    typedef itk::IterativeCavityFillingImageFilter<CScalarVoImageType>                                                  CavityFillingImageFilter;
    typedef itk::NearestNeighborExtrapolateImageFunction<CScalarVoImageType, CoordRepType>                              NearestNeighborExtrapolatorImageFunctionType;
    typedef itk::NearestNeighborInterpolateImageFunction<CScalarVoImageType, CoordRepType>                              NearestNeighborInterpolatorImageFunctionType;
    typedef itk::LabelMapToLabelImageFilter<ImageToShapeLabelMapFilterType::OutputImageType, CScalarVoImageType>        LabelMapToLabelImageFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarVoImageType, CScalarVoImageType, RGBVoImageType>                        LabelOverlayImageFilterType;
    typedef itk::MedianImageFilter<CScalarVoImageType, CScalarVoImageType >                                             MedianImageFilterType;
    typedef itk::MultiChannelBinaryThresholdImageFilter<Vector2DVoImageType, CScalarVoImageType>                        Threshold2ChannelFilterType;
    typedef itk::OtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                       OtsuThresholdImageFilterType;
    typedef itk::OrImageFilter<CScalarVoImageType, CScalarVoImageType, CScalarVoImageType>                              OrImageFilterType;
    typedef itk::ResampleImageFilter<CScalarVoImageType, CScalarVoImageType>                                            ResampleImageFilterType;
    typedef itk::ShapeOpeningLabelMapFilter<ImageToShapeLabelMapFilterType::OutputImageType>                            ShapeOpeningLabelMapFilterType;
    typedef itk::SubtractImageFilter<CScalarVoImageType, CScalarVoImageType>                                            SubstractImageFilterType;
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter<CScalarVoImageType>                                        HoleFillingImageFilterType;

public:
    SegmentSinusoidalNetworkOnTwoChannels60x();
    virtual ~SegmentSinusoidalNetworkOnTwoChannels60x();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteDataSetSummary();

	QString m_fullFilenameChannel1;
    QFileInfo m_infoFullFilenameChannel1;
    std::string m_pathChannel1;
    std::string m_filenameChannel1;
    std::string m_fileExtensionChannel1;

	QString m_fullFilenameChannel2;
    QFileInfo m_infoFullFilenameChannel2;
    std::string m_pathChannel2;
    std::string m_filenameChannel2;
    std::string m_fileExtensionChannel2;

    QString m_fullFilenameNecroticRegion;
    QFileInfo m_infoFullFilenameNecroticRegion;
    std::string m_pathNecroticRegion;
    std::string m_filenameNecroticRegion;
    std::string m_fileExtensionNecroticRegion;

    QString m_fullFilenameSegCV;
    QFileInfo m_infoFullFilenameSegCV;
    QString m_fullFilenameSegPV;
    QFileInfo m_infoFullFilenameSegPV;

    bool m_withCVMask;
    bool m_withPVMask;

    CScalarVoImageType::SpacingType m_spacing;

    bool m_hasNR;

    unsigned int m_entryPoint;

    int m_medianRadius;
    int m_greyscaleClosingRadius;

    int m_thresholdingMode[2];                //0 - adaptive otsu, 1 - normal otsu, 2 - manual
    CScalarVoImageType::SizeType m_adapOtsuRadius[2];
    int m_adapOtsuSamplePoints[2];
    CScalarPixelType m_otsuThreshold[2];
    CVector2DPixelType m_lowerThreshold;
    CVector2DPixelType m_upperThreshold;

    CScalarVoImageType::SizeType m_inverseHoleFillingNeighborhoodRadius;
    unsigned int m_inverseHoleFillingMajThreshold;

    bool m_cavWithRescaling;
    int m_cavityRadius;
    double m_cavityMinFrac;
    double m_cavityMaxFrac;

    StructuringElementType m_closingStructuringElement;
    StructuringElementType m_openingStructuringElement;

    CScalarVoImageType::SizeType m_maskVeinRadius;

    unsigned int m_minimalSinusoidSize;

    float m_overlayOpacity;

    bool m_saveEverything;                   //1 - save everything, 0 - save only essentials
    std::string m_logFilenameSave;
    std::string m_filenameSave;
    std::string m_finalSaveSuffixes[3];
};

#endif /* SEGMENTSINUSOIDALNETWORKONTWOCHANNELS60X_H_ */

