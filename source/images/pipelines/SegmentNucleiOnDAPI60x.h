///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentNucleiOnDAPI60x.h                                             //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-07                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SEGMENTNUCLEIONDAPI60X_H_
#define SEGMENTNUCLEIONDAPI60X_H_

#include "BasePipeline.h"

#include <math.h>

#include "itkRGBPixel.h"

#include "itkLabelMap.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkMath.h"

#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"
#include "itkNearestNeighborExtrapolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"


#include "../filters/imageFilters/itkAdaptiveOtsuThresholdingImageFilter.h"
#include "../filters/imageFilters/IterativeCavityFillingImageFilter.h"


template< unsigned int VImageDimension > class SegmentNucleiOnDAPI60x : public BasePipeline<VImageDimension>
{
public:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

protected:
    typedef double                                                      CoordRepType;
    typedef typename BasePipeline<VImageDimension>::CScalarPixelType    CScalarPixelType;
    typedef typename BasePipeline<VImageDimension>::FScalarPixelType    FScalarPixelType;
    typedef unsigned short                                              SScalarPixelType;
    typedef unsigned int                                                IScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType>                             RGBPixelType;
    typedef itk::ShapeLabelObject<IScalarPixelType, VImageDimension>    ShapeLabelObjectType;

    typedef typename BasePipeline<VImageDimension>::CScalarVoImageType  CScalarImageType;
    typedef typename BasePipeline<VImageDimension>::FScalarVoImageType  FScalarImageType;
    typedef itk::Image<SScalarPixelType, VImageDimension>               SScalarImageType;
    typedef itk::Image<IScalarPixelType, VImageDimension>               IScalarImageType;
    typedef itk::Image<RGBPixelType, VImageDimension>                   RGBImageType;
    typedef itk::LabelMap<ShapeLabelObjectType>                         ShapeLabelMapType;


    typedef typename BasePipeline<VImageDimension>::ScalarVoReaderType  CScalarImageReaderType;
    typedef typename BasePipeline<VImageDimension>::ScalarVoWriterType  CScalarImageWriterType;
    typedef itk::ImageFileWriter<SScalarImageType>                      SScalarWriterType;
    typedef itk::ImageFileWriter<RGBImageType>                          RGBWriterType;


    typedef itk::BinaryBallStructuringElement<CScalarPixelType, VImageDimension>                                    StructuringElementType;

    typedef itk::AdaptiveOtsuThresholdImageFilter<CScalarImageType, CScalarImageType>                               AdaptiveOtsuThresholdImageFilterType;
    typedef itk::OrImageFilter<CScalarImageType>                                                                    OrImageFilterType;
    typedef itk::BinaryErodeImageFilter<CScalarImageType, CScalarImageType, StructuringElementType>                 ErodeImageFilterType;
    typedef itk::BinaryMorphologicalClosingImageFilter<CScalarImageType, CScalarImageType, StructuringElementType>  ClosingImageFilterType;
    typedef itk::BinaryMorphologicalOpeningImageFilter<CScalarImageType, CScalarImageType, StructuringElementType>  OpeningImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarImageType>                                                 ImageToShapeLabelMapFilterType;
    typedef itk::BinaryThresholdImageFilter<CScalarImageType, CScalarImageType>                                     ThresholdFilterType;
    typedef itk::SignedMaurerDistanceMapImageFilter<CScalarImageType, FScalarImageType>                             MaurerDistanceMapImageFilterType;
    typedef itk::InvertIntensityImageFilter<CScalarImageType>                                                       InvertIntensityImageFilterType;
    typedef itk::InvertIntensityImageFilter<FScalarImageType>                                                       FInvertIntensityImageFilterType;
    typedef itk::MaskImageFilter<IScalarImageType, CScalarImageType, IScalarImageType>                              MaskImageFilterType;
    typedef itk::MorphologicalWatershedImageFilter<FScalarImageType, IScalarImageType>                              MorphoWatershedImageFilterType;
    typedef itk::NearestNeighborExtrapolateImageFunction<CScalarImageType, CoordRepType>                            NearestNeighborExtrapolatorImageFunctionType;
    typedef itk::NearestNeighborInterpolateImageFunction<CScalarImageType, CoordRepType>                            NearestNeighborInterpolatorImageFunctionType;
    typedef itk::LabelImageToShapeLabelMapFilter<IScalarImageType>                                                  LabelImageToShapeLabelMapFilterType;
    typedef itk::LabelMapToLabelImageFilter<ShapeLabelMapType, CScalarImageType>                                    LabelMapToLabelImageFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarImageType, CScalarImageType, RGBImageType>                          LabelOverlayImageFilterType;
    typedef itk::OtsuThresholdImageFilter<CScalarImageType, CScalarImageType>                                       OtsuThresholdImageFilterType;
    typedef itk::ResampleImageFilter<CScalarImageType, CScalarImageType>                                            ResampleImageFilterType;
    typedef itk::RescaleIntensityImageFilter<FScalarImageType, SScalarImageType>                                    RescaleImageFilterType;
    typedef itk::ShapeOpeningLabelMapFilter<ShapeLabelMapType>                                                      ShapeOpeningLabelMapFilterType;
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter<CScalarImageType>                                      HoleFillingImageFilterType;

    typedef itk::GrayscaleDilateImageFilter<CScalarImageType, CScalarImageType, StructuringElementType>             GrayscaleDilateImageFilterType;
    typedef itk::GrayscaleErodeImageFilter<CScalarImageType, CScalarImageType, StructuringElementType>              GrayscaleErodeImageFilterType;
    typedef itk::IterativeCavityFillingImageFilter<CScalarImageType>                                                CavityFillingImageFilter;
    typedef itk::MedianImageFilter<CScalarImageType, CScalarImageType >                                             MedianImageFilterType;
    typedef itk::MedianImageFilter<FScalarImageType, CScalarImageType >                                             MedianImageFilterType2;


public:
    SegmentNucleiOnDAPI60x();
    virtual ~SegmentNucleiOnDAPI60x();

    void Update();

private:
    void ClassifyNuclei2D(itk::SmartPointer<ShapeLabelMapType> nucleiLabelMap, std::vector<unsigned long> &labelsToRemoveFromHepNucleiLabelMap,
            std::vector<unsigned long> &labelsToRemoveFromNonHepNucleiLabelMap, bool withNucleiAnalysisOutput);
    void ClassifyNuclei3D(itk::SmartPointer<ShapeLabelMapType> nucleiLabelMap, itk::SmartPointer<FScalarImageType> distanceMap, std::vector<unsigned long> &labelsToRemoveFromHepNucleiLabelMap,
            std::vector<unsigned long> &labelsToRemoveFromNonHepNucleiLabelMap, bool withNucleiAnalysisOutput);

    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteDataSetSummary();

	QString m_fullFilenameChannelDAPI;
    QFileInfo m_infoFullFilenameChannelDAPI;
    std::string m_pathChannelDAPI;
    std::string m_filenameChannelDAPI;
    std::string m_fileExtensionChannelDAPI;

    QString m_fullFilenameSegCV;
    QFileInfo m_infoFullFilenameSegCV;
    QString m_fullFilenameSegPV;
    QFileInfo m_infoFullFilenameSegPV;

    bool m_withCVMask;
    bool m_withPVMask;

    typename CScalarImageType::SpacingType m_spacing;

    int m_medianRadius;
    int m_greyscaleOpeningRadius;

    int m_thresholdingMode;                //0 - adaptive otsu, 1 - normal otsu, 2 - manual
    typename CScalarImageType::SizeType m_adapOtsuRadius;
    int m_adapOtsuSamplePoints;
    CScalarPixelType m_otsuThreshold;
    CScalarPixelType m_lowerThreshold;
    CScalarPixelType m_upperThreshold;

    typename CScalarImageType::SizeType m_invHoleFillingNeighborhoodRadius;
    unsigned int m_invHoleFillingMajThreshold;

    bool m_cavWithRescaling;
    int m_cavitiyFillingRadius;
    double m_cavitiyFillingMinThreshold;
    double m_cavitiyFillingMaxThreshold;

    StructuringElementType m_closingStructuringElement;
    StructuringElementType m_openingStructuringElement;

    double m_minNonHepRadius;
    double m_minNonHepVolume;
    double m_maxNonHepRadius;
    double m_maxNonHepVolume;
    double m_minHepRadius;
    double m_minHepVolume;
    double m_maxHepRadius;
    double m_maxHepVolume;
    double m_roundnessThreshold;

    float m_floodLevel;

    float m_overlayOpacity;

    bool m_saveEverything;                   //1 - save everything, 0 - save only essentials
    std::string m_logFilenameSave;
    std::string m_filenameSave;
    std::string m_saveSuffixesForFinals[4];
};

#include "SegmentNucleiOnDAPI60x.tpp"

#endif /* SEGMENTNUCLEIONDAPI60X_H_ */
