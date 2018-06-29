///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentAndClassifyStellateCells.h                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-01-14                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SEGMENTANDCLASSIFYSTELLATECELLS_H_
#define SEGMENTANDCLASSIFYSTELLATECELLS_H_

#include "BasePipeline.h"

#include "itkRGBPixel.h"

#include "itkLabelMap.h"

#include "itkBinaryBallStructuringElement.h"

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkHConvexImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include "../filters/imageFilters/itkAdaptiveOtsuThresholdingImageFilter.h"



class SegmentAndClassifyStellateCells : public BasePipeline<3>
{
    typedef itk::RGBPixel<CScalarPixelType>             RGBPixelType;

    typedef itk::Image<RGBPixelType,3>                  RGBVoImageType;

    typedef itk::ImageFileWriter<RGBVoImageType>        RGBVoWriterType;

    typedef itk::BinaryBallStructuringElement<CScalarVoImageType::PixelType, 3>     StructuringElementType;

    typedef itk::AdaptiveOtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                               AdaptiveOtsuThresholdImageFilterType;
    typedef itk::BinaryMorphologicalOpeningImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>  OpeningImageFilterType;
    typedef itk::GrayscaleDilateImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>             GreyscaleDilateImageFilterType;
    typedef itk::GrayscaleErodeImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>              GreyscaleErodeImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                                                   ImageToShapeLabelMapFilterType;
    typedef itk::BinaryThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                     ThresholdFilterType;
    typedef itk::HConvexImageFilter<CScalarVoImageType, CScalarVoImageType>                                             ConvexImageFilterType;
    typedef itk::InvertIntensityImageFilter<CScalarVoImageType>                                                         InvertIntensityImageFilterType;
    typedef itk::MaskImageFilter<CScalarVoImageType, CScalarVoImageType>                                                MaskImageFilterType;
    typedef itk::MedianImageFilter<CScalarVoImageType, CScalarVoImageType >                                             MedianImageFilterType;
    typedef itk::LabelMapToLabelImageFilter<ImageToShapeLabelMapFilterType::OutputImageType, CScalarVoImageType>        LabelMapToLabelImageFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarVoImageType, CScalarVoImageType, RGBVoImageType>                        LabelOverlayImageFilterType;
    typedef itk::OtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                       OtsuThresholdImageFilterType;
    typedef itk::ShapeOpeningLabelMapFilter<ImageToShapeLabelMapFilterType::OutputImageType>                            ShapeOpeningLabelMapFilterType;
    typedef itk::SubtractImageFilter<CScalarVoImageType, CScalarVoImageType>                                            SubtractImageFilterType;
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter<CScalarVoImageType>                                        HoleFillingImageFilterType;

public:
    SegmentAndClassifyStellateCells();
    virtual ~SegmentAndClassifyStellateCells();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteDataSetSummary();

    std::string m_fullFilenameDesminChannel;
    std::string m_pathDesminChannel;
    std::string m_filenameDesminChannel;
    std::string m_fileExtensionDesminChannel;

    std::string m_fullFilenameDAPIChannel;
    std::string m_pathDAPIChannel;
    std::string m_filenameDAPIChannel;
    std::string m_fileExtensionDAPIChannel;

    std::string m_fullFilenameNonHepNucleiSeg;
    std::string m_pathNonHepNucleiSeg;
    std::string m_filenameNonHepNucleiSeg;
    std::string m_fileExtensionNonHepNucleiSeg;

    std::string m_fullFilenameHepNucleiSeg;
    std::string m_pathHepNucleiSeg;
    std::string m_filenameHepNucleiSeg;
    std::string m_fileExtensionHepNucleiSeg;

    double m_spacing[3];

    int m_medianRadius;
    int m_greyscaleOpeningRadius;
    bool m_withConvexFilterPreprocessing;
    int m_convexFilterHeight;

    int m_thresholdingMode;                //0 - adaptive otsu, 1 - normal otsu, 2 - manual
    CScalarVoImageType::SizeType m_adapOtsuRadius;
    int m_adapOtsuSamplePoints;
    CScalarPixelType m_otsuThreshold;
    CScalarPixelType m_lowerThreshold;
    CScalarPixelType m_upperThreshold;

    CScalarVoImageType::SizeType m_inverseHoleFillingNeighborhoodRadius;
    unsigned int m_inverseHoleFillingMajThreshold;

    unsigned int m_minimalStellateCellSize;

    double m_nucleiBodyOverlap;

    float m_overlayOpacity;

    bool m_saveEverything;                   //1 - save everything, 0 - save only essentials
    std::string m_logFilenameSave;
    std::string m_filenameSave;
    std::string m_finalSaveSuffixes[3];
};

#endif /* SEGMENTANDCLASSIFYSTELLATECELLS_H_ */
