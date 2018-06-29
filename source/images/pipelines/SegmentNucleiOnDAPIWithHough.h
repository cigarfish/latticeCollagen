///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentNucleiOnDAPIWithHough.h                                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-23-04                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SEGMENTNUCLEIONDAPIWITHHOUGH_H_
#define SEGMENTNUCLEIONDAPIWITHHOUGH_H_

#include "BasePipeline.h"

#include <math.h>

#include "itkRGBPixel.h"

#include "itkLabelMap.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkMath.h"

#include "itkEllipseSpatialObject.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGroupSpatialObject.h"
#include "itkHoughTransformRadialVotingImageFilter.h"
#include "itkMaskedSpatialObjectToImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkPointSet.h"
#include "itkSpatialObjectToImageFilter.h"


class SegmentNucleiOnDAPIWithHough : public BasePipeline<3>
{
    static const unsigned int Dimension = 3;

    typedef float                                       FScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType>             RGBPixelType;

    typedef itk::Image<FScalarPixelType, Dimension>     FScalarVoImageType;
    typedef itk::Image<RGBPixelType, Dimension>         RGBVoImageType;

    typedef itk::ImageFileWriter<RGBVoImageType>        RGBVoWriterType;

    typedef itk::BinaryBallStructuringElement<CScalarVoImageType::PixelType, Dimension>     StructuringElementType;
    typedef itk::EllipseSpatialObject<Dimension>                                            SphereType;
    typedef itk::GroupSpatialObject<Dimension>                                              GroupType;
    typedef GroupType::TransformType                                                        TransformType;

    typedef itk::GrayscaleDilateImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>     GrayscaleDilateImageFilterType;
    typedef itk::GrayscaleErodeImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>      GrayscaleErodeImageFilterType;
    typedef itk::HoughTransformRadialVotingImageFilter<CScalarVoImageType, FScalarVoImageType>                  HoughTransformFilterType;
    typedef HoughTransformFilterType::SpheresListType                                                           SpheresListType;
    typedef itk::MaskedSpatialObjectToImageFilter<GroupType, CScalarVoImageType>                                SpatialObjectToImageFilterType;
    typedef itk::MedianImageFilter<CScalarVoImageType, CScalarVoImageType >                                     MedianImageFilterType;
    typedef itk::LabelImageToShapeLabelMapFilter<CScalarVoImageType>                                            LabelImageToShapeLabelMapFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarVoImageType, CScalarVoImageType, RGBVoImageType>                LabelOverlayImageFilterType;

    typedef itk::PointSet<LabelImageToShapeLabelMapFilterType::OutputImageType::PixelType, Dimension>   PointSetType;
    typedef PointSetType::PointType                                                                     PointSetPointType;

public:
    SegmentNucleiOnDAPIWithHough();
    virtual ~SegmentNucleiOnDAPIWithHough();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteDataSetSummary();

    std::string m_fullFilenameChannelDAPI;
    std::string m_pathChannelDAPI;
    std::string m_filenameChannelDAPI;
    std::string m_fileExtensionChannelDAPI;

    CScalarVoImageType::SpacingType m_spacing;

    int m_medianRadius;
    int m_greyscaleOpeningRadius;

    int m_houghNumSpheres;
    int m_houghThreshold;
    int m_houghGradientThreshold;
    double m_houghOutputThreshold;
    double m_houghMinRadius;
    double m_houghMaxRadius;
    double m_houghSigmaGradient;
    double m_houghVariance;
    double m_houghSphereRadiusRatio;
    double m_houghVotingRadiusRatio;
    double m_houghSamplingRatio;
    int m_sphereToImageResampleFactor;
    int m_sphereToImageDilationRadius;

    double m_minNonHepRadius;           //Dummies: Not implemented yet
    double m_minHepRadius;
    double m_maxHepRadius;

    float m_overlayOpacity;

    bool m_saveEverything;                   //1 - save everything, 0 - save only essentials
    std::string m_logFilenameSave;
    std::string m_filenameSave;
    std::string m_saveSuffixesForFinals[4];
};

#endif /* SEGMENTNUCLEIONDAPIWITHHOUGH_H_ */
