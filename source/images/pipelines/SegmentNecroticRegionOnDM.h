///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentNecroticRegionOnDM.h                                          //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-05-14                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SEGMENTNECROTICREGIONONDM_H_
#define SEGMENTNECROTICREGIONONDM_H_

#include "BasePipeline.h"

#include "itkRGBPixel.h"

#include "itkLabelMap.h"

#include "itkBinaryBallStructuringElement.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"

#include "../imageFilters/itkAdaptiveOtsuThresholdingImageFilter.h"


class SegmentNecroticRegionOnDM : public BasePipeline<3>
{
    typedef itk::RGBPixel<CScalarPixelType>             RGBPixelType;

    typedef itk::Image<RGBPixelType,3>                  RGBVoImageType;

    typedef itk::ImageFileWriter<RGBVoImageType>        RGBVoWriterType;

    typedef itk::BinaryBallStructuringElement<CScalarVoImageType::PixelType, 3>     StructuringElementType;

    typedef itk::AdaptiveOtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                           AdaptiveOtsuThresholdImageFilterType;
    typedef itk::BinaryDilateImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>            DilateImageFilterType;
    typedef itk::BinaryErodeImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>             ErodeImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                                               ImageToShapeLabelMapFilterType;
    typedef itk::BinaryThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                 ThresholdFilterType;
	typedef itk::InvertIntensityImageFilter<CScalarVoImageType>                                                     InvertIntensityImageFilterType;
    typedef itk::LabelMapToLabelImageFilter<ImageToShapeLabelMapFilterType::OutputImageType, CScalarVoImageType>    LabelMapToLabelImageFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarVoImageType, CScalarVoImageType, RGBVoImageType>                    LabelOverlayImageFilterType;
	typedef itk::OtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                   OtsuThresholdImageFilterType;
    typedef itk::ShapeOpeningLabelMapFilter<ImageToShapeLabelMapFilterType::OutputImageType>                        ShapeOpeningLabelMapFilterType;

public:
    SegmentNecroticRegionOnDM();
    virtual ~SegmentNecroticRegionOnDM();

    void Update();

private:
    void ParseParameterContext();
    void WriteDataSetSummary();
    void WriteLogFile(std::string timeStamp);

	QString m_fullFilenameChannelDPPIV;
    QFileInfo m_infoFullFilenameChannelDPPIV;
    std::string m_pathChannelDPPIV;
    std::string m_filenameChannelDPPIV;
    std::string m_fileExtensionChannelDPPIV;

    QString m_fullFilenameChannelDM;
    QFileInfo m_infoFullFilenameChannelDM;
    std::string m_pathChannelDM;
    std::string m_filenameChannelDM;
    std::string m_fileExtensionChannelDM;

    int m_thresholdingMode;                //0 - adaptive otsu, 1 - normal otsu, 2 - manual
    CScalarVoImageType::SizeType m_adapOtsuRadius;
    int m_adapOtsuSamplePoints;
	CScalarPixelType m_otsuThreshold;
	CScalarPixelType m_lowerThreshold;
    CScalarPixelType m_upperThreshold;

    StructuringElementType m_erode1StructuringElement;
    StructuringElementType m_dilate1StructuringElement;
    StructuringElementType m_erode2StructuringElement;
    StructuringElementType m_dilate2StructuringElement;

    unsigned int m_minimalNecroticRegionSize;

    float m_overlayOpacity;

    bool m_saveEverything;                   //1 - save everything, 0 - save only essentials
    std::string m_logFilenameSave;
    std::string m_filenameSave;
    std::string m_finalSaveSuffixes[3];
};

#endif /* SEGMENTNECROTICREGIONONDM_H_ */
