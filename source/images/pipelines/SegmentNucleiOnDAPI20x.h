///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentNucleiOnDAPI20x.h				                             //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-04-18                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SEGMENTNUCLEIONDAPI20X_H_
#define SEGMENTNUCLEIONDAPI20X_H_

#include "BasePipeline.h"

#include <math.h>

#include "itkRGBPixel.h"

#include "itkLabelMap.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkMath.h"

#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include "../filters/imageFilters/itkAdaptiveOtsuThresholdingImageFilter.h"


class SegmentNucleiOnDAPI20x : public BasePipeline<3>			// dervied from BasePipeline: all pipelines should inherit from this class
{
	//Type definitions of all Pixel Types, Image Types and Filter Types that will be deployed by the pipeline
    typedef float                                       FScalarPixelType;
    typedef unsigned short                              SScalarPixelType;
    typedef unsigned int                                IScalarPixelType;
	typedef itk::RGBPixel<CScalarPixelType>             RGBPixelType;

	typedef itk::Image<FScalarPixelType, 3>             FScalarVoImageType;
	typedef itk::Image<SScalarPixelType, 3>             SScalarVoImageType;
	typedef itk::Image<IScalarPixelType, 3>             IScalarVoImageType;
	typedef itk::Image<RGBPixelType,3>                  RGBVoImageType;

	typedef itk::ImageFileWriter<SScalarVoImageType>    SScalarVoWriterType;
	typedef itk::ImageFileWriter<RGBVoImageType>        RGBVoWriterType;

	typedef itk::BinaryBallStructuringElement<CScalarVoImageType::PixelType, 3>     StructuringElementType;

	typedef itk::AdaptiveOtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                   AdaptiveOtsuThresholdImageFilterType;
	typedef itk::OrImageFilter<CScalarVoImageType>                                                                          OrImageFilterType;
	typedef itk::BinaryMorphologicalOpeningImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>		OpeningImageFilterType;
	typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                                               		ImageToShapeLabelMapFilterType;
	typedef itk::BinaryThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                 		ThresholdFilterType;
	typedef itk::SignedMaurerDistanceMapImageFilter<CScalarVoImageType, FScalarVoImageType>                                 MaurerDistanceMapImageFilterType;
	typedef itk::InvertIntensityImageFilter<CScalarVoImageType>                                                             InvertIntensityImageFilterType;
	typedef itk::InvertIntensityImageFilter<IScalarVoImageType>                                                             IInvertIntensityImageFilterType;
	typedef itk::MaskImageFilter<IScalarVoImageType, CScalarVoImageType, IScalarVoImageType>                                MaskImageFilterType;
	typedef itk::MedianImageFilter<CScalarVoImageType, CScalarVoImageType >                                                 MedianImageFilterType;
	typedef itk::MorphologicalWatershedImageFilter<FScalarVoImageType, IScalarVoImageType>                                  MorphoWatershedImageFilterType;
	typedef itk::LabelImageToShapeLabelMapFilter<MorphoWatershedImageFilterType::OutputImageType>                           LabelImageToShapeLabelMapFilterType;
	typedef itk::LabelMapToLabelImageFilter<ImageToShapeLabelMapFilterType::OutputImageType, CScalarVoImageType>   			LabelMapToLabelImageFilterType1;
	typedef itk::LabelMapToLabelImageFilter<LabelImageToShapeLabelMapFilterType::OutputImageType, CScalarVoImageType>       LabelMapToLabelImageFilterType2;
	typedef itk::LabelOverlayImageFilter<CScalarVoImageType, CScalarVoImageType, RGBVoImageType>                    		LabelOverlayImageFilterType;
	typedef itk::OtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>											OtsuThresholdImageFilterType;
	typedef itk::RescaleIntensityImageFilter<FScalarVoImageType, SScalarVoImageType>                                        RescaleImageFilterType;
	typedef itk::ShapeOpeningLabelMapFilter<ImageToShapeLabelMapFilterType::OutputImageType>                        		ShapeOpeningLabelMapFilterType1;
	typedef itk::ShapeOpeningLabelMapFilter<LabelImageToShapeLabelMapFilterType::OutputImageType>                           ShapeOpeningLabelMapFilterType2;
	typedef itk::ShapeOpeningLabelMapFilter<ShapeOpeningLabelMapFilterType2::OutputImageType>                               ShapeOpeningLabelMapFilterType3;
	typedef itk::SubtractImageFilter<CScalarVoImageType, CScalarVoImageType>                                                SubtractImageFilterType;
	typedef itk::VotingBinaryIterativeHoleFillingImageFilter<CScalarVoImageType>                                    		HoleFillingImageFilterType;


public:
	SegmentNucleiOnDAPI20x();
	virtual ~SegmentNucleiOnDAPI20x();

	//Call the Update method to invoke pipeline execution
	void Update();	

private:
	//This method is used to extract all necessary parameters from the parameter tree
	void ParseParameterContext();
	//Writes a log file of used parameters, using the parameter tree
    void WriteLogFile(std::string timeStamp);
	//Writes location of pipeline's end results (segmentation image, overlay image) to a data set summary file, that can be used by subsequent pipelines to locate these results 
    void WriteDataSetSummary();

	//Paramters used by the pipeline
	QString m_fullFilenameChannelDAPI;
    QFileInfo m_infoFullFilenameChannelDAPI;
    std::string m_pathChannelDAPI;
    std::string m_filenameChannelDAPI;
    std::string m_fileExtensionChannelDAPI;

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

    int m_medianRadius;

	int m_thresholdingMode;                //0 - adaptive otsu, 1 - normal otsu, 2 - manual
	CScalarVoImageType::SizeType m_adapOtsuRadius;
	int m_adapOtsuSamplePoints;
	CScalarPixelType m_otsuThreshold;
    CScalarPixelType m_lowerThreshold;
    CScalarPixelType m_upperThreshold;

    CScalarVoImageType::SizeType m_holeFillingNeighborhoodRadius;
    unsigned int m_holeFillingMajThreshold;

    CScalarVoImageType::SizeType m_openingNeighborhoodRadius;

    double m_minimalNucleiRadius;
    double m_minimalNucleiVolume;

    //TODO: maximal radius/volume filter not implemented so far

    double m_maximalNucleiRadius;
    double m_maximalNucleiVolume;

    float m_floodLevel;

    float m_overlayOpacity;

    bool m_saveEverything;                   //1 - save everything, 0 - save only essentials
    std::string m_logFilenameSave;
    std::string m_filenameSave;
    std::string m_saveSuffixes[7];
};

#endif /* SEGMENTNUCLEIONDAPI20X_H_ */
