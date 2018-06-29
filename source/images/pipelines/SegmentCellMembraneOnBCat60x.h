///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentCellMembraneOnBCat60x.h                                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-09-02                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef SEGMENTCELLMEMBRANEONBCAT60X_H_
#define SEGMENTCELLMEMBRANEONBCAT60X_H_

#include "BasePipeline.h"

#include <math.h>

#include "itkRGBPixel.h"

//#include "itkLabelMap.h"

#include "itkBinaryBallStructuringElement.h"
//#include "itkMath.h"

#include "itkAddImageFilter.h"
//#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryThresholdImageFilter.h"
//#include "itkGrayscaleDilateImageFilter.h"
//#include "itkGrayscaleErodeImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkMaskNegatedImageFilter.h"
//#include "itkMedianImageFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"
//#include "itkNearestNeighborExtrapolateImageFunction.h"
//#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
//#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"


#include "../filters/imageFilters/itkAdaptiveOtsuThresholdingImageFilter.h"


class SegmentCellMembraneOnBCat60x : public BasePipeline<3>
{
    typedef float                                       FScalarPixelType;
    typedef unsigned long                               LScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType>             RGBPixelType;

    typedef itk::Image<FScalarPixelType, 3>             FScalarVoImageType;
    typedef itk::Image<LScalarPixelType, 3>             LScalarVoImageType;
    typedef itk::Image<RGBPixelType,3>                  RGBVoImageType;

    typedef itk::ImageFileWriter<RGBVoImageType>        RGBVoWriterType;

    typedef itk::BinaryBallStructuringElement<CScalarVoImageType::PixelType, 3>     StructuringElementType;

    typedef itk::AdaptiveOtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                   AdaptiveOtsuThresholdImageFilterType;
    typedef itk::AddImageFilter<CScalarVoImageType>                                                                         AddImageFilterType;
    typedef itk::OtsuThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                           OtsuThresholdImageFilterType;
    typedef itk::InvertIntensityImageFilter<CScalarVoImageType>                                                             InvertIntensityImageFilterType;
    typedef itk::BinaryThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                         ThresholdFilterType;
    typedef itk::BinaryThresholdImageFilter<LScalarVoImageType, CScalarVoImageType>                                         ThresholdLScalarFilterType;
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter<CScalarVoImageType>                                            HoleFillingImageFilterType;
    typedef itk::BinaryMorphologicalClosingImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>      ClosingImageFilterType;
    typedef itk::BinaryMorphologicalOpeningImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>      OpeningImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                                                       ImageToShapeLabelMapFilterType;
    typedef itk::ShapeOpeningLabelMapFilter<ImageToShapeLabelMapFilterType::OutputImageType>                                ShapeOpeningLabelMapFilterType;
    typedef itk::LabelMapToLabelImageFilter<ShapeOpeningLabelMapFilterType::OutputImageType, LScalarVoImageType>            LabelMapToLabelImageFilterType;
    typedef itk::SignedMaurerDistanceMapImageFilter<CScalarVoImageType, FScalarVoImageType>                                 SignedMaurerDistanceMapImageFilterType;
    typedef itk::RescaleIntensityImageFilter<FScalarVoImageType, CScalarVoImageType>                                        RescaleImageFilterType;
    typedef itk::MorphologicalWatershedImageFilter<FScalarVoImageType, LScalarVoImageType>                                  MorphoWatershedImageFilterType;
    typedef itk::MaskNegatedImageFilter<LScalarVoImageType, CScalarVoImageType, LScalarVoImageType>                         MaskImageFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarVoImageType, LScalarVoImageType, RGBVoImageType>                            LabelOverlayImageFilterType;
//    typedef itk::LabelImageToShapeLabelMapFilter<LScalarVoImageType>                                                        LabelImageToShapeLabelMapFilterType;

//    typedef itk::GrayscaleDilateImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>                 GrayscaleDilateImageFilterType;
//    typedef itk::GrayscaleErodeImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>                  GrayscaleErodeImageFilterType;
//    typedef itk::LabelMapToLabelImageFilter<LabelImageToShapeLabelMapFilterType::OutputImageType, CScalarVoImageType>       LabelMapToLabelImageFilterType2;

public:
    SegmentCellMembraneOnBCat60x();
    virtual ~SegmentCellMembraneOnBCat60x();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteDataSetSummary();

    void BuildCellNucleiAlignmentMaps(LScalarVoImageType::Pointer nucleiLabelImage, LScalarVoImageType::Pointer cellLabelImage,
            itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> nucleiLabelMap, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap);
    void RemoveCellsAtDatasetBorder(LScalarVoImageType::Pointer cellLabelImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithOneNucleusLabelMap,
            itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithTwoNucleiLabelMap, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithMoreNucleiLabelMap);
    void WriteNumNucleiOutlineFiles(itk::SmartPointer<CScalarVoImageType> dppivImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithOneNucleusLabelMap,
            itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithTwoNucleiLabelMap, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithMoreNucleiLabelMap);

    QString m_fullFilenameChannelBCat;
    QFileInfo m_infoFullFilenameChannelDAPI;
    std::string m_pathChannelBCat;
    std::string m_filenameChannelBCat;
    std::string m_fileExtensionChannelBCat;

    QString m_fullFilenameSegmentationSin;
    QString m_fullFilenameSegmentationNuc;

    QString m_fullFilenameSegmentationCV;
    QString m_fullFilenameSegmentationPV;

    CScalarVoImageType::SpacingType m_spacing;

    int m_greyscaleOpeningRadius;

    int m_thresholdingMode;                //0 - adaptive otsu, 1 - normal otsu, 2 - manual
    CScalarVoImageType::SizeType m_adapOtsuRadius;
    int m_adapOtsuSamplePoints;
    CScalarPixelType m_otsuThreshold;
    CScalarPixelType m_lowerThreshold;
    CScalarPixelType m_upperThreshold;

    CScalarVoImageType::SizeType m_holeFillingNeighborhoodRadius;
    unsigned int m_holeFillingMajThreshold;

    CScalarVoImageType::SizeType m_invHoleFillingNeighborhoodRadius;
    unsigned int m_invHoleFillingMajThreshold;

    StructuringElementType m_closingStructuringElement;
    StructuringElementType m_openingStructuringElement1;
    StructuringElementType m_openingStructuringElement2;

    double m_minimalCellDiameter;
    double m_minimalCellVolume;

    float m_floodLevel;

    std::multimap<unsigned long int, std::pair<unsigned long int, double> > m_cellToNuclei;
    std::vector<unsigned long> m_labelsWithZeroNuclei, m_labelsWithOneNucleus, m_labelsWithTwoNuclei, m_labelsWithMoreNuclei;

    float m_overlayOpacity;

    bool m_saveEverything;                   //1 - save everything, 0 - save only essentials
    std::string m_logFilenameSave;
    std::string m_filenameSave;
    std::string m_saveSuffixesForFinals[6];
};

#endif /* SEGMENTCELLMEMBRANEONBCAT60X_H_ */
