///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  EstimateCellShape.h                                                  //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-06-04                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ESTIMATECELLSHAPE_H_
#define ESTIMATECELLSHAPE_H_

#include "BasePipeline.h"

#include "itkAddImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkIndex.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkLabelContourImageFilter.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMaskNegatedImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkRegionalMaximaImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"

#include <vtkSmartPointer.h>
#include <vtkUndirectedGraph.h>
#include <vtkVertexListIterator.h>
#include <vtkGraphWriter.h>

#include "../filters/convertFilters/SkeletonImageToGraphFilter.h"
#include "../filters/graphFilters/ResampleUndirectedGraphFilter.h"
#include "../filters/imageFilters/LabelShapeKeepNObjectsNextToPosImageFilter.h"
#include "../tools/GraphAnnotationHelper.h"


class EstimateCellShape : public BasePipeline<3>
{
    typedef bool                                        BScalarPixelType;
    typedef long                                        LScalarPixelType;
    typedef unsigned int                                IScalarPixelType;
    typedef float                                       FScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType>             RGBPixelType;

    typedef itk::Image<BScalarPixelType, 3>             BScalarVoImageType;
    typedef itk::Image<LScalarPixelType, 3>             LScalarVoImageType;
    typedef itk::Image<IScalarPixelType, 3>             IScalarVoImageType;
    typedef itk::Image<FScalarPixelType, 3>             FScalarVoImageType;
    typedef itk::Image<RGBPixelType, 3>                 RGBVoImageType;

    typedef itk::BinaryBallStructuringElement<CScalarVoImageType::PixelType, 3> StructuringElementType;
    typedef itk::ImageRegion<3>                                                 ImageRegionType;

    typedef itk::ImageFileReader<BScalarVoImageType>    BScalarVoReaderType;
    typedef itk::ImageFileWriter<IScalarVoImageType>    IScalarVoWriterType;
    typedef itk::ImageFileWriter<RGBVoImageType>        RGBVoWriterType;

    typedef itk::AddImageFilter<CScalarVoImageType, CScalarVoImageType, CScalarVoImageType>                             AddCImageFilterType;
    typedef itk::AddImageFilter<FScalarVoImageType, FScalarVoImageType, FScalarVoImageType>                             AddFImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                                                   ImageToShapeLabelMapFilterType;
    typedef itk::BinaryMorphologicalOpeningImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>  OpeningImageFilterType;
    typedef itk::BinaryThresholdImageFilter<LScalarVoImageType, CScalarVoImageType>                                     ThresholdFilterType;
    typedef itk::SignedMaurerDistanceMapImageFilter<CScalarVoImageType, FScalarVoImageType>                             SignedMaurerDistanceMapImageFilterType;
    typedef itk::MinimumMaximumImageCalculator<FScalarVoImageType>                                                      ImageCalculatorFilterType;
    typedef itk::IntensityWindowingImageFilter<FScalarVoImageType, FScalarVoImageType>                                  IntensityWindowingImageFilter;
    typedef itk::DivideImageFilter<FScalarVoImageType, FScalarVoImageType, FScalarVoImageType>                          DivideImageFilterType;
    typedef itk::MaskImageFilter<LScalarVoImageType, CScalarVoImageType, LScalarVoImageType>                            MaskImageFilterType;
    typedef itk::MaskNegatedImageFilter<LScalarVoImageType, CScalarVoImageType, LScalarVoImageType>                     MaskNegatedImageFilterType;
    typedef itk::MorphologicalWatershedImageFilter<FScalarVoImageType, LScalarVoImageType>                              MorphoWatershedImageFilterType;
    typedef itk::MultiplyImageFilter<FScalarVoImageType, FScalarVoImageType, FScalarVoImageType>                        MultiplyImageFilterType;
    typedef itk::LabelImageToShapeLabelMapFilter<MorphoWatershedImageFilterType::OutputImageType>                       LabelImageToShapeLabelMapFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarVoImageType, LScalarVoImageType, RGBVoImageType>                        LabelOverlayImageFilterType;
    typedef itk::PasteImageFilter<CScalarVoImageType, CScalarVoImageType>                                               PasteImageFilterType;
    typedef itk::RegionalMaximaImageFilter<FScalarVoImageType, CScalarVoImageType>                                      RegionalMaximaImageFilterType;
    typedef itk::RelabelComponentImageFilter<MorphoWatershedImageFilterType::OutputImageType, LScalarVoImageType>       RelabelImageFilterType;
    typedef itk::RescaleIntensityImageFilter<LScalarVoImageType, CScalarVoImageType>                                    RescaleLImageFilterType;
    typedef itk::RescaleIntensityImageFilter<LScalarVoImageType, IScalarVoImageType>                                    RescaleLIImageFilterType;
    typedef itk::RescaleIntensityImageFilter<FScalarVoImageType, CScalarVoImageType>                                    RescaleImageFilterType;
    typedef itk::RescaleIntensityImageFilter<CScalarVoImageType, FScalarVoImageType>                                    RescaleBackImageFilterType;
    typedef itk::ShapeOpeningLabelMapFilter<ImageToShapeLabelMapFilterType::OutputImageType>                            ShapeOpeningLabelMapFilterType;
    typedef itk::LabelMapToLabelImageFilter<ShapeOpeningLabelMapFilterType::OutputImageType, LScalarVoImageType>        LabelMapToLabelImageFilterType;
    typedef itk::LabelContourImageFilter<LabelMapToLabelImageFilterType::OutputImageType, LScalarVoImageType>           LabelContourImageFilterType;
    typedef itk::LabelShapeKeepNObjectsNextToPosImageFilter<ShapeOpeningLabelMapFilterType::OutputImageType>            KeepObjectsNearMiddleFilterType;

    typedef itk::ShapeOpeningLabelMapFilter<ImageToShapeLabelMapFilterType::OutputImageType>                            ShapeOpeningLabelMapFilterType2;
    typedef itk::LabelMapToLabelImageFilter<ShapeOpeningLabelMapFilterType2::OutputImageType, LScalarVoImageType>       LabelMapToLabelImageFilterType2;
    typedef itk::LabelContourImageFilter<LabelMapToLabelImageFilterType2::OutputImageType, LScalarVoImageType>          LabelContourImageFilterType2;
    typedef itk::LabelImageToLabelMapFilter<LabelContourImageFilterType2::OutputImageType>                              LabelImageToShapeLabelMapFilterType2;

public:
    EstimateCellShape();
    virtual ~EstimateCellShape();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteDataSetSummary();

    void SetupBorderRegions(CScalarVoImageType::Pointer &image);
    void SetBorderPixel(CScalarVoImageType::Pointer &image, CScalarPixelType value);
    void SetBorderPixelAtLocalMaxima(CScalarVoImageType::Pointer &image, FScalarVoImageType::Pointer distMap, CScalarPixelType value);
    void SetBorderPixelAtLocalMaximaExceedingThreshold(CScalarVoImageType::Pointer &image, FScalarVoImageType::Pointer distMap, FScalarPixelType threshold, CScalarPixelType value);
    void RemoveNucleiInContactWithSinusoids(CScalarVoImageType::Pointer &nucImage, CScalarVoImageType::Pointer sinImage);
    void BuildCellNucleiAlignmentMaps(LScalarVoImageType::Pointer nucleiLabelImage, LScalarVoImageType::Pointer cellLabelImage,
            itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> nucleiLabelMap, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap);
    void RemoveCellsAtDatasetBorder(LScalarVoImageType::Pointer cellLabelImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithOneNucleusLabelMap,
            itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithTwoNucleiLabelMap, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithMoreNucleiLabelMap,
            double allowedBorderPixel);
    void WriteNumNucleiOutlineFiles(itk::SmartPointer<CScalarVoImageType> dppivImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithOneNucleusLabelMap,
            itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithTwoNucleiLabelMap, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithMoreNucleiLabelMap);
    void WriteSinusoidGraphToFile();
    void WriteSiunsoidOutlineToFile();

	QFileInfo m_infoFile;
	QString m_filenameDPPIVChannel;
	QString m_filenameNucleiBin;
    QString m_filenameBileBin;
    QString m_filenameSinusoidBin;
    QString m_filenameCVBin;
    QString m_filenamePVBin;
    QString m_filenameSinusoidSkeleton;
    QString m_fullFilenameNecroticRegion;
    
	std::string m_path;
    std::string m_filenameExtension;

	bool m_hasCV;
	bool m_hasPV;
    bool m_hasNR;

    LScalarVoImageType::SpacingType m_spacing;

    double m_XMax;
    double m_YMax;
    double m_ZMax;

    double m_bileSinWeight;
    double m_nucleiWeight;

    double m_minimalCellDiameter;
    double m_minimalCellVolume;

    double m_watershedFloodLevel;

    ImageRegionType mDataSetFaces[6];

    std::multimap<unsigned long int, std::pair<unsigned long int, double> > m_cellToNuclei;
    std::vector<unsigned long> m_labelsWithZeroNuclei, m_labelsWithOneNucleus, m_labelsWithTwoNuclei, m_labelsWithMoreNuclei;

    bool m_withPositionFilter;
    int m_numberCells;
    itk::Point<double, 3> m_pos;

    double m_overlayOpacity;

    bool m_saveEverything;                   //1 - save everything, 0 - save only essentials
    std::string m_logFilenameSave;
    std::string m_filenameSave;
    std::string m_saveSuffixesForFinals[7];
    bool m_writeCellOutlineToFile;
    bool m_writeSinusoidOutlineToFile;
    bool m_writeSinusoidGraphToFile;
};

#endif /* ESTIMATECELLSHAPE_H_ */
