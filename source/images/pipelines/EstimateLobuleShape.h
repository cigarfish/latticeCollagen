///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  EstimateLobuleShape.h                                                //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-08-05                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ESTIMATELOBULESHAPE_H_
#define ESTIMATELOBULESHAPE_H_

#include "BasePipeline.h"

#include "itkIndex.h"

#include "itkAddImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkLabelContourImageFilter.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkMaskNegatedImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"

#include <vtkSmartPointer.h>
#include <vtkUndirectedGraph.h>
#include <vtkVertexListIterator.h>
#include <vtkGraphWriter.h>

#include "../filters/convertFilters/LabelImageToGraphFilter.h"
#include "../filters/convertFilters/SkeletonImageToGraphFilter.h"
#include "../filters/graphFilters/ResampleUndirectedGraphFilter.h"
#include "../filters/imageFilters/LabelShapeKeepNObjectsNextToPosImageFilter.h"
#include "../tools/GraphAnnotationHelper.h"


struct LobuleAnalysisContainer
{
    LobuleAnalysisContainer() {
        mLabel = -1;
        mVolume = -1;
        mRoundness = -1;
        mElongation = -1;
        mNumPxlOnBorder = -1;

        mPos[0] = 0;
        mPos[1] = 0;
        mPos[2] = 0;

        mNumberCV = -1;
        mNumberPV = -1;
    }

    unsigned long int mLabel;
    double mVolume;
    double mRoundness;
    double mElongation;
    int mNumPxlOnBorder;

    int mPos[3];

    int mNumberCV;
    int mNumberPV;
};


template< unsigned int VImageDimension > class EstimateLobuleShape : public BasePipeline<VImageDimension>
{
public:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

private:
    typedef typename BasePipeline<VImageDimension>::CScalarPixelType    CScalarPixelType;
    typedef typename BasePipeline<VImageDimension>::FScalarPixelType    FScalarPixelType;
    typedef unsigned long                                               LScalarPixelType;
    typedef unsigned int                                                IScalarPixelType;
    typedef unsigned short                                              SScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType>                             RGBPixelType;

    typedef typename BasePipeline<VImageDimension>::CScalarVoImageType  CScalarImageType;
    typedef typename BasePipeline<VImageDimension>::FScalarVoImageType  FScalarImageType;
    typedef itk::Image<LScalarPixelType, ImageDimension>                LScalarImageType;
    typedef itk::Image<IScalarPixelType, ImageDimension>                IScalarImageType;
    typedef itk::Image<SScalarPixelType, ImageDimension>                SScalarImageType;
    typedef itk::Image<RGBPixelType, ImageDimension>                    RGBImageType;

    typedef typename BasePipeline<VImageDimension>::ScalarVoReaderType  CScalarImageReaderType;
    typedef typename BasePipeline<VImageDimension>::ScalarVoWriterType  CScalarImageWriterType;
    typedef itk::ImageFileReader<SScalarImageType>                      SScalarImageReaderType;
    typedef itk::ImageFileWriter<SScalarImageType>                      SScalarImageWriterType;
    typedef itk::ImageFileWriter<RGBImageType>                          RGBImageWriterType;

    typedef itk::AddImageFilter<CScalarImageType, CScalarImageType, CScalarImageType>           AddCImageFilterType;
    typedef itk::AddImageFilter<FScalarImageType, FScalarImageType, FScalarImageType>           AddFImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarImageType>                             ImageToShapeLabelMapFilterType;
    typedef typename ImageToShapeLabelMapFilterType::OutputImageType                            LabelMapType;
    typedef typename LabelMapType::LabelObjectType                                              LabelObjectType;
    typedef itk::BinaryThresholdImageFilter<CScalarImageType, CScalarImageType>                 ThresholdCCFilterType;
    typedef itk::BinaryThresholdImageFilter<LScalarImageType, CScalarImageType>                 ThresholdLCFilterType;
    typedef itk::DivideImageFilter<FScalarImageType, FScalarImageType, FScalarImageType>        DivideImageFilterType;
    typedef itk::IntensityWindowingImageFilter<FScalarImageType, FScalarImageType>              IntensityWindowingImageFilterType;
    typedef itk::LabelMapToLabelImageFilter<LabelMapType, LScalarImageType>                     LabelMapToLabelImageFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarImageType, LScalarImageType, RGBImageType>      LabelOverlayImageFilterType;
    typedef itk::MaskImageFilter<CScalarImageType, CScalarImageType, CScalarImageType>          MaskImageFilterType;
    typedef itk::MaskNegatedImageFilter<CScalarImageType, CScalarImageType, CScalarImageType>   MaskNegatedImageFilterType;
    typedef itk::MinimumMaximumImageCalculator<FScalarImageType>                                ImageCalculatorFilterType;
    typedef itk::MorphologicalWatershedImageFilter<FScalarImageType, LScalarImageType>          MorphoWatershedImageFilterType;
    typedef itk::MultiplyImageFilter<FScalarImageType, FScalarImageType, FScalarImageType>      MultiplyImageFilterType;
    typedef itk::RescaleIntensityImageFilter<FScalarImageType, CScalarImageType>                RescaleFCImageFilterType;
    typedef itk::SignedMaurerDistanceMapImageFilter<CScalarImageType, FScalarImageType>         SignedMaurerDistanceMapImageFilterType;

//    typedef LabelImageToGraphFilter<VImageDimension>                                            LabelImageToGraphFilterType;

public:
    EstimateLobuleShape();
    virtual ~EstimateLobuleShape();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteDataSetSummary();

    void OrganizeLobulesAndVeins();
    void ClassifyLobules();

    void AnalyzeLobules();
    void WriteAnalysisFile();


    std::string mDataSetID;
    std::string mPath;
    std::string mFilenameExtension;

    QFileInfo mInfoFullFilename;

    std::string mFilenameRawData;
    std::string mFilenameCVBin;
    std::string mFilenamePVBin;
    std::string mFilenameTissueBin;

    int mCentralVeinThreshold;
    int mPortalVeinThreshold;
    int mTissueThreshold;
    int mVoidThreshold;

    bool mTissueFileExists;

    typename CScalarImageType::SpacingType mSpacing;

    double mPortalVeinWeight;
    double mCentralVeinWeight;

    double mMinimalLobuleDiameter;
    double mMinimalLobuleVolume;
    double mMaximalLobuleDiameter;
    double mMaximalLobuleVolume;

    itk::SmartPointer<LabelMapType> mpLobuleLabelMap;
    vtkSmartPointer<vtkMutableUndirectedGraph> mLobuleGraph;

    double mWatershedFloodLevel;

    double mOverlayOpacity;

    std::map<unsigned long int, LobuleAnalysisContainer> mLobules;

    bool mWithAnalysis;

    bool mSaveEverything;                   //1 - save everything, 0 - save only essentials
    std::string mLogFilenameSave;
    std::string mFilenameSave;
    std::string mSaveSuffixesForFinals[4];
};

#include "EstimateLobuleShape.tpp"

#endif /* ESTIMATELOBULESHAPE_H_ */
