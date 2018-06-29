///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeLobuleFilter.h                                                //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-03-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ANALYZELOBULEFILTER_H_
#define ANALYZELOBULEFILTER_H_

#include "BasePipeline.h"
#include "Lobule2DHandlingFilter.h"

#include <itkAddImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkDivideImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkLabelImageToShapeLabelMapFilter.h>
#include <itkLineIterator.h>
#include <itkMaskImageFilter.h>
#include <itkMaskNegatedImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMorphologicalWatershedImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkTileImageFilter.h>

#include <vtkMutableUndirectedGraph.h>
#include <vtkSmartPointer.h>

#include <string>
#include <vector>


//struct LobuleAnalysisContainer
//{
//    LobuleAnalysisContainer() {
//        label = 0;
//        volume = -1;
//        numCV = 0;
//        numPV = 0;
//    }
//    unsigned long int label;
//    double volume;
//    unsigned int numCV;
//    unsigned int numPV;
//};


class AnalyzeLobuleFilter : public BasePipeline<3>, public Lobule2DHandlingFilter
{
public:
    const static int ImageDimension = 3;

protected:
//    typedef itk::Vector<double, ImageDimension>                                        VectorType;
    typedef typename BasePipeline<ImageDimension>::CScalarPixelType                    CScalarPixelType;
    typedef typename BasePipeline<ImageDimension>::FScalarPixelType                    FScalarPixelType;
    typedef typename BasePipeline<ImageDimension>::LScalarPixelType                    LScalarPixelType;

    typedef typename BasePipeline<ImageDimension>::CScalarVoImageType                  CScalarImageType;
    typedef typename BasePipeline<ImageDimension>::FScalarVoImageType                  FScalarImageType;
    typedef typename BasePipeline<ImageDimension>::LScalarVoImageType                  LScalarImageType;
    typedef itk::LabelMap< itk::ShapeLabelObject<CScalarPixelType, ImageDimension> >   LabelMapType;

    typedef typename BasePipeline<ImageDimension>::CScalarImagePointerType             CScalarImagePointerType;
    typedef typename BasePipeline<ImageDimension>::FScalarImagePointerType             FScalarImagePointerType;
    typedef typename BasePipeline<ImageDimension>::LScalarImagePointerType             LScalarImagePointerType;
    typedef typename LabelMapType::Pointer                                             LabelMapPointerType;

    typedef typename CScalarImageType::IndexType                                       IndexType;
    typedef typename IndexType::IndexValueType                                         IndexValueType;
    typedef typename CScalarImageType::SizeType                                        SizeType;
    typedef typename CScalarImageType::RegionType                                      RegionType;
    typedef typename CScalarImageType::SpacingType                                     SpacingType;
    typedef typename LabelMapType::LabelObjectType                                     LabelObjectType;
    typedef itk::LineIterator<CScalarImageType>                                        CLineIteratorType;
    typedef itk::LineIterator<FScalarImageType>                                        FLineIteratorType;

    typedef itk::BinaryBallStructuringElement<CScalarPixelType, ImageDimension>        StructuringElementType;

    typedef typename BasePipeline<ImageDimension>::ScalarVoReaderType                  CScalarImageReaderType;
    typedef typename BasePipeline<ImageDimension>::ScalarVoWriterType                  CScalarImageWriterType;

    typedef itk::AddImageFilter<CScalarImageType, CScalarImageType, CScalarImageType>                   AddCImageFilterType;
    typedef itk::AddImageFilter<FScalarImageType, FScalarImageType, FScalarImageType>                   AddFImageFilterType;
    typedef itk::BinaryDilateImageFilter<CScalarImageType, CScalarImageType, StructuringElementType>    DilateImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarImageType, LabelMapType>                       BinaryImageToShapeLabelMapFilterType;
    typedef itk::SignedMaurerDistanceMapImageFilter<CScalarImageType, FScalarImageType>                 SignedMaurerDistanceMapImageFilterType;
    typedef itk::DivideImageFilter<FScalarImageType, FScalarImageType, FScalarImageType>                DivideImageFilterType;
    typedef itk::LabelImageToShapeLabelMapFilter<CScalarImageType, LabelMapType>                        LabelImageToShapeLabelMapFilterType;
    typedef itk::LabelMapToLabelImageFilter<LabelMapType, CScalarImageType>                             LabelMapToLabelImageFilterType;
    typedef itk::MaskImageFilter<FScalarImageType, CScalarImageType>                                    MaskImageFilterType1;
    typedef itk::MaskImageFilter<CScalarImageType, CScalarImageType>                                    MaskImageFilterType2;
    typedef itk::MaskNegatedImageFilter<FScalarImageType, CScalarImageType>                             MaskNegImageFilterType1;
    typedef itk::MaskNegatedImageFilter<CScalarImageType, CScalarImageType>                             MaskNegImageFilterType2;
    typedef itk::MinimumMaximumImageCalculator<CScalarImageType>                                        CImageCalculatorFilterType;
    typedef itk::MinimumMaximumImageCalculator<FScalarImageType>                                        FImageCalculatorFilterType;
    typedef itk::MorphologicalWatershedImageFilter<FScalarImageType, CScalarImageType>                  MorphoWatershedImageFilterType;
    typedef itk::RescaleIntensityImageFilter<FScalarImageType, CScalarImageType>                        RescaleImageFilterType;
    typedef itk::TileImageFilter<CScalar2DImageType, CScalarImageType>                                  TileImageFilterType;
    typedef itk::BinaryThresholdImageFilter<CScalarImageType, CScalarImageType>                         ThresholdFilterType1;
    typedef itk::BinaryThresholdImageFilter<FScalarImageType, CScalarImageType>                         ThresholdFilterType2;
    typedef itk::ExtractImageFilter<CScalarImageType, CScalar2DImageType>                               ExtractImageFilterType;

public:
    AnalyzeLobuleFilter();
    virtual ~AnalyzeLobuleFilter();

    void Update();

protected:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);

    void CollectBasicImageInformation();
    void CollectBasicLobuleInformation();

    void WriteBasicInformationFile();
//    void WriteLobuleInformationFile(std::string filenamePostfix, std::vector<LobuleAnalysisContainer> lobules);
    void WriteGraph(vtkSmartPointer<vtkMutableUndirectedGraph> graph, std::string filename);

    void ComputeRelativeLobularPositions(FScalarImagePointerType cvDistMap, FScalarImagePointerType lobuleDistMap);
    void Decompose3DImageInto2DImages(CScalarImagePointerType inImage, std::vector<int> slicesToCollect, std::map<int, CScalar2DImagePointerType> &outSlices);

    void PrepareVeinData(std::string loadFilename, std::string path, std::string saveFilename, std::string ext, CScalarImagePointerType &bin, LabelMapPointerType &labelMap,
            std::map<int, std::map<int, CScalar2DIndexType> > &veinCentersPerSlice, std::map<vtkIdType, int> &vertexIdToVeinState, int veinState);
    void ComputeVeinCatchmentBasins(std::string path, std::string suffix, std::string ext, CScalarImagePointerType veinBin, LabelMapPointerType veinLabelMap,
            CScalarImagePointerType &veinCatchLabelImage, FScalarImagePointerType &veinDistmap);
    void SyncVeinAndCatchmentLabels(LabelMapPointerType veinLabelMap, LabelMapPointerType &catchmentLabelMap);
    void CorrectVeinCatchmentBasin(FScalarImagePointerType &distanceMap);

    void EstablishVeinAxes();
    bool IsAxisValid(int axisType, vtkIdType start, vtkIdType end, std::map<int, CScalarImagePointerType> cvMasks, std::map<int, CScalarImagePointerType> pvMasks);

    void PrepareVeinAxesData(std::map<int, CScalar2DImagePointerType> &veinAxesImages, std::string path, std::string ext);
    void ComputeVeinAxesCatchmentBasins(CScalarImagePointerType veinAxesImage, std::string path, std::string ext);
    void SyncVeinAxesAndCatchmentLabels(LabelMapPointerType &catchmentLabelMap);
    void CorrectVeinAxisCatchmentBasin(FScalarImagePointerType &distanceMap, vtkIdType veinStartId, vtkIdType veinEndId);

    void ComputeVeinCenterPoints(int zStepWidth, CScalarImagePointerType veinImage, CScalarImagePointerType veinLabelImage, vtkSmartPointer<vtkMutableUndirectedGraph> &veinGraph,
            std::map<int, std::map<int, CScalar2DIndexType> > &veinCentersPerSlice, std::map<vtkIdType, int> &vertexIdToVeinState, int veinState);

    std::string mDataSetID;
    std::string mPath;
    std::string mFilename;
    std::string mFilenameExtension;

    std::string mDataSetFullFilenameCV;
    std::string mDataSetFullFilenamePV;

    QFileInfo mInfoFullFilename;

    bool mCVFileExists;
    bool mPVFileExists;

    std::map<int, std::map<int, CScalar2DIndexType> > mCVeinCentersPerSlice;       //slice nr // vein id // central index
    std::map<int, std::map<int, CScalar2DIndexType> > mPVeinCentersPerSlice;       //slice nr // vein id // central index

    vtkSmartPointer<vtkMutableUndirectedGraph> mVeinSkeletonGraph;

    CScalarImagePointerType mCVBin;
    CScalarImagePointerType mPVBin;

    LabelMapPointerType mCVLabelMap;
    LabelMapPointerType mPVLabelMap;

    CScalarImagePointerType mCVCatchLabelImage;
    CScalarImagePointerType mPVCatchLabelImage;

    FScalarImagePointerType mCVDistmap;
    FScalarImagePointerType mPVDistmap;

    vtkSmartPointer<vtkMutableUndirectedGraph> mVeinAxesGraph;
    std::map<vtkIdType, vtkIdType> mAxisIdToVeinStartId;
    std::map<vtkIdType, vtkIdType> mAxisIdToVeinEndId;


    FScalarImagePointerType mRelativeLobularPositionMap;

    double mMaxCorrectionDistanceToBorder;
    double mMinAxisDistToVein;
    double mCenterRUSCRadius;

    typename CScalarImageType::SpacingType mSpacing;
    typename CScalarImageType::SizeType mSize;

    long double mVoxelVolume;
    long double mDatasetVolume;
    long double mCVVolume;
    long double mPVVolume;

    bool mSaveEverything;
};

#endif /* ANALYZELOBULEFILTER_H_ */
