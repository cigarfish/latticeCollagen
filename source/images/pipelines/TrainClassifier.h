///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  TrainClassifier.h                                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-08-01                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef TRAINCLASSIFIER_H_
#define TRAINCLASSIFIER_H_

#include "BasePipeline.h"

#include "ObjectBasedSegmentation.h"

#include "itkListSample.h"
#include "itkMembershipSample.h"
#include "itkVector.h"

#include <vtkUndirectedGraph.h>
#include <vtkSmartPointer.h>


class SampleTypePointer;


template< unsigned int VImageDimension > class TrainClassifier : public BasePipeline<VImageDimension>
{
public:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

private:
    static const unsigned int MaxNumberFeaturesPerPair = 100;

    typedef typename BasePipeline<VImageDimension>::CScalarPixelType    CScalarPixelType;
    typedef typename BasePipeline<VImageDimension>::LScalarPixelType    LScalarPixelType;
    typedef unsigned int                                                IScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType>                             CRGBPixelType;

    typedef typename BasePipeline<VImageDimension>::CScalarVoImageType  CScalarImageType;
    typedef typename BasePipeline<VImageDimension>::LScalarVoImageType  LScalarImageType;
    typedef itk::Image<IScalarPixelType, ImageDimension>                IScalarImageType;
    typedef itk::Image<CRGBPixelType, ImageDimension>                   CRGBImageType;

    typedef typename BasePipeline<VImageDimension>::ScalarVoReaderType  CScalarReaderType;
    typedef itk::ImageFileReader<CRGBImageType>                         CRGBReaderType;
    typedef itk::ImageFileReader<IScalarImageType>                      IScalarReaderType;

    typedef itk::CastImageFilter<IScalarImageType, LScalarImageType>                        CastILImageFilterType;
    typedef itk::MinimumMaximumImageCalculator<LScalarImageType>                            LMinMaxCalculatorType;
    typedef itk::MinimumMaximumImageCalculator<IScalarImageType>                            IMinMaxCalculatorType;
    typedef itk::LabelImageToLabelMapGraphFilter<CRGBImageType, LScalarImageType>           LabelImageToLabelMapGraphFilterType;
    typedef typename LabelImageToLabelMapGraphFilterType::OutputImageType                   LabelMapGraphType;
    typedef typename LabelMapGraphType::LabelObjectType                                     FeatureLabelObjectType;
    typedef typename LabelMapGraphType::Pointer                                             LabelMapGraphPointerType;
    typedef ObjectBasedSegmentation<ImageDimension>                                         ObjectBasedSegmentationType;

    typedef typename LabelMapGraphType::UnaryFeatures                                       UnaryFeatureType;
    typedef typename LabelMapGraphType::BinaryFeatures                                      BinaryFeatureType;
    typedef typename ObjectBasedSegmentation<ImageDimension>::Class                         ClassType;

    typedef itk::Vector<double, MaxNumberFeaturesPerPair>       MeasurementVectorType;
    typedef itk::Statistics::ListSample<MeasurementVectorType>  SampleType;
    typedef SampleType::Pointer                                 SampleTypePointer;

public:
    TrainClassifier();
    virtual ~TrainClassifier();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);

    void Init();
    void RetrieveVerticesAndEdgesFromMaskImage(std::string filename, std::set<unsigned long> &labelObjects, std::vector< std::pair<unsigned long,unsigned long> > &labelObjectPairs);

    void ExtractTrainingData(std::set<unsigned long> trainingLabelObjects, std::vector< std::pair<unsigned long,unsigned long> > trainingLabelObjectPairs,
            SampleTypePointer &unarySamples, SampleTypePointer &binarySamples);

    void SaveToTrainingDatabase(std::set<unsigned long> trainingLabelObjects, std::vector< std::pair<unsigned long,unsigned long> > trainingLabelObjectPairs,
            SampleTypePointer unarySamples, SampleTypePointer binarySamples);

    void LoadFromTrainingDatabase(std::map<int, ClassType> &idToClass, SampleTypePointer &sample);

    void AddTrainingDataToTrainingDatabase();


    //Paramters used by the pipeline
    std::string mPath;

    std::string mRawImageFullFilename;
    std::string mRawImagePath;
    std::string mRawImageFilename;
    std::string mRawImageFileExtension;

    std::string mObjectImageFullFilename;
    std::string mObjectImagePath;
    std::string mObjectImageFilename;
    std::string mObjectImageFileExtension;

    std::string mObjectGraphFullFilename;

    std::string mTrainingDatasetFilename;
    std::string mTrainingDatabaseFilename;

    typename CScalarImageType::SpacingType mVoxelSpacing;

    typename LabelMapGraphType::Pointer mLabelMapGraph;
    itk::SmartPointer<CRGBImageType> mOriginalImage;
    itk::SmartPointer<LScalarImageType> mObjectImage;
    itk::SmartPointer<LScalarImageType> mObjectImage2;
    vtkSmartPointer<vtkUndirectedGraph> mObjectGraph;

    bool mWithIncrementalTraining;
    ClassType mTrainingClass;

    bool* mWithFeature;
};

#include "TrainClassifier.tpp"

#endif /* TRAINCLASSIFIER_H_ */
