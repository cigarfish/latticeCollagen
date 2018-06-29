///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CompareSegmentations.h                                               //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-04-16                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef COMPARESEGMENTATIONS_H_
#define COMPARESEGMENTATIONS_H_

#include <string>
#include <vector>

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelMap.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkSmartPointer.h"


class CSParameterContext;

struct ObjectPairs
{
    unsigned long mGoldLabel;
    unsigned long mEvalLabel;

    int mOverlapVolume;
    int mGoldExcess;
    int mEvalExcess;

	bool mGoldInner;
	bool mEvalInner;
};

//very basic segmentation comparison:
//pairs are given by: goldstandard label object -> evaluation label object @ goldstandard label object centroid
//computes for these pairs overlap region + excess of goldstandard + evaluation label objects
class CompareSegmentations
{
protected:
    typedef unsigned char                               CScalarPixelType;
    typedef float                                       FScalarPixelType;
    typedef long                                        LScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType>             RGBPixelType;

    typedef itk::Image<CScalarPixelType, 3>             CScalarVoImageType;
    typedef itk::Image<FScalarPixelType, 3>             FScalarVoImageType;
    typedef itk::Image<LScalarPixelType, 3>             LScalarVoImageType;
    typedef itk::Image<RGBPixelType, 3>                 RGBVoImageType;

    typedef itk::ImageFileReader<CScalarVoImageType>    ScalarVoReaderType;
    typedef itk::ImageFileWriter<CScalarVoImageType>    CScalarVoWriterType;
    typedef itk::ImageFileWriter<RGBVoImageType>        RGBVoWriterType;

    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                                               ImageToShapeLabelMapFilterType;
    typedef itk::LabelMapToLabelImageFilter<ImageToShapeLabelMapFilterType::OutputImageType, LScalarVoImageType>    LabelMapToLabelImageFilterType;

public:
    CompareSegmentations();
    virtual ~CompareSegmentations();

    void SetParameterContext(CSParameterContext *paramContext) { mpParamContext = paramContext; };

    void Update();
protected:
    void ParseParameterContext();
    void WriteDataFile();

    CSParameterContext *mpParamContext;

    int mDatasetID;

    std::string mGoldstandardDatasetFullFilename;
    std::string mGoldstandardDatasetPath;
    std::string mGoldstandardDatasetName;
    std::string mGoldstandardDatasetExtension;

    std::string mEvaluationDatasetFullFilename;
    std::string mEvaluationDatasetPath;
    std::string mEvaluationDatasetName;
    std::string mEvaluationDatasetFile;

    std::map<unsigned long, ObjectPairs> mObjectPairs;
};

#endif /* COMPARESEGMENTATIONS_H_ */
