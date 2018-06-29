///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphBaseClassifier.h                                                //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-03-23                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHBASECLASSIFIER_H_
#define GRAPHBASECLASSIFIER_H_

#include "../../tools/parameters/CSParameterContext.h"
#include "../../tools/LabelMapGraph.h"
#include "../../tools/FeatureLabelObject.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkRGBPixel.h>

#include <vtkUndirectedGraph.h>
#include <vtkSmartPointer.h>

#include <set>

#include "../../pipelines/ObjectBasedSegmentation.h"


//encapsulates base classifier working on LabelMapGraph class (e.g. hard coded feature thresholds, NN, SVM, ..)
template< unsigned int VImageDimension > class GraphBaseClassifier
{
protected:
    typedef unsigned char                                                                   CScalarPixelType;
    typedef unsigned long int                                                               LScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType>                                                 CRGBPixelType;

    typedef itk::Image<CScalarPixelType, VImageDimension>                                   CScalarImageType;
    typedef itk::Image<LScalarPixelType, VImageDimension>                                   LScalarImageType;
    typedef itk::Image<CRGBPixelType, VImageDimension>                                      CRGBImageType;

    typedef itk::FeatureLabelObject<LScalarPixelType, VImageDimension>                      FeatureLabelObjectType;
    typedef itk::LabelMapGraph<CRGBImageType, LScalarImageType, FeatureLabelObjectType >    LabelMapGraphType;
    typedef typename LabelMapGraphType::Pointer                                             LabelMapGraphPointerType;

    typedef typename ObjectBasedSegmentation<VImageDimension>::Class                        ClassType;

    typedef itk::ImageFileReader<CScalarImageType>                                          ScalarReaderType;

public:
    GraphBaseClassifier();
    virtual ~GraphBaseClassifier();

    virtual void SetPath(std::string path) { mPath = path; };
    virtual void SetGraph(vtkSmartPointer<vtkUndirectedGraph> graph) { mGraph = graph; };       //or do we need the LabelMap also? faster?
    virtual void SetLabelMapGraph(LabelMapGraphPointerType labelMapGraph) { mLabelMapGraph = labelMapGraph; mGraph = mLabelMapGraph->GetLabelMapGraph(); };
    virtual void SetTrainingDatabase(std::string filename) { mTrainingDatabaseFilename = filename; };

    virtual void SetTargetClass(ClassType targetClass) { mTargetClass = targetClass; };

    virtual void Init();
    virtual void SetupClassifier() {};

protected:

    std::string mPath;
    std::string mTrainingDatabaseFilename;

    ClassType mTargetClass;

    unsigned int mNumActiveUnaryFeatures;
    unsigned int mNumActiveUnaryFeatureComponents;
    unsigned int mNumActiveBinaryFeatures;

    vtkSmartPointer<vtkUndirectedGraph> mGraph;
    LabelMapGraphPointerType mLabelMapGraph;
};

#include "GraphBaseClassifier.tpp"

#endif /* GRAPHBASECLASSIFIER_H_ */
