///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphEdgeKNNClassifier.h                                             //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-03-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHEDGEKNNCLASSIFIER_H_
#define GRAPHEDGEKNNCLASSIFIER_H_

#include "GraphEdgeBaseClassifier.h"

#include "itkEuclideanDistanceMetric.h"
#include "itkKdTree.h"
#include "itkListSample.h"
#include "itkMembershipSample.h"
#include "itkVector.h"
#include "itkWeightedCentroidKdTreeGenerator.h"


class SampleTypePointer;

template< unsigned int VImageDimension > class GraphEdgeKNNClassifier : public GraphEdgeBaseClassifier<VImageDimension>
{
    typedef typename GraphBaseClassifier<VImageDimension>::ClassType            ClassType;
    typedef typename GraphBaseClassifier<VImageDimension>::LabelMapGraphType    LabelMapGraphType;

    typedef itk::Vector<double, 3>                                              MeasurementVectorType;
    typedef itk::Statistics::ListSample<MeasurementVectorType>                  SampleType;
    typedef SampleType::Pointer                                                 SampleTypePointer;
    typedef itk::Statistics::MembershipSample<SampleType>                       SampleClassType;

    typedef itk::Statistics::WeightedCentroidKdTreeGenerator<SampleType>        TreeGeneratorType;
    typedef TreeGeneratorType::KdTreeType                                       TreeType;

    typedef itk::Statistics::EuclideanDistanceMetric<MeasurementVectorType>     DistanceMetricType;

public:
    GraphEdgeKNNClassifier();
    virtual ~GraphEdgeKNNClassifier();

    std::vector< std::pair<vtkIdType,vtkIdType> > EvaluateAllLabelObjectPairs();
    std::vector< std::pair<vtkIdType,vtkIdType> > EvaluateModifiedLabelObjectPairs(std::set<vtkIdType> &modifiedInLastIteration);

    void Init();
    void SetupClassifier();

protected:
    TreeType::Pointer mDecisionTree;
    SampleClassType::Pointer mClassifiedTrainingData;

    unsigned int mNumberOfClasses;
    unsigned int mNumberOfNeighbors;
};

#include "GraphEdgeKNNClassifier.tpp"

#endif /* GRAPHEDGEKNNCLASSIFIER_H_ */
