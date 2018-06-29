///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphEdgeBaseClassifier.h                                            //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-02-19                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHEDGETESTCLASSIFIER_H_
#define GRAPHEDGETESTCLASSIFIER_H_

#include "GraphEdgeBaseClassifier.h"


template< unsigned int VImageDimension > class GraphEdgeTestClassifier : public GraphEdgeBaseClassifier<VImageDimension>
{
    typedef typename GraphEdgeBaseClassifier<VImageDimension>::ClassType            ClassType;
    typedef typename GraphEdgeBaseClassifier<VImageDimension>::LabelMapGraphType    LabelMapGraphType;

    typedef ObjectBasedSegmentation<VImageDimension>                                ObjectBasedSegmentationType;

public:
    GraphEdgeTestClassifier();
    virtual ~GraphEdgeTestClassifier();

    std::vector< std::pair<vtkIdType,vtkIdType> > EvaluateAllLabelObjectPairs();
    std::vector< std::pair<vtkIdType,vtkIdType> > EvaluateModifiedLabelObjectPairs(std::set<vtkIdType> &modifiedInLastIteration);
    void Init();
    void TrainClassifier();

protected:
    //all my classifier-specific variables
    double mMinAverageRedIntensity;
    double mMinAverageGreenIntensity;
    double mMinAverageBlueIntensity;
    double mMinAverageRedIntensitySimilarity;
    double mMinAverageGreenIntensitySimilarity;
    double mMinAverageBlueIntensitySimilarity;
};

#include "GraphEdgeTestClassifier.tpp"

#endif /* GRAPHEDGETESTCLASSIFIER_H_ */
