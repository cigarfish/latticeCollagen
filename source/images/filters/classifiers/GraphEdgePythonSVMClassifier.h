///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphEdgePythonSVMClassifier.h                                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-08-08                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHEDGEPYTHONSVMCLASSIFIER_H_
#define GRAPHEDGEPYTHONSVMCLASSIFIER_H_

#ifdef CS_BUILD_PYTHON
#include <Python.h>
#endif

#include "GraphEdgeBaseClassifier.h"


template< unsigned int VImageDimension > class GraphEdgePythonSVMClassifier : public GraphEdgeBaseClassifier<VImageDimension>
{
private:
    typedef typename GraphBaseClassifier<VImageDimension>::ClassType            ClassType;
    typedef typename GraphBaseClassifier<VImageDimension>::LabelMapGraphType    LabelMapGraphType;

    typedef typename LabelMapGraphType::UnaryFeatures                           UnaryFeatureType;
    typedef typename LabelMapGraphType::BinaryFeatures                          BinaryFeatureType;

public:
    GraphEdgePythonSVMClassifier();
    virtual ~GraphEdgePythonSVMClassifier();

    std::vector< std::pair<vtkIdType,vtkIdType> > EvaluateAllLabelObjectPairs();
    std::vector< std::pair<vtkIdType,vtkIdType> > EvaluateModifiedLabelObjectPairs(std::set<vtkIdType> &modifiedInLastIteration);

    void Init();
    void SetupClassifier();

private:

#ifdef CS_BUILD_PYTHON
    PyObject *mpTrainDBHandlerInstance;
    PyObject *mpSVMHandlingInstance;
#endif
};

#include "GraphEdgePythonSVMClassifier.tpp"

#endif /* GRAPHEDGEPYTHONSVMCLASSIFIER_H_ */
