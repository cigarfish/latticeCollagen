///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphVertexPythonSVMClassifier.h                                     //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-03-23                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHVERTEXPYTHONSVMCLASSIFIER_H_
#define GRAPHVERTEXPYTHONSVMCLASSIFIER_H_

#ifdef CS_BUILD_PYTHON
#include <Python.h>
#endif

#include "GraphVertexBaseClassifier.h"


template< unsigned int VImageDimension > class GraphVertexPythonSVMClassifier : public GraphVertexBaseClassifier<VImageDimension>
{
private:
    typedef typename GraphBaseClassifier<VImageDimension>::ClassType            ClassType;
    typedef typename GraphBaseClassifier<VImageDimension>::LabelMapGraphType    LabelMapGraphType;

    typedef typename LabelMapGraphType::UnaryFeatures                           UnaryFeatureType;
    typedef typename LabelMapGraphType::BinaryFeatures                          BinaryFeatureType;

public:
    GraphVertexPythonSVMClassifier();
    virtual ~GraphVertexPythonSVMClassifier();

    std::set<vtkIdType> EvaluateAllLabelObjects();
    std::set<vtkIdType> EvaluateModifiedLabelObjects(std::set<vtkIdType> &modifiedInLastIteration);

    void Init();
    void SetupClassifier();

private:

#ifdef CS_BUILD_PYTHON
    PyObject *mpTrainDBHandlerInstance;
    PyObject *mpSVMHandlingInstance;
#endif
};

#include "GraphVertexPythonSVMClassifier.tpp"

#endif /* GRAPHVERTEXPYTHONSVMCLASSIFIER_H_ */
