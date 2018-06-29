///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphEdgeBaseClassifier.cpp                                          //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-02-19                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHEDGEBASECLASSIFIER_H_
#define GRAPHEDGEBASECLASSIFIER_H_

#include "GraphBaseClassifier.h"


//encapsulates classifier working on edges of LabelMapGraph class (e.g. hard coded feature thresholds, NN, SVM, ..)
template< unsigned int VImageDimension > class GraphEdgeBaseClassifier : public GraphBaseClassifier<VImageDimension>
{
public:
    GraphEdgeBaseClassifier();
    virtual ~GraphEdgeBaseClassifier();

    //returns vector of vertex pairs given by their pedigreeId
    virtual std::vector< std::pair<vtkIdType,vtkIdType> > EvaluateAllLabelObjectPairs() = 0;
    virtual std::vector< std::pair<vtkIdType,vtkIdType> > EvaluateModifiedLabelObjectPairs(std::set<vtkIdType> &modifiedInLastIteration) = 0;

    virtual void SetupClassifier() {};

protected:

    std::vector< std::pair<vtkIdType,vtkIdType> > mLabelObjectPairsOfTargetClass;
};

#include "GraphEdgeBaseClassifier.tpp"

#endif /* GraphEdgeBaseClassifier_H_ */
