///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphVertexBaseClassifier.h                                          //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-03-23                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHVERTEXBASECLASSIFIER_H_
#define GRAPHVERTEXBASECLASSIFIER_H_

#include "GraphBaseClassifier.h"


//encapsulates classifier working on vertices of LabelMapGraph class (e.g. hard coded feature thresholds, NN, SVM, ..)
template< unsigned int VImageDimension > class GraphVertexBaseClassifier : public GraphBaseClassifier<VImageDimension>
{
public:
    GraphVertexBaseClassifier();
    virtual ~GraphVertexBaseClassifier();

    std::map<std::string, std::map<vtkIdType, float> > GetLabelObjectClassProbabilities() { return mClassProbabilites; };

    //returns vector of vertices given by their pedigreeId
    virtual std::set<vtkIdType> EvaluateAllLabelObjects() = 0;
    virtual std::set<vtkIdType> EvaluateModifiedLabelObjects(std::set<vtkIdType> &modifiedInLastIteration) = 0;

    virtual void SetupClassifier() {};

protected:
    std::set<vtkIdType> mLabelObjectsOfTargetClass;
    std::map<std::string, std::map<vtkIdType, float> > mClassProbabilites;
};

#include "GraphVertexBaseClassifier.tpp"

#endif /* GRAPHVERTEXBASECLASSIFIER_H_ */
