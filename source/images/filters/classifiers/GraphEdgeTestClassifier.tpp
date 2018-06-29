///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphEdgeTestClassifier.tpp                                          //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-02-19                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "GraphEdgeTestClassifier.h"

#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkEdgeListIterator.h>
#include <vtkFloatArray.h>
#include <vtkOutEdgeIterator.h>
#include <vtkUnsignedLongArray.h>
#include <vtkVariant.h>



template< unsigned int VImageDimension > GraphEdgeTestClassifier< VImageDimension >::GraphEdgeTestClassifier()
    : GraphEdgeBaseClassifier<VImageDimension>()
{
    // TODO Auto-generated constructor stub
    mMinAverageRedIntensity = 60;
    mMinAverageGreenIntensity = 50;
    mMinAverageBlueIntensity = 110;
    mMinAverageRedIntensitySimilarity = 0.6;
    mMinAverageGreenIntensitySimilarity = 0.5;
    mMinAverageBlueIntensitySimilarity = 0.7;
}


template< unsigned int VImageDimension > GraphEdgeTestClassifier< VImageDimension >::~GraphEdgeTestClassifier()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void GraphEdgeTestClassifier< VImageDimension >::Init()
{
    GraphEdgeBaseClassifier<VImageDimension>::Init();
    TrainClassifier();
}


template< unsigned int VImageDimension > void GraphEdgeTestClassifier< VImageDimension >::TrainClassifier()
{
    //Nothing to do here for this simple classifier example
}


//TODO: Once parameter handling is dynamic (user specified) abandon hardcoded class->feature relation
template< unsigned int VImageDimension > std::vector< std::pair<vtkIdType,vtkIdType> > GraphEdgeTestClassifier< VImageDimension >::EvaluateAllLabelObjectPairs()
{
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIdArray = vtkUnsignedLongArray::SafeDownCast(this->mGraph->GetVertexData()->GetPedigreeIds());
    vtkSmartPointer<vtkFloatArray> avgRedIntArray, avgGreenIntArray, avgBlueIntArray, avgRedIntSimArray, avgGreenIntSimArray, avgBlueIntSimArray;
    if(this->mTargetClass == ObjectBasedSegmentationType::NUCLEUS) {
        avgBlueIntArray = vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(LabelMapGraphType::UnaryFeatureName[LabelMapGraphType::MEAN_BLUE_INTENSITY].c_str()));
        avgBlueIntSimArray = vtkFloatArray::SafeDownCast(this->mGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[LabelMapGraphType::MEAN_BLUE_SIMILARITY].c_str()));
    }
    else if(this->mTargetClass == ObjectBasedSegmentationType::BILE || this->mTargetClass == ObjectBasedSegmentationType::SINUSOID) {
        avgRedIntArray = vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(LabelMapGraphType::UnaryFeatureName[LabelMapGraphType::MEAN_RED_INTENSITY].c_str()));
        avgGreenIntArray = vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(LabelMapGraphType::UnaryFeatureName[LabelMapGraphType::MEAN_GREEN_INTENSITY].c_str()));
        avgRedIntSimArray = vtkFloatArray::SafeDownCast(this->mGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[LabelMapGraphType::MEAN_RED_SIMILARITY].c_str()));
        avgGreenIntSimArray = vtkFloatArray::SafeDownCast(this->mGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[LabelMapGraphType::MEAN_GREEN_SIMILARITY].c_str()));
    }

    this->mPairsOfSameClass.clear();

    vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
    this->mGraph->GetEdges(it);
    while (it->HasNext()) {
        vtkEdgeType e = it->Next();

        //Evaluate pair
        float avgRedInt, avgGreenInt, avgBlueInt, avgRedIntSim, avgGreenIntSim, avgBlueIntSim;
        if(this->mTargetClass == ObjectBasedSegmentationType::NUCLEUS) {
            avgBlueInt = ( avgBlueIntArray->GetValue(e.Source) + avgBlueIntArray->GetValue(e.Target) ) / 2.;
            avgBlueIntSim = avgBlueIntSimArray->GetValue(e.Id);
        }
        else if(this->mTargetClass == ObjectBasedSegmentationType::BILE || this->mTargetClass == ObjectBasedSegmentationType::SINUSOID) {
            avgRedInt = ( avgRedIntArray->GetValue(e.Source) + avgRedIntArray->GetValue(e.Target) ) / 2.;
            avgGreenInt = ( avgGreenIntArray->GetValue(e.Source) + avgGreenIntArray->GetValue(e.Target) ) / 2.;
            avgRedIntSim = avgRedIntSimArray->GetValue(e.Id);
            avgGreenIntSim = avgGreenIntSimArray->GetValue(e.Id);
        }

        //Decision: pair to merge?
        if(this->mTargetClass == ObjectBasedSegmentationType::NUCLEUS) {
            if(avgBlueIntSim > mMinAverageBlueIntensitySimilarity && avgBlueInt > mMinAverageBlueIntensity)
                this->mPairsOfSameClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(e.Source), pedigreeIdArray->GetValue(e.Target)));
        }
        else if(this->mTargetClass == ObjectBasedSegmentationType::BILE) {
            if(avgGreenIntSim > mMinAverageGreenIntensitySimilarity && avgGreenInt > mMinAverageGreenIntensity && avgRedInt <= mMinAverageRedIntensity)
                this->mPairsOfSameClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(e.Source), pedigreeIdArray->GetValue(e.Target)));
        }
        else if(this->mTargetClass == ObjectBasedSegmentationType::SINUSOID) {
            if(avgGreenIntSim > mMinAverageGreenIntensitySimilarity && avgGreenInt > mMinAverageGreenIntensity && avgRedInt > mMinAverageRedIntensity && avgRedIntSim > mMinAverageRedIntensitySimilarity)
                this->mPairsOfSameClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(e.Source), pedigreeIdArray->GetValue(e.Target)));
        }
    }
    return this->mPairsOfSameClass;
}


//TODO: Once parameter handling is dynamic (user specified) abandon hardcoded class->feature relation
template< unsigned int VImageDimension > std::vector< std::pair<vtkIdType,vtkIdType> > GraphEdgeTestClassifier< VImageDimension >::EvaluateModifiedLabelObjectPairs(std::set<vtkIdType> &modifiedInLastIteration)
{
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIdArray = vtkUnsignedLongArray::SafeDownCast(this->mGraph->GetVertexData()->GetPedigreeIds());
    vtkSmartPointer<vtkFloatArray> avgRedIntArray, avgGreenIntArray, avgBlueIntArray, avgRedIntSimArray, avgGreenIntSimArray, avgBlueIntSimArray;
    if(this->mTargetClass == ObjectBasedSegmentationType::NUCLEUS) {
        avgBlueIntArray = vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(LabelMapGraphType::UnaryFeatureName[LabelMapGraphType::MEAN_BLUE_INTENSITY].c_str()));
        avgBlueIntSimArray = vtkFloatArray::SafeDownCast(this->mGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[LabelMapGraphType::MEAN_BLUE_SIMILARITY].c_str()));
    }
    else if(this->mTargetClass == ObjectBasedSegmentationType::BILE || this->mTargetClass == ObjectBasedSegmentationType::SINUSOID) {
        avgRedIntArray = vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(LabelMapGraphType::UnaryFeatureName[LabelMapGraphType::MEAN_RED_INTENSITY].c_str()));
        avgGreenIntArray = vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(LabelMapGraphType::UnaryFeatureName[LabelMapGraphType::MEAN_GREEN_INTENSITY].c_str()));
        avgRedIntSimArray = vtkFloatArray::SafeDownCast(this->mGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[LabelMapGraphType::MEAN_RED_SIMILARITY].c_str()));
        avgGreenIntSimArray = vtkFloatArray::SafeDownCast(this->mGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[LabelMapGraphType::MEAN_GREEN_SIMILARITY].c_str()));
    }

    std::set<vtkIdType> edgesAlreadyVisited;
    this->mPairsOfSameClass.clear();

    for(std::set<vtkIdType>::iterator it=modifiedInLastIteration.begin(); it!=modifiedInLastIteration.end(); ++it) {
        vtkIdType source = this->mGraph->FindVertex(*it);

        if(source != -1) {
            vtkSmartPointer<vtkOutEdgeIterator> outIt = vtkSmartPointer<vtkOutEdgeIterator>::New();
            this->mGraph->GetOutEdges(source, outIt);

            while (outIt->HasNext()) {
                vtkOutEdgeType e = outIt->Next();

                if(edgesAlreadyVisited.count(e.Id)==0) {
                    vtkIdType target = e.Target;

                    //Evaluate pair
                    float avgRedInt, avgGreenInt, avgBlueInt, avgRedIntSim, avgGreenIntSim, avgBlueIntSim;
                    if(this->mTargetClass == ObjectBasedSegmentationType::NUCLEUS) {
                        avgBlueInt = ( avgBlueIntArray->GetValue(source) + avgBlueIntArray->GetValue(target) ) / 2.;
                        avgBlueIntSim = avgBlueIntSimArray->GetValue(e.Id);
                    }
                    else if(this->mTargetClass == ObjectBasedSegmentationType::BILE || this->mTargetClass == ObjectBasedSegmentationType::SINUSOID) {
                        avgRedInt = ( avgRedIntArray->GetValue(source) + avgRedIntArray->GetValue(target) ) / 2.;
                        avgGreenInt = ( avgGreenIntArray->GetValue(source) + avgGreenIntArray->GetValue(target) ) / 2.;
                        avgRedIntSim = avgRedIntSimArray->GetValue(e.Id);
                        avgGreenIntSim = avgGreenIntSimArray->GetValue(e.Id);
                    }

                    //Decision: pair to merge?
                    if(this->mTargetClass == ObjectBasedSegmentationType::NUCLEUS) {
                        if(avgBlueIntSim > mMinAverageBlueIntensitySimilarity && avgBlueInt > mMinAverageBlueIntensity)
                            this->mPairsOfSameClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(source), pedigreeIdArray->GetValue(target)));
                    }
                    else if(this->mTargetClass == ObjectBasedSegmentationType::BILE) {
                        if(avgGreenIntSim > mMinAverageGreenIntensitySimilarity && avgGreenInt > mMinAverageGreenIntensity && avgRedInt <= mMinAverageRedIntensity)
                            this->mPairsOfSameClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(source), pedigreeIdArray->GetValue(target)));
                    }
                    else if(this->mTargetClass == ObjectBasedSegmentationType::SINUSOID) {
                        if(avgGreenIntSim > mMinAverageGreenIntensitySimilarity && avgGreenInt > mMinAverageGreenIntensity && avgRedInt > mMinAverageRedIntensity && avgRedIntSim > mMinAverageRedIntensitySimilarity)
                            this->mPairsOfSameClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(source), pedigreeIdArray->GetValue(target)));
                    }
                    edgesAlreadyVisited.insert(e.Id);
                }
            }
        }
    }
    return this->mPairsOfSameClass;
}
