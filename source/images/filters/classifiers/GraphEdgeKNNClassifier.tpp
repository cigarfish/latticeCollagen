///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphEdgeKNNClassifier.tpp                                           //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-03-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "GraphEdgeKNNClassifier.h"


#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkEdgeListIterator.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkOutEdgeIterator.h>
#include <vtkUnsignedLongArray.h>


template< unsigned int VImageDimension > GraphEdgeKNNClassifier< VImageDimension >::GraphEdgeKNNClassifier()
    : GraphEdgeBaseClassifier<VImageDimension>()
{
    // TODO Auto-generated constructor stub
    mNumberOfNeighbors = 50;
}


template< unsigned int VImageDimension > GraphEdgeKNNClassifier< VImageDimension >::~GraphEdgeKNNClassifier()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void GraphEdgeKNNClassifier< VImageDimension >::Init()
{
    GraphEdgeBaseClassifier<VImageDimension>::Init();
    SetupClassifier();
}


template< unsigned int VImageDimension > void GraphEdgeKNNClassifier< VImageDimension >::SetupClassifier()
{
    std::map< int, ClassType > idToClass;
    SampleTypePointer sample = SampleType::New();

    //Use Load method provided by TrainClassifier class
//    LoadTrainingDatabase(idToClass, sample);

    mClassifiedTrainingData = SampleClassType::New();
    mClassifiedTrainingData->SetSample(sample);
    mClassifiedTrainingData->SetNumberOfClasses(mNumberOfClasses);

    for(unsigned int i=0; i<sample->Size(); i++)
        mClassifiedTrainingData->AddInstance(idToClass[i], i);

    TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
    treeGenerator->SetSample(sample);
    treeGenerator->SetBucketSize(16);
    treeGenerator->Update();

    mDecisionTree = treeGenerator->GetOutput();

//    Test code************************************************************************************************************************
//    typedef TreeType::NearestNeighbors NeighborsType;
//    typedef TreeType::KdTreeNodeType NodeType;
//
//    NodeType* root = tree->GetRoot();
//    if(root->IsTerminal())
//        std::cout << "Root node is a terminal node." << std::endl;
//    else
//        std::cout << "Root node is not a terminal node." << std::endl;
//
//    unsigned int partitionDimension;
//    double partitionValue;
//    root->GetParameters(partitionDimension, partitionValue);
//    std::cout << "Dimension chosen to split the space = " << partitionDimension << std::endl;
//    std::cout << "Split point on the partition dimension = " << partitionValue << std::endl;
//    std::cout << "Address of the left chile of the root node = " << root->Left() << std::endl;
//    std::cout << "Address of the right chile of the root node = " << root->Right() << std::endl;
//    std::cout << "Number of the measurement vectors under the root node" << " in the tree hierarchy = " << root->Size() << std::endl;
//    NodeType::CentroidType centroid;
//    root->GetWeightedCentroid( centroid );
//    std::cout << "Sum of the measurement vectors under the root node = " << centroid << std::endl;
//    std::cout << "Number of the measurement vectors under the left child" << " of the root node = " << root->Left()->Size() << std::endl;
}


//TODO: Dynamic evaluation of parameters
template< unsigned int VImageDimension > std::vector< std::pair<vtkIdType,vtkIdType> > GraphEdgeKNNClassifier< VImageDimension >::EvaluateAllLabelObjectPairs()
{
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIdArray = vtkUnsignedLongArray::SafeDownCast(this->mGraph->GetVertexData()->GetPedigreeIds());
    vtkSmartPointer<vtkFloatArray> avgBlueIntArray = vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(LabelMapGraphType::UnaryFeatureName[LabelMapGraphType::MEAN_BLUE_INTENSITY].c_str()));
    vtkSmartPointer<vtkFloatArray> avgBlueIntSimArray = vtkFloatArray::SafeDownCast(this->mGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[LabelMapGraphType::MEAN_BLUE_SIMILARITY].c_str()));

    this->mLabelObjectPairsOfTargetClass.clear();

    vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
    this->mGraph->GetEdges(it);
    while (it->HasNext()) {
        vtkEdgeType e = it->Next();

        //Evaluate pair
        MeasurementVectorType queryPoint;
        queryPoint[0] = avgBlueIntArray->GetValue(e.Source) / 255.;
        queryPoint[1] = avgBlueIntArray->GetValue(e.Target) / 255.;
        queryPoint[2] = avgBlueIntSimArray->GetValue(e.Id);

        //Decision: pair to merge?
        DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
        DistanceMetricType::OriginType origin(3);

//        std::cout << "sample->GetMeasurementVectorSize() = " << mDecisionTree->GetSample()->GetMeasurementVectorSize() << std::endl;
        for(unsigned int i=0; i<mDecisionTree->GetSample()->GetMeasurementVectorSize(); ++i)
            origin[i] = queryPoint[i];

        distanceMetric->SetOrigin(origin);

        TreeType::InstanceIdentifierVectorType neighbors;
        mDecisionTree->Search(queryPoint, mNumberOfNeighbors, neighbors);

//        std::cout << "weighted centroid kd-tree knn search result:" << std::endl
//                << "query point = [" << queryPoint << "]" << std::endl
//                << "k = " << mNumberOfNeighbors << std::endl;
//        std::cout << "class : measurement vector : distance" << std::endl;
//
//        for(unsigned int i=0; i<mNumberOfNeighbors; ++i) {
//            std::cout << mClassifiedTrainingData->GetClassLabel(neighbors[i]);
//            std::cout << " : [" << mDecisionTree->GetMeasurementVector(neighbors[i]) << "] : " << distanceMetric->Evaluate(mDecisionTree->GetMeasurementVector( neighbors[i])) << std::endl;
//        }

        if(pedigreeIdArray->GetValue(e.Source)==10845 || pedigreeIdArray->GetValue(e.Target)==10845) {
            std::cout << "Encountered LabelObjects " << pedigreeIdArray->GetValue(e.Source) << " and " << pedigreeIdArray->GetValue(e.Target) << std::endl;
            std::cout << "with blue intensity " << avgBlueIntArray->GetValue(e.Source) << " and " << avgBlueIntArray->GetValue(e.Target) << std::endl;
            std::cout << "and similarity " << avgBlueIntSimArray->GetValue(e.Id) << std::endl;
        }

        int classNuclei = 0, classBackground = 0;
        for(unsigned int i=0; i<mNumberOfNeighbors; ++i) {
//            std::cout << "neighbors[" << i << "] = " << neighbors[i] << std::endl;

            if(mClassifiedTrainingData->GetClassLabel(neighbors[i]) == 1) {       //it is a nucleus superpixel pair!!!
                classNuclei++;
//                std::cout << "nucleus match" << std::endl;
            }
            else {
                classBackground++;
//                std::cout << "background match" << std::endl;
            }
        }

        if(pedigreeIdArray->GetValue(e.Source)==10845 || pedigreeIdArray->GetValue(e.Target)==10845)
            std::cout << "overall nucleus vs background match = " << classNuclei << ", " << classBackground << std::endl;

        if(classNuclei > classBackground)
            this->mLabelObjectPairsOfTargetClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(e.Source), pedigreeIdArray->GetValue(e.Target)));
    }

    return this->mLabelObjectPairsOfTargetClass;
}


//TODO: Dynamic evaluation of parameters
template< unsigned int VImageDimension > std::vector< std::pair<vtkIdType,vtkIdType> > GraphEdgeKNNClassifier< VImageDimension >::EvaluateModifiedLabelObjectPairs(std::set<vtkIdType> &modifiedInLastIteration)
{
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIdArray = vtkUnsignedLongArray::SafeDownCast(this->mGraph->GetVertexData()->GetPedigreeIds());
    vtkSmartPointer<vtkFloatArray> avgBlueIntArray = vtkFloatArray::SafeDownCast(this->mGraph->GetVertexData()->GetArray(LabelMapGraphType::UnaryFeatureName[LabelMapGraphType::MEAN_BLUE_INTENSITY].c_str()));
    vtkSmartPointer<vtkFloatArray> avgBlueIntSimArray = vtkFloatArray::SafeDownCast(this->mGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[LabelMapGraphType::MEAN_BLUE_SIMILARITY].c_str()));

    std::set<vtkIdType> edgesAlreadyVisited;
    this->mLabelObjectPairsOfTargetClass.clear();

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
                    MeasurementVectorType queryPoint;
                    queryPoint[0] = avgBlueIntArray->GetValue(source) / 255.;
                    queryPoint[1] = avgBlueIntArray->GetValue(target) / 255.;
                    queryPoint[2] = avgBlueIntSimArray->GetValue(e.Id);

                    //Decision: pair to merge?
                    DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
                    DistanceMetricType::OriginType origin(3);

//                    std::cout << "sample->GetMeasurementVectorSize() = " << mDecisionTree->GetSample()->GetMeasurementVectorSize() << std::endl;
                    for(unsigned int i=0; i<mDecisionTree->GetSample()->GetMeasurementVectorSize(); ++i)
                        origin[i] = queryPoint[i];

                    distanceMetric->SetOrigin(origin);

                    TreeType::InstanceIdentifierVectorType neighbors;
                    mDecisionTree->Search(queryPoint, mNumberOfNeighbors, neighbors);

//                    std::cout << "weighted centroid kd-tree knn search result:" << std::endl
//                            << "query point = [" << queryPoint << "]" << std::endl
//                            << "k = " << mNumberOfNeighbors << std::endl;
//                    std::cout << "class : measurement vector : distance" << std::endl;
//
//                    for(unsigned int i=0; i<mNumberOfNeighbors; ++i) {
//                        std::cout << mClassifiedTrainingData->GetClassLabel(neighbors[i]);
//                        std::cout << " : [" << mDecisionTree->GetMeasurementVector(neighbors[i]) << "] : " << distanceMetric->Evaluate(mDecisionTree->GetMeasurementVector( neighbors[i])) << std::endl;
//                    }

                    int classNuclei = 0, classBackground = 0;
                    for(unsigned int i=0; i<mNumberOfNeighbors; ++i) {
                        if(mClassifiedTrainingData->GetClassLabel(neighbors[i]) == 1)       //it is a nucleus superpixel pair!!!
                            classNuclei++;
                        else
                            classBackground++;
                    }

                    if(classNuclei > classBackground)
                        this->mLabelObjectPairsOfTargetClass.push_back(std::pair<vtkIdType,vtkIdType>(pedigreeIdArray->GetValue(source), pedigreeIdArray->GetValue(target)));
                    edgesAlreadyVisited.insert(e.Id);
                }
            }
        }
    }
    return this->mLabelObjectPairsOfTargetClass;
}
