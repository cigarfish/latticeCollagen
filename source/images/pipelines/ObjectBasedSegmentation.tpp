///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ObjectBasedSegmentation.tpp                                          //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-07-31                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "ObjectBasedSegmentation.h"

#ifdef CS_BUILD_PYTHON
#include <Python.h>
#endif

#include <limits>
#include <iostream>

#include <QDateTime>

#include "../../../Core.h"

#include "../filters/classifiers/GraphEdgeBaseClassifier.h"
#include "../filters/classifiers/GraphEdgeKNNClassifier.h"
#include "../filters/classifiers/GraphEdgePythonSVMClassifier.h"
#include "../filters/classifiers/GraphEdgeTestClassifier.h"
#include "../filters/classifiers/GraphVertexPythonSVMClassifier.h"
#include "../tools/GraphAnnotationHelper.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



template< unsigned int VImageDimension > const std::string ObjectBasedSegmentation< VImageDimension >::ClassName[] =
{
        "Background",
        "Nuclei",
        "Sinusoids",
        "Bile",
        "Veins",
        "Necrotic region",
        "" // be sure this empty string is the last entry
};

template< unsigned int VImageDimension > const std::string ObjectBasedSegmentation< VImageDimension >::SegmentationModeName[] =
{
        "Label-Object-based",
        "Label-Object-pair-based",
        "" // be sure this empty string is the last entry

};

template< unsigned int VImageDimension > const std::string ObjectBasedSegmentation< VImageDimension >::ClassifierName[] =
{
        "WithoutSegmentation",
        "KNNBasedMerging",
        "PythonSVMBasedMerging",
        "" // be sure this empty string is the last entry

};


template< unsigned int VImageDimension > ObjectBasedSegmentation< VImageDimension >::ObjectBasedSegmentation()
{
}


template< unsigned int VImageDimension > ObjectBasedSegmentation< VImageDimension >::~ObjectBasedSegmentation()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void ObjectBasedSegmentation< VImageDimension >::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    this->m_paramContext->dump(parameter);

    file.open((this->mRawImagePath + mLogFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-object-based-segmentation-------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-object-based-segmentation---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


template< unsigned int VImageDimension > void ObjectBasedSegmentation< VImageDimension >::WriteDataSetSummary()
{
    //TODO
//    ImageAnalysisSummaryFileIO::AddEntry(HepNucleiSegmentationBin, m_pathChannelDAPI, m_pathChannelDAPI + m_filenameSave + m_saveSuffixes[5] + m_fileExtensionChannelDAPI);
}


template< unsigned int VImageDimension > void ObjectBasedSegmentation< VImageDimension >::ParseParameterContext()
{
    if(this->m_paramContext->findContext("Object-based Segmentation",0)==NULL) {
        std::cout << "Error: ObjectBasedSegmentation: Invalid parameter context" << std::endl;
        return;
    }

    this->mRawImageFullFilename = *(std::string*)(this->m_paramContext->findContext("Superpixel Partitioning Reader")->findParameter("Image filename", 0)->dataPointer());
    this->mObjectImageFullFilename = *(std::string*)(this->m_paramContext->findContext("Superpixel Partitioning Reader")->findParameter("Object dataset filename", 0)->dataPointer());
    this->mObjectGraphFullFilename = *(std::string*)(this->m_paramContext->findContext("Superpixel Partitioning Reader")->findParameter("Object graph filename", 0)->dataPointer());

    this->mVoxelSpacing[0] = *(double*)(this->m_paramContext->findParameter("Image voxel spacing x", 0)->dataPointer());
    this->mVoxelSpacing[1] = *(double*)(this->m_paramContext->findParameter("Image voxel spacing y", 0)->dataPointer());
    if(ImageDimension==3)
        this->mVoxelSpacing[2] = *(double*)(this->m_paramContext->findParameter("Image voxel spacing z", 0)->dataPointer());

    mTrainingDatabaseFilename = *(std::string*)(this->m_paramContext->findContext("Segmentation")->findParameter("Training database filename", 0)->dataPointer());

    std::string segmentationMode = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Segmentation mode", 0)->dataPointer()) )->currentString();
    for(int i=0; i<NumberSegmentationModes; i++) {
        if( segmentationMode.compare(SegmentationModeName[i])==0 ) {
            mSegmentationMode = (SegmentationMode)i;
            break;
        }
    }

    std::string segmentationClassifier = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Classifier", 0)->dataPointer()) )->currentString();
    for(int i=0; i<NumberClassifiers; i++) {
        if( segmentationClassifier.compare(ClassifierName[i])==0 ) {
            mSegmentationClassifier = (Classifier)i;
            break;
        }
    }

    std::string segmentationTargetClass = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Segmentation target class", 0)->dataPointer()) )->currentString();
    for(int i=0; i<NumberClasses; i++) {
        if( segmentationTargetClass.compare(ClassName[i])==0 ) {
            mSegmentationTargetClass = (Class)i;
            break;
        }
    }

    mWithFeature = new bool[LabelMapType::NumberUnaryFeatures+LabelMapType::NumberBinaryFeatures];
    int objI=0;
    for(objI=0; objI<LabelMapType::NumberUnaryFeatures; objI++)
        mWithFeature[objI] = *(bool*)(this->m_paramContext->findParameter("With " + LabelMapType::UnaryFeatureName[objI], 0)->dataPointer());
    for(int j=0; j<LabelMapType::NumberBinaryFeatures; j++)
        mWithFeature[objI+j] = *(bool*)(this->m_paramContext->findParameter("With " + LabelMapType::BinaryFeatureName[j], 0)->dataPointer());

    mKeepOnlyTargetClassLabelObjects = *(bool*)(this->m_paramContext->findParameter("Keep only target class label objects", 0)->dataPointer());

    //TODO: color for superpixel outline

    std::string saveMode = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        mSaveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        mSaveEverything = 0;
    mLogFilenameSave = *(std::string*)(this->m_paramContext->findParameter("Log file name", 0)->dataPointer());
    mFilenameSave = *(std::string*)(this->m_paramContext->findParameter("Save prefix", 0)->dataPointer());
}


template< unsigned int VImageDimension > void ObjectBasedSegmentation< VImageDimension >::ActivateFeaturesAtLabelMap()
{
    int i;
    for(i=0; i<LabelMapType::NumberUnaryFeatures; i++)
        if(mWithFeature[i])
            this->mLabelMap->ActivateUnaryFeature(static_cast<UnaryFeatureType>(i));

    for(int j=0; j<LabelMapType::NumberBinaryFeatures; j++)
        if(mWithFeature[i+j])
            this->mLabelMap->ActivateBinaryFeature(static_cast<BinaryFeatureType>(j));

    this->mLabelMap->InitializeAllLabelObjects(true);
}


template< unsigned int VImageDimension > void ObjectBasedSegmentation< VImageDimension >::ClassifierBasedMerging()
{
    //----------EXPERIMENTAL-LABEL-SELECTION-FILTER----------------------------------------------------------------------------
    std::cout << "Classify stuff" << std::endl;
    std::cout << "Setup classifier" << std::endl;

#ifdef CS_BUILD_PYTHON
    Core *core;

    PyEval_AcquireLock();
    PyInterpreterState * mainInterpreterState = core->pMainThreadState->interp;
    PyThreadState * myThreadState = PyThreadState_New(mainInterpreterState);

    PyThreadState_Swap(myThreadState);
#endif

    GraphEdgeBaseClassifier<ImageDimension> *classifier;
    if(mSegmentationClassifier == KNNBasedMerging)
        classifier = new GraphEdgeKNNClassifier<ImageDimension>();
    else if(mSegmentationClassifier == PythonSVMBasedMerging) {
        classifier = new GraphEdgePythonSVMClassifier<ImageDimension>();
    }
    else
        std::cout << "Unknown segmentation mode in ClassifierBasedMerging(). We stop here now. Sorry!" << std::endl;

    classifier->SetTargetClass(mSegmentationTargetClass);
    classifier->SetLabelMapGraph(this->mLabelMap);
    classifier->SetTrainingDatabase(mTrainingDatabaseFilename);
    classifier->Init();

    std::set<vtkIdType> alteredLabelObjects;
    bool mergingGoingOn = true;
    int mergingOperations = 0, iteration = 0;

    std::cout << "Start classification/merging loop" << std::endl;
    while(mergingGoingOn) {
        classifier->SetLabelMapGraph(this->mLabelMap);

        std::vector< std::pair<vtkIdType, vtkIdType> > pairsToMerge;
        if(mergingOperations==0)    pairsToMerge = classifier->EvaluateAllLabelObjectPairs();
        else                        pairsToMerge = classifier->EvaluateModifiedLabelObjectPairs(alteredLabelObjects);

        alteredLabelObjects.clear();
        mergingGoingOn = false;

        for(unsigned int i=0; i<pairsToMerge.size(); i++) {
            if(alteredLabelObjects.count(pairsToMerge[i].first)==0 && alteredLabelObjects.count(pairsToMerge[i].second)==0) {
                this->mLabelMap->LazyMergeLabelObjects(pairsToMerge[i].first, pairsToMerge[i].second);
                alteredLabelObjects.insert(pairsToMerge[i].first);
                alteredLabelObjects.insert(pairsToMerge[i].second);
                mergingGoingOn = true;
                mergingOperations++;
            }
        }
        std::cout << "Finished loop iteration: performed merging operations = " << mergingOperations << std::endl;

        if(mSaveEverything) {
            std::stringstream name;
            name << this->mRawImagePath << mFilenameSave << "_TestClassifier" << iteration << this->mRawImageFileExtension ;
            this->mObjectImage = this->mLabelMap->GetLabelImage();
            this->SaveObjectContourOverlayImage(name.str());
        }
        iteration++;
    }

#ifdef CS_BUILD_PYTHON
    PyThreadState_Swap(NULL);
    PyThreadState_Clear(myThreadState);
    PyThreadState_Delete(myThreadState);
    PyEval_ReleaseLock();
#endif

    this->mObjectImage = this->mLabelMap->GetLabelImage();
    this->mObjectGraph = this->mLabelMap->GetLabelMapGraph();
    std::cout << "-----number label objects/vertices = " << this->mObjectGraph->GetNumberOfVertices() << ", number label object pairs/edges = " << this->mObjectGraph->GetNumberOfEdges() << std::endl;
}


template< unsigned int VImageDimension > void ObjectBasedSegmentation< VImageDimension >::ClassifyLabelObjects()
{
    //----------EXPERIMENTAL-LABEL-SELECTION-FILTER----------------------------------------------------------------------------
    std::cout << "Classify stuff" << std::endl;
    std::cout << "Setup classifier" << std::endl;
    GraphAnnotationHelper anno;

#ifdef CS_BUILD_PYTHON
    Core *core;

    PyEval_AcquireLock();
    PyInterpreterState * mainInterpreterState = core->pMainThreadState->interp;
    PyThreadState * myThreadState = PyThreadState_New(mainInterpreterState);

    PyThreadState_Swap(myThreadState);
#endif

    GraphVertexBaseClassifier<ImageDimension> *classifier;
    if(mSegmentationClassifier == PythonSVMBasedMerging) {
        classifier = new GraphVertexPythonSVMClassifier<ImageDimension>();
    }
    else
        std::cout << "Unknown segmentation mode in ClassifyLabelObjects(). We stop here now. Sorry!" << std::endl;

    classifier->SetTargetClass(mSegmentationTargetClass);
    classifier->SetLabelMapGraph(this->mLabelMap);
    classifier->SetTrainingDatabase(mTrainingDatabaseFilename);
    classifier->Init();

    std::set<vtkIdType> labelObjectsOfTargetClass = classifier->EvaluateAllLabelObjects();

    std::map<std::string, std::map<vtkIdType, float> > classProbabilities = classifier->GetLabelObjectClassProbabilities();
    for(std::map<std::string, std::map<vtkIdType, float> >::iterator it = classProbabilities.begin(); it!=classProbabilities.end(); ++it)
        anno.AddCustomVertexAnnotation(this->mLabelMap->GetLabelMapGraph(), it->first, it->second, -1);

    std::map<vtkIdType, int> isTargetClass;
    for(std::set<vtkIdType>::iterator it=labelObjectsOfTargetClass.begin(); it!=labelObjectsOfTargetClass.end(); ++it)
        isTargetClass.insert(std::pair<vtkIdType, int>(this->mLabelMap->GetLabelMapGraph()->FindVertex(*it), 1));
    anno.AddCustomVertexAnnotation(this->mLabelMap->GetLabelMapGraph(), ClassName[mSegmentationTargetClass], isTargetClass, 0);
    this->mObjectGraph = this->mLabelMap->GetLabelMapGraph();
    this->SaveObjectGraph(this->mRawImagePath + mFilenameSave + "_" + SegmentationModeName[mSegmentationMode] + "_" + ClassName[mSegmentationTargetClass] + "_debug.txt");

    std::cout << "LabelObjects to merge: " << labelObjectsOfTargetClass.size() << std::endl;
    this->mLabelMap->LazyCollapseLabelObjects(labelObjectsOfTargetClass);
    if(mKeepOnlyTargetClassLabelObjects)
        this->mLabelMap->KeepLabelsRemoveRest(labelObjectsOfTargetClass);

#ifdef CS_BUILD_PYTHON
    PyThreadState_Swap(NULL);
    PyThreadState_Clear(myThreadState);
    PyThreadState_Delete(myThreadState);
    PyEval_ReleaseLock();
#endif

    this->mObjectImage = this->mLabelMap->GetLabelImage();
    this->mObjectGraph = this->mLabelMap->GetLabelMapGraph();
    std::cout << "-----number label objects/vertices = " << this->mObjectGraph->GetNumberOfVertices() << ", number label object pairs/edges = " << this->mObjectGraph->GetNumberOfEdges() << std::endl;
}


template< unsigned int VImageDimension > void ObjectBasedSegmentation< VImageDimension >::Update()
{
    ParseParameterContext();

    bool f1 = FilenameParser::ParseFilename(this->mRawImageFullFilename, this->mRawImagePath, this->mRawImageFilename, this->mRawImageFileExtension);
    FilenameParser::ParseFilename(this->mObjectImageFullFilename, this->mObjectImagePath, this->mObjectImageFilename, this->mObjectImageFileExtension);

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    if(!f1) {
        std::cout << "Error: Object-based segmentation: Could not execute pipeline, because input is invalid: " << this->mRawImageFullFilename << std::endl;
        return;
    }

    std::cout << "Object-based segmentation - raw image: " << std::endl;
    std::cout << " dir: " << this->mRawImagePath << std::endl;
    std::cout << " file: " << this->mRawImageFilename << std::endl;
    std::cout << " ext: " << this->mRawImageFileExtension << std::endl;
    std::cout << "Object-based segmentation - object image: " << std::endl;
    std::cout << " dir: " << this->mObjectImagePath << std::endl;
    std::cout << " file: " << this->mObjectImageFilename << std::endl;
    std::cout << " ext: " << this->mObjectImageFileExtension << std::endl;
    std::cout << "Object-based segmentation - object graph: " << std::endl;
    std::cout << " dir: " << this->mObjectGraphFullFilename << std::endl;

    std::cout << "Read images and graph" << std::endl;
    this->ReadOriginalImage(this->mRawImagePath + this->mRawImageFilename + this->mRawImageFileExtension);
    this->ReadObjectImage(this->mObjectImagePath + this->mObjectImageFilename + this->mObjectImageFileExtension);
    this->ReadObjectGraph(this->mObjectGraphFullFilename);
    this->InitLabelMap();
    ActivateFeaturesAtLabelMap();

    switch ( mSegmentationMode )
    {
    case LabelObjects:
    {
        switch ( mSegmentationClassifier )
        {
        case WithoutSegmentation:
        {
            std::cout << "Sorry! Segmentation classifier = 0 -> Nothing to do." << std::endl;
        }
        break;
        case KNNBasedMerging:
        {
            std::cout << "Sorry! Segmentation classifier = 1 -> Nothing to do." << std::endl;
        }
        break;
        case PythonSVMBasedMerging:
        {
            ClassifyLabelObjects();

            std::cout << "Save images and graph" << std::endl;
            this->SaveObjectImage(this->mRawImagePath + mFilenameSave + "_" + SegmentationModeName[mSegmentationMode] + "_" + ClassName[mSegmentationTargetClass] + ".nrrd");
            this->SaveObjectContourOverlayImage(this->mRawImagePath + mFilenameSave + "_" + SegmentationModeName[mSegmentationMode] + "_" + ClassName[mSegmentationTargetClass] + this->mRawImageFileExtension);
            this->SaveObjectGraph(this->mRawImagePath + mFilenameSave + "_" + SegmentationModeName[mSegmentationMode] + "_" + ClassName[mSegmentationTargetClass] + ".txt");
        }
        break;
        default:
            std::cout << "Sorry, unknown segmentation mode!" << std::endl;
        }
    }
    break;
    case LabelObjectPairs:
    {
        std::cout << "Sorry! No segmentation strategy for label-object-pair based segmentation implemented yet" << std::endl;
    }
    break;
    default:
        std::cout << "Sorry, unknown segmentation mode!" << std::endl;
    }

    WriteLogFile(timeStamp);
    WriteDataSetSummary();
}

