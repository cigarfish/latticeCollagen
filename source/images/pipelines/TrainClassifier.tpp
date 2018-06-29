///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  TrainClassifier.tpp                                                  //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-08-01                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "TrainClassifier.h"

#include <itkImage.h>
#include <itkNrrdImageIO.h>

#include <vtkDataSetAttributes.h>
#include <vtkEdgeListIterator.h>
#include <vtkGraphReader.h>
#include <vtkUnsignedLongArray.h>
#include <vtkVariant.h>

#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



template< unsigned int VImageDimension > TrainClassifier< VImageDimension >::TrainClassifier()
{
    mObjectGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
}


template< unsigned int VImageDimension > TrainClassifier< VImageDimension >::~TrainClassifier()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void TrainClassifier< VImageDimension >::ParseParameterContext()
{
    if(this->m_paramContext->findContext("Train Classifier",0)==NULL) {
        std::cout << "Error: Train Classifier: Invalid parameter context" << std::endl;
        return;
    }

    std::string classifier = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Classifier", 0)->dataPointer()) )->currentString();

    mRawImageFullFilename = *(std::string*)(this->m_paramContext->findContext("Data")->findParameter("Image filename", 0)->dataPointer());
    mObjectImageFullFilename = *(std::string*)(this->m_paramContext->findContext("Data")->findParameter("Object dataset filename", 0)->dataPointer());
    mObjectGraphFullFilename = *(std::string*)(this->m_paramContext->findContext("Data")->findParameter("Object graph filename", 0)->dataPointer());
    mTrainingDatasetFilename = *(std::string*)(this->m_paramContext->findContext("Data")->findParameter("Training dataset filename", 0)->dataPointer());
    mTrainingDatabaseFilename = *(std::string*)(this->m_paramContext->findContext("Data")->findParameter("Training database filename", 0)->dataPointer());

    QString filenameQ;
    QFileInfo infoFilenameQ;

    filenameQ = QString::fromStdString(mRawImageFullFilename);
    infoFilenameQ.setFile(filenameQ);

    if(!infoFilenameQ.exists())
        throw std::string("Please specify original image");

    mRawImagePath = (infoFilenameQ.path() + QString("/")).toStdString();
    mRawImageFilename = infoFilenameQ.baseName().toStdString();
    mRawImageFileExtension = (QString(".") + infoFilenameQ.suffix()).toStdString();

    filenameQ = QString::fromStdString(mObjectImageFullFilename);
    infoFilenameQ.setFile(filenameQ);

    if(!infoFilenameQ.exists())
        throw std::string("Please specify object image");

    mObjectImagePath = (infoFilenameQ.path() + QString("/")).toStdString();
    mObjectImageFilename = infoFilenameQ.baseName().toStdString();
    mObjectImageFileExtension = (QString(".") + infoFilenameQ.suffix()).toStdString();

    filenameQ = QString::fromStdString(mObjectGraphFullFilename);
    infoFilenameQ.setFile(filenameQ);

    if(!infoFilenameQ.exists())
        throw std::string("Please specify object graph");

    filenameQ = QString::fromStdString(mTrainingDatasetFilename);
    infoFilenameQ.setFile(filenameQ);

    if(!infoFilenameQ.exists())
        throw std::string("Please specify training dataset");

    filenameQ = QString::fromStdString(mTrainingDatabaseFilename);
    infoFilenameQ.setFile(filenameQ);

    mPath = mRawImagePath;
    if(!infoFilenameQ.exists())
        mTrainingDatabaseFilename = mPath + "trainingDB.csv";

    mWithIncrementalTraining = *(bool*)(this->m_paramContext->findParameter("Incremental training", 0)->dataPointer());

    if( classifier.compare("KNN")==0 )
        mTrainingClass = (ClassType)( (CSParameterChoice*)(this->m_paramContext->findParameter("Class selection", 0)->dataPointer()) )->currentIndex();

    mWithFeature = new bool[LabelMapGraphType::NumberUnaryFeatures+LabelMapGraphType::NumberBinaryFeatures];
    int objI=0, numFeatureComponentsPerPair=0;
    for(objI=0; objI<LabelMapGraphType::NumberUnaryFeatures; objI++) {
        mWithFeature[objI] = *(bool*)(this->m_paramContext->findParameter("With " + LabelMapGraphType::UnaryFeatureName[objI], 0)->dataPointer());
        if(mWithFeature[objI])
            numFeatureComponentsPerPair+=2*LabelMapGraphType::UnaryFeatureComponents[objI];
    }
    for(int j=0; j<LabelMapGraphType::NumberBinaryFeatures; j++) {
        mWithFeature[objI+j] = *(bool*)(this->m_paramContext->findParameter("With " + LabelMapGraphType::BinaryFeatureName[j], 0)->dataPointer());
        if(mWithFeature[objI+j])
            numFeatureComponentsPerPair++;
    }
    if(numFeatureComponentsPerPair > MaxNumberFeaturesPerPair)
        throw std::string("TrainClassifier: Error, number of features exceeds maximal number of features.");
}


template< unsigned int VImageDimension > void TrainClassifier< VImageDimension >::WriteLogFile(std::string timeStamp)
{
//    std::fstream file;
//
//    std::stringstream parameter;
////    mpParamContext->dump(parameter);
//
//    file.open(mParamFile.c_str(), std::ios::out | std::ios::app);
//
////    file << "class name " << typeid(this).name << "\n";
//
//    file << "Begin-object-based-segmentation-------------------------------------------------------------------------------------------------------------------\n";
//    file << "-------------------------time stamp start: " << timeStamp << "\n";
//    file << parameter.str();
//    file << "End-object-based-segmentation---------------------------------------------------------------------------------------------------------------------\n\n";
//
//    file.close();
}


template< unsigned int VImageDimension > void TrainClassifier< VImageDimension >::Init()
{
    //Read object image
    typename IScalarReaderType::Pointer objectReader = IScalarReaderType::New();
    objectReader->SetFileName(mObjectImagePath + mObjectImageFilename + mObjectImageFileExtension);
#if (ITK_VERSION_MAJOR >= 4)
    objectReader->SetImageIO( itk::NrrdImageIO::New() );
#endif
    objectReader->Update();

    typename IMinMaxCalculatorType::Pointer minMaxCalc1 = IMinMaxCalculatorType::New();
    minMaxCalc1->SetImage(objectReader->GetOutput());
    minMaxCalc1->Compute();

    typename CastILImageFilterType::Pointer castFilter = CastILImageFilterType::New();
    castFilter->SetInput(objectReader->GetOutput());
    castFilter->Update();

    mObjectImage = castFilter->GetOutput();
    mObjectImage->DisconnectPipeline();

    //TODO: get rid of this duplication
    castFilter->Update();

    mObjectImage2 = castFilter->GetOutput();
    mObjectImage2->DisconnectPipeline();

    typename LMinMaxCalculatorType::Pointer minMaxCalc2 = LMinMaxCalculatorType::New();
    minMaxCalc2->SetImage(mObjectImage);
    minMaxCalc2->Compute();

    if(minMaxCalc1->GetMaximum() != minMaxCalc2->GetMaximum()) {
        std::cout << "Error: lost superpixels during casting operation" << std::endl;
        return;
    }

    //Read orginal image
    typename CRGBReaderType::Pointer reader = CRGBReaderType::New();
    reader->SetFileName(mRawImagePath + mRawImageFilename + mRawImageFileExtension);
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    reader->Update();

    mOriginalImage = reader->GetOutput();
    mOriginalImage->DisconnectPipeline();
    mOriginalImage->SetSpacing(mObjectImage->GetSpacing());

    //Read object graph
    vtkSmartPointer<vtkGraphReader> graphReader = vtkSmartPointer<vtkGraphReader>::New();
    graphReader->SetFileName(mObjectGraphFullFilename.c_str());
    graphReader->Update();
    graphReader->GetOutput()->ToUndirectedGraph(mObjectGraph);

    typename LabelImageToLabelMapGraphFilterType::Pointer labelImageToLabelMapGraph = LabelImageToLabelMapGraphFilterType::New();           //Build label map from slic label image
    labelImageToLabelMapGraph->SetOriginalImage(mOriginalImage);
    labelImageToLabelMapGraph->SetInput(mObjectImage);
    labelImageToLabelMapGraph->SetLabelImage(mObjectImage2);
    labelImageToLabelMapGraph->UsePrecalculatedGraph(mObjectGraph);
    labelImageToLabelMapGraph->SetBackgroundValue(0);
    labelImageToLabelMapGraph->Update();

    mLabelMapGraph = labelImageToLabelMapGraph->GetOutput();
    mLabelMapGraph->DisconnectPipeline();

    std::cout << "Initialize label objects features" << std::endl;

    int i;
    for(i=0; i<LabelMapGraphType::NumberUnaryFeatures; i++)
        if(mWithFeature[i])
            mLabelMapGraph->ActivateUnaryFeature(static_cast<UnaryFeatureType>(i));

    for(int j=0; j<LabelMapGraphType::NumberBinaryFeatures; j++)
        if(mWithFeature[i+j])
            mLabelMapGraph->ActivateBinaryFeature(static_cast<BinaryFeatureType>(j));

    mObjectGraph = mLabelMapGraph->GetLabelMapGraph();
}


template< unsigned int VImageDimension > void TrainClassifier< VImageDimension >::LoadFromTrainingDatabase(std::map<int, ClassType> &idToClass, SampleTypePointer &sample)
{
//    vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader->SetFileName(mTrainingDatabaseFilename.c_str());
//    reader->SetHaveHeaders(true);
//    reader->DetectNumericColumnsOn();
//    reader->SetFieldDelimiterCharacters(" ");
//    reader->Update();
//
//    vtkSmartPointer<vtkTable> readTable = reader->GetOutput();
//
//    for(unsigned int i=1; i<readTable->GetNumberOfColumns(); i++) {
//        MeasurementVectorType mv;
//        mv[0] = vtkFloatArray::SafeDownCast(readTable->GetColumn(i))->GetValue(2);
//        mv[1] = vtkFloatArray::SafeDownCast(readTable->GetColumn(i))->GetValue(3);
//        mv[2] = vtkFloatArray::SafeDownCast(readTable->GetColumn(i))->GetValue(4);
//
//        idToClass.insert( std::pair<int, ObjectBasedSegmentation::Class>( vtkIntArray::SafeDownCast(readTable->GetColumn(i))->GetValue(1), (ObjectBasedSegmentation::Class)(vtkIntArray::SafeDownCast(readTable->GetColumn(i))->GetValue(0)) ) );
//
//        sample->PushBack(mv);
//    }
}


template< unsigned int VImageDimension > void TrainClassifier< VImageDimension >::SaveToTrainingDatabase(std::set<unsigned long> trainingLabelObjects,
        std::vector< std::pair<unsigned long,unsigned long> > trainingLabelObjectPairs, SampleTypePointer unarySamples, SampleTypePointer binarySamples)
{
    std::fstream file1, file2;
    bool hasHeader = false;

    file1.open((mTrainingDatabaseFilename).c_str(), fstream::in);

    std::string line;
    getline(file1, line);

    std::vector<std::string> wordsInFirstLine;
    std::stringstream firstLineStream(line);
    std::string word;

    while(getline(firstLineStream, word,','))
        wordsInFirstLine.push_back(word);

    int numberUnarySamples = 0;
    int numberBinarySamples = 0;
    if(wordsInFirstLine.size()>=6 && wordsInFirstLine[0].compare("header")==0) {
        hasHeader = true;
        numberUnarySamples = atoi(wordsInFirstLine[2].c_str());
        numberBinarySamples = atoi(wordsInFirstLine[4].c_str());
    }

    if(!hasHeader)
        file2.open((mTrainingDatabaseFilename).c_str(), fstream::out);
    else
        file2.open((mPath + "trainingDB_temp.csv").c_str(), fstream::out);

    file2.flags(fstream::left | fstream::scientific);
    if(!hasHeader) {
        file2 << "header" << ", ";
        file2 << "numberUnarySamples" << ", ";
        file2 << unarySamples->Size() << ", ";
        file2 << "numberBinarySamples" << ", ";
        file2 << binarySamples->Size() << ", ";
        file2 << "numberUnaryFeatures" << ", ";
        file2 << mLabelMapGraph->GetNumberActiveUnaryFeatures() << ", ";
        file2 << "numberUnaryFeatureComponents" << ", ";
        file2 << mLabelMapGraph->GetNumberActiveUnaryFeatureComponents() << ", ";
        file2 << "numberBinaryFeatures" << ", ";
        file2 << mLabelMapGraph->GetNumberActiveBinaryFeatures() << ", ";
        file2 << "numberClasses" << ", ";
        file2 << ObjectBasedSegmentationType::NumberClasses << ", ";
        file2 << "unaryFeatureName/numberComponents" << ", ";
        int objI;
        for(objI=0; objI<LabelMapGraphType::NumberUnaryFeatures; objI++) {
            if(mWithFeature[objI]) {
                file2 << LabelMapGraphType::UnaryFeatureName[objI] << ", ";
                file2 << LabelMapGraphType::UnaryFeatureComponents[objI] << ", ";
            }
        }
        file2 << "binaryFeatureName/numberComponents" << ", ";
        for(int j=0; j<LabelMapGraphType::NumberBinaryFeatures; j++) {
            if(mWithFeature[objI+j]) {
                file2 << LabelMapGraphType::BinaryFeatureName[j] << ", ";
                file2 << 1 << ", ";                                             //LabelMapGraphType::BinaryFeatureComponents[objI] not implemented yet
            }
        }
        file2 << "className" << ", ";
        for(int i=0; i<ObjectBasedSegmentationType::NumberClasses; i++) {
            file2 << ObjectBasedSegmentationType::ClassName[i];
            if(i!=ObjectBasedSegmentationType::NumberClasses-1)
                file2<< ", ";
        }
        file2 << "\n";
    }
    else {
        file2 << "header" << ", ";
        file2 << "numberUnarySamples" << ", ";
        file2 << numberUnarySamples + unarySamples->Size() << ", ";
        file2 << "numberBinarySamples" << ", ";
        file2 << numberBinarySamples + binarySamples->Size() << ", ";
        for(unsigned int i=5; i<wordsInFirstLine.size(); i++) {
            file2 << wordsInFirstLine[i];
            if(i!=wordsInFirstLine.size()-1)
                file2<< ", ";
        }
        file2 << "\n";

        int i=0;
        while(getline(file1, line) && i<numberUnarySamples) {
            file2 << line << "\n";
            i++;
        }
    }

    unsigned int i=0;
    for(std::set<unsigned long>::iterator it=trainingLabelObjects.begin(); it!=trainingLabelObjects.end(); ++it) {
        file2 << *it << ", ";
        for(unsigned int j=0; j<mLabelMapGraph->GetNumberActiveUnaryFeatureComponents(); j++)
            file2 << unarySamples->GetMeasurementVector(i)[j] << ", ";
        file2 << mTrainingClass;
        file2 << "\n";
        i++;
    }

    if(hasHeader)
        while(getline(file1, line))
            file2 << line << "\n";

    for(unsigned int i=0; i<trainingLabelObjectPairs.size(); i++) {
        file2 << trainingLabelObjectPairs[i].first << ", " << trainingLabelObjectPairs[i].second << ", ";
        for(unsigned int j=0; j<mLabelMapGraph->GetNumberActiveBinaryFeatures(); j++)
            file2 << binarySamples->GetMeasurementVector(i)[j] << ", ";
        file2 << mTrainingClass;
        file2 << "\n";
    }
    file1.close();
    file2.close();

    if(hasHeader) {
        std::remove((mTrainingDatabaseFilename).c_str());
        std::rename((mPath + "trainingDB_temp.csv").c_str(), (mTrainingDatabaseFilename).c_str());
    }
}


template< unsigned int VImageDimension > void TrainClassifier< VImageDimension >::RetrieveVerticesAndEdgesFromMaskImage(std::string filename, std::set<unsigned long> &labelObjects, std::vector< std::pair<unsigned long,unsigned long> > &labelObjectPairs)
{
    typename CScalarReaderType::Pointer reader = CScalarReaderType::New();
    reader->SetFileName(filename);
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    reader->Update();

    typedef itk::BinaryImageToLabelMapFilter<CScalarImageType> BinaryImageToLabelMapFilterType;
    typename BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
    binaryImageToLabelMapFilter->SetInput(reader->GetOutput());
    binaryImageToLabelMapFilter->Update();

    std::vector< std::set<unsigned long> > labelsInTrainingRegion;
    for(unsigned int i = 0; i < binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
        typename BinaryImageToLabelMapFilterType::OutputImageType::LabelObjectType* trainingRegion = binaryImageToLabelMapFilter->GetOutput()->GetNthLabelObject(i);

        std::set<unsigned long> lITR;
        for(unsigned int pixelId = 0; pixelId < trainingRegion->Size(); pixelId++) {
//            std::cout << "Object " << i << " contains pixel " << trainingRegion->GetIndex(pixelId) << std::endl;

            lITR.insert(mLabelMapGraph->GetPixel(trainingRegion->GetIndex(pixelId)));
        }
//        std::cout << "Labels in training region " << i << ":" << std::endl;
//        for(std::set<unsigned long>::iterator it=lITR.begin(); it!=lITR.end(); ++it)
//            std::cout << ' ' << *it;
//        std::cout << std::endl;

        labelsInTrainingRegion.push_back(lITR);
    }

    for(unsigned int i=0; i<labelsInTrainingRegion.size(); i++) {
        std::vector< std::pair<unsigned long,unsigned long> > pITR;

        for(std::set<unsigned long>::iterator it1=labelsInTrainingRegion[i].begin(); it1!=labelsInTrainingRegion[i].end(); ++it1) {
            if(labelObjects.count(*it1)==0)
                labelObjects.insert(*it1);

            for(std::set<unsigned long>::iterator it2=it1; it2!=labelsInTrainingRegion[i].end(); ++it2) {
                if(mObjectGraph->GetEdgeId(mObjectGraph->FindVertex(*it1), mObjectGraph->FindVertex(*it2)) != -1)
                    labelObjectPairs.push_back(std::pair<unsigned long,unsigned long>(*it1, *it2));
            }
        }
//        std::cout << "Label-pairs in training region " << i << ":" << std::endl;
//        for(unsigned int j=0; j<pITR.size(); j++)
//            std::cout << ' ' << pITR[j].first << "," << pITR[j].second;
//        for(unsigned int j=0; j<pairsInTrainingRegion.size(); j++)
//            std::cout << ' ' << pairsInTrainingRegion[j].first << "," << pairsInTrainingRegion[j].second;
//        std::cout << std::endl;
    }
}


template< unsigned int VImageDimension > void TrainClassifier< VImageDimension >::ExtractTrainingData(std::set<unsigned long> trainingLabelObjects,
        std::vector< std::pair<unsigned long,unsigned long> > trainingLabelObjectPairs, SampleTypePointer &unarySamples, SampleTypePointer &binarySamples)
{
    unsigned int k=0;
    for(std::set<unsigned long>::iterator it=trainingLabelObjects.begin(); it!=trainingLabelObjects.end(); ++it) {
        MeasurementVectorType mv;
        mv.Fill(0);
        vtkIdType v = mObjectGraph->FindVertex(*it);

        int count = 0;
        for(int j=0; j<LabelMapGraphType::NumberUnaryFeatures; ++j) {
            if(mLabelMapGraph->GetUnaryFeatureState(static_cast<UnaryFeatureType>(j)) == LabelMapGraphType::IsActive) {
                for(unsigned int k=0; k<LabelMapGraphType::UnaryFeatureComponents[j]; k++) {
                    std::stringstream featureName;
                    featureName << LabelMapGraphType::UnaryFeatureName[j];
                    if(LabelMapGraphType::UnaryFeatureComponents[j] > 1)
                        featureName << k;

                    mv[count] = vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(featureName.str().c_str()))->GetValue(v);
                    count++;
                }
            }
        }
        unarySamples->PushBack(mv);
        k++;
    }

    for(unsigned int i=0; i<trainingLabelObjectPairs.size(); i++) {
        MeasurementVectorType mv;
        mv.Fill(0);
        vtkIdType v1 = mObjectGraph->FindVertex(trainingLabelObjectPairs[i].first);
        vtkIdType v2 = mObjectGraph->FindVertex(trainingLabelObjectPairs[i].second);
        vtkIdType e  = mObjectGraph->GetEdgeId(v1, v2);

        for(int j=0; j<LabelMapGraphType::NumberBinaryFeatures; ++j) {
            if(mLabelMapGraph->GetBinaryFeatureState(static_cast<BinaryFeatureType>(j)) == LabelMapGraphType::IsActive)
                mv[j] = vtkFloatArray::SafeDownCast(mObjectGraph->GetEdgeData()->GetArray(LabelMapGraphType::BinaryFeatureName[j].c_str()))->GetValue(e);
        }
        binarySamples->PushBack(mv);
    }
}


template< unsigned int VImageDimension > void TrainClassifier< VImageDimension >::AddTrainingDataToTrainingDatabase()
{
    std::set<unsigned long> labelObjects;
    std::vector< std::pair<unsigned long,unsigned long> > labelObjectPairs;
    RetrieveVerticesAndEdgesFromMaskImage(mTrainingDatasetFilename, labelObjects, labelObjectPairs);

    SampleTypePointer unarySamples = SampleType::New();
    SampleTypePointer binarySamples = SampleType::New();
    ExtractTrainingData(labelObjects, labelObjectPairs, unarySamples, binarySamples);

    SaveToTrainingDatabase(labelObjects, labelObjectPairs, unarySamples, binarySamples);
}


template< unsigned int VImageDimension > void TrainClassifier< VImageDimension >::Update()
{
    ParseParameterContext();
    Init();

    AddTrainingDataToTrainingDatabase();

    //For testing and debug purpose
    vtkSmartPointer<vtkGraphWriter> graphWriter = vtkSmartPointer<vtkGraphWriter>::New();
    graphWriter->SetFileName((mPath + "debugTrainingGraph.txt").c_str());
    graphWriter->SetInput(mObjectGraph);
    graphWriter->Update();
}
