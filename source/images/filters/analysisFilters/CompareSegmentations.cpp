///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CompareSegmentations.cpp                                             //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-04-16                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "CompareSegmentations.h"

#include <math.h>
#include <fstream>
#include <string>
#include <sstream>

#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>

#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIOFactory.h>
#include <itkTIFFImageIO.h>
#endif // (ITK_VERSION_MAJOR >= 4)

#include "../../../tools/parameters/CSParameter.h"
#include "../../../tools/parameters/CSParameterContext.h"
#include "../../../tools/input/FilenameParser.h"



CompareSegmentations::CompareSegmentations()
{
    // TODO Auto-generated constructor stub

}


CompareSegmentations::~CompareSegmentations()
{
    // TODO Auto-generated destructor stub
}


void CompareSegmentations::ParseParameterContext()
{
    if(mpParamContext->findContext("Compare Segmentations",0)==NULL) {
        std::cout << "Error: AnalyzeCellsFilter: Invalid parameter context" << std::endl;
        return;
    }

    mDatasetID = *(int*)(mpParamContext->findParameter("Dataset ID", 0)->dataPointer());
    mGoldstandardDatasetFullFilename = *(std::string*)(mpParamContext->findParameter("Goldstandard segmentation filename", 0)->dataPointer());
    mEvaluationDatasetFullFilename = *(std::string*)(mpParamContext->findParameter("Evaluation segmentation filename", 0)->dataPointer());

    std::cout << "dataset id = " << mDatasetID << std::endl;
}


void CompareSegmentations::WriteDataFile()
{
    std::fstream file1, file2, tempfile;
    bool hasHeader = false;

    file1.open((mGoldstandardDatasetPath + "../" + "Segmentation_comparison_file.txt").c_str(), std::fstream::in);

    std::string line;
    getline(file1, line);

    if(line.find("dataSet")!=std::string::npos) {
        hasHeader = true;

        tempfile.open((mGoldstandardDatasetPath + "../" + "Segmentation_comparison_tempfile.txt").c_str(), std::fstream::out);
        tempfile << line;

        while(!file1.eof()) {
            getline(file1, line);
            if(line.find(mDatasetID)!=0) {
                tempfile << std::endl;
                tempfile << line;
            }
        }
        file1.close();
        tempfile.close();

        std::remove((mGoldstandardDatasetPath + "../" + "Segmentation_comparison_file.txt").c_str());
        std::rename((mGoldstandardDatasetPath + "../" + "Segmentation_comparison_tempfile.txt").c_str(), (mGoldstandardDatasetPath + "../" + "Segmentation_comparison_file.txt").c_str());
    }
    else
        file1.close();

    if(!hasHeader)
        file2.open((mGoldstandardDatasetPath + "../" + "Segmentation_comparison_file.txt").c_str(), std::fstream::out);
    else
        file2.open((mGoldstandardDatasetPath + "../" + "Segmentation_comparison_file.txt").c_str(), std::fstream::out | std::fstream::app);

    file2.flags(std::fstream::left | std::fstream::scientific);
    if(!hasHeader) {
        file2.width(40);
        file2 << "dataSet";
        file2.width(25);
        file2 << "goldstandardLabel";
        file2.width(25);
        file2 << "evaluationLabel";
        file2.width(20);
        file2 << "volumeOverlap";
        file2.width(20);
        file2 << "volumeGoldExcess";
        file2.width(20);
        file2 << "volumeEvalExcess";
		file2.width(15);
        file2 << "isGoldInner";
		file2.width(15);
        file2 << "isEvalInner";
        file2 << std::endl;
    }

    for(std::map<unsigned long, ObjectPairs>::iterator it = mObjectPairs.begin(); it != mObjectPairs.end(); ++it) {
        file2.width(40);
        file2 << mDatasetID;
        file2.width(25);
        file2 << it->second.mGoldLabel;
        file2.width(25);
        file2 << it->second.mEvalLabel;
        file2.width(20);
        file2 << it->second.mOverlapVolume;
        file2.width(20);
        file2 << it->second.mGoldExcess;
        file2.width(20);
        file2 << it->second.mEvalExcess;
		file2.width(15);
		file2 << it->second.mGoldInner;
		file2.width(15);
		file2 << it->second.mEvalInner;
        file2 << std::endl;
    }

    file2.close();
}


void CompareSegmentations::Update()
{
    ParseParameterContext();

    bool f1 = FilenameParser::ParseFilename(mGoldstandardDatasetFullFilename, mGoldstandardDatasetPath, mGoldstandardDatasetName, mGoldstandardDatasetExtension);
    bool f2 = FilenameParser::ParseFilename(mEvaluationDatasetFullFilename, mEvaluationDatasetPath, mEvaluationDatasetName, mEvaluationDatasetFile);

    std::string timeStamp;

    if(!f1 || !f2) {
        std::cout << "Error: CompareSegmentations: Could not execute pipeline, because input is invalid! " << std::endl;
        return;
    }

    std::cout << "Start segmentation comparison: " << std::endl;


    itk::SmartPointer<ScalarVoReaderType>               goldstandardReader;
    itk::SmartPointer<ScalarVoReaderType>               evaluationReader;
    itk::SmartPointer<ImageToShapeLabelMapFilterType>   imageToShaLabMapFilter;
    itk::SmartPointer<LabelMapToLabelImageFilterType>   labelMapToLabelImageFilter;


    //----------READER---------------------------------------------------------------------------------------------------------
    goldstandardReader = ScalarVoReaderType::New();
    goldstandardReader->SetFileName(mGoldstandardDatasetPath + mGoldstandardDatasetName + mGoldstandardDatasetExtension);
#if (ITK_VERSION_MAJOR >= 4)
    goldstandardReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    goldstandardReader->Update();

    evaluationReader = ScalarVoReaderType::New();
    evaluationReader->SetFileName(mEvaluationDatasetPath + mEvaluationDatasetName + mEvaluationDatasetFile);
#if (ITK_VERSION_MAJOR >= 4)
    evaluationReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    evaluationReader->Update();


    itk::SmartPointer<CScalarVoImageType> goldstandardImage = goldstandardReader->GetOutput();
    goldstandardImage->DisconnectPipeline();

    itk::SmartPointer<CScalarVoImageType> evaluationImage = evaluationReader->GetOutput();
    evaluationImage->DisconnectPipeline();

    RGBPixelType background;
    background.Fill(0);

    itk::SmartPointer<RGBVoImageType> validationImage1 = RGBVoImageType::New();
    validationImage1->SetRegions(goldstandardImage->GetLargestPossibleRegion());
    validationImage1->Allocate();
    validationImage1->FillBuffer(background);

    itk::SmartPointer<RGBVoImageType> validationImage2 = RGBVoImageType::New();
    validationImage2->SetRegions(goldstandardImage->GetLargestPossibleRegion());
    validationImage2->Allocate();
    validationImage2->FillBuffer(background);
    //-------------------------------------------------------------------------------------------------------------------------

    //----------COMPUTE-CELL-VOLUME-AND-BOUNDARY-CONTACT-AREA------------------------------------------------------------------
    imageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    imageToShaLabMapFilter->SetInput(goldstandardImage);
    imageToShaLabMapFilter->SetFullyConnected(false);
    imageToShaLabMapFilter->Update();

    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> goldstandardLabelMap = imageToShaLabMapFilter->GetOutput();
    goldstandardLabelMap->DisconnectPipeline();
    goldstandardLabelMap->Optimize();

    labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
    labelMapToLabelImageFilter->SetInput(goldstandardLabelMap);
    labelMapToLabelImageFilter->Update();

    itk::SmartPointer<LScalarVoImageType> goldstandardLabelImage = labelMapToLabelImageFilter->GetOutput();
    goldstandardLabelImage->DisconnectPipeline();

    imageToShaLabMapFilter->SetInput(evaluationImage);
    imageToShaLabMapFilter->SetFullyConnected(false);
    imageToShaLabMapFilter->Update();

    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> evaluationLabelMap = imageToShaLabMapFilter->GetOutput();
    evaluationLabelMap->DisconnectPipeline();
    evaluationLabelMap->Optimize();

    labelMapToLabelImageFilter->SetInput(evaluationLabelMap);
    labelMapToLabelImageFilter->Update();

    itk::SmartPointer<LScalarVoImageType> evaluationLabelImage = labelMapToLabelImageFilter->GetOutput();
    evaluationLabelImage->DisconnectPipeline();

    std::map<unsigned long, bool> goldLabelIsInside;
    std::map<unsigned long, bool> evalLabelIsInside;

    std::multimap<unsigned long, unsigned long> evalToGoldAlign;

    //only pairs are considered, where centroid of goldstandard label object is within evaluation label object
    for(unsigned int i=0; i<goldstandardLabelMap->GetNumberOfLabelObjects(); i++) {
        ImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::Pointer goldLabelObject = goldstandardLabelMap->GetNthLabelObject(i);

        itk::Index<3> idx;
        idx[0] = goldLabelObject->GetCentroid()[0];
        idx[1] = goldLabelObject->GetCentroid()[1];
        idx[2] = goldLabelObject->GetCentroid()[2];

        ObjectPairs pair;
        pair.mGoldLabel = goldLabelObject->GetLabel();
        pair.mEvalLabel = evaluationLabelImage->GetPixel(idx);

        if(pair.mGoldLabel!=0 && pair.mEvalLabel!=0) {
            pair.mOverlapVolume = 0;
            pair.mGoldExcess = 0;
            pair.mEvalExcess = 0;
			pair.mGoldInner = true;
			pair.mEvalInner = true;

            mObjectPairs.insert(std::pair<unsigned long, ObjectPairs>(pair.mGoldLabel, pair));

            evalToGoldAlign.insert(std::pair<unsigned long, unsigned long>(pair.mEvalLabel, pair.mGoldLabel));
            goldLabelIsInside.insert(std::pair<unsigned long, bool>(goldLabelObject->GetLabel(), true));
        }
    }

    for(unsigned int i=0; i<evaluationLabelMap->GetNumberOfLabelObjects(); i++)
        if(evaluationLabelMap->GetNthLabelObject(i)->GetLabel()!=0)
            evalLabelIsInside.insert(std::pair<unsigned long, bool>(evaluationLabelMap->GetNthLabelObject(i)->GetLabel(), true));

    itk::ImageRegionConstIterator<LScalarVoImageType> goldIter(goldstandardLabelImage, goldstandardLabelImage->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<LScalarVoImageType> evalIter(evaluationLabelImage, evaluationLabelImage->GetLargestPossibleRegion());

	itk::Index<3> upperCorner = goldstandardLabelImage->GetLargestPossibleRegion().GetUpperIndex();

    while(!goldIter.IsAtEnd()) {
        itk::Index<3> currIdx = goldIter.GetIndex();

        if(mObjectPairs.count(goldIter.Get())!=0)
            if(mObjectPairs[goldIter.Get()].mEvalLabel == evalIter.Get())
                mObjectPairs[goldIter.Get()].mOverlapVolume++;

        if(currIdx[0] == 0 || currIdx[1] == 0 || currIdx[2] == 0 || currIdx[0] == upperCorner[0] || currIdx[1] == upperCorner[1] || currIdx[2] == upperCorner[2]) {
            if(goldIter.Get()!=0)   goldLabelIsInside[goldIter.Get()] = false;
            if(evalIter.Get()!=0)   evalLabelIsInside[evalIter.Get()] = false;
        }
        ++goldIter;
        ++evalIter;
    }

    for(std::map<unsigned long, ObjectPairs>::iterator it = mObjectPairs.begin(); it != mObjectPairs.end(); ++it) {
        if(!goldLabelIsInside[(*it).second.mGoldLabel])
            (*it).second.mGoldInner = false;

        if(!evalLabelIsInside[(*it).second.mEvalLabel])
            (*it).second.mEvalInner = false;
    }

    RGBPixelType overlapInside;
    overlapInside.SetRed(255); overlapInside.SetGreen(255); overlapInside.SetBlue(255);
    RGBPixelType overlapOutside;
    overlapOutside.SetRed(80); overlapOutside.SetGreen(80); overlapOutside.SetBlue(80);
    RGBPixelType goldExcess;
    goldExcess.SetRed(255); goldExcess.SetGreen(255); goldExcess.SetBlue(0);
    RGBPixelType evalExcess;
    evalExcess.SetRed(255); evalExcess.SetGreen(0); evalExcess.SetBlue(0);

    itk::ImageRegionIterator<RGBVoImageType> val1Iter(validationImage1, validationImage1->GetLargestPossibleRegion());
    itk::ImageRegionIterator<RGBVoImageType> val2Iter(validationImage2, validationImage2->GetLargestPossibleRegion());

    goldIter.GoToBegin();
    evalIter.GoToBegin();
    while(!goldIter.IsAtEnd()) {
        if(mObjectPairs.count(goldIter.Get())!=0) {
            if(mObjectPairs[goldIter.Get()].mEvalLabel == evalIter.Get()) {
                if(goldLabelIsInside[goldIter.Get()] && evalLabelIsInside[evalIter.Get()]) {
                    val1Iter.Set(overlapInside);
                    val2Iter.Set(overlapInside);
                }
                else {
                    val1Iter.Set(overlapOutside);
                    val2Iter.Set(overlapOutside);
                }
            }
        }
        ++goldIter;
        ++evalIter;
        ++val1Iter;
        ++val2Iter;
    }
    itk::SmartPointer<RGBVoWriterType> overlapWriter = RGBVoWriterType::New();
    overlapWriter->SetInput(validationImage1);
    overlapWriter->SetFileName(mGoldstandardDatasetPath + "comparison_step1" + mGoldstandardDatasetExtension);
#if (ITK_VERSION_MAJOR >= 4)
    overlapWriter->SetImageIO( itk::TIFFImageIO::New() );
#endif
    overlapWriter->Update();

    goldIter.GoToBegin();
    evalIter.GoToBegin();
    val1Iter.GoToBegin();
    val2Iter.GoToBegin();
    while(!goldIter.IsAtEnd()) {
        if(mObjectPairs.count(goldIter.Get())!=0)
            if(mObjectPairs[goldIter.Get()].mEvalLabel != evalIter.Get())
                if(val1Iter.Get()!=overlapOutside)
                    val1Iter.Set(goldExcess);

        if(evalToGoldAlign.count(evalIter.Get())!=0) {
            for(std::multimap<unsigned long, unsigned long>::iterator it=evalToGoldAlign.equal_range(evalIter.Get()).first; it!=evalToGoldAlign.equal_range(evalIter.Get()).second; ++it) {
                if((*it).second != goldIter.Get()) {
                    if(val2Iter.Get()!=overlapOutside) {
                        val2Iter.Set(evalExcess);
                    }
                }
            }
        }

        ++goldIter;
        ++evalIter;
        ++val1Iter;
        ++val2Iter;
    }

    for(std::map<unsigned long, ObjectPairs>::iterator it = mObjectPairs.begin(); it != mObjectPairs.end(); ++it) {
        ImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::Pointer goldLabelObject = goldstandardLabelMap->GetLabelObject((*it).second.mGoldLabel);
        ImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::Pointer evalLabelObject = evaluationLabelMap->GetLabelObject((*it).second.mEvalLabel);

        (*it).second.mGoldExcess = goldLabelObject->Size() - (*it).second.mOverlapVolume;
        (*it).second.mEvalExcess = evalLabelObject->Size() - (*it).second.mOverlapVolume;
    }
    std::cout << "checkpoint 5" << std::endl;

    overlapWriter->SetInput(validationImage1);
    overlapWriter->SetFileName(mGoldstandardDatasetPath + "comparison_validation1" + mGoldstandardDatasetExtension);
#if (ITK_VERSION_MAJOR >= 4)
    overlapWriter->SetImageIO( itk::TIFFImageIO::New() );
#endif
    overlapWriter->Update();

    overlapWriter->SetInput(validationImage2);
    overlapWriter->SetFileName(mGoldstandardDatasetPath + "comparison_validation2" + mGoldstandardDatasetExtension);
#if (ITK_VERSION_MAJOR >= 4)
    overlapWriter->SetImageIO( itk::TIFFImageIO::New() );
#endif
    overlapWriter->Update();

    WriteDataFile();
}
