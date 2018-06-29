///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SuperpixelPreparation.tpp                                            //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-03-31                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "SuperpixelPreparation.h"

#include <limits>
#include <iostream>

#include <QDateTime>

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



template< unsigned int VImageDimension > SuperpixelPreparation< VImageDimension >::SuperpixelPreparation()
{
}


template< unsigned int VImageDimension > SuperpixelPreparation< VImageDimension >::~SuperpixelPreparation()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void SuperpixelPreparation< VImageDimension >::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    this->m_paramContext->dump(parameter);

    file.open((this->mRawImagePath + mLogFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-superpixel-preparation----------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-superpixel-preparation------------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


template< unsigned int VImageDimension > void SuperpixelPreparation< VImageDimension >::ParseParameterContext()
{
    if(this->m_paramContext->findContext("Superpixel Preparation",0)==NULL) {
        std::cout << "Error: SuperpixelPreparation: Invalid parameter context" << std::endl;
        return;
    }

    this->mRawImageFullFilename = *(std::string*)(this->m_paramContext->findContext("Reader")->findParameter("Image filename", 0)->dataPointer());

    this->mVoxelSpacing[0] = *(double*)(this->m_paramContext->findParameter("Image voxel spacing x", 0)->dataPointer());
    this->mVoxelSpacing[1] = *(double*)(this->m_paramContext->findParameter("Image voxel spacing y", 0)->dataPointer());
    if(ImageDimension==3)
        this->mVoxelSpacing[2] = *(double*)(this->m_paramContext->findParameter("Image voxel spacing z", 0)->dataPointer());


    std::string inputMode = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Start configuration", 0)->dataPointer()) )->currentString();

    if( inputMode.compare("With Superpixel-Partitioning")==0 ) {
        mEntrypoint = 0;

        std::string superpixelMode = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Superpixel mode", 0)->dataPointer()) )->currentString();
        if( superpixelMode.compare("Specify size")==0 )
            mSuperpixelBySize = true;
        else if( superpixelMode.compare("Specify number")==0 )
            mSuperpixelBySize = false;
        else {
            std::cout << "Error: SuperpixelPreparation: Unknown superpixel mode: " << superpixelMode << std::endl;
            return;
        }

        if(mSuperpixelBySize) {
            mSuperpixelSpacing[0] = *(double*)(this->m_paramContext->findParameter("Superpixel spacing x", 0)->dataPointer());
            mSuperpixelSpacing[1] = *(double*)(this->m_paramContext->findParameter("Superpixel spacing y", 0)->dataPointer());
            if(ImageDimension==3)
                mSuperpixelSpacing[2] = *(double*)(this->m_paramContext->findParameter("Superpixel spacing z", 0)->dataPointer());
        }
        else
            mSuperpixelNumber = *(unsigned int*)(this->m_paramContext->findParameter("Superpixel number", 0)->dataPointer());
    }
    else if( inputMode.compare("Add Features only")==0 ) {
        mEntrypoint = 1;

        this->mObjectImageFullFilename = *(std::string*)(this->m_paramContext->findContext("Superpixel Partitioning Reader")->findParameter("Object dataset filename", 0)->dataPointer());
        this->mObjectGraphFullFilename = *(std::string*)(this->m_paramContext->findContext("Superpixel Partitioning Reader")->findParameter("Object graph filename", 0)->dataPointer());
    }
    else {
        std::cout << "Error: SuperpixelPreparation: Unknown start configuration mode: " << inputMode << std::endl;
        return;
    }

    mWithFeature = new bool[LabelMapType::NumberUnaryFeatures+LabelMapType::NumberBinaryFeatures];
    int objI=0;
    for(objI=0; objI<LabelMapType::NumberUnaryFeatures; objI++)
        mWithFeature[objI] = *(bool*)(this->m_paramContext->findParameter("With " + LabelMapType::UnaryFeatureName[objI], 0)->dataPointer());
    for(int j=0; j<LabelMapType::NumberBinaryFeatures; j++)
        mWithFeature[objI+j] = *(bool*)(this->m_paramContext->findParameter("With " + LabelMapType::BinaryFeatureName[j], 0)->dataPointer());

    //TODO: color for superpixel outline

    mLogFilenameSave = *(std::string*)(this->m_paramContext->findParameter("Log file name", 0)->dataPointer());
    mFilenameSave = *(std::string*)(this->m_paramContext->findParameter("Save prefix", 0)->dataPointer());
}


template< unsigned int VImageDimension > void SuperpixelPreparation< VImageDimension >::AddFeaturesToLabelMap()
{
    std::cout << "Set label object feature status" << std::endl;

    int i;
    for(i=0; i<LabelMapType::NumberUnaryFeatures; i++)
        if(mWithFeature[i])
            this->mLabelMap->AttachUnaryFeature(static_cast<UnaryFeatureType>(i));

    for(int j=0; j<LabelMapType::NumberBinaryFeatures; j++)
        if(mWithFeature[i+j])
            this->mLabelMap->AttachBinaryFeature(static_cast<BinaryFeatureType>(j));

    this->mLabelMap->InitializeAllLabelObjects(false);
}


template< unsigned int VImageDimension > void SuperpixelPreparation< VImageDimension >::Update()
{
    ParseParameterContext();

    bool f1 = FilenameParser::ParseFilename(this->mRawImageFullFilename, this->mRawImagePath, this->mRawImageFilename, this->mRawImageFileExtension);
    FilenameParser::ParseFilename(this->mObjectImageFullFilename, this->mObjectImagePath, this->mObjectImageFilename, this->mObjectImageFileExtension);

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    if(!f1) {
        std::cout << "Error: SuperpixelPreparation: Could not execute pipeline, because input is invalid: " << this->mRawImageFullFilename << std::endl;
        return;
    }

    std::cout << "SuperpixelPreparation - raw image: " << std::endl;
    std::cout << " dir: " << this->mRawImagePath << std::endl;
    std::cout << " file: " << this->mRawImageFilename << std::endl;
    std::cout << " ext: " << this->mRawImageFileExtension << std::endl;
    std::cout << "SuperpixelPreparation - object image: " << std::endl;
    std::cout << " dir: " << this->mObjectImagePath << std::endl;
    std::cout << " file: " << this->mObjectImageFilename << std::endl;
    std::cout << " ext: " << this->mObjectImageFileExtension << std::endl;
    std::cout << "SuperpixelPreparation - object graph: " << std::endl;
    std::cout << " dir: " << this->mObjectGraphFullFilename << std::endl;


    this->ReadOriginalImage(this->mRawImagePath + this->mRawImageFilename + this->mRawImageFileExtension);

    if(mEntrypoint==0) {
        typename SLICFilterType::Pointer slicFilter = SLICFilterType::New();
        slicFilter->SetReleaseDataFlag(false);
        slicFilter->SetInput(this->mOriginalImage);
        if(mSuperpixelBySize)   slicFilter->SetSuperpixelSpacing(this->mVoxelSpacing, mSuperpixelSpacing);
        else                    slicFilter->SetNumberSuperPixel(mSuperpixelNumber);
        slicFilter->Update();

        this->mObjectImage = slicFilter->GetOutput();
        this->mObjectImage->DisconnectPipeline();
        this->mObjectImage->SetSpacing(this->mVoxelSpacing);

        this->SaveObjectImage(this->mRawImagePath + mFilenameSave + "_labelObjectImage" + ".nrrd");     //a bit ugly, but needed atm, cause we need two instances of mObjectImage (second instance is mObjectImage2)
        this->ReadObjectImage(this->mRawImagePath + mFilenameSave + "_labelObjectImage" + ".nrrd");

        this->mObjectGraph = LabelImageToGraphFilterType::LabelImageToGraph(this->mObjectImage);
    }
    if(mEntrypoint==1) {
        this->ReadObjectImage(this->mObjectImagePath + this->mObjectImageFilename + this->mObjectImageFileExtension);
        this->ReadObjectGraph(this->mObjectGraphFullFilename);
    }

    this->InitLabelMap();
    AddFeaturesToLabelMap();

    this->mObjectImage = this->mLabelMap->GetLabelImage();
    this->mObjectGraph = this->mLabelMap->GetLabelMapGraph();

    this->SaveObjectGraph(this->mRawImagePath + mFilenameSave + "_labelObjectGraph" + ".txt");
    this->SaveObjectImage(this->mRawImagePath + mFilenameSave + "_labelObjectImage" + ".nrrd");
    this->SaveObjectContourOverlayImage(this->mRawImagePath + mFilenameSave + "_labelObjectOverlay" + this->mRawImageFileExtension);

    WriteLogFile(timeStamp);
}
