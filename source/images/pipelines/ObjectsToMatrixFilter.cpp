///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ObjectsToMatrixFilter.cpp                                            //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-04-14                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "ObjectsToMatrixFilter.h"

#include <iostream>
#include <fstream>

#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIO.h>
#endif

#include <QDateTime>

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"


ObjectsToMatrixFilter::ObjectsToMatrixFilter()
{
    mRow = 5;
    mCol = 5;

    mDimX = 500;
    mDimY = 500;
    mDimZ = 50;

    mFullyConnected = false;

    mLogFilenameSave = "log";
}


ObjectsToMatrixFilter::~ObjectsToMatrixFilter()
{
    // TODO Auto-generated destructor stub
}


void ObjectsToMatrixFilter::ParseParameterContext()
{
    if(m_paramContext->findContext("Objects to Matrix",0)==NULL) {
        std::cout << "Error: Objects to Matrix Filter: Invalid parameter context" << std::endl;
        return;
    }

    mFullFilename = *(std::string*)(m_paramContext->findParameter("Image filename", 0)->dataPointer());
    mOverlayFullFilename = *(std::string*)(m_paramContext->findParameter("Image overlay filename", 0)->dataPointer());

    mSpacing[0] = *(double*)(m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    mSpacing[1] = *(double*)(m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    mSpacing[2] = *(double*)(m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    mVoxelVolume = mSpacing[0]*mSpacing[1]*mSpacing[2];

    mRow = *(int*)(m_paramContext->findParameter("Matrix rows", 0)->dataPointer());
    mCol = *(int*)(m_paramContext->findParameter("Matrix columns", 0)->dataPointer());
    mDimX = *(int*)(m_paramContext->findParameter("Matrix dimension x", 0)->dataPointer());
    mDimY = *(int*)(m_paramContext->findParameter("Matrix dimension y", 0)->dataPointer());
    mDimZ = *(int*)(m_paramContext->findParameter("Matrix dimension z", 0)->dataPointer());
    mFullyConnected = *(bool*)(m_paramContext->findParameter("Objects fully connected?", 0)->dataPointer());

    mWithFiltering = *(bool*)(m_paramContext->findParameter("Apply filter to objects", 0)->dataPointer());
    mFilterBySlide = *(int*)(m_paramContext->findParameter("Filter by slide", 0)->dataPointer());
    mWriteMetaInfoFile = *(bool*)(m_paramContext->findParameter("Write object identification file", 0)->dataPointer());
}


void ObjectsToMatrixFilter::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

    file.open((mPath + mLogFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-objects-to-matrix-filter----------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-objects-to-matrix-filter------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


void ObjectsToMatrixFilter::WriteMetaFile()
{
    std::fstream file;

    std::stringstream parameter;

    file.open((mPath + "metaInfo.txt").c_str(), std::ios::out);

    file.width(10);
    file << "Column";
    file.width(10);
    file << "Row";
    file.width(10);
    file << "R";
    file.width(10);
    file << "G";
    file.width(10);
    file << "B";
    file.width(10);
    file << "x";
    file.width(10);
    file << "y";
    file.width(10);
    file << "z";
    file.width(20);
    file << "volume";
    file << std::endl;

    for(unsigned int i=0; i<mpMatrixMap->GetNumberOfLabelObjects(); i++) {
        unsigned long label = mpMatrixMap->GetNthLabelObject(i)->GetLabel();

        file.width(10);
        file << (int)mLabelObjectToMatrixPosition[label][0];
        file.width(10);
        file << (int)mLabelObjectToMatrixPosition[label][1];
        file.width(10);
        file << (int)mLabelObjectToColor[label][0];
        file.width(10);
        file << (int)mLabelObjectToColor[label][1];
        file.width(10);
        file << (int)mLabelObjectToColor[label][2];
        file.width(10);
        file << mLabelObjectToCentroid[label][0];
        file.width(10);
        file << mLabelObjectToCentroid[label][1];
        file.width(10);
        file << mLabelObjectToCentroid[label][2];
        file.width(20);
        file << mLabelObjectToVolume[label];
        file << std::endl;
    }

    file.close();
}


void ObjectsToMatrixFilter::SetupMatrixImage(itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> labelMap, itk::SmartPointer<CRGBVoImageType> rgbImage)
{
    itk::Index<3> idx;
    idx[0] = 0;
    idx[1] = 0;
    idx[2] = 0;

    itk::Size<3> size;
    size[0] = mDimX;
    size[1] = mDimY;
    size[2] = mDimZ;

    typedef itk::ImageRegion<3> ImageRegionType;
    ImageRegionType matrixRegion;
    matrixRegion.SetIndex(idx);
    matrixRegion.SetSize(size);

    mpMatrixMap = ShapeLabelMapType::New();
    mpMatrixMap->SetRegions(matrixRegion);
    mpMatrixMap->Allocate();
    mpMatrixMap->SetBackgroundValue(0);

    int colWidth = size[0] / mCol;
    int rowWidth = size[1] / mRow;
    int depth = size[2] / 2;

    itk::Index<3>** objectCentroids;
	objectCentroids = new itk::Index<3>*[mCol];
	for(unsigned int i=0; i<mCol; i++)
		objectCentroids[i] = new itk::Index<3>[mRow];

    for(int i=0; i<mCol; i++) {
        for(int j=0; j<mRow; j++) {
            objectCentroids[i][j][0] = colWidth/2 + i*colWidth;
            objectCentroids[i][j][1] = rowWidth/2 + j*rowWidth;
            objectCentroids[i][j][2] = depth;

//            std::cout << "centroid " << j+i*mRow << "[" << i << "," << j << "] = " << objectCentroids[i][j] << std::endl;
        }
    }

    int cell = mCol*mRow;

    unsigned int stop = cell;
    if(labelMap->GetNumberOfLabelObjects()<stop)
        stop = labelMap->GetNumberOfLabelObjects();

    for(unsigned int i=0; i<stop; i++) {
        ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(i);
        itk::Index<3> centroid;
        centroid[0] = labelObject->GetCentroid()[0];
        centroid[1] = labelObject->GetCentroid()[1];
        centroid[2] = labelObject->GetCentroid()[2];

        int* matrixPos = new int[2];
        matrixPos[0] = i/mRow;
        matrixPos[1] = i%mRow;

        if(mWithFiltering) {
            mLabelObjectToColor.insert(std::pair<unsigned long, CRGBPixelType>(labelObject->GetLabel(), rgbImage->GetPixel(centroid)));
            mLabelObjectToMatrixPosition.insert(std::pair<unsigned long, int*>(labelObject->GetLabel(), matrixPos));
            mLabelObjectToCentroid.insert(std::pair<unsigned long, itk::Index<3> >(labelObject->GetLabel(), centroid));
            mLabelObjectToVolume.insert(std::pair<unsigned long, double >(labelObject->GetLabel(), labelObject->GetNumberOfPixels()*mVoxelVolume));
        }

        itk::Offset<3> offset = objectCentroids[matrixPos[0]][matrixPos[1]] - centroid;
        labelObject->Shift(offset);

//        std::cout << "Shift LabelObject " << i << " from " <<  centroid << " to " << objectCentroids[i/mRow][i%mRow] << " by " << offset << std::endl;

        mpMatrixMap->AddLabelObject(labelObject);
//        std::cout << "New centroid is " << mpMatrixMap->GetNthLabelObject(i)->GetCentroid() << std::endl;
    }
}


void ObjectsToMatrixFilter::SetupBorderRegions(CScalarVoImageType::Pointer &image)
{
    itk::Index<3> indices[6];
    indices[0][0] = 0;  indices[1][0] = 0;
    indices[0][1] = 0;  indices[1][1] = 0;
    indices[0][2] = 0;  indices[1][2] = image->GetLargestPossibleRegion().GetSize()[2]-1;

    indices[2][0] = 0;  indices[3][0] = 0;
    indices[2][1] = 0;  indices[3][1] = image->GetLargestPossibleRegion().GetSize()[1]-1;
    indices[2][2] = 1;  indices[3][2] = 1;

    indices[4][0] = 0;  indices[5][0] = image->GetLargestPossibleRegion().GetSize()[0]-1;
    indices[4][1] = 1;  indices[5][1] = 1;
    indices[4][2] = 1;  indices[5][2] = 1;

    itk::Size<3> sizes[6];
    sizes[0][0] = image->GetLargestPossibleRegion().GetSize()[0];   sizes[1][0] = image->GetLargestPossibleRegion().GetSize()[0];
    sizes[0][1] = image->GetLargestPossibleRegion().GetSize()[1];   sizes[1][1] = image->GetLargestPossibleRegion().GetSize()[1];
    sizes[0][2] = 1;                                                sizes[1][2] = 1;

    sizes[2][0] = image->GetLargestPossibleRegion().GetSize()[0];   sizes[3][0] = image->GetLargestPossibleRegion().GetSize()[0];
    sizes[2][1] = 1;                                                sizes[3][1] = 1;
    sizes[2][2] = image->GetLargestPossibleRegion().GetSize()[2]-2; sizes[3][2] = image->GetLargestPossibleRegion().GetSize()[2]-2;

    sizes[4][0] = 1;                                                sizes[5][0] = 1;
    sizes[4][1] = image->GetLargestPossibleRegion().GetSize()[1]-2; sizes[5][1] = image->GetLargestPossibleRegion().GetSize()[1]-2;
    sizes[4][2] = image->GetLargestPossibleRegion().GetSize()[2]-2; sizes[5][2] = image->GetLargestPossibleRegion().GetSize()[2]-2;


    mDataSetFaces[0].SetIndex(indices[0]);   mDataSetFaces[0].SetSize(sizes[0]);
    mDataSetFaces[1].SetIndex(indices[1]);   mDataSetFaces[1].SetSize(sizes[1]);
    mDataSetFaces[2].SetIndex(indices[2]);   mDataSetFaces[2].SetSize(sizes[2]);
    mDataSetFaces[3].SetIndex(indices[3]);   mDataSetFaces[3].SetSize(sizes[3]);
    mDataSetFaces[4].SetIndex(indices[4]);   mDataSetFaces[4].SetSize(sizes[4]);
    mDataSetFaces[5].SetIndex(indices[5]);   mDataSetFaces[5].SetSize(sizes[5]);

    if(mWithFiltering && mFilterBySlide < (image->GetLargestPossibleRegion().GetSize()[2]-1)) {
        itk::Index<3> index;
        index[0] = 0;
        index[1] = 0;
        index[2] = mFilterBySlide;

        itk::Size<3> size;
        size[0] = image->GetLargestPossibleRegion().GetSize()[0];
        size[1] = image->GetLargestPossibleRegion().GetSize()[1];
        size[2] = 1;

        mSlide.SetIndex(index); mSlide.SetSize(size);
    }
}


void ObjectsToMatrixFilter::RemoveLabelObjectsAtDatasetBorders(LScalarVoImageType::Pointer cellLabelImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap)
{
    std::map<unsigned long, int> cellsToNumberBorderPixel;

    for(unsigned int i=0; i<6; i++) {
        itk::ImageRegionConstIterator<LScalarVoImageType> iter(cellLabelImage, mDataSetFaces[i]);

        while(!iter.IsAtEnd()) {
            if(iter.Value() != 0) {
                if(cellsToNumberBorderPixel.count(iter.Value()) == 0)
                    cellsToNumberBorderPixel.insert( std::pair<unsigned long, int>(iter.Value(), 1) );
                else
                    cellsToNumberBorderPixel[iter.Value()]++;
            }
            ++iter;
        }
    }

    for(std::map<unsigned long, int>::iterator it = cellsToNumberBorderPixel.begin(); it != cellsToNumberBorderPixel.end(); ++it) {
        cellLabelMap->RemoveLabel(it->first);
    }
}


void ObjectsToMatrixFilter::FilterLabelObjectsBySlide(LScalarVoImageType::Pointer cellLabelImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap)
{
    std::map<unsigned long, int> cellsToNumberSlidePixel;

    itk::ImageRegionConstIterator<LScalarVoImageType> iter(cellLabelImage, mSlide);

    while(!iter.IsAtEnd()) {
        if(iter.Value() != 0) {
            if(cellsToNumberSlidePixel.count(iter.Value()) == 0)
                cellsToNumberSlidePixel.insert( std::pair<unsigned long, int>(iter.Value(), 1) );
            else
                cellsToNumberSlidePixel[iter.Value()]++;
        }
        ++iter;
    }

    std::vector<unsigned long> cellsToRemove;
    for(unsigned int i=0; i<cellLabelMap->GetNumberOfLabelObjects(); i++)
        if(cellsToNumberSlidePixel.count(cellLabelMap->GetNthLabelObject(i)->GetLabel()) == 0 || cellsToNumberSlidePixel[cellLabelMap->GetNthLabelObject(i)->GetLabel()] < 5000)
            cellsToRemove.push_back(cellLabelMap->GetNthLabelObject(i)->GetLabel());

    for(std::vector<unsigned long>::iterator it = cellsToRemove.begin(); it != cellsToRemove.end(); ++it) {
        cellLabelMap->RemoveLabel((*it));
    }
}


void ObjectsToMatrixFilter::Update()
{
    ParseParameterContext();

    bool f = FilenameParser::ParseFilename(mFullFilename, mPath, mFilename, mExtension);

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    std::cout << "Objects to Matrix Filter: " << std::endl;
    std::cout << " dir: " << mPath << std::endl;
    std::cout << " file: " << mFilename << std::endl;
    std::cout << " ext: " << mExtension << std::endl;

    if(!f) {
        std::cout << "Error: Objects to Matrix Filter: Could not execute pipeline, because input is invalid" << std::endl;
        return;
    }

    ScalarVoReaderType::Pointer                     reader;
    ImageToShapeLabelMapFilterType::Pointer         imageToShapeLabelMap;
    LabelMapToLabelImageFilterType::Pointer         labelMapToLabelImage;
    ThresholdFilterType::Pointer                    labelImageToBin;
    ScalarVoWriterType::Pointer                     writer;

    //----------READER--------------------------------------------------------------------------------------------------------
    reader = ScalarVoReaderType::New();
    reader->SetFileName(mPath + mFilename + mExtension);
    reader->ReleaseDataFlagOn();
    reader->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    reader->Update();

    itk::SmartPointer<CScalarVoImageType> image = reader->GetOutput();
    image->DisconnectPipeline();

    itk::SmartPointer<CRGBVoImageType> mRGBImage = CRGBVoImageType::New();
    if(mWriteMetaInfoFile) {
        CRGBVoReaderType::Pointer rgbReader = CRGBVoReaderType::New();
        rgbReader->SetFileName(mOverlayFullFilename.c_str());
#if (ITK_VERSION_MAJOR >= 4)
        rgbReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
        rgbReader->Update();

        mRGBImage = rgbReader->GetOutput();
        mRGBImage->DisconnectPipeline();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER---CREATE-LABEL-MAP--------------------------------------------------------------------------------------
    imageToShapeLabelMap = ImageToShapeLabelMapFilterType::New();
    imageToShapeLabelMap->SetInput(image);
    imageToShapeLabelMap->SetFullyConnected(mFullyConnected);
    imageToShapeLabelMap->Update();

    itk::SmartPointer<ShapeLabelMapType> labelMap = imageToShapeLabelMap->GetOutput();
    labelMap->DisconnectPipeline();

    labelMapToLabelImage = LabelMapToLabelImageFilterType::New();
    labelMapToLabelImage->SetInput(labelMap);
    labelMapToLabelImage->Update();

    itk::SmartPointer<LScalarVoImageType> labelImage = labelMapToLabelImage->GetOutput();
    labelImage->DisconnectPipeline();
    //--------------------------------------------------------------------------------------------------------------------

    //----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------
    SetupBorderRegions(image);
    RemoveLabelObjectsAtDatasetBorders(labelImage, labelMap);
    if(mWithFiltering)  FilterLabelObjectsBySlide(labelImage, labelMap);

    SetupMatrixImage(labelMap, mRGBImage);

    if(mWithFiltering && mWriteMetaInfoFile)
        WriteMetaFile();

    //----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
    labelMapToLabelImage = LabelMapToLabelImageFilterType::New();
    labelMapToLabelImage->SetInput(mpMatrixMap);
    labelMapToLabelImage->Update();

    itk::SmartPointer<LScalarVoImageType> matrixLabelImage = labelMapToLabelImage->GetOutput();
    matrixLabelImage->DisconnectPipeline();

    labelImageToBin = ThresholdFilterType::New();
    labelImageToBin->SetOutsideValue(0);
    labelImageToBin->SetInsideValue(itk::NumericTraits<CScalarPixelType>::max());
    labelImageToBin->SetLowerThreshold(1);
    labelImageToBin->SetUpperThreshold(itk::NumericTraits<LabelMapToLabelImageFilterType::OutputImagePixelType>::max());
    labelImageToBin->SetInput(matrixLabelImage);
    //-------------------------------------------------------------------------------------------------------------------------

    //----------WRITER--------------------------------------------------------------------------------------------------------
    writer = ScalarVoWriterType::New();
    writer->SetFileName(mPath + mFilename + "_matrix" + mExtension);
    writer->SetInput(labelImageToBin->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();
    //-------------------------------------------------------------------------------------------------------------------------
}

