///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeLobuleFilter.cpp                                              //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-03-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "AnalyzeLobuleFilter.h"

#include <math.h>
#include <fstream>
#include <string>
#include <sstream>

#include <QDateTime>

#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIOFactory.h>
#include <itkTIFFImageIO.h>
#endif // (ITK_VERSION_MAJOR >= 4)

#include <vtkDataSetAttributes.h>
#include <vtkDelimitedTextReader.h>
#include <vtkDelimitedTextWriter.h>
#include <vtkEdgeListIterator.h>
#include <vtkGraphWriter.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkTable.h>
#include <vtkUnsignedLongArray.h>
#include <vtkVariantArray.h>

#include "../../tools/GraphAnnotationHelper.h"
#include "../../../tools/parameters/CSParameter.h"
#include "../../../tools/parameters/CSParameterContext.h"


#define MMToMicron              1000



AnalyzeLobuleFilter::AnalyzeLobuleFilter()
{
    mDataSetID = "testDataSet";

    mVoxelVolume = 1;
    mDatasetVolume = 1;
    mCVVolume = 0;                  //TODO input method for cv/pv similar to analyzeBileNetworkFilter
    mPVVolume = 0;

    mSpacing[0] = 1.;
    mSpacing[1] = 1.;
    mSpacing[2] = 1.;

    mSize[0] = 1024;
    mSize[1] = 1024;
    mSize[2] = 100;

    mMaxCorrectionDistanceToBorder = 4.*23.;    //within this distance to dataset border no lobule catchment basin correction happens (i.e. distance to closest vein behind dataset border)
    mMinAxisDistToVein = 2.*23.;
    mCenterRUSCRadius = 4.*23.;                 //radius of vein-axes ellipses: 4 * cell_diameter

    mSaveEverything = false;

#if (ITK_VERSION_MAJOR >= 4)
    itk::TIFFImageIOFactory::RegisterOneFactory();
#endif // (ITK_VERSION_MAJOR >= 4)
}


AnalyzeLobuleFilter::~AnalyzeLobuleFilter()
{
    // TODO Auto-generated destructor stub
}


void AnalyzeLobuleFilter::WriteLogFile(std::string timeStamp)
{
    //Log not implemented yet
}


void AnalyzeLobuleFilter::ParseParameterContext()
{
    if(this->m_paramContext->findContext("Analyze Lobule",0)==NULL) {
        std::cout << "Error: AnalyzeLobuleFilter: Invalid parameter context" << std::endl;
        return;
    }

    mDataSetFullFilenameCV = *(std::string*)(this->m_paramContext->findParameter("Central vein segmentation", 0)->dataPointer());
    mDataSetFullFilenamePV = *(std::string*)(this->m_paramContext->findParameter("Portal vein segmentation", 0)->dataPointer());

    mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenameCV) );
    mCVFileExists = mInfoFullFilename.exists();

    mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenamePV) );
    mPVFileExists = mInfoFullFilename.exists();

    if(!mPVFileExists && mCVFileExists)
        mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenameCV) );

    mPath = (mInfoFullFilename.path() + QString("/")).toStdString();
    mFilenameExtension = (QString(".") + mInfoFullFilename.suffix()).toStdString();

    mSpacing[0] = *(double*)(this->m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    mSpacing[1] = *(double*)(this->m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    mSpacing[2] = *(double*)(this->m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    std::string saveMode = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        mSaveEverything = true;
    else if( saveMode.compare("Save only essentials")==0 )
        mSaveEverything = false;
}


void AnalyzeLobuleFilter::ComputeRelativeLobularPositions(FScalarImagePointerType cvDistMap, FScalarImagePointerType lobuleDistMap)
{
    typename AddFImageFilterType::Pointer addDist;
    typename DivideImageFilterType::Pointer divide;

    addDist = AddFImageFilterType::New();
    addDist->ReleaseDataFlagOn();
    addDist->SetInput1(cvDistMap);
    addDist->SetInput2(lobuleDistMap);

    divide = DivideImageFilterType::New();
    divide->ReleaseDataFlagOn();
    divide->SetInput1(cvDistMap);
    divide->SetInput2(addDist->GetOutput());
    divide->Update();

    mRelativeLobularPositionMap = divide->GetOutput();
    mRelativeLobularPositionMap->DisconnectPipeline();
}


void AnalyzeLobuleFilter::CollectBasicImageInformation()
{
    std::string pathToVeinImage = mDataSetFullFilenameCV;
    if(!mCVFileExists)
        pathToVeinImage = mDataSetFullFilenamePV;

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO((pathToVeinImage).c_str(), itk::ImageIOFactory::ReadMode);
    if(imageIO->CanReadFile((pathToVeinImage).c_str())) {
        imageIO->ReadImageInformation();

        mVoxelVolume = 1;
        mDatasetVolume = 1;

        for(unsigned int i=0; i<ImageDimension; i++) {
            mSize[i] = imageIO->GetDimensions(i);
            mVoxelVolume *= mSpacing[i];
            mDatasetVolume *= mSize[i]*mSpacing[i];
        }
        int dim = ImageDimension;
        mVoxelVolume /= std::pow(MMToMicron, dim);
        mDatasetVolume /= std::pow(MMToMicron, dim);
    }
//    mEffectiveDatasetVolume = mDatasetVolume;
}


void AnalyzeLobuleFilter::CollectBasicLobuleInformation()
{
//    typename ImageToShapeLabelMapFilterType::Pointer lobuleImageToShaLabMapFilter;
//    typename LabelMapToLabelImageFilterType::Pointer lobuleShaLabMapToLabelImageFilter;
//
//
//    typename LScalarImageType::Pointer lobuleImage;
//    if(mLobulesFileExists) {
//        lobuleImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
//        lobuleImageToShaLabMapFilter->FullyConnectedOff();
//        lobuleImageToShaLabMapFilter->SetInput(mLobuleBin);
//        lobuleImageToShaLabMapFilter->Update();
//
//        lobuleShaLabMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
//        lobuleShaLabMapToLabelImageFilter->SetInput(lobuleImageToShaLabMapFilter->GetOutput());
//        lobuleShaLabMapToLabelImageFilter->Update();
//
//        lobuleImage = lobuleShaLabMapToLabelImageFilter->GetOutput();
//        lobuleImage->DisconnectPipeline();
//
//        for(unsigned int i=0; i<mNumberLobules; i++) {
//            std::vector<double> sliceVolume = lobuleToSliceVolume[lobuleImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetLabel()];
//
//            LobulesAnalysisContainer n;
//            n.label = i;
//            n.volume = lobuleImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetPhysicalSize();
//
//            mLobules.push_back(n);
//        }
//
//    }
}


void AnalyzeLobuleFilter::WriteGraph(vtkSmartPointer<vtkMutableUndirectedGraph> graph, std::string filename)
{
    vtkSmartPointer<vtkGraphWriter> writer = vtkSmartPointer<vtkGraphWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(graph);
    writer->Update();
}


void AnalyzeLobuleFilter::WriteBasicInformationFile()
{
//    std::string suffixVolORArea = "Area";
//    if(ImageDimension == 3) suffixVolORArea = "Volume";
//
//    vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader->SetFileName((mPath + "../" + "Nuclei_Output_file_1.txt").c_str());
//    reader->SetHaveHeaders(true);
//    reader->DetectNumericColumnsOn();
//    reader->SetFieldDelimiterCharacters(" ");
//    reader->Update();
//
//    vtkSmartPointer<vtkTable> readTable = reader->GetOutput();
//
//    vtkSmartPointer<vtkVariantArray> head =  vtkVariantArray::SafeDownCast(readTable->GetColumnByName("header"));
//    if(head==NULL) {
//        vtkSmartPointer<vtkVariantArray> header = vtkSmartPointer<vtkVariantArray>::New();
//        header->InsertNextValue( vtkVariant("numDim") );
//        header->InsertNextValue( vtkVariant("dimX") );
//        header->InsertNextValue( vtkVariant("dimY") );
//        if(ImageDimension == 3) header->InsertNextValue( vtkVariant("dimZ") );
//        header->InsertNextValue( vtkVariant("dataSet"+suffixVolORArea) );
//        header->InsertNextValue( vtkVariant("effective"+suffixVolORArea) );
//        header->InsertNextValue( vtkVariant("tissue"+suffixVolORArea) );
//        header->InsertNextValue( vtkVariant("void"+suffixVolORArea) );
//        header->InsertNextValue( vtkVariant("centralVein"+suffixVolORArea) );
//        header->InsertNextValue( vtkVariant("portalVein"+suffixVolORArea) );
//        header->InsertNextValue( vtkVariant("necReg"+suffixVolORArea) );
//        header->InsertNextValue( vtkVariant("analysed"+suffixVolORArea) );
//        header->InsertNextValue( vtkVariant("analysedLobules") );
//        header->InsertNextValue( vtkVariant("numAnalysedNonProlifHepaticNuclei") );
//        header->InsertNextValue( vtkVariant("numAnalysedProlifHepaticNuclei") );
//        header->InsertNextValue( vtkVariant("numAnalysedNonProlifNonHepaticNuclei") );
//        header->InsertNextValue( vtkVariant("numAnalysedProlifNonHepaticNuclei") );
//        header->SetName("header");
//
//        readTable->AddColumn(header);
//    }
//
//    vtkSmartPointer<vtkVariantArray> dataSet = vtkSmartPointer<vtkVariantArray>::New();
//    dataSet->InsertNextValue( vtkVariant(ImageDimension) );
//    dataSet->InsertNextValue( vtkVariant(mSize[0]) );
//    dataSet->InsertNextValue( vtkVariant(mSize[1]) );
//    if(ImageDimension == 3) dataSet->InsertNextValue( vtkVariant(mSize[2]) );
//    dataSet->InsertNextValue( vtkVariant((double)mDatasetVolume) );
//    dataSet->InsertNextValue( vtkVariant((double)mEffectiveDatasetVolume) );
//    dataSet->InsertNextValue( vtkVariant((double)(mTissueVolume)) );
//    dataSet->InsertNextValue( vtkVariant((double)(mVoidVolume)) );
//    dataSet->InsertNextValue( vtkVariant((double)(mCVVolume)) );
//    dataSet->InsertNextValue( vtkVariant((double)(mPVVolume)) );
//    dataSet->InsertNextValue( vtkVariant((double)(mNecRegVolume)) );
//    dataSet->InsertNextValue( vtkVariant((double)(mAnalysedDatasetVolume)) );
//    dataSet->InsertNextValue( vtkVariant((double)(mNumberLobules)) );
//    dataSet->InsertNextValue( vtkVariant(mNumAnalysedNonProlifHepaticNuclei) );
//    dataSet->InsertNextValue( vtkVariant(mNumAnalysedProlifHepaticNuclei) );
//    dataSet->InsertNextValue( vtkVariant(mNumAnalysedNonProlifNonHepaticNuclei) );
//    dataSet->InsertNextValue( vtkVariant(mNumAnalysedProlifNonHepaticNuclei) );
//    dataSet->SetName(mDataSetID.c_str());
//
//    readTable->AddColumn(dataSet);
//
//    readTable->Dump(25);
//
//    vtkSmartPointer<vtkDelimitedTextWriter> writer = vtkSmartPointer<vtkDelimitedTextWriter>::New();
//    writer->SetFileName((mPath + "../" + "Nuclei_Output_file_1.txt").c_str());
//    writer->SetFieldDelimiter(" ");
//    writer->SetInput(readTable);
//    writer->Update();
}


//void AnalyzeLobuleFilter::WriteLobuleInformationFile(std::string filenamePostfix, std::vector<LobuleAnalysisContainer> lobules)
//{
//    std::string suffixVolORArea = "area";
//    if(ImageDimension == 3) suffixVolORArea = "volume";
//
//    std::fstream file1, file2, tempfile;
//    bool hasHeader = false;
//
//    file1.open((mPath + "../Lobule_Output_file" + filenamePostfix + ".txt").c_str(), fstream::in);
//
//    std::string line;
//    getline(file1, line);
//
//    if(line.find("dataset")!=std::string::npos) {
//        hasHeader = true;
//
//        tempfile.open((mPath + "../Lobule_Output_file_tempfile" + filenamePostfix + ".txt").c_str(), fstream::out);
//        tempfile << line;
//
//        while(!file1.eof()) {
//            getline(file1, line);
//            if(line.find(mDataSetID)!=0) {
//                tempfile << std::endl;
//                tempfile << line;
//            }
//        }
//        file1.close();
//        tempfile.close();
//
//        std::remove((mPath + "../Lobule_Output_file" + filenamePostfix + ".txt").c_str());
//        std::rename((mPath + "../Lobule_Output_file_tempfile" + filenamePostfix + ".txt").c_str(), (mPath + "../Lobule_Output_file" + filenamePostfix + ".txt").c_str());
//    }
//    else
//        file1.close();
//
//    if(!hasHeader)
//        file2.open((mPath + "../Lobule_Output_file" + filenamePostfix + ".txt").c_str(), fstream::out);
//    else
//        file2.open((mPath + "../Lobule_Output_file" + filenamePostfix + ".txt").c_str(), fstream::out | fstream::app);
//
//    file2.flags(fstream::left | fstream::scientific);
//    if(!hasHeader) {
//        file2.width(40);
//        file2 << "dataset";
//        file2.width(10);
//        file2 << "label";
//        file2.width(20);
//        file2 << suffixVolORArea;
//        file2.width(30);
//        file2 << "numCV";
//        file2.width(30);
//        file2 << "numPV" << std::endl;
//    }
//
//    for(unsigned int i=0; i<lobules.size(); i++) {
//        file2.width(40);
//        file2 << mDataSetID;
//        file2.width(10);
//        file2 << lobules[i].label;
//        file2.width(20);
//        file2 << lobules[i].volume;
//        file2.width(30);
//        file2 << lobules[i].numCV;
//        file2.width(30);
//        file2 << lobules[i].numPV << std::endl;
//    }
//
//    file2.close();
//}


void AnalyzeLobuleFilter::Decompose3DImageInto2DImages(CScalarImagePointerType inImage,
        std::vector<int> slicesToCollect, std::map<int, CScalar2DImagePointerType> &outSlices)
{
    CScalar2DIndexType destIndex;
    destIndex.Fill(0);

    RegionType sourceRegion = inImage->GetLargestPossibleRegion();
    sourceRegion.SetSize(2, 0);

    CScalar2DRegionType destRegion;
    CScalar2DSizeType destSize;
    destSize[0] = sourceRegion.GetSize()[0];
    destSize[1] = sourceRegion.GetSize()[1];

    destRegion.SetSize(destSize);
    destRegion.SetIndex(destIndex);

    for(unsigned int i=0; i<slicesToCollect.size(); i++) {
        sourceRegion.SetIndex(2, slicesToCollect[i]);

        typename ExtractImageFilterType::Pointer extractFilter = ExtractImageFilterType::New();
        extractFilter->SetNumberOfThreads(4);
        extractFilter->SetInput(inImage);
        extractFilter->SetExtractionRegion(sourceRegion);
#if ITK_VERSION_MAJOR >= 4
        extractFilter->SetDirectionCollapseToIdentity(); // This is required.
#endif
        extractFilter->Update();

        outSlices.insert(std::pair<int, CScalar2DImagePointerType>(slicesToCollect[i], extractFilter->GetOutput()));
    }
}


void AnalyzeLobuleFilter::ComputeVeinCenterPoints(int zStepWidth, CScalarImagePointerType veinImage, CScalarImagePointerType veinLabelImage, vtkSmartPointer<vtkMutableUndirectedGraph> &veinGraph,
        std::map<int, std::map<int, CScalar2DIndexType> > &veinCentersPerSlice, std::map<vtkIdType, int> &vertexIdToVeinState, int veinState)
{
    vtkSmartPointer<vtkPoints> points = veinGraph->GetPoints();
    vtkSmartPointer<vtkIntArray> labelId = vtkIntArray::SafeDownCast(veinGraph->GetVertexData()->GetArray("label id"));

    if(labelId==NULL) {
        labelId = vtkSmartPointer<vtkIntArray>::New();
        labelId->SetNumberOfComponents(1);
        labelId->SetName("label id");
    }

    std::map<int, CScalar2DImagePointerType> slicesToVisit;
    std::vector<int> slicesToCollect;

    double zSize = veinImage->GetLargestPossibleRegion().GetSize(2);
    int stepWidth = zSize / (double)ceil(zSize / zStepWidth);               //correct step width to fit uniformly in number of z-slices
    std::cout << "corrected stepWidth = " << stepWidth << std::endl;

    slicesToCollect.push_back(0);
    int j=1;
    while(slicesToCollect[slicesToCollect.size()-1] <= zSize-1) {
        slicesToCollect.push_back((stepWidth*j)-1);
        j++;
    }
    slicesToCollect.pop_back();
    if(slicesToCollect[slicesToCollect.size()-1] != zSize-1)
        slicesToCollect.push_back(zSize-1);

    Decompose3DImageInto2DImages(veinLabelImage, slicesToCollect, slicesToVisit);

    std::map<int, vtkIdType> labelIdToVertexFromPrevSlice;
    for(std::map<int, CScalar2DImagePointerType>::iterator jt=slicesToVisit.begin(); jt!=slicesToVisit.end(); ++jt) {
        std::map<int, CScalar2DIndexType> veinIdToCenterPoints;

        ComputeVeinCenterPointsIn2DImage(jt->second, veinIdToCenterPoints);

        for(std::map<int, CScalar2DIndexType>::iterator it=veinIdToCenterPoints.begin(); it!=veinIdToCenterPoints.end(); ++it) {
            vtkIdType newVertId = veinGraph->AddVertex();
            points->InsertNextPoint(it->second[0], it->second[1], jt->first);
            labelId->InsertValue(newVertId, it->first);
            vertexIdToVeinState.insert(std::pair<vtkIdType, int>(newVertId, veinState));

            if(labelIdToVertexFromPrevSlice.count(it->first)>0) {
                vtkIdType prevVert = labelIdToVertexFromPrevSlice[it->first];
                veinGraph->AddEdge(newVertId, labelIdToVertexFromPrevSlice[it->first]);
            }
            labelIdToVertexFromPrevSlice[it->first] = newVertId;
        }
        veinCentersPerSlice.insert(std::pair<int, std::map<int, CScalar2DIndexType> >(jt->first, veinIdToCenterPoints));
    }
    veinGraph->SetPoints(points);
    veinGraph->GetVertexData()->AddArray(labelId);
}


void AnalyzeLobuleFilter::PrepareVeinData(std::string loadFilename, std::string path, std::string saveFilename, std::string ext, CScalarImagePointerType &bin,
        LabelMapPointerType &labelMap, std::map<int, std::map<int, CScalar2DIndexType> > &veinCentersPerSlice, std::map<vtkIdType, int> &vertexIdToVeinState, int veinState)
{
    typename BinaryImageToShapeLabelMapFilterType::Pointer      imageToShapeLabelMap;
    typename LabelMapToLabelImageFilterType::Pointer            labelMapToImage;
    typename CScalarImageType::Pointer                          labelImage;
    typename CScalarImageWriterType::Pointer                    writer;

    std::stringstream savename;

    this->ReadImage(loadFilename, bin, mSpacing);

    imageToShapeLabelMap = BinaryImageToShapeLabelMapFilterType::New();
    imageToShapeLabelMap->SetNumberOfThreads(4);
    imageToShapeLabelMap->ReleaseDataFlagOn();
    imageToShapeLabelMap->FullyConnectedOff();
    imageToShapeLabelMap->SetInput(bin);
    imageToShapeLabelMap->Update();

    labelMap = imageToShapeLabelMap->GetOutput();
    labelMap->DisconnectPipeline();
    std::cout << "number " << loadFilename << " = " << labelMap->GetNumberOfLabelObjects() << std::endl;

    labelMapToImage = LabelMapToLabelImageFilterType::New();
    labelMapToImage->SetInput(labelMap);
    labelMapToImage->Update();

    labelImage = labelMapToImage->GetOutput();
    labelImage->DisconnectPipeline();

    ComputeVeinCenterPoints(1, bin, labelImage, mVeinSkeletonGraph, veinCentersPerSlice, vertexIdToVeinState, veinState);

    if(mSaveEverything) {
        savename << path << saveFilename << ext;

        writer = CScalarImageWriterType::New();
        writer->SetFileName(savename.str());
        writer->SetInput(labelImage);
#if (ITK_VERSION_MAJOR >= 4)
        writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer->Update();
    }
}


void AnalyzeLobuleFilter::SyncVeinAndCatchmentLabels(LabelMapPointerType veinLabelMap, LabelMapPointerType &catchmentLabelMap)
{
    for(unsigned int i=0; i<veinLabelMap->GetNumberOfLabelObjects(); i++)
        catchmentLabelMap->GetLabelObject( veinLabelMap->GetNthLabelObject(i)->GetIndex(0) )->SetLabel( veinLabelMap->GetNthLabelObject(i)->GetLabel() );
}


void AnalyzeLobuleFilter::CorrectVeinCatchmentBasin(FScalarImagePointerType &distanceMap)
{
    typename CScalarImageType::Pointer artBoundVeins = CScalarImageType::New();
    artBoundVeins->SetRegions(distanceMap->GetLargestPossibleRegion());
    artBoundVeins->Allocate();
    artBoundVeins->FillBuffer(0);
    artBoundVeins->SetSpacing(mSpacing);

    typename itk::ImageRegionIterator<CScalarImageType> iter(artBoundVeins, artBoundVeins->GetLargestPossibleRegion());

    while(!iter.IsAtEnd()) {
        if(iter.GetIndex()[0] == 0 || iter.GetIndex()[0] == artBoundVeins->GetLargestPossibleRegion().GetSize()[0]-1)
            if(distanceMap->GetPixel(iter.GetIndex()) >= mMaxCorrectionDistanceToBorder)
                iter.Set(255);
        if(iter.GetIndex()[1] == 0 || iter.GetIndex()[1] == artBoundVeins->GetLargestPossibleRegion().GetSize()[1]-1)
            if(distanceMap->GetPixel(iter.GetIndex()) >= mMaxCorrectionDistanceToBorder)
                iter.Set(255);
        ++iter;
    }

    typename SignedMaurerDistanceMapImageFilterType::Pointer distanceMapFilter = SignedMaurerDistanceMapImageFilterType::New();
    distanceMapFilter->SetNumberOfThreads(8);
    distanceMapFilter->UseImageSpacingOn();
    distanceMapFilter->SetInput(artBoundVeins);
    distanceMapFilter->InsideIsPositiveOff();
    distanceMapFilter->Update();

    FScalarImagePointerType distToBoundaryVeins = distanceMapFilter->GetOutput();
    distToBoundaryVeins->DisconnectPipeline();

    typename itk::ImageRegionIterator<FScalarImageType> iter2(distanceMap, distanceMap->GetLargestPossibleRegion());
    typename itk::ImageRegionIterator<FScalarImageType> iter3(distToBoundaryVeins, distToBoundaryVeins->GetLargestPossibleRegion());

    while(!iter2.IsAtEnd() || !iter3.IsAtEnd()) {
        if(iter3.Value() < iter2.Value())
            iter2.Set(-1.);

        ++iter2;
        ++iter3;
    }
}


void AnalyzeLobuleFilter::ComputeVeinCatchmentBasins(std::string path, std::string suffix, std::string ext,
        CScalarImagePointerType veinBin, LabelMapPointerType veinLabelMap, CScalarImagePointerType &veinCatchLabelImage, FScalarImagePointerType &veinDistmap)
{
    typename SignedMaurerDistanceMapImageFilterType::Pointer    distanceMapFilter;
    typename MaskNegImageFilterType1::Pointer                   maskNegFilter;
    typename MorphoWatershedImageFilterType::Pointer            morphWatershed;
    typename LabelImageToShapeLabelMapFilterType::Pointer       labelImageToLabelMap;
    typename LabelMapToLabelImageFilterType::Pointer            labelMapToLabelImage;
    typename FImageCalculatorFilterType::Pointer                imageCalculatorFilter;
    typename MaskImageFilterType2::Pointer                      maskFilter;
    typename ThresholdFilterType2::Pointer                      threshold;
    typename RescaleImageFilterType::Pointer                    rescaler;
    typename CScalarImageWriterType::Pointer                    writer;

    distanceMapFilter = SignedMaurerDistanceMapImageFilterType::New();
    distanceMapFilter->SetNumberOfThreads(8);
    distanceMapFilter->UseImageSpacingOn();
    distanceMapFilter->SetInput(veinBin);
    distanceMapFilter->InsideIsPositiveOff();

    maskNegFilter = MaskNegImageFilterType1::New();
    maskNegFilter->SetInput(distanceMapFilter->GetOutput());
    maskNegFilter->SetMaskImage(veinBin);
    maskNegFilter->SetOutsideValue(0.0);
    maskNegFilter->Update();

    typename FScalarImageType::Pointer distanceMap = maskNegFilter->GetOutput();
    distanceMap->DisconnectPipeline();

    morphWatershed = MorphoWatershedImageFilterType::New();
    morphWatershed->ReleaseDataFlagOn();
    morphWatershed->SetLevel(1);
    morphWatershed->SetMarkWatershedLine(false);
    morphWatershed->FullyConnectedOn();
    morphWatershed->SetInput(distanceMap);
    morphWatershed->Update();

    typename CScalarImageType::Pointer catchmentLabelImage = morphWatershed->GetOutput();
    catchmentLabelImage->DisconnectPipeline();

    labelImageToLabelMap = LabelImageToShapeLabelMapFilterType::New();
    labelImageToLabelMap->SetInput(catchmentLabelImage);
    labelImageToLabelMap->Update();

    LabelMapPointerType catchmentLabelMap = labelImageToLabelMap->GetOutput();
    catchmentLabelMap->DisconnectPipeline();

    std::cout << suffix << " number of catchment basins = " << catchmentLabelMap->GetNumberOfLabelObjects() << std::endl;

    SyncVeinAndCatchmentLabels(veinLabelMap, catchmentLabelMap);

    labelMapToLabelImage = LabelMapToLabelImageFilterType::New();
    labelMapToLabelImage->SetInput(catchmentLabelMap);
    labelMapToLabelImage->Update();

    catchmentLabelImage = labelMapToLabelImage->GetOutput();
    catchmentLabelImage->DisconnectPipeline();

    if(mSaveEverything) {
        std::stringstream s1;
        s1 << path << suffix << "catchLabelImage" << ext;

        writer = CScalarImageWriterType::New();
        writer->SetFileName(s1.str());
        writer->SetInput(catchmentLabelImage);
#if (ITK_VERSION_MAJOR >= 4)
        writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer->Update();
    }

    typename FScalarImageType::Pointer distanceMapCorr = distanceMap;

    CorrectVeinCatchmentBasin(distanceMapCorr);

    veinDistmap = distanceMapCorr;

    imageCalculatorFilter = FImageCalculatorFilterType::New ();
    imageCalculatorFilter->SetImage(distanceMapCorr);
    imageCalculatorFilter->Compute();

    threshold = ThresholdFilterType2::New();
    threshold->SetOutsideValue(0);
    threshold->SetInsideValue(255);
    threshold->SetLowerThreshold(1);
    threshold->SetUpperThreshold(imageCalculatorFilter->GetMaximum());
    threshold->SetInput(distanceMapCorr);

    maskFilter = MaskImageFilterType2::New();
    maskFilter->SetInput(catchmentLabelImage);
    maskFilter->SetMaskImage(threshold->GetOutput());
    maskFilter->SetOutsideValue(0.0);
    maskFilter->Update();

    veinCatchLabelImage = maskFilter->GetOutput();
    veinCatchLabelImage->DisconnectPipeline();

    if(mSaveEverything) {
        std::stringstream s2;
        s2 << path << suffix << "distMapBorderCorrection" << ext;

        rescaler = RescaleImageFilterType::New();
        rescaler->SetInput(distanceMapCorr);

        writer = CScalarImageWriterType::New();
        writer->SetFileName(s2.str());
        writer->SetInput(rescaler->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer->Update();

        std::stringstream s3;
        s3 << path << suffix << "corrCatchLabelImage" << ext;

        writer->SetFileName(s3.str());
        writer->SetInput(veinCatchLabelImage);
#if (ITK_VERSION_MAJOR >= 4)
        writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer->Update();
    }
}


//bool AnalyzeLobuleFilter::IsAxisValid(int axisType, vtkIdType start, vtkIdType end,
//        std::map<int, CScalarImagePointerType> cvMasks, std::map<int, CScalarImagePointerType> pvMasks)
//{
//    vtkSmartPointer<vtkIntArray> labelId = vtkIntArray::SafeDownCast(mVeinAxesGraph->GetVertexData()->GetArray("label id"));
//
//    StructuringElementType structuringElement;
//    structuringElement.SetRadius(1);
//    structuringElement.CreateStructuringElement();
//
//    typename DilateImageFilterType::Pointer dilateFilter = DilateImageFilterType::New();
//    if(axisType==0)         dilateFilter->SetInput(cvMasks[labelId->GetValue(start)]);
//    else if(axisType==1)    dilateFilter->SetInput(pvMasks[labelId->GetValue(start)]);
//    dilateFilter->SetKernel(structuringElement);
//
//    typename MaskImageFilterType2::Pointer maskFilter = MaskImageFilterType2::New();
//    if(axisType<2)  maskFilter->SetInput(dilateFilter->GetOutput());
//    else            maskFilter->SetInput(cvMasks[labelId->GetValue(start)]);
//    if(axisType==0) maskFilter->SetMaskImage(cvMasks[labelId->GetValue(end)]);
//    else            maskFilter->SetMaskImage(pvMasks[labelId->GetValue(end)]);
//    maskFilter->SetOutsideValue(0.0);
//    maskFilter->Update();
//
//    typename CImageCalculatorFilterType::Pointer imageCalculatorFilter = CImageCalculatorFilterType::New ();
//    imageCalculatorFilter->SetImage(maskFilter->GetOutput());
//    imageCalculatorFilter->Compute();
//
//    if(imageCalculatorFilter->GetMaximum()>0)
//        return true;
//
//    return false;
//}


bool AnalyzeLobuleFilter::IsAxisValid(int axisType, vtkIdType start, vtkIdType end,
        std::map<int, CScalarImagePointerType> cvMasks, std::map<int, CScalarImagePointerType> pvMasks)
{
    bool isValid = false;

    vtkSmartPointer<vtkIntArray> labelId = vtkIntArray::SafeDownCast(mVeinAxesGraph->GetVertexData()->GetArray("label id"));
    vtkSmartPointer<vtkPoints> points = mVeinAxesGraph->GetPoints();

    StructuringElementType structuringElement;
    structuringElement.SetRadius(1);
    structuringElement.CreateStructuringElement();

    typename DilateImageFilterType::Pointer dilateFilter = DilateImageFilterType::New();
    if(axisType==0)         dilateFilter->SetInput(cvMasks[labelId->GetValue(start)]);
    else if(axisType==1)    dilateFilter->SetInput(pvMasks[labelId->GetValue(start)]);
    dilateFilter->SetKernel(structuringElement);

    typename MaskImageFilterType2::Pointer maskFilter = MaskImageFilterType2::New();
    if(axisType<2)  maskFilter->SetInput(dilateFilter->GetOutput());
    else            maskFilter->SetInput(cvMasks[labelId->GetValue(start)]);
    if(axisType==0) maskFilter->SetMaskImage(cvMasks[labelId->GetValue(end)]);
    else            maskFilter->SetMaskImage(pvMasks[labelId->GetValue(end)]);
    maskFilter->SetOutsideValue(0.0);
    maskFilter->Update();

    typename CImageCalculatorFilterType::Pointer imageCalculatorFilter = CImageCalculatorFilterType::New ();
    imageCalculatorFilter->SetImage(maskFilter->GetOutput());
    imageCalculatorFilter->Compute();

    if(imageCalculatorFilter->GetMaximum()>0)
        isValid = true;

    if(isValid) {
        IndexType startIdx, endIdx;
        startIdx[0] = points->GetPoint(start)[0]; startIdx[1] = points->GetPoint(start)[1]; startIdx[2] = points->GetPoint(start)[2];
        endIdx[0] = points->GetPoint(end)[0]; endIdx[1] = points->GetPoint(end)[1]; endIdx[2] = points->GetPoint(end)[2];

        std::cout << "test axis (axistype=" << axisType << ") from " << startIdx << " to " << endIdx << std::endl;

        FLineIteratorType it1(mCVDistmap, startIdx, endIdx);
        CLineIteratorType it2(mCVCatchLabelImage, startIdx, endIdx);
        FLineIteratorType it3(mPVDistmap, startIdx, endIdx);
        CLineIteratorType it4(mPVCatchLabelImage, startIdx, endIdx);
        while(!it1.IsAtEnd() || !it2.IsAtEnd() || !it3.IsAtEnd() || !it4.IsAtEnd()) {
            if(axisType==0) {
                std::cout << "start labelId " << labelId->GetValue(start) << " end labelId " << labelId->GetValue(end) << " it2.Value() " << (int)it2.Value() << std::endl;
                std::cout << "it1.Value() " << (float)it1.Value() << " it3.Value() " << (float)it3.Value() << std::endl;

                if(it2.Value() != 0 && it2.Value() != labelId->GetValue(start) && it2.Value() != labelId->GetValue(end))
                    if(it1.Value() < mMinAxisDistToVein)
                        isValid = false;
                if(it3.Value() != -1 && it3.Value() < mMinAxisDistToVein)
                    isValid = false;
            }
            else if(axisType==1) {
                if(it4.Value() != 0 && it4.Value() != labelId->GetValue(start) && it4.Value() != labelId->GetValue(end))
                    if(it3.Value() < mMinAxisDistToVein)
                        isValid = false;
                if(it1.Value() != -1 && it1.Value() < mMinAxisDistToVein)
                    isValid = false;
            }
            else if(axisType==2) {
                if(it2.Value() != 0 && it2.Value() != labelId->GetValue(start))
                    if(it1.Value() < mMinAxisDistToVein)
                        isValid = false;
                if(it4.Value() != 0 && it4.Value() != labelId->GetValue(end))
                    if(it3.Value() < mMinAxisDistToVein)
                        isValid = false;
            }
            if(!isValid)
                break;

            ++it1; ++it2; ++it3; ++it4;
        }
    }

    return isValid;
}


void AnalyzeLobuleFilter::EstablishVeinAxes()
{
    mVeinAxesGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    std::map<vtkIdType, int> vertIdToState;
    std::map<vtkIdType, int> vertIdToLabelId;
    std::map<int, vtkIdType> cvLabelIdToVertId, pvLabelIdToVertId;
    std::multimap<vtkIdType, vtkIdType> cvcvAxes, pvpvAxes, cvpvAxes;
    std::map<vtkIdType, int> edgeState;
    std::map<vtkIdType, int> edgeId;

    GraphAnnotationHelper *anno = new GraphAnnotationHelper();

    for(unsigned int i=0; i<mCVLabelMap->GetNumberOfLabelObjects(); i++) {
        vtkIdType newVertId = mVeinAxesGraph->AddVertex();
        IndexType idx;
        mCVLabelMap->TransformPhysicalPointToIndex(mCVLabelMap->GetNthLabelObject(i)->GetCentroid(), idx);
        points->InsertNextPoint(idx[0], idx[1], idx[2]);
        vertIdToState.insert(std::pair<vtkIdType, int>(newVertId, 0));
        vertIdToLabelId.insert(std::pair<vtkIdType, int>(newVertId, mCVLabelMap->GetNthLabelObject(i)->GetLabel()));
        cvLabelIdToVertId.insert(std::pair<int, vtkIdType>(mCVLabelMap->GetNthLabelObject(i)->GetLabel(), newVertId));
    }
    for(unsigned int i=0; i<mPVLabelMap->GetNumberOfLabelObjects(); i++) {
        vtkIdType newVertId = mVeinAxesGraph->AddVertex();
        IndexType idx;
        mPVLabelMap->TransformPhysicalPointToIndex(mPVLabelMap->GetNthLabelObject(i)->GetCentroid(), idx);
        points->InsertNextPoint(idx[0], idx[1], idx[2]);
        vertIdToState.insert(std::pair<vtkIdType, int>(newVertId, 1));
        vertIdToLabelId.insert(std::pair<vtkIdType, int>(newVertId, mPVLabelMap->GetNthLabelObject(i)->GetLabel()));
        pvLabelIdToVertId.insert(std::pair<int, vtkIdType>(mPVLabelMap->GetNthLabelObject(i)->GetLabel(), newVertId));
    }
    mVeinAxesGraph->SetPoints(points);
    anno->EnableVertexIDAnnotation();
    anno->AddPredefinedAnnotations(mVeinAxesGraph);
    anno->AddCustomVertexAnnotation(mVeinAxesGraph, "vein state", vertIdToState, -100);
    anno->AddCustomVertexAnnotation(mVeinAxesGraph, "label id", vertIdToLabelId, -100);

    for(unsigned int i=0; i<mCVLabelMap->GetNumberOfLabelObjects(); i++)
        for(unsigned int j=i+1; j<mCVLabelMap->GetNumberOfLabelObjects(); j++)
            cvcvAxes.insert(std::pair<vtkIdType, vtkIdType>(cvLabelIdToVertId[mCVLabelMap->GetNthLabelObject(i)->GetLabel()], cvLabelIdToVertId[mCVLabelMap->GetNthLabelObject(j)->GetLabel()]));
    std::cout << "cvcvAxes.size() = " << cvcvAxes.size() << std::endl;

    for(unsigned int i=0; i<mPVLabelMap->GetNumberOfLabelObjects(); i++)
        for(unsigned int j=i+1; j<mPVLabelMap->GetNumberOfLabelObjects(); j++)
            pvpvAxes.insert(std::pair<vtkIdType, vtkIdType>(pvLabelIdToVertId[mPVLabelMap->GetNthLabelObject(i)->GetLabel()], pvLabelIdToVertId[mPVLabelMap->GetNthLabelObject(j)->GetLabel()]));
    std::cout << "pvpvAxes.size() = " << pvpvAxes.size() << std::endl;

    for(unsigned int i=0; i<mCVLabelMap->GetNumberOfLabelObjects(); i++)
        for(unsigned int j=0; j<mPVLabelMap->GetNumberOfLabelObjects(); j++)
            cvpvAxes.insert(std::pair<vtkIdType, vtkIdType>(cvLabelIdToVertId[mCVLabelMap->GetNthLabelObject(i)->GetLabel()], pvLabelIdToVertId[mPVLabelMap->GetNthLabelObject(j)->GetLabel()]));
    std::cout << "cvpvAxes.size() = " << cvpvAxes.size() << std::endl;

    std::map<int, CScalarImagePointerType> cvMasks, pvMasks;
    for(unsigned int i=0; i<mCVLabelMap->GetNumberOfLabelObjects(); i++) {
        CScalarImagePointerType cvMask;
        this->LabelToBinImage(mCVCatchLabelImage, mCVLabelMap->GetNthLabelObject(i)->GetLabel(), cvMask);
        cvMasks.insert(std::pair<int, CScalarImagePointerType>(mCVLabelMap->GetNthLabelObject(i)->GetLabel(), cvMask));
    }
    for(unsigned int i=0; i<mPVLabelMap->GetNumberOfLabelObjects(); i++) {
        CScalarImagePointerType pvMask;
        this->LabelToBinImage(mPVCatchLabelImage, mPVLabelMap->GetNthLabelObject(i)->GetLabel(), pvMask);
        pvMasks.insert(std::pair<int, CScalarImagePointerType>(mPVLabelMap->GetNthLabelObject(i)->GetLabel(), pvMask));
    }

    int i=0;
//    std::cout << "cv-cv-axes tests" << std::endl;
    for(std::multimap<vtkIdType, vtkIdType>::iterator it=cvcvAxes.begin(); it!=cvcvAxes.end(); ++it) {
        if(IsAxisValid(0, it->first, it->second, cvMasks, pvMasks)) {
            vtkEdgeType e = mVeinAxesGraph->AddEdge(it->first, it->second);
            edgeState.insert(std::pair<vtkIdType, int>(e.Id, 0));
            edgeId.insert(std::pair<vtkIdType, int>(e.Id, i++));
        }
    }
//    std::cout << "pv-pv-axes tests" << std::endl;
    for(std::multimap<vtkIdType, vtkIdType>::iterator it=pvpvAxes.begin(); it!=pvpvAxes.end(); ++it) {
        if(IsAxisValid(1, it->first, it->second, cvMasks, pvMasks)) {
            vtkEdgeType e = mVeinAxesGraph->AddEdge(it->first, it->second);
            edgeState.insert(std::pair<vtkIdType, int>(e.Id, 1));
            edgeId.insert(std::pair<vtkIdType, int>(e.Id, i++));
        }
    }
//    std::cout << "cv-pv-axes tests" << std::endl;
    for(std::multimap<vtkIdType, vtkIdType>::iterator it=cvpvAxes.begin(); it!=cvpvAxes.end(); ++it) {
        if(IsAxisValid(2, it->first, it->second, cvMasks, pvMasks)) {
            vtkEdgeType e = mVeinAxesGraph->AddEdge(it->first, it->second);
            edgeState.insert(std::pair<vtkIdType, int>(e.Id, 2));
            edgeId.insert(std::pair<vtkIdType, int>(e.Id, i++));
        }
    }
    anno->AddCustomEdgeAnnotation(mVeinAxesGraph, "axis state", edgeState, -100);
    anno->AddCustomEdgeAnnotation(mVeinAxesGraph, "axis id", edgeId, -100);
}


void AnalyzeLobuleFilter::PrepareVeinAxesData(std::map<int, CScalar2DImagePointerType> &veinAxesImages, std::string path, std::string ext)
{
    vtkSmartPointer<vtkIntArray> labelId = vtkIntArray::SafeDownCast(mVeinAxesGraph->GetVertexData()->GetArray("label id"));
    vtkSmartPointer<vtkIntArray> veinState = vtkIntArray::SafeDownCast(mVeinAxesGraph->GetVertexData()->GetArray("vein state"));
    vtkSmartPointer<vtkIntArray> axisState = vtkIntArray::SafeDownCast(mVeinAxesGraph->GetEdgeData()->GetArray("axis state"));
    vtkSmartPointer<vtkPoints> points = mVeinAxesGraph->GetPoints();

    typename TileImageFilterType::Pointer       tiler;
    typename MaskNegImageFilterType2::Pointer   maskFilter1;
    typename MaskNegImageFilterType2::Pointer   maskFilter2;
    typename CScalarImageWriterType::Pointer    writer;

    CScalar2DIndexType index2d;
    CScalar2DSizeType size2d;
    CScalar2DRegionType region2d;
    CScalar2DSpacingType spacing2d;

    index2d[0] = mCVBin->GetLargestPossibleRegion().GetIndex()[0];
    index2d[1] = mCVBin->GetLargestPossibleRegion().GetIndex()[1];
    size2d[0] = mCVBin->GetLargestPossibleRegion().GetSize()[0];
    size2d[1] = mCVBin->GetLargestPossibleRegion().GetSize()[1];

    region2d.SetIndex(index2d);
    region2d.SetSize(size2d);

    spacing2d[0] = mSpacing[0];
    spacing2d[1] = mSpacing[1];

    itk::FixedArray< unsigned int, 3> layout;
    layout[0] = 1;
    layout[1] = 1;
    layout[2] = 0;

    tiler = TileImageFilterType::New();
    tiler->SetLayout(layout);

    vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();

    std::map<int, std::map<int, CScalar2DIndexType> >::iterator jtCV=mCVeinCentersPerSlice.begin();
    std::map<int, std::map<int, CScalar2DIndexType> >::iterator jtPV=mPVeinCentersPerSlice.begin();
    for(; jtCV!=mCVeinCentersPerSlice.end() && jtPV!=mPVeinCentersPerSlice.end(); ++jtCV, ++jtPV) {
        typename CScalar2DImageType::Pointer imageTile = CScalar2DImageType::New();
        imageTile->SetRegions(region2d);
        imageTile->Allocate();
        imageTile->FillBuffer(0);

        mVeinAxesGraph->GetEdges(it);
        while(it->HasNext()) {
            vtkEdgeType e = it->Next();
            int state = axisState->GetValue(e.Id);

//            std::cout << "axis " << e.Id << " of state " << state << std::endl;
            if(state==0) {                                                          //cv-cv axis
                CScalar2DIndexType p1, p2, p3;
                p1 = jtCV->second[labelId->GetValue(e.Source)];
                p2 = jtCV->second[labelId->GetValue(e.Target)];
                p3[0] = (p1[0]+p2[0]) / 2;
                p3[1] = (p1[1]+p2[1]) / 2;

                DrawEllipse(imageTile, p1, p3, mCenterRUSCRadius/spacing2d[0]);
                DrawEllipse(imageTile, p2, p3, mCenterRUSCRadius/spacing2d[0]);
            }
            else if(state==2) {                                                      //cv-pv axis
                CScalar2DIndexType p1, p2;

                if(veinState->GetValue(e.Source)==0) {
                    p1 = jtCV->second[labelId->GetValue(e.Source)];
                    p2 = jtPV->second[labelId->GetValue(e.Target)];
                }
                else {
                    p1 = jtPV->second[labelId->GetValue(e.Source)];
                    p2 = jtCV->second[labelId->GetValue(e.Target)];
                }
                DrawEllipse(imageTile, p1, p2, mCenterRUSCRadius/spacing2d[0]);
            }
        }
        veinAxesImages.insert(std::pair<int, CScalar2DImagePointerType>(jtCV->first, imageTile));
        tiler->SetInput(jtCV->first, imageTile);
    }
    tiler->Update();

    typename CScalarImageType::Pointer veinAxesImage = tiler->GetOutput();
    veinAxesImage->DisconnectPipeline();
    veinAxesImage->SetSpacing(mSpacing);

    maskFilter1 = MaskNegImageFilterType2::New();
    maskFilter1->SetInput(veinAxesImage);
    maskFilter1->SetMaskImage(mCVBin);
    maskFilter1->SetOutsideValue(0.0);

    maskFilter2 = MaskNegImageFilterType2::New();
    maskFilter2->SetInput(maskFilter1->GetOutput());
    maskFilter2->SetMaskImage(mPVBin);
    maskFilter2->SetOutsideValue(0.0);
    maskFilter2->Update();

    veinAxesImage = maskFilter2->GetOutput();
    veinAxesImage->DisconnectPipeline();

    std::stringstream s;
    s << path << "vein_axes_area_bin" << ext;

    writer = CScalarImageWriterType::New();
    writer->SetFileName(s.str());
    writer->SetInput(veinAxesImage);
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();
}


void AnalyzeLobuleFilter::Update()
{
    ParseParameterContext();

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    std::cout << "Analyze lobule pipeline: " << std::endl;
    std::cout << " path: " << mPath << std::endl;
    std::cout << " extension: " << mFilenameExtension << std::endl;


    mVeinSkeletonGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    GraphAnnotationHelper *anno = new GraphAnnotationHelper();
    std::map<vtkIdType, int> vertexIdToVeinState;

    PrepareVeinData(mDataSetFullFilenameCV, mPath, "vein_central_labelImage", mFilenameExtension, mCVBin, mCVLabelMap, mCVeinCentersPerSlice, vertexIdToVeinState, 0);
    PrepareVeinData(mDataSetFullFilenamePV, mPath, "vein_portal_labelImage", mFilenameExtension, mPVBin, mPVLabelMap, mPVeinCentersPerSlice, vertexIdToVeinState, 1);

    anno->EnableVertexIDAnnotation();
    anno->AddPredefinedAnnotations(mVeinSkeletonGraph);
    anno->AddCustomVertexAnnotation(mVeinSkeletonGraph, "vein state", vertexIdToVeinState, -1);

    ComputeVeinCatchmentBasins(mPath, "vein_central_", mFilenameExtension, mCVBin, mCVLabelMap, mCVCatchLabelImage, mCVDistmap);
    ComputeVeinCatchmentBasins(mPath, "vein_portal_", mFilenameExtension, mPVBin, mPVLabelMap, mPVCatchLabelImage, mPVDistmap);

    EstablishVeinAxes();

    std::map<int, CScalar2DImagePointerType> veinAxesImages;
    PrepareVeinAxesData(veinAxesImages, mPath, mFilenameExtension);

    WriteGraph(mVeinSkeletonGraph, mPath+"veinSkeletonGraph.txt");
    WriteGraph(mVeinAxesGraph, mPath+"veinAxisGraph.txt");


    WriteLogFile(timeStamp);

    mCVBin->ReleaseData();
    mCVBin = NULL;

    mCVLabelMap->ReleaseData();
    mCVLabelMap = NULL;

    mCVCatchLabelImage->ReleaseData();
    mCVCatchLabelImage = NULL;

    mPVBin->ReleaseData();
    mPVBin = NULL;

    mPVLabelMap->ReleaseData();
    mPVLabelMap = NULL;

    mPVCatchLabelImage->ReleaseData();
    mPVCatchLabelImage = NULL;

//    mRelativeLobularPositionMap->ReleaseData();
//    mRelativeLobularPositionMap = NULL;
}

