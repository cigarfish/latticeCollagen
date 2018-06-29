///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  EstimateLobuleShape.tpp                                              //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-08-05                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "EstimateLobuleShape.h"

#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIO.h>
#endif

#include <vtkFloatArray.h>
#include <vtkGraphWriter.h>
#include <vtkIntArray.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkUndirectedGraph.h>
#include <vtkSmartPointer.h>

#include <QDateTime>

#include <iostream>
#include <fstream>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



template< unsigned int VImageDimension > EstimateLobuleShape< VImageDimension >::EstimateLobuleShape()
{
    mSpacing[0] = 7.275;        //From Ole/Michael Schwier data
    mSpacing[1] = 7.275;
    if(ImageDimension == 3) mSpacing[2] = 7.275;

    mWatershedFloodLevel = 0.15;

    mOverlayOpacity = 0.5;

    mDataSetID = "myDataset";

    mFilenameSave = "lobule";
    mSaveEverything = true;

    mSaveSuffixesForFinals[0] = "_step3_lobLabels";
    mSaveSuffixesForFinals[1] = "_step3_lobOverlay";
    mSaveSuffixesForFinals[2] = "_step4_postproLobLabels";
    mSaveSuffixesForFinals[3] = "_step4_postproLobOverlay";
}


template< unsigned int VImageDimension > EstimateLobuleShape< VImageDimension >::~EstimateLobuleShape()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void EstimateLobuleShape< VImageDimension >::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    this->m_paramContext->dump(parameter);

    file.open((mPath + mLogFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-lobule-shape-approximation-------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-lobule-shape-approximation---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


template< unsigned int VImageDimension > void EstimateLobuleShape< VImageDimension >::WriteDataSetSummary()
{
    ImageAnalysisSummaryFileIO::AddEntry(LobuleShapeBin, mPath, mPath + mFilenameSave + mSaveSuffixesForFinals[0] + mFilenameExtension);
    ImageAnalysisSummaryFileIO::AddEntry(LobuleShapeOverlay, mPath, mPath + mFilenameSave + mSaveSuffixesForFinals[1] + mFilenameExtension);
}


template< unsigned int VImageDimension > void EstimateLobuleShape< VImageDimension >::ParseParameterContext()
{
    if(this->m_paramContext->findContext("Approximate Lobule Shape",0)==NULL) {
        std::cout << "Error: EstimateLobuleShape: Invalid parameter context" << std::endl;
        return;
    }

    mFilenameRawData = *(std::string*)(this->m_paramContext->findParameter("Raw data", 0)->dataPointer());
    mFilenameCVBin = *(std::string*)(this->m_paramContext->findParameter("Central vein segmentation", 0)->dataPointer());
    mFilenamePVBin = *(std::string*)(this->m_paramContext->findParameter("Portal vein segmentation", 0)->dataPointer());
    mFilenameTissueBin = *(std::string*)(this->m_paramContext->findParameter("Tissue segmentation", 0)->dataPointer());

    mCentralVeinThreshold = *(int*)(this->m_paramContext->findParameter("Central vein intensity", 0)->dataPointer());
    mPortalVeinThreshold = *(int*)(this->m_paramContext->findParameter("Portal vein intensity", 0)->dataPointer());
    mTissueThreshold = *(int*)(this->m_paramContext->findParameter("Tissue intensity", 0)->dataPointer());
    mVoidThreshold = *(int*)(this->m_paramContext->findParameter("Void intensity", 0)->dataPointer());

    mInfoFullFilename.setFile(QString::fromStdString( mFilenameCVBin) );
    bool cvFileExists = mInfoFullFilename.exists();
    mInfoFullFilename.setFile(QString::fromStdString( mFilenamePVBin) );
    bool pvFileExists = mInfoFullFilename.exists();

    if(!cvFileExists || !pvFileExists)
        throw std::string("Please specify both vein segmentation files.");

    mInfoFullFilename.setFile(QString::fromStdString( mFilenameTissueBin) );
    mTissueFileExists = mInfoFullFilename.exists();

    mPath = (mInfoFullFilename.path() + QString("/")).toStdString();
    mFilenameExtension = (QString(".") + mInfoFullFilename.suffix()).toStdString();

    mSpacing[0] = *(double*)(this->m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    mSpacing[1] = *(double*)(this->m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    if(ImageDimension == 3) mSpacing[2] = *(double*)(this->m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    mPortalVeinWeight = *(double*)(this->m_paramContext->findParameter("Portal vein weight [0,1]", 0)->dataPointer());
    mCentralVeinWeight = 1. - mPortalVeinWeight;

    mMinimalLobuleDiameter = *(double*)(this->m_paramContext->findParameter("Minimal lobule diameter", 0)->dataPointer());
    if(ImageDimension == 2)
        mMinimalLobuleVolume = 1./4. * itk::Math::pi * pow(mMinimalLobuleDiameter, 2);
    else if(ImageDimension == 3)
        mMinimalLobuleVolume = 1./6. * itk::Math::pi * pow(mMinimalLobuleDiameter, 3);

    mMaximalLobuleDiameter = *(double*)(this->m_paramContext->findParameter("Maximal lobule diameter", 0)->dataPointer());
    if(ImageDimension == 2)
        mMaximalLobuleVolume = 1./4. * itk::Math::pi * pow(mMaximalLobuleDiameter, 2);
    else if(ImageDimension == 3)
        mMaximalLobuleVolume = 1./6. * itk::Math::pi * pow(mMaximalLobuleDiameter, 3);

    mWatershedFloodLevel = *(double*)(this->m_paramContext->findParameter("Alpha", 0)->dataPointer());

    std::string analysisMode = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Analysis mode", 0)->dataPointer()) )->currentString();
    if( analysisMode.compare("Without Analysis")==0 )
        mWithAnalysis = false;
    else if( analysisMode.compare("With Analysis")==0 )
        mWithAnalysis = true;

    mDataSetID = *(std::string*)(this->m_paramContext->findParameter("Specify data set name", 0)->dataPointer());

    std::string saveMode = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        mSaveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        mSaveEverything = 0;
    mLogFilenameSave = *(std::string*)(this->m_paramContext->findParameter("Log file name", 0)->dataPointer());
    mFilenameSave = *(std::string*)(this->m_paramContext->findParameter("Save prefix", 0)->dataPointer());
}


template< unsigned int VImageDimension > void EstimateLobuleShape< VImageDimension >::OrganizeLobulesAndVeins()
{
    typename CScalarImageReaderType::Pointer            readerCV, readerPV;
    typename ThresholdCCFilterType::Pointer             thresholdCV, thresholdPV;
    typename ImageToShapeLabelMapFilterType::Pointer    CVShaLabMapFilter, PVShaLabMapFilter;
    typename LabelMapToLabelImageFilterType::Pointer    labelMapToImageFilter;

    readerCV = CScalarImageReaderType::New();
    readerCV->SetFileName(mFilenameCVBin);
    readerCV->ReleaseDataFlagOn();
    readerCV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerCV->SetImageIO( itk::TIFFImageIO::New() );
#endif

    CVShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    CVShaLabMapFilter->ReleaseDataFlagOn();
    CVShaLabMapFilter->SetInput(readerCV->GetOutput());
    CVShaLabMapFilter->SetInputForegroundValue(mCentralVeinThreshold);
    CVShaLabMapFilter->SetFullyConnected(true);
    CVShaLabMapFilter->Update();

    itk::SmartPointer<LabelMapType> CVLabelMap = CVShaLabMapFilter->GetOutput();
    CVLabelMap->DisconnectPipeline();
    CVLabelMap->SetSpacing(mSpacing);

    readerPV = CScalarImageReaderType::New();
    readerPV->SetFileName(mFilenamePVBin);
    readerPV->ReleaseDataFlagOn();
    readerPV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerPV->SetImageIO( itk::TIFFImageIO::New() );
#endif

    PVShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    PVShaLabMapFilter->ReleaseDataFlagOn();
    PVShaLabMapFilter->SetInput(readerPV->GetOutput());
    PVShaLabMapFilter->SetInputForegroundValue(mPortalVeinThreshold);
    PVShaLabMapFilter->SetFullyConnected(true);
    PVShaLabMapFilter->Update();

    itk::SmartPointer<LabelMapType> PVLabelMap = PVShaLabMapFilter->GetOutput();
    PVLabelMap->DisconnectPipeline();
    PVLabelMap->SetSpacing(mSpacing);

    itk::SmartPointer<LabelMapType> objectLabelMap = LabelMapType::New();
    objectLabelMap->SetRegions(mpLobuleLabelMap->GetLargestPossibleRegion());
    objectLabelMap->Allocate();
    objectLabelMap->SetBackgroundValue(0);
    objectLabelMap->SetSpacing(mSpacing);

    std::map<unsigned long,int> pedIdToTissueState;
    std::map<vtkIdType,int> vertIdToTissueState;

    int nextLabel = 1;
    for(unsigned int i=0; i<mpLobuleLabelMap->GetNumberOfLabelObjects(); i++) {
        typename LabelObjectType::Pointer labelObject = mpLobuleLabelMap->GetNthLabelObject(i);
        vtkIdType v = labelObject->GetLabel();
        labelObject->SetLabel(nextLabel);
        objectLabelMap->AddLabelObject(labelObject);
        pedIdToTissueState.insert(std::pair<vtkIdType,int>(nextLabel, 0));
        nextLabel++;
    }

    for(unsigned int i=0; i<CVLabelMap->GetNumberOfLabelObjects(); i++) {
        typename LabelObjectType::Pointer labelObject = CVLabelMap->GetNthLabelObject(i);
        vtkIdType v = labelObject->GetLabel();
        labelObject->SetLabel(nextLabel);
        objectLabelMap->AddLabelObject(labelObject);
        pedIdToTissueState.insert(std::pair<vtkIdType,int>(nextLabel, 1));
        nextLabel++;
    }

    for(unsigned int i=0; i<PVLabelMap->GetNumberOfLabelObjects(); i++) {
        typename LabelObjectType::Pointer labelObject = PVLabelMap->GetNthLabelObject(i);
        vtkIdType v = labelObject->GetLabel();
        labelObject->SetLabel(nextLabel);
        objectLabelMap->AddLabelObject(labelObject);
        pedIdToTissueState.insert(std::pair<vtkIdType,int>(nextLabel, 2));
        nextLabel++;
    }

    labelMapToImageFilter = LabelMapToLabelImageFilterType::New();
    labelMapToImageFilter->ReleaseDataFlagOn();
    labelMapToImageFilter->SetInput(objectLabelMap);
    labelMapToImageFilter->Update();

    itk::SmartPointer<LScalarImageType> labelImage = labelMapToImageFilter->GetOutput();
    labelImage->DisconnectPipeline();

    mLobuleGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    mLobuleGraph->DeepCopy( LabelImageToGraphFilter<VImageDimension>::LabelImageToGraph(labelImage, true) );

    std::map<vtkIdType,float> vertIdToVolume;
    std::map<vtkIdType,float> vertIdToElong;
    std::map<vtkIdType,int> vertIdToNumCV;
    std::map<vtkIdType,int> vertIdToNumPV;

    for(std::map<unsigned long,int>::iterator it=pedIdToTissueState.begin(); it!=pedIdToTissueState.end(); ++it) {
        unsigned long pedId = it->first;
        vtkIdType vertId = mLobuleGraph->FindVertex(pedId);

        vertIdToTissueState.insert(std::pair<vtkIdType,int>(vertId, it->second));
    }

    for(std::map<unsigned long,int>::iterator it=pedIdToTissueState.begin(); it!=pedIdToTissueState.end(); ++it) {
        unsigned long pedId = it->first;
        vtkIdType vertId = mLobuleGraph->FindVertex(pedId);

        if(it->second==0) {
            vertIdToVolume.insert(std::pair<vtkIdType,float>(vertId, mpLobuleLabelMap->GetLabelObject(pedId)->GetPhysicalSize()));
            vertIdToElong.insert(std::pair<vtkIdType,float>(vertId, mpLobuleLabelMap->GetLabelObject(pedId)->GetElongation()));

            int numCVs = 0;
            int numPVs = 0;
            for(unsigned int i=0; i<mLobuleGraph->GetDegree(vertId); i++) {
                vtkOutEdgeType edge = mLobuleGraph->GetOutEdge(vertId, i);

                if(vertIdToTissueState[edge.Target]==1)
                    numCVs++;
                if(vertIdToTissueState[edge.Target]==2)
                    numPVs++;
            }
            vertIdToNumCV.insert(std::pair<vtkIdType,int>(vertId, numCVs));
            vertIdToNumPV.insert(std::pair<vtkIdType,int>(vertId, numPVs));
        }
    }

    GraphAnnotationHelper* anno = new GraphAnnotationHelper();
    anno->AddCustomVertexAnnotation(mLobuleGraph, "tissue state", vertIdToTissueState, -1);
    anno->AddCustomVertexAnnotation(mLobuleGraph, "lobule volume", vertIdToVolume, -1);
    anno->AddCustomVertexAnnotation(mLobuleGraph, "lobule elongation", vertIdToElong, -1);
    anno->AddCustomVertexAnnotation(mLobuleGraph, "lobule's CVs", vertIdToNumCV, -1);
    anno->AddCustomVertexAnnotation(mLobuleGraph, "lobule's PVs", vertIdToNumPV, -1);
}


template< unsigned int VImageDimension > void EstimateLobuleShape< VImageDimension >::ClassifyLobules()
{
    vtkUnsignedLongArray* pedigreeIdArray = vtkUnsignedLongArray::SafeDownCast( mLobuleGraph->GetVertexData()->GetPedigreeIds() );
    vtkIntArray* tissueStateArray = vtkIntArray::SafeDownCast( mLobuleGraph->GetVertexData()->GetArray("tissue state") );
    vtkFloatArray* volumeArray = vtkFloatArray::SafeDownCast( mLobuleGraph->GetVertexData()->GetArray("lobule volume") );
    vtkFloatArray* elongationArray = vtkFloatArray::SafeDownCast( mLobuleGraph->GetVertexData()->GetArray("lobule elongation") );
    vtkIntArray* cvNumArray = vtkIntArray::SafeDownCast( mLobuleGraph->GetVertexData()->GetArray("lobule's CVs") );
    vtkIntArray* pvNumArray = vtkIntArray::SafeDownCast( mLobuleGraph->GetVertexData()->GetArray("lobule's PVs") );

    std::vector<unsigned long> labelsToRemove;
    for(int i=0; i<mLobuleGraph->GetNumberOfVertices(); i++) {
        if(tissueStateArray->GetValue(i)==0) {
            if(volumeArray->GetValue(i) < mMinimalLobuleVolume || volumeArray->GetValue(i) > mMaximalLobuleVolume)
                labelsToRemove.push_back(pedigreeIdArray->GetValue(i));
            else if(elongationArray->GetValue(i) > 1.75)
                labelsToRemove.push_back(pedigreeIdArray->GetValue(i));
            else if(cvNumArray->GetValue(i) < 1 || cvNumArray->GetValue(i) > 3)
                labelsToRemove.push_back(pedigreeIdArray->GetValue(i));
            else if(pvNumArray->GetValue(i) < 2 || cvNumArray->GetValue(i) > 5)
                labelsToRemove.push_back(pedigreeIdArray->GetValue(i));
        }
    }

    for(unsigned int i=0; i<labelsToRemove.size(); ++i) {
        mpLobuleLabelMap->RemoveLabel(labelsToRemove[i]);
        mLobuleGraph->RemoveVertex(mLobuleGraph->FindVertex(labelsToRemove[i]));
    }
}


//TODO: own class? at least option for analysis only; filename handling!!
template< unsigned int VImageDimension > void EstimateLobuleShape< VImageDimension >::AnalyzeLobules()
{
    for(unsigned int i=0; i<mpLobuleLabelMap->GetNumberOfLabelObjects(); ++i) {
        LobuleAnalysisContainer l;
        l.mLabel = mpLobuleLabelMap->GetNthLabelObject(i)->GetLabel();
        l.mVolume = mpLobuleLabelMap->GetNthLabelObject(i)->GetPhysicalSize();
        l.mRoundness = mpLobuleLabelMap->GetNthLabelObject(i)->GetRoundness();
        l.mElongation = mpLobuleLabelMap->GetNthLabelObject(i)->GetElongation();
        l.mNumPxlOnBorder = mpLobuleLabelMap->GetNthLabelObject(i)->GetNumberOfPixelsOnBorder();
        l.mPos[0] = mpLobuleLabelMap->GetNthLabelObject(i)->GetCentroid()[0];
        l.mPos[1] = mpLobuleLabelMap->GetNthLabelObject(i)->GetCentroid()[1];
        l.mPos[2] = mpLobuleLabelMap->GetNthLabelObject(i)->GetCentroid()[2];

        if(l.mLabel != 0)                                                                                //Label #0 is background label
            mLobules.insert(std::pair<unsigned long int, LobuleAnalysisContainer>(l.mLabel, l) );
    }
}


template< unsigned int VImageDimension > void EstimateLobuleShape< VImageDimension >::WriteAnalysisFile()
{
    std::fstream file1, file2, tempfile;
    bool hasHeader = false;

    file1.open((mPath + "../" + "Lobule_Output_file.txt").c_str(), fstream::in);

    std::string line;
    getline(file1, line);

    if(line.find("dataSet")!=std::string::npos) {
        hasHeader = true;

        tempfile.open((mPath + "../" + "Lobule_tempfile.txt").c_str(), fstream::out);
        tempfile << line;

        while(!file1.eof()) {
            getline(file1, line);
            if(line.find(mDataSetID)!=0) {
                tempfile << std::endl;
                tempfile << line;
            }
        }
        file1.close();
        tempfile.close();

        std::remove((mPath + "../" + "Lobule_Output_file.txt").c_str());
        std::rename((mPath + "../" + "Lobule_tempfile.txt").c_str(), (mPath + "../" + "Lobule_Output_file.txt").c_str());
    }
    else
        file1.close();

    if(!hasHeader)
        file2.open((mPath + "../" + "Lobule_Output_file.txt").c_str(), fstream::out);
    else
        file2.open((mPath + "../" + "Lobule_Output_file.txt").c_str(), fstream::out | fstream::app);

    file2.flags(fstream::left | fstream::scientific);
    if(!hasHeader) {
        file2.width(40);
        file2 << "dataSet";
        file2.width(15);
        file2 << "label";
        file2.width(15);
        file2 << "volume";
        file2.width(15);
        file2 << "roundness";
        file2.width(15);
        file2 << "elongation";
        file2.width(15);
        file2 << "numPxlOnBorder";
        file2.width(15);
        file2 << "xCentroid";
        file2.width(15);
        file2 << "yCentroid";
        file2.width(15);
        file2 << "zCentroid" << std::endl;
    }

    for(std::map<unsigned long int, LobuleAnalysisContainer>::iterator lobuleIt = mLobules.begin(); lobuleIt != mLobules.end(); ++lobuleIt) {
        file2.width(40);
        file2 << mDataSetID;
        file2.width(15);
        file2 << lobuleIt->second.mLabel;
        file2.width(15);
        file2 << lobuleIt->second.mVolume;
        file2.width(15);
        file2 << lobuleIt->second.mRoundness;
        file2.width(15);
        file2 << lobuleIt->second.mElongation;
        file2.width(15);
        file2 << lobuleIt->second.mNumPxlOnBorder;
        file2.width(15);
        file2 << lobuleIt->second.mPos[0];
        file2.width(15);
        file2 << lobuleIt->second.mPos[1];
        file2.width(15);
        file2 << lobuleIt->second.mPos[2] << std::endl;
    }

    file2.close();
}


template< unsigned int VImageDimension > void EstimateLobuleShape< VImageDimension >::Update()
{
    ParseParameterContext();

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    std::cout << "Approximate lobule shape: " << std::endl;
    std::cout << " path: " << mPath << std::endl;
    std::cout << " extension: " << mFilenameExtension << std::endl;

    std::cout << " mFilenameRawData: " << mFilenameRawData << std::endl;
    std::cout << " mFilenameCVBin: " << mFilenameCVBin << std::endl;
    std::cout << " mFilenamePVBin: " << mFilenamePVBin << std::endl;
    std::cout << " mFilenameTissueBin: " << mFilenameTissueBin << std::endl;

    typename CScalarImageReaderType::Pointer                    readerRawData, readerCV, readerPV, readerTissue;
    typename ThresholdCCFilterType::Pointer                     thresholdCV, thresholdPV, thresholdTissue, thresholdVoid;
    typename SignedMaurerDistanceMapImageFilterType::Pointer    distanceMapCV, distanceMapPV;
    typename IntensityWindowingImageFilterType::Pointer         distCVCapped, distPVCapped;
    typename RescaleFCImageFilterType::Pointer                  rescaler1, rescaler2, rescaler3, rescaler4;
    typename MultiplyImageFilterType::Pointer                   multCV, multPV;
    typename AddFImageFilterType::Pointer                       addDist;
    typename ImageCalculatorFilterType::Pointer                 addDistCalc, divideCalc;
    typename IntensityWindowingImageFilterType::Pointer         addDistCapped, divideCapped;
    typename DivideImageFilterType::Pointer                     divide;
    typename MorphoWatershedImageFilterType::Pointer            morphWatershed;
    typename MaskImageFilterType::Pointer                       rmMaskImageFilter;
    typename LabelMapToLabelImageFilterType::Pointer            lobuleLabelMapToImageFilter;
    typename ImageToShapeLabelMapFilterType::Pointer            lobuleBinImageToShapeLabelMap;
    typename ThresholdLCFilterType::Pointer                     lobuleImageToBin, thresholdImageFilter;
    typename LabelOverlayImageFilterType::Pointer               overlayImage;
    typename CScalarImageWriterType::Pointer                    writer1, writer2, writer3, writer4, writer5, writer6;
    typename RGBImageWriterType::Pointer                        writerRGB;

    //----------READER--------------------------------------------------------------------------------------------------------
    itk::SmartPointer<CScalarImageType> rawImage, CVImage, PVImage, tissueImage;

    readerRawData = CScalarImageReaderType::New();
    readerRawData->SetFileName(mFilenameRawData);
    readerRawData->ReleaseDataFlagOn();
    readerRawData->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerRawData->SetImageIO( itk::TIFFImageIO::New() );
#endif
    readerRawData->Update();

    rawImage = readerRawData->GetOutput();
    rawImage->DisconnectPipeline();
    rawImage->SetSpacing(mSpacing);

    readerCV = CScalarImageReaderType::New();
    readerCV->SetFileName(mFilenameCVBin);
    readerCV->ReleaseDataFlagOn();
    readerCV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerCV->SetImageIO( itk::TIFFImageIO::New() );
#endif

    thresholdCV = ThresholdCCFilterType::New();
    thresholdCV->ReleaseDataFlagOn();
    thresholdCV->SetOutsideValue(0);
    thresholdCV->SetInsideValue(255);
    thresholdCV->SetLowerThreshold(mCentralVeinThreshold);
    thresholdCV->SetUpperThreshold(mCentralVeinThreshold);
    thresholdCV->SetInput(readerCV->GetOutput());
    thresholdCV->Update();

    CVImage = thresholdCV->GetOutput();
    CVImage->DisconnectPipeline();
    CVImage->SetSpacing(mSpacing);

    readerPV = CScalarImageReaderType::New();
    readerPV->SetFileName(mFilenamePVBin);
    readerPV->ReleaseDataFlagOn();
    readerPV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerPV->SetImageIO( itk::TIFFImageIO::New() );
#endif

    thresholdPV = ThresholdCCFilterType::New();
    thresholdPV->ReleaseDataFlagOn();
    thresholdPV->SetOutsideValue(0);
    thresholdPV->SetInsideValue(255);
    thresholdPV->SetLowerThreshold(mPortalVeinThreshold);
    thresholdPV->SetUpperThreshold(mPortalVeinThreshold);
    thresholdPV->SetInput(readerPV->GetOutput());
    thresholdPV->Update();

    PVImage = thresholdPV->GetOutput();
    PVImage->DisconnectPipeline();
    PVImage->SetSpacing(mSpacing);

    if(mTissueFileExists) {
        readerTissue = CScalarImageReaderType::New();
        readerTissue->SetFileName(mFilenameTissueBin);
        readerTissue->ReleaseDataFlagOn();
        readerTissue->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerTissue->SetImageIO( itk::TIFFImageIO::New() );
#endif
    readerTissue->Update();

    thresholdTissue = ThresholdCCFilterType::New();
    thresholdTissue->ReleaseDataFlagOn();
    thresholdTissue->SetOutsideValue(0);
    thresholdTissue->SetInsideValue(255);
    thresholdTissue->SetLowerThreshold(mTissueThreshold);
    thresholdTissue->SetUpperThreshold(mTissueThreshold);
    thresholdTissue->SetInput(readerTissue->GetOutput());
    thresholdTissue->Update();

    tissueImage = thresholdTissue->GetOutput();
    tissueImage->DisconnectPipeline();
    tissueImage->SetSpacing(mSpacing);

    //not used atm
    thresholdVoid = ThresholdCCFilterType::New();
    thresholdVoid->ReleaseDataFlagOn();
    thresholdVoid->SetOutsideValue(0);
    thresholdVoid->SetInsideValue(255);
    thresholdVoid->SetLowerThreshold(mVoidThreshold);
    thresholdVoid->SetUpperThreshold(mVoidThreshold);
    thresholdVoid->SetInput(readerTissue->GetOutput());
    }
    else {
        tissueImage = CScalarImageType::New();
        this->CreateImage(tissueImage, CVImage->GetLargestPossibleRegion(), mTissueThreshold, mSpacing);

        thresholdTissue = ThresholdCCFilterType::New();
        thresholdTissue->ReleaseDataFlagOn();
        thresholdTissue->SetOutsideValue(0);
        thresholdTissue->SetInsideValue(255);
        thresholdTissue->SetLowerThreshold(mTissueThreshold);
        thresholdTissue->SetUpperThreshold(mTissueThreshold);
        thresholdTissue->SetInput(tissueImage);
        thresholdTissue->Update();

        tissueImage = thresholdTissue->GetOutput();
        tissueImage->DisconnectPipeline();
        tissueImage->SetSpacing(mSpacing);

        //not used atm
        thresholdVoid = ThresholdCCFilterType::New();
        thresholdVoid->ReleaseDataFlagOn();
        thresholdVoid->SetOutsideValue(0);
        thresholdVoid->SetInsideValue(255);
        thresholdVoid->SetLowerThreshold(mVoidThreshold);
        thresholdVoid->SetUpperThreshold(mVoidThreshold);
        thresholdVoid->SetInput(tissueImage);
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------DISTANCE-MAP---------------------------------------------------------------------------------------------------
    distanceMapCV = SignedMaurerDistanceMapImageFilterType::New();
    distanceMapCV->ReleaseDataFlagOn();
    distanceMapCV->UseImageSpacingOn();
    distanceMapCV->SquaredDistanceOff();
    distanceMapCV->SetInput(CVImage);

    distCVCapped = IntensityWindowingImageFilterType::New();        //sets all values >mMaximalLobuleDiameter to mMaximalLobuleDiameter
    distCVCapped->ReleaseDataFlagOn();
    distCVCapped->SetWindowMinimum(0.);
    distCVCapped->SetWindowMaximum(mMaximalLobuleDiameter/2.);
    distCVCapped->SetOutputMinimum(0.);
    distCVCapped->SetOutputMaximum(mMaximalLobuleDiameter/2.);
    distCVCapped->SetInput(distanceMapCV->GetOutput());
    distCVCapped->Update();

    distanceMapPV = SignedMaurerDistanceMapImageFilterType::New();
    distanceMapPV->ReleaseDataFlagOn();
    distanceMapPV->UseImageSpacingOn();
    distanceMapPV->SquaredDistanceOff();
    distanceMapPV->SetInput(PVImage);

    distPVCapped = IntensityWindowingImageFilterType::New();        //sets all values >mMaximalLobuleDiameter to mMaximalLobuleDiameter
    distPVCapped->ReleaseDataFlagOn();
    distPVCapped->SetWindowMinimum(0.);
    distPVCapped->SetWindowMaximum(mMaximalLobuleDiameter/2.);
    distPVCapped->SetOutputMinimum(0.);
    distPVCapped->SetOutputMaximum(mMaximalLobuleDiameter/2.);
    distPVCapped->SetInput(distanceMapPV->GetOutput());
    distPVCapped->Update();

    itk::SmartPointer<FScalarImageType> distMapCVImage = distCVCapped->GetOutput();
    distMapCVImage->DisconnectPipeline();

    itk::SmartPointer<FScalarImageType> distMapPVImage = distPVCapped->GetOutput();
    distMapPVImage->DisconnectPipeline();

    CVImage->ReleaseData();
    CVImage = NULL;
    PVImage->ReleaseData();
    PVImage = NULL;

    std::cout << "Write distMap intermediate results " << mSaveEverything << std::endl;
    if(mSaveEverything) {
        rescaler1 = RescaleFCImageFilterType::New();
        rescaler1->ReleaseDataFlagOn();
        rescaler1->SetInput(distMapCVImage);

        rescaler2 = RescaleFCImageFilterType::New();
        rescaler2->ReleaseDataFlagOn();
        rescaler2->SetInput(distMapPVImage);

        writer2 = CScalarImageWriterType::New();
        writer2->SetFileName(mPath + mFilenameSave + "_step1_distMapCV.tif");
        writer2->ReleaseDataFlagOn();
        writer2->SetInput(rescaler1->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer2->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer2->Update();

        writer3 = CScalarImageWriterType::New();
        writer3->SetFileName(mPath + mFilenameSave + "_step1_distMapPV.tif");
        writer3->ReleaseDataFlagOn();
        writer3->SetInput(rescaler2->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer3->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer3->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------MULTIPLY-WITH-WEIGHTS------------------------------------------------------------------------------------------
    multCV = MultiplyImageFilterType::New();
    multCV->SetInput(distMapCVImage);
    multCV->SetConstant(mCentralVeinWeight);

    multPV = MultiplyImageFilterType::New();
    multPV->ReleaseDataFlagOn();
    multPV->SetInput(distMapPVImage);
    multPV->SetConstant(mPortalVeinWeight);
    //------------------------------------------------------------------------------------------------------------------------

    //----------CALCULATE-REFERENCE-IMAGE-------------------------------------------------------------------------------------
    addDist = AddFImageFilterType::New();
    addDist->ReleaseDataFlagOn();
    addDist->SetInput1(multCV->GetOutput());
    addDist->SetInput2(multPV->GetOutput());
    addDist->Update();

    addDistCalc = ImageCalculatorFilterType::New();
    addDistCalc->SetImage(addDist->GetOutput());
    addDistCalc->Compute();

    std::cout << "addDistCalc->GetMinimum() of divide = " << addDistCalc->GetMinimum() << std::endl;
    std::cout << "addDistCalc->GetMaximum() of divide = " << addDistCalc->GetMaximum() << std::endl;

    addDistCapped = IntensityWindowingImageFilterType::New();           //sets all values <0.0000001 to 0.0000001
    addDistCapped->SetWindowMinimum(0.0000001);
    addDistCapped->SetWindowMaximum(addDistCalc->GetMaximum());
    addDistCapped->SetOutputMinimum(0.0000001);
    addDistCapped->SetOutputMaximum(addDistCalc->GetMaximum());
    addDistCapped->SetInput(addDist->GetOutput());
    addDistCapped->Update();

    itk::SmartPointer<FScalarImageType> intermediateGradientImage = addDistCapped->GetOutput();
    intermediateGradientImage->DisconnectPipeline();

    divide = DivideImageFilterType::New();
    divide->ReleaseDataFlagOn();
    divide->SetInput1(multCV->GetOutput());
    divide->SetInput2(intermediateGradientImage);
    divide->Update();

    divideCalc = ImageCalculatorFilterType::New();
    divideCalc->SetImage(divide->GetOutput());
    divideCalc->Compute();

    std::cout << "divideCalc->GetMinimum() of divide = " << divideCalc->GetMinimum() << std::endl;
    std::cout << "divideCalc->GetMaximum() of divide = " << divideCalc->GetMaximum() << std::endl;

    divideCapped = IntensityWindowingImageFilterType::New();        //scales everything to 0..100 for comparable watershed flood level values
    divideCapped->ReleaseDataFlagOn();
    divideCapped->SetWindowMinimum(divideCalc->GetMinimum());
    divideCapped->SetWindowMaximum(divideCalc->GetMaximum());
    divideCapped->SetOutputMinimum(0.);
    divideCapped->SetOutputMaximum(100.);
    divideCapped->SetInput(divide->GetOutput());
    divideCapped->Update();

    itk::SmartPointer<FScalarImageType> gradientImage = divideCapped->GetOutput();
    gradientImage->DisconnectPipeline();

    gradientImage->SetSpacing(mSpacing);

    distMapCVImage->ReleaseData();
    distMapCVImage = NULL;
    distMapPVImage->ReleaseData();
    distMapPVImage = NULL;
    std::cout << "EstimateLobuleShape: released distMapCVImage & distMapPVImage" << std::endl;

    divideCalc->SetImage(gradientImage);
    divideCalc->Compute();

    std::cout << "divideCalc->GetMinimum() of divideCapped = " << divideCalc->GetMinimum() << std::endl;
    std::cout << "divideCalc->GetMaximum() of divideCapped = " << divideCalc->GetMaximum() << std::endl;

    std::cout << "Write formula intermediate results " << mSaveEverything << std::endl;
    if(mSaveEverything) {
        rescaler3 = RescaleFCImageFilterType::New();
        rescaler3->ReleaseDataFlagOn();
        rescaler3->SetInput(intermediateGradientImage);

        rescaler4 = RescaleFCImageFilterType::New();
        rescaler4->ReleaseDataFlagOn();
        rescaler4->SetInput(gradientImage);

        writer4 = CScalarImageWriterType::New();
        writer4->SetFileName(mPath + mFilenameSave + "_step2_add.tif");
        writer4->ReleaseDataFlagOn();
        writer4->SetInput(rescaler3->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer4->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer4->Update();

        writer5 = CScalarImageWriterType::New();
        writer5->SetFileName(mPath + mFilenameSave + "_step2_divide.tif");
        writer5->ReleaseDataFlagOn();
        writer5->SetInput(rescaler4->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer5->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer5->Update();
    }
    intermediateGradientImage->ReleaseData();
    intermediateGradientImage = NULL;
    std::cout << "EstimateLobuleShape: released intermediateGradientImage" << std::endl;
    //------------------------------------------------------------------------------------------------------------------------

    //----------WATERSHED-----------------------------------------------------------------------------------------------------
    morphWatershed = MorphoWatershedImageFilterType::New();
    morphWatershed->ReleaseDataFlagOn();
    morphWatershed->SetLevel(mWatershedFloodLevel);
    morphWatershed->FullyConnectedOff();
    morphWatershed->SetInput(gradientImage);
    morphWatershed->Update();

    itk::SmartPointer<LScalarImageType> lobuleImage = morphWatershed->GetOutput();
    lobuleImage->DisconnectPipeline();
    std::cout << "lobuleImage->GetSpacing() = " << lobuleImage->GetSpacing() << std::endl;
    std::cout << "rawImage->GetSpacing() = " << rawImage->GetSpacing() << std::endl;

    gradientImage->ReleaseData();
    gradientImage = NULL;
    std::cout << "EstimateLobuleShape: released gradientImage" << std::endl;
    //------------------------------------------------------------------------------------------------------------------------

    //----------MASK-AND-RELABEL-WATERSHED-LABEL-OBJECTS-----------------------------------------------------------------------
    lobuleImageToBin = ThresholdLCFilterType::New();
    lobuleImageToBin->ReleaseDataFlagOn();
    lobuleImageToBin->SetOutsideValue(0);
    lobuleImageToBin->SetInsideValue(itk::NumericTraits<CScalarPixelType>::max());
    lobuleImageToBin->SetLowerThreshold(1);
    lobuleImageToBin->SetUpperThreshold(itk::NumericTraits<LScalarPixelType>::max());
    lobuleImageToBin->SetInput(lobuleImage);

    rmMaskImageFilter = MaskImageFilterType::New();
    rmMaskImageFilter->ReleaseDataFlagOn();
    rmMaskImageFilter->SetInput(lobuleImageToBin->GetOutput());
    rmMaskImageFilter->SetMaskImage(tissueImage);
    rmMaskImageFilter->Update();

    lobuleImage->ReleaseData();
    lobuleImage = NULL;
    std::cout << "EstimateLobuleShape: released lobuleImage" << std::endl;

    itk::SmartPointer<CScalarImageType> maskedLobuleImage = rmMaskImageFilter->GetOutput();
    maskedLobuleImage->DisconnectPipeline();

    lobuleBinImageToShapeLabelMap = ImageToShapeLabelMapFilterType::New();
    lobuleBinImageToShapeLabelMap->ReleaseDataFlagOn();
    lobuleBinImageToShapeLabelMap->SetFullyConnected(false);
    lobuleBinImageToShapeLabelMap->SetInput(maskedLobuleImage);
    lobuleBinImageToShapeLabelMap->Update();

    mpLobuleLabelMap = lobuleBinImageToShapeLabelMap->GetOutput();
    mpLobuleLabelMap->DisconnectPipeline();

    lobuleLabelMapToImageFilter = LabelMapToLabelImageFilterType::New();
    lobuleLabelMapToImageFilter->ReleaseDataFlagOn();
    lobuleLabelMapToImageFilter->SetInput(mpLobuleLabelMap);
    lobuleLabelMapToImageFilter->Update();

    maskedLobuleImage->ReleaseData();
    maskedLobuleImage = NULL;
    std::cout << "EstimateLobuleShape: released maskedLobuleImage" << std::endl;

    lobuleImage = lobuleLabelMapToImageFilter->GetOutput();
    lobuleImage->DisconnectPipeline();
    //-------------------------------------------------------------------------------------------------------------------------

    overlayImage = LabelOverlayImageFilterType::New();
    overlayImage->ReleaseDataFlagOn();
    overlayImage->SetOpacity(mOverlayOpacity);
    overlayImage->SetInput(rawImage);
    overlayImage->SetLabelImage(lobuleImage);

    writerRGB = RGBImageWriterType::New();
    writerRGB->ReleaseDataFlagOn();
    writerRGB->SetFileName(mPath + mFilenameSave + mSaveSuffixesForFinals[1] + mFilenameExtension);
    writerRGB->SetInput(overlayImage->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writerRGB->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerRGB->Update();

    thresholdImageFilter = ThresholdLCFilterType::New();
    thresholdImageFilter->ReleaseDataFlagOn();
    thresholdImageFilter->SetOutsideValue(0);
    thresholdImageFilter->SetInsideValue(itk::NumericTraits<CScalarPixelType>::max());
    thresholdImageFilter->SetLowerThreshold(1);
    thresholdImageFilter->SetUpperThreshold(itk::NumericTraits<LScalarPixelType>::max());
    thresholdImageFilter->SetInput(lobuleImage);

    writer6 = CScalarImageWriterType::New();
    writer6->ReleaseDataFlagOn();
    writer6->SetFileName(mPath + mFilenameSave + mSaveSuffixesForFinals[0] + mFilenameExtension);
    writer6->SetInput(thresholdImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer6->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer6->Update();

    lobuleImage->ReleaseData();
    lobuleImage = NULL;
    std::cout << "EstimateLobuleShape: released lobuleImage" << std::endl;

    OrganizeLobulesAndVeins();

    vtkSmartPointer<vtkGraphWriter> graphWriter = vtkSmartPointer<vtkGraphWriter>::New();
    graphWriter->SetFileName((mPath + "lobeArchitectureGraph.txt").c_str());
    graphWriter->SetInput(mLobuleGraph);
    graphWriter->Update();

    ClassifyLobules();

    graphWriter->SetFileName((mPath + "lobeArchitectureGraphPostprocessed.txt").c_str());
    graphWriter->SetInput(mLobuleGraph);
    graphWriter->Update();

    lobuleLabelMapToImageFilter->SetInput(mpLobuleLabelMap);
    lobuleLabelMapToImageFilter->Update();

    lobuleImage = lobuleLabelMapToImageFilter->GetOutput();
    lobuleImage->DisconnectPipeline();

    overlayImage->SetOpacity(mOverlayOpacity);
    overlayImage->SetInput(rawImage);
    overlayImage->SetLabelImage(lobuleImage);

    writerRGB->SetInput(overlayImage->GetOutput());
    writerRGB->SetFileName(mPath + mFilenameSave + mSaveSuffixesForFinals[3] + mFilenameExtension);
    writerRGB->Update();

    rawImage->ReleaseData();
    rawImage = NULL;
    std::cout << "EstimateLobuleShape: released rawImage" << std::endl;

    thresholdImageFilter->SetInput(lobuleImage);

    writer6->SetFileName(mPath + mFilenameSave + mSaveSuffixesForFinals[2] + mFilenameExtension);
    writer6->SetInput(thresholdImageFilter->GetOutput());
    writer6->Update();

    if(mWithAnalysis) {
        AnalyzeLobules();
        WriteAnalysisFile();
    }

    WriteLogFile(timeStamp);
    WriteDataSetSummary();
}
