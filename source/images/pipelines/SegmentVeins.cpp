///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentVeins.cpp                                                     //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-03-16                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "SegmentVeins.h"

#if (ITK_VERSION_MAJOR < 4)
#define NUMBER_OF_PIXELS SIZE
#endif

#include <QFileInfo>
#include <QDateTime>
#include <QString>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



SegmentVeins::SegmentVeins()
{
    m_openingStrucElem.SetRadius(6);
    m_openingStrucElem.CreateStructuringElement();

    m_closingStrucElem.SetRadius(6);
    m_closingStrucElem.CreateStructuringElement();

//    m_maskDilateStrucElem.SetRadius(50);
//    m_maskDilateStrucElem.CreateStructuringElement();

    m_minimalVeinSize = 5000;

    m_overlayOpacity = 0.5;

    m_cvSaveBin = "vein_central_bin";
    m_cvSaveOverlay = "vein_central_overlay";
    m_pvSaveBin = "vein_portal_bin";
    m_pvSaveOverlay = "vein_portal_overlay";
//    m_cvMaskSaveBin = "vein_central_mask_bin";
//    m_cvMaskSaveOverlay = "vein_central_mask_overlay";
//    m_pvMaskSaveBin = "vein_portal_mask_bin";
//    m_pvMaskSaveOverlay = "vein_portal_mask_overlay";
}


SegmentVeins::~SegmentVeins()
{
    // TODO Auto-generated destructor stub
}

void SegmentVeins::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

    file.open((m_pathChannel1 + "log_veins" + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-segment-veins-------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-segment-veins---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


void SegmentVeins::WriteDataSetSummary()
{
    if(m_withSegmentation) {
        ImageAnalysisSummaryFileIO::AddEntry(CentralVeinSegmentationBin, m_pathChannel1, m_pathChannel1 + m_cvSaveBin + m_fileExtensionChannel1);
        ImageAnalysisSummaryFileIO::AddEntry(CentralVeinSegmentationOverlay, m_pathChannel1, m_pathChannel1 + m_cvSaveOverlay + m_fileExtensionChannel1);

        ImageAnalysisSummaryFileIO::AddEntry(PortalVeinSegmentationBin, m_pathChannel1, m_pathChannel1 + m_pvSaveBin + m_fileExtensionChannel1);
        ImageAnalysisSummaryFileIO::AddEntry(PortalVeinSegmentationOverlay, m_pathChannel1, m_pathChannel1 + m_pvSaveOverlay + m_fileExtensionChannel1);

//        ImageAnalysisSummaryFileIO::AddEntry(CentralVeinMaskSegmentationBin, m_pathChannel1, m_pathChannel1 + m_cvMaskSaveBin + m_fileExtensionChannel1);
//        ImageAnalysisSummaryFileIO::AddEntry(CentralVeinMaskSegmentationOverlay, m_pathChannel1, m_pathChannel1 + m_cvMaskSaveOverlay + m_fileExtensionChannel1);
//
//        ImageAnalysisSummaryFileIO::AddEntry(PortalVeinMaskSegmentationBin, m_pathChannel1, m_pathChannel1 + m_pvMaskSaveBin + m_fileExtensionChannel1);
//        ImageAnalysisSummaryFileIO::AddEntry(PortalVeinMaskSegmentationOverlay, m_pathChannel1, m_pathChannel1 + m_pvMaskSaveOverlay + m_fileExtensionChannel1);
    }
}


void SegmentVeins::ParseParameterContext()
{
    if(m_paramContext->findContext("Segment Veins",0)==NULL) {
        std::cout << "Error: SegmentVeins: Invalid parameter context: " << std::endl;
        return;
    }

    m_fullFilenameChannel0 = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("DAPI channel", 0)->dataPointer()) );
	m_infoFullFilenameChannel0.setFile(m_fullFilenameChannel0);
	m_hasCh0 = m_infoFullFilenameChannel0.exists();
	if(m_hasCh0) {
		m_pathChannel0 = (m_infoFullFilenameChannel0.path() + QString("/")).toStdString();
		m_filenameChannel0 = m_infoFullFilenameChannel0.baseName().toStdString();
		m_fileExtensionChannel0 = (QString(".") + m_infoFullFilenameChannel0.suffix()).toStdString();
	}

    m_fullFilenameChannel1 = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("DPPIV channel", 0)->dataPointer()) );
	m_infoFullFilenameChannel1.setFile(m_fullFilenameChannel1);
	if(!m_infoFullFilenameChannel1.exists())
		throw std::string("Please specify DPPIV channel");
	m_pathChannel1 = (m_infoFullFilenameChannel1.path() + QString("/")).toStdString();
    m_filenameChannel1 = m_infoFullFilenameChannel1.baseName().toStdString();
    m_fileExtensionChannel1 = (QString(".") + m_infoFullFilenameChannel1.suffix()).toStdString();

    m_fullFilenameChannel2 = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("DMs channel", 0)->dataPointer()) );
	m_infoFullFilenameChannel2.setFile(m_fullFilenameChannel2);
	if(!m_infoFullFilenameChannel2.exists())
		throw std::string("Please specify DMs channel");
	m_pathChannel2 = (m_infoFullFilenameChannel2.path() + QString("/")).toStdString();
    m_filenameChannel2 = m_infoFullFilenameChannel2.baseName().toStdString();
    m_fileExtensionChannel2 = (QString(".") + m_infoFullFilenameChannel2.suffix()).toStdString();

    m_thresholdDAPI = *(int*)(m_paramContext->findParameter("DAPI threshold", 0)->dataPointer());
    m_thresholdDPPIV = *(int*)(m_paramContext->findParameter("DPPIV threshold", 0)->dataPointer());
    m_thresholdDM = *(int*)(m_paramContext->findParameter("DMs threshold", 0)->dataPointer());

    int oRad = *(int*)(m_paramContext->findParameter("Opening kernel radius", 0)->dataPointer());
    int cRad = *(int*)(m_paramContext->findParameter("Closing kernel radius", 0)->dataPointer());

    m_openingStrucElem.SetRadius(oRad);
    m_openingStrucElem.CreateStructuringElement();

    m_closingStrucElem.SetRadius(cRad);
    m_closingStrucElem.CreateStructuringElement();

    m_withResampling = *(bool*)(m_paramContext->findParameter("With resampling", 0)->dataPointer());
    m_withSegmentation = *(bool*)(m_paramContext->findParameter("With segmentation", 0)->dataPointer());
    m_withLobulePrediction = *(bool*)(m_paramContext->findParameter("With lobule prediction", 0)->dataPointer());
}


void SegmentVeins::Update()
{
    ParseParameterContext();

    std::string timeStamp, segmentationCVResult, segmentationPVResult, segmentationMaskCVResult, segmentationMaskPVResult;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    std::cout << "Segment veins sources: " << std::endl;
    std::cout << " dir0: " << m_pathChannel0 << std::endl;
    std::cout << " file0: " << m_filenameChannel0 << std::endl;
    std::cout << " ext0: " << m_fileExtensionChannel0 << std::endl;

    std::cout << " dir1: " << m_pathChannel1 << std::endl;
    std::cout << " file1: " << m_filenameChannel1 << std::endl;
    std::cout << " ext1: " << m_fileExtensionChannel1 << std::endl;

    std::cout << " dir2: " << m_pathChannel2 << std::endl;
    std::cout << " file2: " << m_filenameChannel2 << std::endl;
    std::cout << " ext2: " << m_fileExtensionChannel2 << std::endl;


    ScalarVoReaderType::Pointer             readerDAPI;
    ScalarVoReaderType::Pointer             readerDPPIV;
    ScalarVoReaderType::Pointer             readerDM;
    ResampleImageFilterType::Pointer        resampleDAPI;
    ResampleImageFilterType::Pointer        resampleDPPIV;
    ResampleImageFilterType::Pointer        resampleDM;
    CScalarVoImageType::SizeType            inputSize;
    CScalarVoImageType::SpacingType         inputSpacing;
    Compose2DVectorImageFilterType::Pointer composeDPPIVDM;
    Compose3DVectorImageFilterType::Pointer composeDAPIDPPIVDM;
    Connected3DThresholdFilterType::Pointer connected3DThresholdCV;
    Connected3DThresholdFilterType::Pointer connected3DThresholdPV;
    Connected2DThresholdFilterType::Pointer connected2DThresholdCV;
    Connected2DThresholdFilterType::Pointer connected2DThresholdPV;
    OpeningImageFilterType::Pointer         openingCV;
    OpeningImageFilterType::Pointer         openingPV;
    ClosingImageFilterType::Pointer         closingCV;
    ClosingImageFilterType::Pointer         closingPV;
    ImageToShapeLabelMapFilterType::Pointer imageToShapeLabelMapCV;
    ImageToShapeLabelMapFilterType::Pointer imageToShapeLabelMapPV;
    ShapeOpeningLabelMapFilterType::Pointer shapeOpeningCV;
    ShapeOpeningLabelMapFilterType::Pointer shapeOpeningPV;
    LabelMapToLabelImageFilterType::Pointer labelMapToImageCV;
    LabelMapToLabelImageFilterType::Pointer labelMapToImagePV;
    ThresholdFilterType::Pointer            thresholdCV;
    ThresholdFilterType::Pointer            thresholdPV;
    ResampleImageFilterType::Pointer        resampleCV1;
    ResampleImageFilterType::Pointer        resamplePV1;
    ScalarVoWriterType::Pointer             writerBinCV;
    ScalarVoWriterType::Pointer             writerBinPV;
    LabelOverlayImageFilterType::Pointer    overlayImageCV;
    LabelOverlayImageFilterType::Pointer    overlayImagePV;
    RGBVoWriterType::Pointer                writerOverlayCV;
    RGBVoWriterType::Pointer                writerOverlayPV;
    DilateImageFilterType::Pointer          maskDilateCV;
    DilateImageFilterType::Pointer          maskDilatePV;
    ImageToShapeLabelMapFilterType::Pointer maskImageToShapeLabelMapCV;
    ImageToShapeLabelMapFilterType::Pointer maskImageToShapeLabelMapPV;
    LabelMapToLabelImageFilterType::Pointer maskLabelMapToImageCV;
    LabelMapToLabelImageFilterType::Pointer maskLabelMapToImagePV;
    LabelOverlayImageFilterType::Pointer    maskOverlayImageCV;
    LabelOverlayImageFilterType::Pointer    maskOverlayImagePV;
    RGBVoWriterType::Pointer                maskWriterOverlayCV;
    RGBVoWriterType::Pointer                maskWriterOverlayPV;

    ScalarVoWriterType::Pointer             writerTest1;

    ScalarVoReaderType::Pointer                     readerCV;
    ScalarVoReaderType::Pointer                     readerPV;
    ResampleImageFilterType::Pointer                resampleCV2;
    ResampleImageFilterType::Pointer                resamplePV2;
    DanielssonDistanceMapImageFilterType::Pointer   distanceMapCV;
    DanielssonDistanceMapImageFilterType::Pointer   distanceMapPV;
    AddImageFilterType::Pointer                     add;
    DivideImageFilterType::Pointer                  divide;
    MorphoWatershedImageFilterType::Pointer         morphWatershed;
    LabelImageToLabelMapFilterType::Pointer         labelImageToLabelMap;
    LabelMapToLabelImageFilterTypeW::Pointer        labelMapToImageW;
    LabelContourImageFilterType::Pointer            labelContour;
    RescaleImageFilterType::Pointer                 rescaler1, rescaler2, rescaler3;
    ResampleImageFilterType::Pointer                resample1, resample2, resample3, resample4;
    ScalarVoWriterType::Pointer                     writerW1, writerW2, writerW3, writerW4;
    LabelOverlayImageFilterType::Pointer            overlayImageW;
    RGBVoWriterType::Pointer                        writerWRGB1;

    std::cout << "jetzt gehts lohos"<< std::endl;

	if(m_withSegmentation && m_CVSeedPoints.size()==0 && m_PVSeedPoints.size()==0)
		throw std::string("Please specify seed points, using the button below the parameter table");

    if(m_withSegmentation)
    {
        //----------READER--------------------------------------------------------------------------------------------------------
        if(m_hasCh0) {
            if(GetNumberOfDimensions(m_pathChannel0 + m_filenameChannel0 + m_fileExtensionChannel0) != 3)
                throw std::string("Please specify a DAPI channel file with three dimensions.");

            readerDAPI = ScalarVoReaderType::New();
            readerDAPI->SetFileName(m_pathChannel0 + m_filenameChannel0 + m_fileExtensionChannel0);
            readerDAPI->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
            readerDAPI->SetImageIO( itk::TIFFImageIO::New() );
#endif
            readerDAPI->Update();
        }

        if(GetNumberOfDimensions(m_pathChannel1 + m_filenameChannel1 + m_fileExtensionChannel1) != 3)
            throw std::string("Please specify a DPPIV channel file with three dimensions.");

        readerDPPIV = ScalarVoReaderType::New();
        readerDPPIV->SetFileName(m_pathChannel1 + m_filenameChannel1 + m_fileExtensionChannel1);
        readerDPPIV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
        readerDPPIV->SetImageIO( itk::TIFFImageIO::New() );
#endif
        readerDPPIV->Update();

        if(GetNumberOfDimensions(m_pathChannel2 + m_filenameChannel2 + m_fileExtensionChannel2) != 3)
            throw std::string("Please specify a DMs channel file with three dimensions.");

        readerDM = ScalarVoReaderType::New();
        readerDM->SetFileName(m_pathChannel2 + m_filenameChannel2 + m_fileExtensionChannel2);
        readerDM->ReleaseDataBeforeUpdateFlagOn();
        readerDM->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
        readerDM->SetImageIO( itk::TIFFImageIO::New() );
#endif
        //-------------------------------------------------------------------------------------------------------------------------

        //----------RESAMPLE-IMAGES------------------------------------------------------------------------------------------------
        NearestNeighborExtrapolatorImageFunctionType::Pointer extrapolator = NearestNeighborExtrapolatorImageFunctionType::New();
        NearestNeighborInterpolatorImageFunctionType::Pointer interpolator = NearestNeighborInterpolatorImageFunctionType::New();

        if(m_withResampling)
        {
            inputSize = readerDPPIV->GetOutput()->GetLargestPossibleRegion().GetSize();
            inputSpacing = readerDPPIV->GetOutput()->GetSpacing();

            std::cout << "Input size: " << inputSize << std::endl;
            std::cout << "Input spacing: " << inputSpacing << std::endl;

            CScalarVoImageType::SizeType outputSize;
            outputSize[0] = std::ceil((double)inputSize[0]/2.);
            outputSize[1] = std::ceil((double)inputSize[1]/2.);
            outputSize[2] = inputSize[2];

            CScalarVoImageType::SpacingType outputSpacing;
            outputSpacing[0] = inputSpacing[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
            outputSpacing[1] = inputSpacing[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
            outputSpacing[2] = inputSpacing[2];

            std::cout << "Output size: " << outputSize << std::endl;
            std::cout << "Output spacing: " << outputSpacing << std::endl;

            if(m_hasCh0) {
                resampleDAPI = ResampleImageFilterType::New();
                resampleDAPI->SetSize(outputSize);
                resampleDAPI->SetOutputSpacing(outputSpacing);
                resampleDAPI->SetExtrapolator(extrapolator);
                resampleDAPI->SetInterpolator(interpolator);
                resampleDAPI->ReleaseDataFlagOn();
                resampleDAPI->SetInput(readerDAPI->GetOutput());
            }

            resampleDPPIV = ResampleImageFilterType::New();
            resampleDPPIV->SetSize(outputSize);
            resampleDPPIV->SetOutputSpacing(outputSpacing);
            resampleDPPIV->SetExtrapolator(extrapolator);
            resampleDPPIV->SetInterpolator(interpolator);
            resampleDPPIV->ReleaseDataFlagOn();
            resampleDPPIV->SetInput(readerDPPIV->GetOutput());

            resampleDM = ResampleImageFilterType::New();
            resampleDM->SetSize(outputSize);
            resampleDM->SetOutputSpacing(outputSpacing);
            resampleDM->SetExtrapolator(extrapolator);
            resampleDM->SetInterpolator(interpolator);
            resampleDM->ReleaseDataFlagOn();
            resampleDM->SetInput(readerDM->GetOutput());

            double scaling[3];
            scaling[0] = (static_cast<double>(outputSize[0]) / static_cast<double>(inputSize[0]));
            scaling[1] = (static_cast<double>(outputSize[1]) / static_cast<double>(inputSize[1]));
            scaling[2] = (static_cast<double>(outputSize[2]) / static_cast<double>(inputSize[2]));

            for(unsigned int i=0; i<m_CVSeedPoints.size(); i++) {
                m_CVSeedPoints[i][0] *= scaling[0];
                m_CVSeedPoints[i][1] *= scaling[1];
                m_CVSeedPoints[i][2] *= scaling[2];
            }

            for(unsigned int i=0; i<m_PVSeedPoints.size(); i++) {
                m_PVSeedPoints[i][0] *= scaling[0];
                m_PVSeedPoints[i][1] *= scaling[1];
                m_PVSeedPoints[i][2] *= scaling[2];
            }
        }
        //-------------------------------------------------------------------------------------------------------------------------

        //----------FILTER-COMPOSE-VECTOR-IMAGE-OF-TWO-CHANNELS-------------------------------------------------------------------
        if(m_hasCh0) {
            composeDAPIDPPIVDM = Compose3DVectorImageFilterType::New();
            composeDAPIDPPIVDM->SetNumberOfThreads(4);
            composeDAPIDPPIVDM->ReleaseDataFlagOn();
            if(m_withResampling) {
                composeDAPIDPPIVDM->SetInput1(resampleDPPIV->GetOutput());
                composeDAPIDPPIVDM->SetInput2(resampleDM->GetOutput());
                composeDAPIDPPIVDM->SetInput3(resampleDAPI->GetOutput());
            }
            else {
                composeDAPIDPPIVDM->SetInput1(readerDPPIV->GetOutput());
                composeDAPIDPPIVDM->SetInput2(readerDM->GetOutput());
                composeDAPIDPPIVDM->SetInput3(readerDAPI->GetOutput());
            }
        }
        else {
            composeDPPIVDM = Compose2DVectorImageFilterType::New();
            composeDPPIVDM->SetNumberOfThreads(4);
            composeDPPIVDM->ReleaseDataFlagOn();
            if(m_withResampling) {
                composeDPPIVDM->SetInput1(resampleDPPIV->GetOutput());
                composeDPPIVDM->SetInput2(resampleDM->GetOutput());
            }
            else {
                composeDPPIVDM->SetInput1(readerDPPIV->GetOutput());
                composeDPPIVDM->SetInput2(readerDM->GetOutput());
            }
        }
        //-------------------------------------------------------------------------------------------------------------------------

        //----------FILTER-REGION-GROWING-ON-VECTOR-IMAGE--------------------------------------------------------------------------
        if(m_hasCh0) {
            itk::Vector<unsigned char, 3> lower;    lower[0] = 0;   lower[1] = 0;   lower[2] = 0;
            itk::Vector<unsigned char, 3> upper;    upper[0] = m_thresholdDPPIV; upper[1] = m_thresholdDM; upper[2] = m_thresholdDAPI;

            connected3DThresholdCV = Connected3DThresholdFilterType::New();
            for(unsigned int i=0; i<m_CVSeedPoints.size(); i++)
                connected3DThresholdCV->AddSeed(m_CVSeedPoints[i]);
            connected3DThresholdCV->SetLower(lower);
            connected3DThresholdCV->SetUpper(upper);
            connected3DThresholdCV->SetReplaceValue(255);
            connected3DThresholdCV->SetNumberOfThreads(4);
            connected3DThresholdCV->ReleaseDataFlagOn();
            connected3DThresholdCV->SetInput(composeDAPIDPPIVDM->GetOutput());

            connected3DThresholdPV = Connected3DThresholdFilterType::New();
            for(unsigned int i=0; i<m_PVSeedPoints.size(); i++)
                connected3DThresholdPV->AddSeed(m_PVSeedPoints[i]);
            connected3DThresholdPV->SetLower(lower);
            connected3DThresholdPV->SetUpper(upper);
            connected3DThresholdPV->SetReplaceValue(255);
            connected3DThresholdPV->SetNumberOfThreads(4);
            connected3DThresholdPV->ReleaseDataFlagOn();
            connected3DThresholdPV->SetInput(composeDAPIDPPIVDM->GetOutput());
        }
        else {
            itk::Vector<unsigned char, 2> lower;    lower[0] = 0;   lower[1] = 0;
            itk::Vector<unsigned char, 2> upper;    upper[0] = m_thresholdDPPIV; upper[1] = m_thresholdDM;

            connected2DThresholdCV = Connected2DThresholdFilterType::New();
            for(unsigned int i=0; i<m_CVSeedPoints.size(); i++)
                connected2DThresholdCV->AddSeed(m_CVSeedPoints[i]);
            connected2DThresholdCV->SetLower(lower);
            connected2DThresholdCV->SetUpper(upper);
            connected2DThresholdCV->SetReplaceValue(255);
            connected2DThresholdCV->SetNumberOfThreads(4);
            connected2DThresholdCV->ReleaseDataFlagOn();
            connected2DThresholdCV->SetInput(composeDPPIVDM->GetOutput());

            connected2DThresholdPV = Connected2DThresholdFilterType::New();
            for(unsigned int i=0; i<m_PVSeedPoints.size(); i++)
                connected2DThresholdPV->AddSeed(m_PVSeedPoints[i]);
            connected2DThresholdPV->SetLower(lower);
            connected2DThresholdPV->SetUpper(upper);
            connected2DThresholdPV->SetReplaceValue(255);
            connected2DThresholdPV->SetNumberOfThreads(4);
            connected2DThresholdPV->ReleaseDataFlagOn();
            connected2DThresholdPV->SetInput(composeDPPIVDM->GetOutput());
        }
        //-------------------------------------------------------------------------------------------------------------------------

        //----------FILTER-ERODE---------------------------------------------------------------------------------------------------
        if(m_withResampling)
            m_openingStrucElem.SetRadius(m_openingStrucElem.GetRadius(0)/2);

        openingCV = OpeningImageFilterType::New();
        openingCV->SetKernel(m_openingStrucElem);
        openingCV->SetNumberOfThreads(4);
        openingCV->ReleaseDataFlagOn();
        if(m_hasCh0)	openingCV->SetInput(connected3DThresholdCV->GetOutput());
        else            openingCV->SetInput(connected2DThresholdCV->GetOutput());

        openingPV = OpeningImageFilterType::New();
        openingPV->SetKernel(m_openingStrucElem);
        openingPV->SetNumberOfThreads(4);
        openingPV->ReleaseDataFlagOn();
        if(m_hasCh0)	openingPV->SetInput(connected3DThresholdPV->GetOutput());
        else            openingPV->SetInput(connected2DThresholdPV->GetOutput());
        //-------------------------------------------------------------------------------------------------------------------------

        //----------FILTER-DILATE--------------------------------------------------------------------------------------------------
        if(m_withResampling)
            m_closingStrucElem.SetRadius(m_closingStrucElem.GetRadius(0)/2);

        closingCV = ClosingImageFilterType::New();
        closingCV->SetKernel(m_closingStrucElem);
        closingCV->SetNumberOfThreads(4);
        closingCV->ReleaseDataFlagOn();
        closingCV->SetInput(openingCV->GetOutput());

        closingPV = ClosingImageFilterType::New();
        closingPV->SetKernel(m_closingStrucElem);
        closingPV->SetNumberOfThreads(4);
        closingPV->ReleaseDataFlagOn();
        closingPV->SetInput(openingPV->GetOutput());
        //-------------------------------------------------------------------------------------------------------------------------

        //----------RESAMPLE-IMAGES------------------------------------------------------------------------------------------------
        if(m_withResampling)
        {
            resampleCV1 = ResampleImageFilterType::New();
            resampleCV1->SetSize(inputSize);
            resampleCV1->SetOutputSpacing(inputSpacing);
            resampleCV1->SetOutputSpacing(inputSpacing);
            resampleCV1->SetExtrapolator(extrapolator);
            resampleCV1->SetInterpolator(interpolator);
            resampleCV1->ReleaseDataFlagOn();
            resampleCV1->SetInput(closingCV->GetOutput());

            resamplePV1 = ResampleImageFilterType::New();
            resamplePV1->SetSize(inputSize);
            resamplePV1->SetOutputSpacing(inputSpacing);
            resamplePV1->SetExtrapolator(extrapolator);
            resamplePV1->SetInterpolator(interpolator);
            resamplePV1->ReleaseDataFlagOn();
            resamplePV1->ReleaseDataFlagOn();
            resamplePV1->SetInput(closingPV->GetOutput());
        }
        //-------------------------------------------------------------------------------------------------------------------------

        //----------FILTER---CREATE-LABEL-MAP--------------------------------------------------------------------------------------
        imageToShapeLabelMapCV = ImageToShapeLabelMapFilterType::New();
        imageToShapeLabelMapCV->SetNumberOfThreads(4);
        imageToShapeLabelMapCV->ReleaseDataFlagOn();
        imageToShapeLabelMapCV->FullyConnectedOff();
        if(m_withResampling)    imageToShapeLabelMapCV->SetInput(resampleCV1->GetOutput());
        else                    imageToShapeLabelMapCV->SetInput(closingCV->GetOutput());
        //        imageToShapeLabelMapCV->Update();

        imageToShapeLabelMapPV = ImageToShapeLabelMapFilterType::New();
        imageToShapeLabelMapPV->SetNumberOfThreads(4);
        imageToShapeLabelMapPV->ReleaseDataFlagOn();
        imageToShapeLabelMapPV->FullyConnectedOff();
        if(m_withResampling)    imageToShapeLabelMapPV->SetInput(resamplePV1->GetOutput());
        else                    imageToShapeLabelMapPV->SetInput(closingPV->GetOutput());
        //        imageToShapeLabelMapPV->Update();

        //        std::cout << "number label objects CV = " << imageToShapeLabelMapCV->GetOutput()->GetNumberOfLabelObjects() << std::endl;
        //        std::cout << "number label objects PV = " << imageToShapeLabelMapPV->GetOutput()->GetNumberOfLabelObjects() << std::endl;
        //-------------------------------------------------------------------------------------------------------------------------

        //----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------
        shapeOpeningCV = ShapeOpeningLabelMapFilterType::New();
        shapeOpeningCV->SetLambda(m_minimalVeinSize);                             //attribute value
        shapeOpeningCV->ReverseOrderingOff();                                     //removes objects with attribute smaller than lambda
        shapeOpeningCV->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::NUMBER_OF_PIXELS);
        shapeOpeningCV->SetNumberOfThreads(4);
        shapeOpeningCV->ReleaseDataFlagOn();
        shapeOpeningCV->SetInput(imageToShapeLabelMapCV->GetOutput());
        shapeOpeningCV->Update();

        shapeOpeningPV = ShapeOpeningLabelMapFilterType::New();
        shapeOpeningPV->SetLambda(m_minimalVeinSize);                             //attribute value
        shapeOpeningPV->ReverseOrderingOff();                                     //removes objects with attribute smaller than lambda
        shapeOpeningPV->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::NUMBER_OF_PIXELS);
        shapeOpeningPV->SetNumberOfThreads(4);
        shapeOpeningPV->ReleaseDataFlagOn();
        shapeOpeningPV->SetInput(imageToShapeLabelMapPV->GetOutput());
        shapeOpeningPV->Update();

        std::vector<unsigned long> labelsToRemove;

        ShapeOpeningLabelMapFilterType::InputImageType::Pointer cvLabelMap = shapeOpeningCV->GetOutput();
        cvLabelMap->DisconnectPipeline();
        std::cout << "number label objects CV after shape opening= " << cvLabelMap->GetNumberOfLabelObjects() << std::endl;
        for(unsigned int i=0; i<cvLabelMap->GetNumberOfLabelObjects(); i++) {
            cvLabelMap->GetNthLabelObject(i)->Optimize();
            std::cout << "label object " << i << " has size " << cvLabelMap->GetNthLabelObject(i)->Size() << std::endl;

            bool isLabelObjectToOneOfTheSeedPoints = false;
            for(unsigned int j=0; j<m_CVSeedPoints.size(); j++)
                if(cvLabelMap->GetNthLabelObject(i)->HasIndex(m_CVSeedPoints[j])) {
                    isLabelObjectToOneOfTheSeedPoints = true;
                    break;
                }
            if(!isLabelObjectToOneOfTheSeedPoints)
                labelsToRemove.push_back(cvLabelMap->GetNthLabelObject(i)->GetLabel());
        }
        for(unsigned int i=0; i<labelsToRemove.size(); i++)
            cvLabelMap->RemoveLabel(labelsToRemove[i]);
        std::cout << "number label objects CV after removal of to seed points unconnected label objects = " << cvLabelMap->GetNumberOfLabelObjects() << std::endl;

        labelsToRemove.clear();

        ShapeOpeningLabelMapFilterType::InputImageType::Pointer pvLabelMap = shapeOpeningPV->GetOutput();
        pvLabelMap->DisconnectPipeline();
        std::cout << "number label objects PV after shape opening= " << pvLabelMap->GetNumberOfLabelObjects() << std::endl;
        for(unsigned int i=0; i<pvLabelMap->GetNumberOfLabelObjects(); i++) {
            pvLabelMap->GetNthLabelObject(i)->Optimize();
            std::cout << "label object " << i << " has size " << pvLabelMap->GetNthLabelObject(i)->Size() << std::endl;

            bool isLabelObjectToOneOfTheSeedPoints = false;
            for(unsigned int j=0; j<m_PVSeedPoints.size(); j++)
                if(pvLabelMap->GetNthLabelObject(i)->HasIndex(m_PVSeedPoints[j])) {
                    isLabelObjectToOneOfTheSeedPoints = true;
                    break;
                }
            if(!isLabelObjectToOneOfTheSeedPoints)
                labelsToRemove.push_back(pvLabelMap->GetNthLabelObject(i)->GetLabel());
        }
        for(unsigned int i=0; i<labelsToRemove.size(); i++)
            pvLabelMap->RemoveLabel(labelsToRemove[i]);
        std::cout << "number label objects PV after removal of to seed points unconnected label objects = " << pvLabelMap->GetNumberOfLabelObjects() << std::endl;
        //-------------------------------------------------------------------------------------------------------------------------

        labelMapToImageCV = LabelMapToLabelImageFilterType::New();
        labelMapToImageCV->SetNumberOfThreads(4);
        labelMapToImageCV->SetInput(cvLabelMap);

        labelMapToImagePV = LabelMapToLabelImageFilterType::New();
        labelMapToImagePV->SetNumberOfThreads(4);
        labelMapToImagePV->SetInput(pvLabelMap);

        //----------VISULALIZATION-ONLY-FILTER---NOT-NEEDED-BY-OTHER-PIPELINE-FILTERS----------------------------------------------
        overlayImageCV = LabelOverlayImageFilterType::New();
        overlayImageCV->SetLabelImage(labelMapToImageCV->GetOutput());
        overlayImageCV->SetOpacity(m_overlayOpacity);
        overlayImageCV->ReleaseDataFlagOn();
        overlayImageCV->SetInput1(readerDPPIV->GetOutput());
        overlayImageCV->Update();

        overlayImagePV = LabelOverlayImageFilterType::New();
        overlayImagePV->SetLabelImage(labelMapToImagePV->GetOutput());
        overlayImagePV->SetOpacity(m_overlayOpacity);
        overlayImagePV->ReleaseDataFlagOn();
        overlayImagePV->SetInput1(readerDPPIV->GetOutput());
        overlayImagePV->Update();
        //-------------------------------------------------------------------------------------------------------------------------

        //----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
        thresholdCV = ThresholdFilterType::New();
        thresholdCV->SetOutsideValue(0);
        thresholdCV->SetInsideValue(255);
        thresholdCV->SetLowerThreshold(1);
        thresholdCV->SetUpperThreshold(255);
        thresholdCV->SetNumberOfThreads(4);
        thresholdCV->ReleaseDataFlagOn();
        thresholdCV->SetInput(labelMapToImageCV->GetOutput());

        thresholdPV = ThresholdFilterType::New();
        thresholdPV->SetOutsideValue(0);
        thresholdPV->SetInsideValue(255);
        thresholdPV->SetLowerThreshold(1);
        thresholdPV->SetUpperThreshold(255);
        thresholdPV->SetNumberOfThreads(4);
        thresholdPV->ReleaseDataFlagOn();
        thresholdPV->SetInput(labelMapToImagePV->GetOutput());
        //-------------------------------------------------------------------------------------------------------------------------

        //----------WRITER--------------------------------------------------------------------------------------------------------
        segmentationCVResult = m_pathChannel1 + m_cvSaveBin + m_fileExtensionChannel1;
        segmentationPVResult = m_pathChannel1 + m_pvSaveBin + m_fileExtensionChannel1;

        writerBinCV = ScalarVoWriterType::New();
        writerBinCV->SetFileName(segmentationCVResult);
        writerBinCV->ReleaseDataFlagOn();
        if(m_withResampling)    writerBinCV->SetInput(resampleCV1->GetOutput());
        else                    writerBinCV->SetInput(thresholdCV->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerBinCV->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerBinCV->Update();

        writerBinPV = ScalarVoWriterType::New();
        writerBinPV->SetFileName(segmentationPVResult);
        writerBinPV->ReleaseDataFlagOn();
        if(m_withResampling)    writerBinPV->SetInput(resamplePV1->GetOutput());
        else                    writerBinPV->SetInput(thresholdPV->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerBinPV->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerBinPV->Update();
        //-------------------------------------------------------------------------------------------------------------------------

        //----------WRITER--------------------------------------------------------------------------------------------------------
        writerOverlayCV = RGBVoWriterType::New();
        writerOverlayCV->ReleaseDataFlagOn();
        writerOverlayCV->SetFileName(m_pathChannel1 + m_cvSaveOverlay + m_fileExtensionChannel1);
        writerOverlayCV->SetInput(overlayImageCV->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerOverlayCV->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerOverlayCV->Update();

        writerOverlayPV = RGBVoWriterType::New();
        writerOverlayPV->ReleaseDataFlagOn();
        writerOverlayPV->SetFileName(m_pathChannel1 + m_pvSaveOverlay + m_fileExtensionChannel1);
        writerOverlayPV->SetInput(overlayImagePV->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerOverlayPV->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerOverlayPV->Update();
        //-------------------------------------------------------------------------------------------------------------------------


//        RGBVoWriterType::Pointer                maskWriterOverlayCV;
//        RGBVoWriterType::Pointer                maskWriterOverlayPV;
//
//        DilateImageFilterType::Pointer maskDilateCV = DilateImageFilterType::New();
//        maskDilateCV->ReleaseDataFlagOn();
//        maskDilateCV->SetKernel(m_maskDilateStrucElem);
//        maskDilateCV->SetInput(thresholdCV->GetOutput());
//
//        ImageToShapeLabelMapFilterType::Pointer maskImageToShapeLabelMapCV = ImageToShapeLabelMapFilterType::New();
//        maskImageToShapeLabelMapCV->SetInput(maskDilateCV->GetOutput());
//        maskImageToShapeLabelMapCV->Update();
//
//        LabelMapToLabelImageFilterType::Pointer maskLabelMapToImageCV = LabelMapToLabelImageFilterType::New();
//        maskLabelMapToImageCV->SetInput(maskImageToShapeLabelMapCV->GetOutput());
//
//        LabelOverlayImageFilterType::Pointer maskOverlayImageCV = LabelOverlayImageFilterType::New();
//        maskOverlayImageCV->SetLabelImage(maskLabelMapToImageCV->GetOutput());
//        maskOverlayImageCV->SetOpacity(0.5);
//        maskOverlayImageCV->ReleaseDataFlagOn();
//        maskOverlayImageCV->SetInput(readerDPPIV->GetOutput());
//        maskOverlayImageCV->Update();
//
//        DilateImageFilterType::Pointer maskDilatePV = DilateImageFilterType::New();
//        maskDilatePV->ReleaseDataFlagOn();
//        maskDilatePV->SetKernel(m_maskDilateStrucElem);
//        maskDilatePV->SetInput(thresholdPV->GetOutput());
//
//        ImageToShapeLabelMapFilterType::Pointer maskImageToShapeLabelMapPV = ImageToShapeLabelMapFilterType::New();
//        maskImageToShapeLabelMapPV->SetInput(maskDilatePV->GetOutput());
//        maskImageToShapeLabelMapPV->Update();
//
//        LabelMapToLabelImageFilterType::Pointer maskLabelMapToImagePV = LabelMapToLabelImageFilterType::New();
//        maskLabelMapToImagePV->SetInput(maskImageToShapeLabelMapPV->GetOutput());
//
//        LabelOverlayImageFilterType::Pointer maskOverlayImagePV = LabelOverlayImageFilterType::New();
//        maskOverlayImagePV->SetLabelImage(maskLabelMapToImagePV->GetOutput());
//        maskOverlayImagePV->SetOpacity(0.5);
//        maskOverlayImagePV->ReleaseDataFlagOn();
//        maskOverlayImagePV->SetInput(readerDPPIV->GetOutput());
//        maskOverlayImagePV->Update();
//
//        segmentationMaskCVResult = m_pathChannel1 + m_cvMaskSaveBin + m_fileExtensionChannel1;
//        segmentationMaskPVResult = m_pathChannel1 + m_pvMaskSaveBin + m_fileExtensionChannel1;
//
//        ScalarVoWriterType::Pointer writerBinCVMask = ScalarVoWriterType::New();
//        writerBinCVMask->SetFileName(segmentationMaskCVResult);
//        writerBinCVMask->SetInput(maskDilateCV->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//        writerBinCVMask->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//        writerBinCVMask->Update();
//
//        RGBVoWriterType::Pointer writerOverlayCVMask = RGBVoWriterType::New();
//        writerOverlayCVMask->ReleaseDataFlagOn();
//        writerOverlayCVMask->SetFileName(m_pathChannel1 + m_cvMaskSaveOverlay + m_fileExtensionChannel1);
//        writerOverlayCVMask->SetInput(maskOverlayImageCV->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//        writerOverlayCVMask->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//        writerOverlayCVMask->Update();
//
//        ScalarVoWriterType::Pointer writerBinPVMask = ScalarVoWriterType::New();
//        writerBinPVMask->SetFileName(segmentationMaskPVResult);
//        writerBinPVMask->SetInput(maskDilatePV->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//        writerBinPVMask->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//        writerBinPVMask->Update();
//
//        RGBVoWriterType::Pointer writerOverlayPVMask = RGBVoWriterType::New();
//        writerOverlayPVMask->ReleaseDataFlagOn();
//        writerOverlayPVMask->SetFileName(m_pathChannel1 + m_pvMaskSaveOverlay + m_fileExtensionChannel1);
//        writerOverlayPVMask->SetInput(maskOverlayImagePV->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//        writerOverlayPVMask->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//        writerOverlayPVMask->Update();
    }

    if(m_withLobulePrediction)
    {
        //----------READER--------------------------------------------------------------------------------------------------------
        readerCV = ScalarVoReaderType::New();
        if(m_withSegmentation)  readerCV->SetFileName(segmentationCVResult);
        else                    readerCV->SetFileName(m_pathChannel1 + m_filenameChannel1 + m_fileExtensionChannel1);
        readerCV->ReleaseDataBeforeUpdateFlagOn();
        readerCV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
        readerCV->SetImageIO( itk::TIFFImageIO::New() );
#endif
        readerCV->Update();

        readerPV = ScalarVoReaderType::New();
        if(m_withSegmentation)  readerPV->SetFileName(segmentationPVResult);
        else                    readerPV->SetFileName(m_pathChannel2 + m_filenameChannel2 + m_fileExtensionChannel2);
        readerPV->ReleaseDataBeforeUpdateFlagOn();
        readerPV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
        readerPV->SetImageIO( itk::TIFFImageIO::New() );
#endif
        //-------------------------------------------------------------------------------------------------------------------------

        //----------RESAMPLE-IMAGES------------------------------------------------------------------------------------------------
        if(m_withResampling)
        {
            inputSize = readerCV->GetOutput()->GetLargestPossibleRegion().GetSize();
            inputSpacing = readerCV->GetOutput()->GetSpacing();

            std::cout << "Input size: " << inputSize << std::endl;
            std::cout << "Input spacing: " << inputSpacing << std::endl;

            CScalarVoImageType::SizeType outputSize;
            outputSize[0] = 256;
            outputSize[1] = 256;
            outputSize[2] = inputSize[2];

            CScalarVoImageType::SpacingType outputSpacing;
            outputSpacing[0] = inputSpacing[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
            outputSpacing[1] = inputSpacing[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
            outputSpacing[2] = inputSpacing[2];

            std::cout << "Output size: " << outputSize << std::endl;
            std::cout << "Output spacing: " << outputSpacing << std::endl;

            resampleCV2 = ResampleImageFilterType::New();
            resampleCV2->SetSize(outputSize);
            resampleCV2->SetOutputSpacing(outputSpacing);
            resampleCV2->ReleaseDataFlagOn();
            resampleCV2->SetInput(readerCV->GetOutput());

            resamplePV2 = ResampleImageFilterType::New();
            resamplePV2->SetSize(outputSize);
            resamplePV2->SetOutputSpacing(outputSpacing);
            resamplePV2->ReleaseDataFlagOn();
            resamplePV2->SetInput(readerPV->GetOutput());
        }
        //-------------------------------------------------------------------------------------------------------------------------

        //----------DISTANCE-MAP---------------------------------------------------------------------------------------------------
        distanceMapCV = DanielssonDistanceMapImageFilterType::New();
        //    distanceMapCV->ReleaseDataFlagOn();
        if(m_withResampling)    distanceMapCV->SetInput(resampleCV2->GetOutput());
        else                    distanceMapCV->SetInput(readerCV->GetOutput());

        distanceMapPV = DanielssonDistanceMapImageFilterType::New();
        //    distanceMapPV->ReleaseDataFlagOn();
        if(m_withResampling)    distanceMapPV->SetInput(resamplePV2->GetOutput());
        else                    distanceMapPV->SetInput(readerPV->GetOutput());
        //-------------------------------------------------------------------------------------------------------------------------

        //----------CALCULATE-REFERENCE-IMAGE--------------------------------------------------------------------------------------
        add = AddImageFilterType::New();
        //    add->InplaceOn();
        add->ReleaseDataFlagOn();
        add->SetInput1(distanceMapCV->GetOutput());
        add->SetInput2(distanceMapPV->GetOutput());

        divide = DivideImageFilterType::New();
        //    divide->InplaceOn();
        //    divide->ReleaseDataFlagOn();
        divide->SetInput1(distanceMapCV->GetOutput());
        divide->SetInput2(add->GetOutput());
        //-------------------------------------------------------------------------------------------------------------------------

        //----------WATERSHED-----------------------------------------------------------------------------------------------------
        morphWatershed = MorphoWatershedImageFilterType::New();
        morphWatershed->SetLevel(0.1);
        morphWatershed->FullyConnectedOff();
        morphWatershed->ReleaseDataFlagOn();
        morphWatershed->SetInput(divide->GetOutput());
        //-------------------------------------------------------------------------------------------------------------------------

        //----------TO-LABEL-MAP-AND-BACK-----------------------------------------------------------------------------------------
        labelImageToLabelMap = LabelImageToLabelMapFilterType::New ();
        labelImageToLabelMap->SetInput(morphWatershed->GetOutput());

        labelMapToImageW = LabelMapToLabelImageFilterTypeW::New();
        labelMapToImageW->SetInput(labelImageToLabelMap->GetOutput());
        //-------------------------------------------------------------------------------------------------------------------------

        //    //----------CONTOUR-WATERSHED-REGIONS-------------------------------------------------------------------------------------
        //    labelContour = LabelContourImageFilterType::New();
        //    labelContour->ReleaseDataFlagOn();
        //    labelContour->SetInput(morphWatershed->GetOutput());
        //    //-------------------------------------------------------------------------------------------------------------------------

        //----------RESCALE-IMAGES------------------------------------------------------------------------------------------------
        rescaler1 = RescaleImageFilterType::New();
        rescaler1->ReleaseDataFlagOn();
        rescaler1->SetInput(distanceMapCV->GetOutput());

        rescaler2 = RescaleImageFilterType::New();
        rescaler2->ReleaseDataFlagOn();
        rescaler2->SetInput(distanceMapPV->GetOutput());

        rescaler3 = RescaleImageFilterType::New();
        rescaler3->ReleaseDataFlagOn();
        //    rescaler3->SetInplaceOn();
        rescaler3->SetInput(divide->GetOutput());
        //-------------------------------------------------------------------------------------------------------------------------


        //----------RESAMPLE-IMAGES------------------------------------------------------------------------------------------------
        if(m_withResampling)
        {
            resample1 = ResampleImageFilterType::New();
            resample1->SetSize(inputSize);
            resample1->SetOutputSpacing(inputSpacing);
            resample1->ReleaseDataFlagOn();
            resample1->SetInput(rescaler1->GetOutput());

            resample2 = ResampleImageFilterType::New();
            resample2->SetSize(inputSize);
            resample2->SetOutputSpacing(inputSpacing);
            resample2->ReleaseDataFlagOn();
            resample2->SetInput(rescaler2->GetOutput());

            resample3 = ResampleImageFilterType::New();
            resample3->SetSize(inputSize);
            resample3->SetOutputSpacing(inputSpacing);
            //        resample3->ReleaseDataFlagOn();
            resample3->SetInput(rescaler3->GetOutput());

            resample4 = ResampleImageFilterType::New();
            resample4->SetSize(inputSize);
            resample4->SetOutputSpacing(inputSpacing);
            resample4->ReleaseDataFlagOn();
            resample4->SetInput(labelMapToImageW->GetOutput());
        }
        //-------------------------------------------------------------------------------------------------------------------------

        //----------VISULALIZATION-ONLY-FILTER---NOT-NEEDED-BY-OTHER-PIPELINE-FILTERS----------------------------------------------
        overlayImageW = LabelOverlayImageFilterType::New();
        overlayImageW->SetOpacity(m_overlayOpacity);
        overlayImageW->ReleaseDataFlagOn();
        if(m_withResampling) {
            overlayImageW->SetInput(resample3->GetOutput());
            overlayImageW->SetLabelImage(resample4->GetOutput());
        }
        else {
            overlayImageW->SetInput(rescaler3->GetOutput());
            overlayImageW->SetLabelImage(labelMapToImageW->GetOutput());
        }
        //-------------------------------------------------------------------------------------------------------------------------

        //----------WRITER--------------------------------------------------------------------------------------------------------
        writerW1 = ScalarVoWriterType::New();
        writerW1->ReleaseDataFlagOn();
        writerW1->SetFileName(m_pathChannel1 + "bin_distMap_CV" + m_fileExtensionChannel1);
        if(m_withResampling)    writerW1->SetInput(resample1->GetOutput());
        else                    writerW1->SetInput(rescaler1->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerW1->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerW1->Update();

        writerW2 = ScalarVoWriterType::New();
        writerW2->ReleaseDataFlagOn();
        writerW2->SetFileName(m_pathChannel1 + "bin_distMap_PV" + m_fileExtensionChannel1);
        if(m_withResampling)    writerW2->SetInput(resample2->GetOutput());
        else                    writerW2->SetInput(rescaler2->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerW2->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerW2->Update();

        writerW3 = ScalarVoWriterType::New();
        writerW3->ReleaseDataFlagOn();
        writerW3->SetFileName(m_pathChannel1 + "bin_refData" + m_fileExtensionChannel1);
        if(m_withResampling)    writerW3->SetInput(resample3->GetOutput());
        else                    writerW3->SetInput(rescaler3->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerW3->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerW3->Update();

        //    writerW4 = ScalarVoWriterType::New();
        //    writerW4->ReleaseDataFlagOn();
        //    writerW4->SetFileName(m_pathChannel1 + "label_contour" + m_fileExtensionChannel1);
        //    writerW4->SetInput(labelContour->GetOutput());
        //#if (ITK_VERSION_MAJOR >= 4)
        //    writerW4->SetImageIO( itk::TIFFImageIO::New() );
        //#endif
        //    writerW4->Update();

        writerWRGB1 = RGBVoWriterType::New();
        writerWRGB1->ReleaseDataFlagOn();
        writerWRGB1->SetFileName(m_pathChannel1 + "overlay_watershed" + m_fileExtensionChannel1);
        writerWRGB1->SetInput(overlayImageW->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerWRGB1->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerWRGB1->Update();
        //-------------------------------------------------------------------------------------------------------------------------
    }

    WriteLogFile(timeStamp);
    WriteDataSetSummary();
}

