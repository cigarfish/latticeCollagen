///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentNecroticRegionOnDM.cpp                                        //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-05-14                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "SegmentNecroticRegionOnDM.h"

#include <QDateTime>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"


SegmentNecroticRegionOnDM::SegmentNecroticRegionOnDM()
{
    // TODO Auto-generated destructor stub
    m_overlayOpacity = 0.6;

    m_finalSaveSuffixes[0] = "_step7_bin";
    m_finalSaveSuffixes[1] = "_step7_DM_overlay";
    m_finalSaveSuffixes[2] = "_step7_DPPIV_overlay";
}


SegmentNecroticRegionOnDM::~SegmentNecroticRegionOnDM()
{
    // TODO Auto-generated destructor stub
}


void SegmentNecroticRegionOnDM::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

    file.open((m_pathChannelDM + m_logFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-segment-necrotic-region---------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-segment-necrotic-region---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


void SegmentNecroticRegionOnDM::WriteDataSetSummary()
{
    ImageAnalysisSummaryFileIO::AddEntry(NecroticRegionSegmentationBin, m_pathChannelDM, m_pathChannelDM + m_filenameSave + m_finalSaveSuffixes[0] + m_fileExtensionChannelDM);
    ImageAnalysisSummaryFileIO::AddEntry(NecroticRegionSegmentationDMsOverlay, m_pathChannelDM, m_pathChannelDM + m_filenameSave + m_finalSaveSuffixes[1] + m_fileExtensionChannelDM);
    ImageAnalysisSummaryFileIO::AddEntry(NecroticRegionSegmentationDPPIVOverlay, m_pathChannelDM, m_pathChannelDM + m_filenameSave + m_finalSaveSuffixes[2] + m_fileExtensionChannelDM);
}


void SegmentNecroticRegionOnDM::ParseParameterContext(){
    if(m_paramContext->findContext("Segment Necrotic Region",0)==NULL) {
        std::cout << "Error: SegmentNecroticRegionOnDMs: Invalid parameter context" << std::endl;
        return;
    }

	m_fullFilenameChannelDPPIV = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("DPPIV channel", 0)->dataPointer()) );
	m_infoFullFilenameChannelDPPIV.setFile(m_fullFilenameChannelDPPIV);

	if(!m_infoFullFilenameChannelDPPIV.exists())
		throw std::string("Please specify DPPIV channel");

	m_fullFilenameChannelDM = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("DMs channel", 0)->dataPointer()) );
	m_infoFullFilenameChannelDM.setFile(m_fullFilenameChannelDM);

	if(!m_infoFullFilenameChannelDM.exists())
		throw std::string("Please specify DMs channel");

	m_pathChannelDPPIV = (m_infoFullFilenameChannelDPPIV.path() + QString("/")).toStdString();
    m_filenameChannelDPPIV = m_infoFullFilenameChannelDPPIV.baseName().toStdString();
    m_fileExtensionChannelDPPIV = (QString(".") + m_infoFullFilenameChannelDPPIV.suffix()).toStdString();

    m_pathChannelDM = (m_infoFullFilenameChannelDM.path() + QString("/")).toStdString();
    m_filenameChannelDM = m_infoFullFilenameChannelDM.baseName().toStdString();
    m_fileExtensionChannelDM = (QString(".") + m_infoFullFilenameChannelDM.suffix()).toStdString();


    std::string thresMode = ( (CSParameterChoice*)(m_paramContext->findParameter("Threshold mode", 0)->dataPointer()) )->currentString();
    if( thresMode.compare("Adaptive Otsu-Threshold")==0 )
        m_thresholdingMode = 0;
    else if( thresMode.compare("Otsu-Threshold")==0 )
        m_thresholdingMode = 1;
    else if( thresMode.compare("Manual Threshold")==0 )
        m_thresholdingMode = 2;

    m_adapOtsuRadius[0] = *(int*)(m_paramContext->findParameter("Sample region size x", 0)->dataPointer());
    m_adapOtsuRadius[1] = *(int*)(m_paramContext->findParameter("Sample region size y", 0)->dataPointer());
    m_adapOtsuRadius[2] = *(int*)(m_paramContext->findParameter("Sample region size z", 0)->dataPointer());
    m_adapOtsuSamplePoints = *(int*)(m_paramContext->findParameter("Number samples", 0)->dataPointer());
    m_lowerThreshold = *(int*)(m_paramContext->findParameter("DMs manual threshold", 0)->dataPointer());
    //TODO max;
    m_upperThreshold = 255;

    m_erode1StructuringElement.SetRadius( *(int*)(m_paramContext->findContext("1.2) Erode A on 1.1",0)->findParameter("Kernel radius", 0)->dataPointer()) );
    m_erode1StructuringElement.CreateStructuringElement();

    m_dilate1StructuringElement.SetRadius( *(int*)(m_paramContext->findContext("1.3) Dilate A on 1.2",0)->findParameter("Kernel radius", 0)->dataPointer()) );
    m_dilate1StructuringElement.CreateStructuringElement();

    m_erode2StructuringElement.SetRadius( *(int*)(m_paramContext->findContext("1.4) Erode B on 1.3",0)->findParameter("Kernel radius", 0)->dataPointer()) );
    m_erode2StructuringElement.CreateStructuringElement();

    m_dilate2StructuringElement.SetRadius( *(int*)(m_paramContext->findContext("1.5) Dilate B on 1.4",0)->findParameter("Kernel radius", 0)->dataPointer()) );
    m_dilate2StructuringElement.CreateStructuringElement();

    m_minimalNecroticRegionSize = *(int*)(m_paramContext->findParameter("Remove objects with volume smaller than", 0)->dataPointer());

    std::string saveMode = ( (CSParameterChoice*)(m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        m_saveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        m_saveEverything = 0;
    m_logFilenameSave = *(std::string*)(m_paramContext->findParameter("Log file name", 0)->dataPointer());
    m_filenameSave = *(std::string*)(m_paramContext->findParameter("Save prefix", 0)->dataPointer());
}


void SegmentNecroticRegionOnDM::Update()
{
    ParseParameterContext();

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();


    std::cout << "Segment necrotic region: " << std::endl;
    std::cout << " dir DPPIV: " << m_pathChannelDPPIV << std::endl;
    std::cout << " file DPPIV: " << m_filenameChannelDPPIV << std::endl;
    std::cout << " ext DPPIV: " << m_fileExtensionChannelDPPIV << std::endl;

    std::cout << " dir DM: " << m_pathChannelDM << std::endl;
    std::cout << " file DM: " << m_filenameChannelDM << std::endl;
    std::cout << " ext DM: " << m_fileExtensionChannelDM << std::endl;


    ScalarVoReaderType::Pointer                     readerDPPIV;
    ScalarVoReaderType::Pointer                     readerDM;
    AdaptiveOtsuThresholdImageFilterType::Pointer   adapOtsuFilter;
	OtsuThresholdImageFilterType::Pointer		    otsuFilter;
	InvertIntensityImageFilterType::Pointer		    invertFilter;
    ThresholdFilterType::Pointer                    thresDMFilter;
    ErodeImageFilterType::Pointer                   erode1DMFilter;
    DilateImageFilterType::Pointer                  dilate1DMFilter;
    ErodeImageFilterType::Pointer                   erode2DMFilter;
    DilateImageFilterType::Pointer                  dilate2DMFilter;
    ImageToShapeLabelMapFilterType::Pointer         imageToNecRegShaLabMapFilter;
    ShapeOpeningLabelMapFilterType::Pointer         shapeOpeningNecRegLabMapFilter;
    LabelMapToLabelImageFilterType::Pointer         necRegLabMapToImageFilter;
    LabelOverlayImageFilterType::Pointer            necRegDMOverlayImageFilter;
    LabelOverlayImageFilterType::Pointer            necRegDPPIVOverlayImageFilter;
    ThresholdFilterType::Pointer                    necRegImageFilter;


    //----------READER---------------------------------------------------------------------------------------------------------
    if(GetNumberOfDimensions(m_pathChannelDPPIV + m_filenameChannelDPPIV + m_fileExtensionChannelDPPIV) != 3)
        throw std::string("Please specify a DPPIV channel file with three dimensions.");

    readerDPPIV = ScalarVoReaderType::New();
    readerDPPIV->SetFileName(m_pathChannelDPPIV + m_filenameChannelDPPIV + m_fileExtensionChannelDPPIV);
    readerDPPIV->ReleaseDataBeforeUpdateFlagOn();
    readerDPPIV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerDPPIV->SetImageIO( itk::TIFFImageIO::New() );
#endif

    if(GetNumberOfDimensions(m_pathChannelDM + m_filenameChannelDM + m_fileExtensionChannelDM) != 3)
        throw std::string("Please specify a DMs channel file with three dimensions.");

    readerDM = ScalarVoReaderType::New();
    readerDM->SetFileName(m_pathChannelDM + m_filenameChannelDM + m_fileExtensionChannelDM);
    readerDM->ReleaseDataBeforeUpdateFlagOn();
    readerDM->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerDM->SetImageIO( itk::TIFFImageIO::New() );
#endif
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-THRESHOLD-ON-DM-CHANNEL---------------------------------------------------------------------------------
    switch(m_thresholdingMode)
    {
    case 0:
        {
        adapOtsuFilter = AdaptiveOtsuThresholdImageFilterType::New();
        adapOtsuFilter->SetInput(readerDM->GetOutput());
        adapOtsuFilter->SetInsideValue(255);
        adapOtsuFilter->SetOutsideValue(0);
        adapOtsuFilter->SetNumberOfHistogramBins(256);
        adapOtsuFilter->SetSplineOrder(3);
        adapOtsuFilter->SetNumberOfControlPoints(10);
        adapOtsuFilter->SetNumberOfLevels(3);
        adapOtsuFilter->SetNumberOfSamples(m_adapOtsuSamplePoints);
        adapOtsuFilter->SetRadius(m_adapOtsuRadius);
        break;
        }
    case 1:
        {
        otsuFilter = OtsuThresholdImageFilterType::New();
        otsuFilter->SetInput(readerDM->GetOutput());
        otsuFilter->Update();
        m_otsuThreshold = (CScalarPixelType)(otsuFilter->GetThreshold());

        invertFilter = InvertIntensityImageFilterType::New();
        invertFilter->SetInput(otsuFilter->GetOutput());
        invertFilter->SetMaximum(255);
        break;
        }
    default:
        {
        thresDMFilter = ThresholdFilterType::New();
        thresDMFilter->ReleaseDataFlagOn();
        thresDMFilter->SetOutsideValue(0);
        thresDMFilter->SetInsideValue(255);
        thresDMFilter->SetLowerThreshold(m_lowerThreshold);
        thresDMFilter->SetUpperThreshold(m_upperThreshold);
        thresDMFilter->SetInput(readerDM->GetOutput());
        break;
        }
    }

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer1 = ScalarVoWriterType::New();
        writer1->ReleaseDataFlagOn();
        writer1->SetFileName(m_pathChannelDM + m_filenameSave + "_step1_bin" + m_fileExtensionChannelDM);
        switch(m_thresholdingMode)
        {
        case 0:
            writer1->SetInput(adapOtsuFilter->GetOutput());
            break;
        case 1:
            writer1->SetInput(invertFilter->GetOutput());
            break;
        default:
            writer1->SetInput(thresDMFilter->GetOutput());
            break;
        }
#if (ITK_VERSION_MAJOR >= 4)
        writer1->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer1->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-ERODE-ON-DM-CHANNEL-------------------------------------------------------------------------------------
    erode1DMFilter = ErodeImageFilterType::New();
    erode1DMFilter->ReleaseDataFlagOn();
    erode1DMFilter->SetKernel(m_erode1StructuringElement);
    switch(m_thresholdingMode)
    {
    case 0:
        erode1DMFilter->SetInput(adapOtsuFilter->GetOutput());
        break;
    case 1:
        erode1DMFilter->SetInput(invertFilter->GetOutput());
        break;
    default:
        erode1DMFilter->SetInput(thresDMFilter->GetOutput());
        break;
    }

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer2 = ScalarVoWriterType::New();
        writer2->ReleaseDataFlagOn();
        writer2->SetFileName(m_pathChannelDM + m_filenameSave + "_step2_erodeA" + m_fileExtensionChannelDM);
        writer2->SetInput(erode1DMFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer2->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer2->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-DILATE-ON-DM-CHANNEL------------------------------------------------------------------------------------
    dilate1DMFilter = DilateImageFilterType::New();
    dilate1DMFilter->ReleaseDataFlagOn();
    dilate1DMFilter->SetKernel(m_dilate1StructuringElement);
    dilate1DMFilter->SetInput(erode1DMFilter->GetOutput());

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer2 = ScalarVoWriterType::New();
        writer2->ReleaseDataFlagOn();
        writer2->SetFileName(m_pathChannelDM + m_filenameSave + "_step3_dilateA" + m_fileExtensionChannelDM);
        writer2->SetInput(dilate1DMFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer2->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer2->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-ERODE-ON-DM-CHANNEL-------------------------------------------------------------------------------------
    erode2DMFilter = ErodeImageFilterType::New();
    erode2DMFilter->ReleaseDataFlagOn();
    erode2DMFilter->SetKernel(m_erode2StructuringElement);
    erode2DMFilter->SetInput(dilate1DMFilter->GetOutput());

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer4 = ScalarVoWriterType::New();
        writer4->ReleaseDataFlagOn();
        writer4->SetFileName(m_pathChannelDM + m_filenameSave + "_step4_erodeB" + m_fileExtensionChannelDM);
        writer4->SetInput(erode2DMFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer4->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer4->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-DILATE-ON-DM-CHANNEL------------------------------------------------------------------------------------
    dilate2DMFilter = DilateImageFilterType::New();
    dilate2DMFilter->ReleaseDataFlagOn();
    dilate2DMFilter->SetKernel(m_dilate2StructuringElement);
    dilate2DMFilter->SetInput(erode2DMFilter->GetOutput());

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer5 = ScalarVoWriterType::New();
        writer5->ReleaseDataFlagOn();
        writer5->SetFileName(m_pathChannelDM + m_filenameSave + "_step5_dilateB" + m_fileExtensionChannelDM);
        writer5->SetInput(dilate2DMFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer5->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer5->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER---CREATE-LABEL-MAP-ON-DM-CHANNEL-DO-SOME-STUFF-AND-CONVERT-IT-BACK-INTO-IMAGE---------------------------
    imageToNecRegShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    imageToNecRegShaLabMapFilter->SetInput(dilate2DMFilter->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------
    shapeOpeningNecRegLabMapFilter = ShapeOpeningLabelMapFilterType::New();
    shapeOpeningNecRegLabMapFilter->SetLambda(m_minimalNecroticRegionSize);                                         //attribute value
    shapeOpeningNecRegLabMapFilter->ReverseOrderingOff();                                                 //removes objects with attribute smaller than lambda
    shapeOpeningNecRegLabMapFilter->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::NUMBER_OF_PIXELS);
    shapeOpeningNecRegLabMapFilter->SetInput(imageToNecRegShaLabMapFilter->GetOutput());
    //--------------------------------------------------------------------------------------------------------------------

    necRegLabMapToImageFilter = LabelMapToLabelImageFilterType::New();
    necRegLabMapToImageFilter->SetInput(shapeOpeningNecRegLabMapFilter->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER---CREATE-LABEL-MAP-ON-DPPIV-CHANNEL-DO-SOME-STUFF-AND-CONVERT-IT-BACK-INTO-IMAGE------------------------
    necRegDMOverlayImageFilter = LabelOverlayImageFilterType::New();
    necRegDMOverlayImageFilter->SetLabelImage(necRegLabMapToImageFilter->GetOutput());
    necRegDMOverlayImageFilter->SetOpacity(m_overlayOpacity);
    necRegDMOverlayImageFilter->SetInput(readerDM->GetOutput());

    necRegDPPIVOverlayImageFilter = LabelOverlayImageFilterType::New();
    necRegDPPIVOverlayImageFilter->SetLabelImage(necRegLabMapToImageFilter->GetOutput());
    necRegDPPIVOverlayImageFilter->SetOpacity(m_overlayOpacity);
    necRegDPPIVOverlayImageFilter->SetInput(readerDPPIV->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
    necRegImageFilter = ThresholdFilterType::New();
    necRegImageFilter->ReleaseDataFlagOn();
    necRegImageFilter->SetOutsideValue(0);
    necRegImageFilter->SetInsideValue(255);
    necRegImageFilter->SetLowerThreshold(1);
    necRegImageFilter->SetUpperThreshold(255);
    necRegImageFilter->SetInput(necRegLabMapToImageFilter->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    ScalarVoWriterType::Pointer writer6 = ScalarVoWriterType::New();
    writer6->ReleaseDataFlagOn();
    writer6->SetFileName(m_pathChannelDM + m_filenameSave + m_finalSaveSuffixes[0] + m_fileExtensionChannelDM);
    writer6->SetInput(necRegImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer6->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer6->Update();

    RGBVoWriterType::Pointer writer7 = RGBVoWriterType::New();
    writer7->ReleaseDataFlagOn();
    writer7->SetFileName(m_pathChannelDM + m_filenameSave + m_finalSaveSuffixes[1] + m_fileExtensionChannelDM);
    writer7->SetInput(necRegDMOverlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer7->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer7->Update();

    RGBVoWriterType::Pointer writer8 = RGBVoWriterType::New();
    writer8->ReleaseDataFlagOn();
    writer8->SetFileName(m_pathChannelDM + m_filenameSave + m_finalSaveSuffixes[2] + m_fileExtensionChannelDM);
    writer8->SetInput(necRegDPPIVOverlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer8->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer8->Update();

    WriteLogFile(timeStamp);
}
