///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentSinusoidalNetworkOnTwoChannels20x.cpp                         //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-01-10                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "SegmentSinusoidalNetworkOnTwoChannels20x.h"

#if (ITK_VERSION_MAJOR < 4)
#define NUMBER_OF_PIXELS SIZE
#endif

#include <fstream>

#include <QDateTime>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



SegmentSinusoidalNetworkOnTwoChannels20x::SegmentSinusoidalNetworkOnTwoChannels20x()
{
	m_entryPoint = 0;

	m_overlayOpacity = 0.5;

	m_saveSuffixes[0] = "_step1_bin";
	m_saveSuffixes[1] = "_step2_hole";
	m_saveSuffixes[2] = "_step3_closopening";
	m_saveSuffixes[3] = "_step4_overlay";
	m_saveSuffixes[4] = "_step4_bin";
	m_saveSuffixes[5] = "_step4_skeleton";
}


SegmentSinusoidalNetworkOnTwoChannels20x::~SegmentSinusoidalNetworkOnTwoChannels20x()
{
}


void SegmentSinusoidalNetworkOnTwoChannels20x::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

	file.open((m_pathChannel1 + m_logFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

	file << "Begin-segment-sinusoids20x--------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
	file << "End-segment-sinusoids20x----------------------------------------------------------------------------------------------------------------\n";

    file.close();
}


void SegmentSinusoidalNetworkOnTwoChannels20x::WriteDataSetSummary()
{
    ImageAnalysisSummaryFileIO::AddEntry(SinusoidSegmentationBin, m_pathChannel1, m_pathChannel1 + m_filenameSave + m_saveSuffixes[4] + m_fileExtensionChannel1);
    ImageAnalysisSummaryFileIO::AddEntry(SinusoidSegmentationOverlay, m_pathChannel1, m_pathChannel1 + m_filenameSave + m_saveSuffixes[3] + m_fileExtensionChannel1);
    ImageAnalysisSummaryFileIO::AddEntry(SinusoidSkeleton, m_pathChannel1, m_pathChannel1 + m_filenameSave + m_saveSuffixes[5] + m_fileExtensionChannel1);
}


void SegmentSinusoidalNetworkOnTwoChannels20x::ParseParameterContext()
{
    if(m_paramContext->findContext("Segment Sinusoids + Bile Canaliculi in 20x Datasets",0)==NULL) {
        std::cout << "Error: SegmentSinusoidalNetworkOnTwoChannels20x: Invalid parameter context" << std::endl;
        return;
    }

    m_fullFilenameChannel1 = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("DPPIV channel", 0)->dataPointer()) );
    m_infoFullFilenameChannel1.setFile(m_fullFilenameChannel1);

	if(!m_infoFullFilenameChannel1.exists())
		throw std::string("Please specify DPPIV channel");

    m_fullFilenameChannel2 = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("DMs channel", 0)->dataPointer()) );
	m_infoFullFilenameChannel2.setFile(m_fullFilenameChannel2);

	if(!m_infoFullFilenameChannel2.exists())
		throw std::string("Please specify DMs channel");

	m_hasNR = *(bool*)(m_paramContext->findParameter("Is there a necrotic region", 0)->dataPointer());

    m_fullFilenameNecroticRegion = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Necrotic Region", 0)->dataPointer()) );
	m_infoFullFilenameNecroticRegion.setFile(m_fullFilenameNecroticRegion);

	if(m_hasNR && !m_infoFullFilenameNecroticRegion.exists())
		throw std::string("Please specify Necrotic Region binary mask");

	m_pathChannel1 = (m_infoFullFilenameChannel1.path() + QString("/")).toStdString();
    m_filenameChannel1 = m_infoFullFilenameChannel1.baseName().toStdString();
    m_fileExtensionChannel1 = (QString(".") + m_infoFullFilenameChannel1.suffix()).toStdString();

	m_pathChannel2 = (m_infoFullFilenameChannel2.path() + QString("/")).toStdString();
    m_filenameChannel2 = m_infoFullFilenameChannel2.baseName().toStdString();
    m_fileExtensionChannel2 = (QString(".") + m_infoFullFilenameChannel2.suffix()).toStdString();

	m_pathNecroticRegion = (m_infoFullFilenameNecroticRegion.path() + QString("/")).toStdString();
    m_filenameNecroticRegion = m_infoFullFilenameNecroticRegion.baseName().toStdString();
    m_fileExtensionNecroticRegion = (QString(".") + m_infoFullFilenameNecroticRegion.suffix()).toStdString();

    m_fullFilenameSegCV = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Central vein segmentation", 0)->dataPointer()) );
    m_infoFullFilenameSegCV.setFile(m_fullFilenameSegCV);
    m_withCVMask = m_infoFullFilenameSegCV.exists();

    m_fullFilenameSegPV = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Portal vein segmentation", 0)->dataPointer()) );
    m_infoFullFilenameSegPV.setFile(m_fullFilenameSegPV);
    m_withPVMask = m_infoFullFilenameSegPV.exists();

    std::string thresMode = ( (CSParameterChoice*)(m_paramContext->findContext("1.1a) Binary Threshold on DPPIV Channel", 0)->findParameter("DPPIV channel threshold mode", 0)->dataPointer()) )->currentString();
    if( thresMode.compare("Adaptive Otsu-Threshold")==0 )
        m_thresholdingMode[0] = 0;
    else if( thresMode.compare("Otsu-Threshold")==0 )
        m_thresholdingMode[0] = 1;
    else if( thresMode.compare("Manual Threshold")==0 )
        m_thresholdingMode[0] = 2;

    m_adapOtsuRadius[0][0] = *(int*)(m_paramContext->findContext("1.1a) Binary Threshold on DPPIV Channel", 0)->findParameter("Sample region size x", 0)->dataPointer());
    m_adapOtsuRadius[0][1] = *(int*)(m_paramContext->findContext("1.1a) Binary Threshold on DPPIV Channel", 0)->findParameter("Sample region size y", 0)->dataPointer());
    m_adapOtsuRadius[0][2] = *(int*)(m_paramContext->findContext("1.1a) Binary Threshold on DPPIV Channel", 0)->findParameter("Sample region size z", 0)->dataPointer());
    m_adapOtsuSamplePoints[0] = *(int*)(m_paramContext->findContext("1.1a) Binary Threshold on DPPIV Channel", 0)->findParameter("Number samples", 0)->dataPointer());
    m_lowerThreshold[0] = *(int*)(m_paramContext->findContext("1.1a) Binary Threshold on DPPIV Channel", 0)->findParameter("DPPIV manual threshold", 0)->dataPointer());
    //TODO max;
    m_upperThreshold[0] = 255;

    thresMode = ( (CSParameterChoice*)(m_paramContext->findContext("1.1b) Binary Threshold on DMs Channel", 0)->findParameter("DMs channel threshold mode", 0)->dataPointer()) )->currentString();
    if( thresMode.compare("Adaptive Otsu-Threshold")==0 )
        m_thresholdingMode[1] = 0;
    else if( thresMode.compare("Otsu-Threshold")==0 )
        m_thresholdingMode[1] = 1;
    else if( thresMode.compare("Manual Threshold")==0 )
        m_thresholdingMode[1] = 2;

    m_adapOtsuRadius[1][0] = *(int*)(m_paramContext->findContext("1.1b) Binary Threshold on DMs Channel", 0)->findParameter("Sample region size x", 0)->dataPointer());
    m_adapOtsuRadius[1][1] = *(int*)(m_paramContext->findContext("1.1b) Binary Threshold on DMs Channel", 0)->findParameter("Sample region size y", 0)->dataPointer());
    m_adapOtsuRadius[1][2] = *(int*)(m_paramContext->findContext("1.1b) Binary Threshold on DMs Channel", 0)->findParameter("Sample region size z", 0)->dataPointer());
    m_adapOtsuSamplePoints[1] = *(int*)(m_paramContext->findContext("1.1b) Binary Threshold on DMs Channel", 0)->findParameter("Number samples", 0)->dataPointer());
    m_lowerThreshold[1] = *(int*)(m_paramContext->findContext("1.1b) Binary Threshold on DMs Channel", 0)->findParameter("DMs manual threshold", 0)->dataPointer());
    //TODO max;
    m_upperThreshold[1] = 255;

    m_holeFillingNeighborhoodRadius[0] = *(int*)(m_paramContext->findContext("1.2) Hole Filling on 1.1", 0)->findParameter("Radius x", 0)->dataPointer());
    m_holeFillingNeighborhoodRadius[1] = *(int*)(m_paramContext->findContext("1.2) Hole Filling on 1.1", 0)->findParameter("Radius y", 0)->dataPointer());
    m_holeFillingNeighborhoodRadius[2] = *(int*)(m_paramContext->findContext("1.2) Hole Filling on 1.1", 0)->findParameter("Radius z", 0)->dataPointer());
    m_holeFillingMajThreshold = *(int*)(m_paramContext->findContext("1.2) Hole Filling on 1.1", 0)->findParameter("Majority threshold", 0)->dataPointer());

    int rad = ( *(int*)(m_paramContext->findContext("1.3) Closing on 1.2", 0)->findParameter("Kernel radius", 0)->dataPointer()) );
    m_closingStructuringElement.SetRadius(rad);
    m_closingStructuringElement.CreateStructuringElement();

    rad = ( *(int*)(m_paramContext->findContext("1.4) Opening on 1.3", 0)->findParameter("Kernel radius", 0)->dataPointer()) );
    m_openingStructuringElement.SetRadius(rad);
    m_openingStructuringElement.CreateStructuringElement();

    m_maskVeinRadius[0] = *(int*)(m_paramContext->findContext("1.5) Optional vein masking on 1.4", 0)->findParameter("Radius x", 0)->dataPointer());
    m_maskVeinRadius[1] = *(int*)(m_paramContext->findContext("1.5) Optional vein masking on 1.4", 0)->findParameter("Radius y", 0)->dataPointer());
    m_maskVeinRadius[2] = *(int*)(m_paramContext->findContext("1.5) Optional vein masking on 1.4", 0)->findParameter("Radius z", 0)->dataPointer());

    m_minimalSinusoidSize = *(unsigned int*)(m_paramContext->findContext("1.6) Remove objects on 1.5", 0)->findParameter("Smallest allowed object size", 0)->dataPointer());


    std::string saveMode = ( (CSParameterChoice*)(m_paramContext->findContext("Segment Sinusoids on DPPIV + DMs Channel", 0)->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        m_saveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        m_saveEverything = 0;
    m_logFilenameSave = *(std::string*)(m_paramContext->findContext("Segment Sinusoids on DPPIV + DMs Channel", 0)->findParameter("Log file name", 0)->dataPointer());
    m_filenameSave = *(std::string*)(m_paramContext->findContext("Segment Sinusoids on DPPIV + DMs Channel", 0)->findParameter("Save prefix", 0)->dataPointer());
}


void SegmentSinusoidalNetworkOnTwoChannels20x::Update()
{
    ParseParameterContext();

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    std::cout << "Segment sinusoids: " << std::endl;
    std::cout << " dir1: " << m_pathChannel1 << std::endl;
    std::cout << " file1: " << m_filenameChannel1 << std::endl;
    std::cout << " ext1: " << m_fileExtensionChannel1 << std::endl;

    std::cout << " dir2: " << m_pathChannel2 << std::endl;
    std::cout << " file2: " << m_filenameChannel2 << std::endl;
    std::cout << " ext2: " << m_fileExtensionChannel2 << std::endl;

    std::cout << " dir3: " << m_pathNecroticRegion << std::endl;
    std::cout << " file3: " << m_filenameNecroticRegion << std::endl;
    std::cout << " ext3: " << m_fileExtensionNecroticRegion << std::endl;


	int entry_point_counter = m_entryPoint;


	OrImageFilterType::Pointer orMaskImagesFilter;
	AdaptiveOtsuThresholdImageFilterType::Pointer adapOtsuDPPIVFilter;
	AdaptiveOtsuThresholdImageFilterType::Pointer adapOtsuDMFilter;
	OtsuThresholdImageFilterType::Pointer otsuDPPIVFilter;
	OtsuThresholdImageFilterType::Pointer otsuDMFilter;
	ThresholdFilterType::Pointer thresDPPIVFilter;
	ThresholdFilterType::Pointer thresDMFilter;
	CScalarVoImageType::Pointer binDPPIV;
	CScalarVoImageType::Pointer binDM;
	InvertIntensityImageFilterType::Pointer invertDPPIVFilter;
	InvertIntensityImageFilterType::Pointer invertDMFilter;
	AndImageFilterType::Pointer andFilter;
	ScalarVoReaderType::Pointer readerThresDPPIVDMFilter;
	HoleFillingImageFilterType::Pointer holeFillingDPPIVDMFilter;
	ScalarVoReaderType::Pointer readerHoleFillingDPPIVDMFilter;
	ClosingImageFilterType::Pointer closingDPPIVDMFilter;
	OpeningImageFilterType::Pointer openingDPPIVDMFilter;
	ScalarVoReaderType::Pointer readerDilateDPPIVDMFilter;
	SubstractImageFilterType::Pointer substractNecRegFromSinusoidsFilter;
	DilateImageFilterType::Pointer dilateVeinMaskFilter;
	SubstractImageFilterType::Pointer subtractVeinMaskFromSinusoidsFilter;
	ImageToShapeLabelMapFilterType::Pointer imageToSinusoidShaLabMapFilter;
	ShapeOpeningLabelMapFilterType::Pointer shapeOpeningSinusoidLabMapFilter;
	LabelMapToLabelImageFilterType::Pointer sinusoidLabMapToImageFilter;
	LabelOverlayImageFilterType::Pointer sinusoidOverlayImageFilter;
	ThresholdFilterType::Pointer sinusoidImageFilter;
	Thinning3DImageFilterType::Pointer thinningFilter;


	//----------READER--------------------------------------------------------------------------------------------------------
	m_spacing.Fill(1);

	if(GetNumberOfDimensions(m_pathChannel1 + m_filenameChannel1 + m_fileExtensionChannel1) != 3)
	    throw std::string("Please specify a DPPIV channel file with three dimensions.");

	CScalarVoImageType::Pointer imageDPPIV = CScalarVoImageType::New();
	ReadImage(m_pathChannel1 + m_filenameChannel1 + m_fileExtensionChannel1, imageDPPIV, m_spacing);

	CScalarVoImageType::Pointer imageNecReg = CScalarVoImageType::New();
	if(m_hasNR)
	    ReadImage(m_pathNecroticRegion + m_filenameNecroticRegion + m_fileExtensionNecroticRegion, imageNecReg, m_spacing);

	CScalarVoImageType::Pointer segCVImage = CScalarVoImageType::New();
	if(m_withCVMask)
	    ReadImage(m_fullFilenameSegCV.toStdString(), segCVImage, m_spacing);

	CScalarVoImageType::Pointer segPVImage = CScalarVoImageType::New();
	if(m_withPVMask)
	    ReadImage(m_fullFilenameSegPV.toStdString(), segPVImage, m_spacing);

	CScalarVoImageType::Pointer maskVeinImage = CScalarVoImageType::New();
	if(m_withCVMask && m_withPVMask) {
	    orMaskImagesFilter = OrImageFilterType::New();
	    orMaskImagesFilter->SetInput1(segCVImage);
	    orMaskImagesFilter->SetInput2(segPVImage);
	    orMaskImagesFilter->Update();

	    maskVeinImage = orMaskImagesFilter->GetOutput();
	    maskVeinImage->DisconnectPipeline();
	}
	else if(m_withCVMask && !m_withPVMask)
	    maskVeinImage = segCVImage;
	else if(!m_withCVMask && m_withPVMask)
	    maskVeinImage = segPVImage;

	//-------------------------------------------------------------------------------------------------------------------------

	if(entry_point_counter==0) {
		//----------READER--------------------------------------------------------------------------------------------------------
	    if(GetNumberOfDimensions(m_pathChannel2 + m_filenameChannel2 + m_fileExtensionChannel2) != 3)
	        throw std::string("Please specify a DMs channel file with three dimensions.");

	    CScalarVoImageType::Pointer imageDM = CScalarVoImageType::New();
	    ReadImage(m_pathChannel2 + m_filenameChannel2 + m_fileExtensionChannel2, imageDM, m_spacing);
		//-------------------------------------------------------------------------------------------------------------------------

		switch(m_thresholdingMode[0])
		{
		case 0:
		    {
		    adapOtsuDPPIVFilter = AdaptiveOtsuThresholdImageFilterType::New();
		    adapOtsuDPPIVFilter->SetInput(imageDPPIV);
		    adapOtsuDPPIVFilter->SetInsideValue(255);
		    adapOtsuDPPIVFilter->SetOutsideValue(0);
		    adapOtsuDPPIVFilter->SetNumberOfHistogramBins(256);
		    adapOtsuDPPIVFilter->SetSplineOrder(3);
		    adapOtsuDPPIVFilter->SetNumberOfControlPoints(10);
		    adapOtsuDPPIVFilter->SetNumberOfLevels(3);
		    adapOtsuDPPIVFilter->SetNumberOfSamples(m_adapOtsuSamplePoints[0]);
		    adapOtsuDPPIVFilter->SetRadius(m_adapOtsuRadius[0]);
		    if(m_withCVMask || m_withPVMask)
		        adapOtsuDPPIVFilter->SetMaskImage(maskVeinImage);
		    adapOtsuDPPIVFilter->Update();

		    binDPPIV = adapOtsuDPPIVFilter->GetOutput();
		    break;
		    }
		case 1:
		    {
		    otsuDPPIVFilter = OtsuThresholdImageFilterType::New();
		    otsuDPPIVFilter->SetInput(imageDPPIV);
		    otsuDPPIVFilter->Update();
		    m_otsuThreshold[0] = (CScalarPixelType)(otsuDPPIVFilter->GetThreshold());

		    invertDPPIVFilter = InvertIntensityImageFilterType::New();
		    invertDPPIVFilter->SetInput(otsuDPPIVFilter->GetOutput());
		    invertDPPIVFilter->SetMaximum(255);
		    invertDPPIVFilter->Update();

		    binDPPIV = invertDPPIVFilter->GetOutput();
		    break;
		    }
		default:
		    {
		    thresDPPIVFilter = ThresholdFilterType::New();
		    thresDPPIVFilter->ReleaseDataFlagOn();
		    thresDPPIVFilter->SetOutsideValue(0);
		    thresDPPIVFilter->SetInsideValue(255);
		    thresDPPIVFilter->SetLowerThreshold(m_lowerThreshold[0]);
		    thresDPPIVFilter->SetUpperThreshold(m_upperThreshold[0]);
		    thresDPPIVFilter->SetInput(imageDPPIV);
		    thresDPPIVFilter->Update();

		    binDPPIV = thresDPPIVFilter->GetOutput();
		    break;
		    }
		}
		binDPPIV->DisconnectPipeline();

		switch(m_thresholdingMode[1])
		{
		case 0:
		    {
		    adapOtsuDMFilter = AdaptiveOtsuThresholdImageFilterType::New();
		    adapOtsuDMFilter->SetInput(imageDM);
		    adapOtsuDMFilter->SetInsideValue(255);
		    adapOtsuDMFilter->SetOutsideValue(0);
		    adapOtsuDMFilter->SetNumberOfHistogramBins(256);
		    adapOtsuDMFilter->SetSplineOrder(3);
		    adapOtsuDMFilter->SetNumberOfControlPoints(10);
		    adapOtsuDMFilter->SetNumberOfLevels(3);
		    adapOtsuDMFilter->SetNumberOfSamples(m_adapOtsuSamplePoints[1]);
		    adapOtsuDMFilter->SetRadius(m_adapOtsuRadius[1]);
		    if(m_withCVMask || m_withPVMask)
		        adapOtsuDMFilter->SetMaskImage(maskVeinImage);
		    adapOtsuDMFilter->Update();

		    binDM = adapOtsuDMFilter->GetOutput();
		    break;
		    }
		case 1:
		    {
		    otsuDMFilter = OtsuThresholdImageFilterType::New();
		    otsuDMFilter->SetInput(imageDM);
		    otsuDMFilter->Update();
		    m_otsuThreshold[1] = (CScalarPixelType)(otsuDPPIVFilter->GetThreshold());

		    invertDMFilter = InvertIntensityImageFilterType::New();
		    invertDMFilter->SetInput(otsuDMFilter->GetOutput());
		    invertDMFilter->SetMaximum(255);
		    invertDMFilter->Update();

		    binDM = invertDMFilter->GetOutput();
		    break;
		    }
		default:
		    {
		    thresDMFilter = ThresholdFilterType::New();
		    thresDMFilter->ReleaseDataFlagOn();
		    thresDMFilter->SetOutsideValue(0);
		    thresDMFilter->SetInsideValue(255);
		    thresDMFilter->SetLowerThreshold(m_lowerThreshold[1]);
		    thresDMFilter->SetUpperThreshold(m_upperThreshold[1]);
		    thresDMFilter->SetInput(imageDM);
		    thresDMFilter->Update();

		    binDM = thresDMFilter->GetOutput();
		    break;
		    }
		}
		binDM->DisconnectPipeline();

		andFilter = AndImageFilterType::New();
		andFilter->SetInput1(binDPPIV);
		andFilter->SetInput2(binDM);

		if(m_saveEverything) {
			ScalarVoWriterType::Pointer writer1 = ScalarVoWriterType::New();
			writer1->ReleaseDataFlagOn();
			writer1->SetFileName(m_pathChannel1 + m_filenameSave + m_saveSuffixes[0] + m_fileExtensionChannel1);
			writer1->SetInput(andFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
			writer1->SetImageIO( itk::TIFFImageIO::New() );
#endif
			writer1->Update();
		}

		entry_point_counter++;
	}
	if(entry_point_counter==1) {
		//----------FILTER-HOLE-FILLING-ON-COMBINED-CHANNEL------------------------------------------------------------------------
		holeFillingDPPIVDMFilter = HoleFillingImageFilterType::New();
		holeFillingDPPIVDMFilter->ReleaseDataFlagOn();
		holeFillingDPPIVDMFilter->SetRadius(m_holeFillingNeighborhoodRadius);
		holeFillingDPPIVDMFilter->SetBackgroundValue(0);
		holeFillingDPPIVDMFilter->SetForegroundValue(255);
		holeFillingDPPIVDMFilter->SetMajorityThreshold(m_holeFillingMajThreshold);					//number of foreground neighbors should be at least (3x3x3-1)/2 + majority
		holeFillingDPPIVDMFilter->SetMaximumNumberOfIterations(20);

		if(m_entryPoint<1) {
		    holeFillingDPPIVDMFilter->SetInput(andFilter->GetOutput());
		}
		else {
			//----------READER---------------------------------------------------------------------------------------------------------
			readerThresDPPIVDMFilter = ScalarVoReaderType::New();
			readerThresDPPIVDMFilter->SetFileName(m_pathChannel2 + m_filenameChannel2 + m_fileExtensionChannel2);
			readerThresDPPIVDMFilter->ReleaseDataBeforeUpdateFlagOn();
			readerThresDPPIVDMFilter->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
			readerThresDPPIVDMFilter->SetImageIO( itk::TIFFImageIO::New() );
#endif
			//-------------------------------------------------------------------------------------------------------------------------

			holeFillingDPPIVDMFilter->SetInput(readerThresDPPIVDMFilter->GetOutput());
		}
		//-------------------------------------------------------------------------------------------------------------------------

		if(m_saveEverything) {
			ScalarVoWriterType::Pointer writer2 = ScalarVoWriterType::New();
			writer2->ReleaseDataFlagOn();
			writer2->SetFileName(m_pathChannel1 + m_filenameSave + m_saveSuffixes[1] + m_fileExtensionChannel1);
			writer2->SetInput(holeFillingDPPIVDMFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
            writer2->SetImageIO( itk::TIFFImageIO::New() );
#endif
			writer2->Update();
		}

		entry_point_counter++;
	}
	if(entry_point_counter==2) {
		//----------FILTER-ERODE-ON-DPPIV-CHANNEL----------------------------------------------------------------------------------
		closingDPPIVDMFilter = ClosingImageFilterType::New();
		closingDPPIVDMFilter->ReleaseDataFlagOn();
		closingDPPIVDMFilter->SetKernel(m_closingStructuringElement);

		if(m_entryPoint<2)
		    closingDPPIVDMFilter->SetInput(holeFillingDPPIVDMFilter->GetOutput());
		else {
			//----------READER---------------------------------------------------------------------------------------------------------
			readerHoleFillingDPPIVDMFilter = ScalarVoReaderType::New();
			readerHoleFillingDPPIVDMFilter->SetFileName(m_pathChannel2 + m_filenameChannel2 + m_fileExtensionChannel2);
			readerHoleFillingDPPIVDMFilter->ReleaseDataBeforeUpdateFlagOn();
			readerHoleFillingDPPIVDMFilter->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
			readerHoleFillingDPPIVDMFilter->SetImageIO( itk::TIFFImageIO::New() );
#endif
			//-------------------------------------------------------------------------------------------------------------------------

			closingDPPIVDMFilter->SetInput(readerHoleFillingDPPIVDMFilter->GetOutput());
		}
		//-------------------------------------------------------------------------------------------------------------------------

		//----------FILTER-DILATE-ON-DPPIV-CHANNEL---------------------------------------------------------------------------------
		openingDPPIVDMFilter = OpeningImageFilterType::New();
		openingDPPIVDMFilter->ReleaseDataFlagOn();
		openingDPPIVDMFilter->SetKernel(m_openingStructuringElement);
		openingDPPIVDMFilter->SetInput(closingDPPIVDMFilter->GetOutput());
		//-------------------------------------------------------------------------------------------------------------------------

		if(m_saveEverything) {
			ScalarVoWriterType::Pointer writer4 = ScalarVoWriterType::New();
			writer4->ReleaseDataFlagOn();
			writer4->SetFileName(m_pathChannel1 + m_filenameSave + m_saveSuffixes[2] + m_fileExtensionChannel1);
			writer4->SetInput(openingDPPIVDMFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
            writer4->SetImageIO( itk::TIFFImageIO::New() );
#endif
			writer4->Update();
		}

		entry_point_counter++;
	}
	if(entry_point_counter==3) {

	    if(m_hasNR) {
	        //----------FILTER-SUBSTRACT-NECREGION-FROM-DPPIV-CHANNEL------------------------------------------------------------------
	        substractNecRegFromSinusoidsFilter = SubstractImageFilterType::New();
	        substractNecRegFromSinusoidsFilter->ReleaseDataFlagOn();

	        if(m_entryPoint<3)
	            substractNecRegFromSinusoidsFilter->SetInput1(openingDPPIVDMFilter->GetOutput());
	        else {
	            //----------READER---------------------------------------------------------------------------------------------------------
	            readerDilateDPPIVDMFilter = ScalarVoReaderType::New();
	            readerDilateDPPIVDMFilter->SetFileName(m_pathChannel2 + m_filenameChannel2 + m_fileExtensionChannel2);
	            readerDilateDPPIVDMFilter->ReleaseDataBeforeUpdateFlagOn();
	            readerDilateDPPIVDMFilter->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
	            readerDilateDPPIVDMFilter->SetImageIO( itk::TIFFImageIO::New() );
#endif
	            //-------------------------------------------------------------------------------------------------------------------------
	            substractNecRegFromSinusoidsFilter->SetInput1(readerDilateDPPIVDMFilter->GetOutput());
	        }
	        substractNecRegFromSinusoidsFilter->SetInput2(imageNecReg);
	        //--------------------------------------------------------------------------------------------------------------------
	    }
	    if(m_withCVMask || m_withPVMask) {
	        StructuringElementType dilateKernel;
	        dilateKernel.SetRadius(m_maskVeinRadius);
	        dilateKernel.CreateStructuringElement();

	        dilateVeinMaskFilter = DilateImageFilterType::New();
	        dilateVeinMaskFilter->SetKernel(dilateKernel);
	        dilateVeinMaskFilter->SetInput(maskVeinImage);

	        subtractVeinMaskFromSinusoidsFilter = SubstractImageFilterType::New();
	        subtractVeinMaskFromSinusoidsFilter->ReleaseDataFlagOn();
	        if(m_hasNR) subtractVeinMaskFromSinusoidsFilter->SetInput1(substractNecRegFromSinusoidsFilter->GetOutput());   //entry point stuff not supported, but not used atm anyway
	        else        subtractVeinMaskFromSinusoidsFilter->SetInput1(openingDPPIVDMFilter->GetOutput());
	        subtractVeinMaskFromSinusoidsFilter->SetInput2(dilateVeinMaskFilter->GetOutput());
	    }

		//----------FILTER---CREATE-LABEL-MAP-ON-DPPIV-CHANNEL-DO-SOME-STUFF-AND-CONVERT-IT-BACK-INTO-IMAGE------------------------
		imageToSinusoidShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
		if(m_withCVMask || m_withPVMask)    imageToSinusoidShaLabMapFilter->SetInput(subtractVeinMaskFromSinusoidsFilter->GetOutput());
		else if(m_hasNR)                    imageToSinusoidShaLabMapFilter->SetInput(substractNecRegFromSinusoidsFilter->GetOutput());
		else {
		    if(m_entryPoint<3)
		        imageToSinusoidShaLabMapFilter->SetInput(openingDPPIVDMFilter->GetOutput());
		    else {
		        //----------READER---------------------------------------------------------------------------------------------------------
		        readerDilateDPPIVDMFilter = ScalarVoReaderType::New();
		        readerDilateDPPIVDMFilter->SetFileName(m_pathChannel2 + m_filenameChannel2 + m_fileExtensionChannel2);
		        readerDilateDPPIVDMFilter->ReleaseDataBeforeUpdateFlagOn();
		        readerDilateDPPIVDMFilter->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
		        readerDilateDPPIVDMFilter->SetImageIO( itk::TIFFImageIO::New() );
#endif
		        //-------------------------------------------------------------------------------------------------------------------------
		        imageToSinusoidShaLabMapFilter->SetInput(readerDilateDPPIVDMFilter->GetOutput());
		    }
		    //--------------------------------------------------------------------------------------------------------------------
		}
		imageToSinusoidShaLabMapFilter->Update();

		//----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------
		shapeOpeningSinusoidLabMapFilter = ShapeOpeningLabelMapFilterType::New();
		shapeOpeningSinusoidLabMapFilter->SetLambda(m_minimalSinusoidSize);							//attribute value
		shapeOpeningSinusoidLabMapFilter->ReverseOrderingOff();										//removes objects with attribute smaller than lambda
		shapeOpeningSinusoidLabMapFilter->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::NUMBER_OF_PIXELS);
		shapeOpeningSinusoidLabMapFilter->SetInput(imageToSinusoidShaLabMapFilter->GetOutput());
		shapeOpeningSinusoidLabMapFilter->Update();

		sinusoidLabMapToImageFilter = LabelMapToLabelImageFilterType::New();
		sinusoidLabMapToImageFilter->SetInput(shapeOpeningSinusoidLabMapFilter->GetOutput());
		//-------------------------------------------------------------------------------------------------------------------------

		//----------VISULALIZATION-ONLY-FILTER---NOT-NEEDED-BY-OTHER-PIPELINE-FILTERS----------------------------------------------
		sinusoidOverlayImageFilter = LabelOverlayImageFilterType::New();
		sinusoidOverlayImageFilter->SetLabelImage(sinusoidLabMapToImageFilter->GetOutput());
		sinusoidOverlayImageFilter->SetOpacity(m_overlayOpacity);
		sinusoidOverlayImageFilter->ReleaseDataFlagOn();
		sinusoidOverlayImageFilter->SetInput(imageDPPIV);
		sinusoidOverlayImageFilter->Update();
		//-------------------------------------------------------------------------------------------------------------------------

		//----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
		sinusoidImageFilter = ThresholdFilterType::New();
		sinusoidImageFilter->SetOutsideValue(0);
		sinusoidImageFilter->SetInsideValue(255);
		sinusoidImageFilter->SetLowerThreshold(1);
		sinusoidImageFilter->SetUpperThreshold(255);
		sinusoidImageFilter->SetInput(sinusoidLabMapToImageFilter->GetOutput());
		//--------------------------------------------------------------------------------------------------------------------------

		//----------THINNING--------------------------------------------------------------------------------------------------------
		thinningFilter = Thinning3DImageFilterType::New();
		thinningFilter->ReleaseDataFlagOn();
		thinningFilter->SetInput(sinusoidImageFilter->GetOutput());
		//--------------------------------------------------------------------------------------------------------------------------

		//----------WRITER----------------------------------------------------------------------------------------------------------
		RGBVoWriterType::Pointer writer6 = RGBVoWriterType::New();
		writer6->ReleaseDataFlagOn();
		writer6->SetFileName(m_pathChannel1 + m_filenameSave + m_saveSuffixes[3] + m_fileExtensionChannel1);
		writer6->SetInput(sinusoidOverlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
		writer6->SetImageIO( itk::TIFFImageIO::New() );
#endif
		writer6->Update();

		ScalarVoWriterType::Pointer writer7 = ScalarVoWriterType::New();
		writer7->ReleaseDataFlagOn();
		writer7->SetFileName(m_pathChannel1 + m_filenameSave + m_saveSuffixes[4] + m_fileExtensionChannel1);
		writer7->SetInput(sinusoidImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
		writer7->SetImageIO( itk::TIFFImageIO::New() );
#endif
		writer7->Update();

		ScalarVoWriterType::Pointer writer8 = ScalarVoWriterType::New();
		writer8->ReleaseDataFlagOn();
		writer8->SetFileName(m_pathChannel1 + m_filenameSave + m_saveSuffixes[5] + m_fileExtensionChannel1);
		writer8->SetInput(thinningFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
		writer8->SetImageIO( itk::TIFFImageIO::New() );
#endif
		writer8->Update();

		ScalarVoWriterType::Pointer writerTemp = ScalarVoWriterType::New();
		writerTemp->ReleaseDataFlagOn();
		writerTemp->SetFileName(m_pathChannel1 + "_delete_me_temp" + m_fileExtensionChannel1);
		writerTemp->SetInput(sinusoidImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
		writerTemp->SetImageIO( itk::TIFFImageIO::New() );
#endif
		writerTemp->Update();
	}

	WriteLogFile(timeStamp);
	WriteDataSetSummary();
}
