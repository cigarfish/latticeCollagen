///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentBileNetworkOnDPPIV20x.h				                         //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-01-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "SegmentBileNetworkOnDPPIV20x.h"

#include <fstream>

#if (ITK_VERSION_MAJOR < 4)
#define NUMBER_OF_PIXELS SIZE
#endif

#include <QDateTime>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



SegmentBileNetworkOnDPPIV20x::SegmentBileNetworkOnDPPIV20x()
{
	m_entryPoint = 0;

	m_overlayOpacity = 0.5;

	m_saveSuffixes[0] = "_step1_bin";
	m_saveSuffixes[1] = "_step2_hole";
	m_saveSuffixes[2] = "_step3_subtract";
	m_saveSuffixes[3] = "_step4_invHole";
	m_saveSuffixes[4] = "_step5_overlay";
	m_saveSuffixes[5] = "_step5_bin";
	m_saveSuffixes[6] = "_step5_skeleton";
}


SegmentBileNetworkOnDPPIV20x::~SegmentBileNetworkOnDPPIV20x()
{
	// TODO Auto-generated destructor stub
}


void SegmentBileNetworkOnDPPIV20x::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

	file.open((m_pathChannelDPPIV + m_logFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

	file << "Begin-segment-bile20x-------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
	file << "End-segment-bile20x---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


void SegmentBileNetworkOnDPPIV20x::WriteDataSetSummary()
{
    ImageAnalysisSummaryFileIO::AddEntry(BileSegmentationBin, m_pathChannelDPPIV, m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[5] + m_fileExtensionChannelDPPIV);
    ImageAnalysisSummaryFileIO::AddEntry(BileSegmentationOverlay, m_pathChannelDPPIV, m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[4] + m_fileExtensionChannelDPPIV);
    ImageAnalysisSummaryFileIO::AddEntry(BileSkeleton, m_pathChannelDPPIV, m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[6] + m_fileExtensionChannelDPPIV);
}


void SegmentBileNetworkOnDPPIV20x::ParseParameterContext()
{
    if(m_paramContext->findContext("Segment Sinusoids + Bile Canaliculi in 20x Datasets",0)==NULL) {
        std::cout << "Error: SegmentBileCanaliculiNetworkOnDPPIVChannel20x: Invalid parameter context" << std::endl;
        return;
    }

    m_fullFilenameChannelDPPIV = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("DPPIV channel", 0)->dataPointer()) );
	m_infoFullFilenameChannelDPPIV.setFile(m_fullFilenameChannelDPPIV);

	if(!m_infoFullFilenameChannelDPPIV.exists())
		throw std::string("Please specify DPPIV channel");

    m_hasNR = *(bool*)(m_paramContext->findParameter("Is there a necrotic region", 0)->dataPointer());

    m_fullFilenameNecroticRegion = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Necrotic Region", 0)->dataPointer()) );
	m_infoFullFilenameNecroticRegion.setFile(m_fullFilenameNecroticRegion);

	if(m_hasNR && !m_infoFullFilenameNecroticRegion.exists())
		throw std::string("Please specify Necrotic Region binary mask");

	m_pathChannelDPPIV = (m_infoFullFilenameChannelDPPIV.path() + QString("/")).toStdString();
    m_filenameChannelDPPIV = m_infoFullFilenameChannelDPPIV.baseName().toStdString();
    m_fileExtensionChannelDPPIV = (QString(".") + m_infoFullFilenameChannelDPPIV.suffix()).toStdString();

	m_pathNecroticRegion = (m_infoFullFilenameNecroticRegion.path() + QString("/")).toStdString();
    m_filenameNecroticRegion = m_infoFullFilenameNecroticRegion.baseName().toStdString();
    m_fileExtensionNecroticRegion = (QString(".") + m_infoFullFilenameNecroticRegion.suffix()).toStdString();
    
    m_fullFilenameSegCV = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Central vein segmentation", 0)->dataPointer()) );
    m_infoFullFilenameSegCV.setFile(m_fullFilenameSegCV);
    m_withCVMask = m_infoFullFilenameSegCV.exists();

    m_fullFilenameSegPV = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Portal vein segmentation", 0)->dataPointer()) );
    m_infoFullFilenameSegPV.setFile(m_fullFilenameSegPV);
    m_withPVMask = m_infoFullFilenameSegPV.exists();


    std::string thresMode = ( (CSParameterChoice*)(m_paramContext->findContext("2.1) Binary Threshold Filter on DPPIV Channel", 0)->findParameter("DPPIV channel threshold mode", 0)->dataPointer()) )->currentString();
    if( thresMode.compare("Adaptive Otsu-Threshold")==0 )
        m_thresholdingMode = 0;
    else if( thresMode.compare("Otsu-Threshold")==0 )
        m_thresholdingMode = 1;
    else if( thresMode.compare("Manual Threshold")==0 )
        m_thresholdingMode = 2;

    m_adapOtsuRadius[0] = *(int*)(m_paramContext->findContext("2.1) Binary Threshold Filter on DPPIV Channel", 0)->findParameter("Sample region size x", 0)->dataPointer());
    m_adapOtsuRadius[1] = *(int*)(m_paramContext->findContext("2.1) Binary Threshold Filter on DPPIV Channel", 0)->findParameter("Sample region size y", 0)->dataPointer());
    m_adapOtsuRadius[2] = *(int*)(m_paramContext->findContext("2.1) Binary Threshold Filter on DPPIV Channel", 0)->findParameter("Sample region size z", 0)->dataPointer());
    m_adapOtsuSamplePoints = *(int*)(m_paramContext->findContext("2.1) Binary Threshold Filter on DPPIV Channel", 0)->findParameter("Number samples", 0)->dataPointer());
    m_lowerThreshold = *(int*)(m_paramContext->findContext("2.1) Binary Threshold Filter on DPPIV Channel", 0)->findParameter("DPPIV manual threshold", 0)->dataPointer());
    //TODO max;
    m_upperThreshold = 255;

    m_holeFilling1NeighborhoodRadius[0] = *(int*)(m_paramContext->findContext("2.2) Hole Filling on 2.1", 0)->findParameter("Radius x", 0)->dataPointer());
    m_holeFilling1NeighborhoodRadius[1] = *(int*)(m_paramContext->findContext("2.2) Hole Filling on 2.1", 0)->findParameter("Radius y", 0)->dataPointer());
    m_holeFilling1NeighborhoodRadius[2] = *(int*)(m_paramContext->findContext("2.2) Hole Filling on 2.1", 0)->findParameter("Radius z", 0)->dataPointer());
    m_holeFilling1MajThreshold = *(int*)(m_paramContext->findContext("2.2) Hole Filling on 2.1", 0)->findParameter("Majority threshold", 0)->dataPointer());

    m_inverseHoleFillingNeighborhoodRadius[0] = *(int*)(m_paramContext->findContext("2.3) Inverse Hole Filling on 2.2", 0)->findParameter("Radius x", 0)->dataPointer());
    m_inverseHoleFillingNeighborhoodRadius[1] = *(int*)(m_paramContext->findContext("2.3) Inverse Hole Filling on 2.2", 0)->findParameter("Radius y", 0)->dataPointer());
    m_inverseHoleFillingNeighborhoodRadius[2] = *(int*)(m_paramContext->findContext("2.3) Inverse Hole Filling on 2.2", 0)->findParameter("Radius z", 0)->dataPointer());
    m_inverseHoleFillingMajThreshold = *(int*)(m_paramContext->findContext("2.3) Inverse Hole Filling on 2.2", 0)->findParameter("Majority threshold", 0)->dataPointer());

    m_maskVeinRadius[0] = *(int*)(m_paramContext->findContext("2.4) Optional vein masking on 2.3", 0)->findParameter("Radius x", 0)->dataPointer());
    m_maskVeinRadius[1] = *(int*)(m_paramContext->findContext("2.4) Optional vein masking on 2.3", 0)->findParameter("Radius y", 0)->dataPointer());
    m_maskVeinRadius[2] = *(int*)(m_paramContext->findContext("2.4) Optional vein masking on 2.3", 0)->findParameter("Radius z", 0)->dataPointer());

    m_minimalBileSize = *(unsigned int*)(m_paramContext->findContext("2.5) Remove objects on 2.4", 0)->findParameter("Smallest allowed object size", 0)->dataPointer());


    std::string saveMode = ( (CSParameterChoice*)(m_paramContext->findContext("Segment Bile Canaliculi on DPPIV Channel", 0)->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        m_saveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        m_saveEverything = 0;
    m_logFilenameSave = *(std::string*)(m_paramContext->findContext("Segment Bile Canaliculi on DPPIV Channel", 0)->findParameter("Log file name", 0)->dataPointer());
    m_filenameSave = *(std::string*)(m_paramContext->findContext("Segment Bile Canaliculi on DPPIV Channel", 0)->findParameter("Save prefix", 0)->dataPointer());
}


void SegmentBileNetworkOnDPPIV20x::Update()
{
    ParseParameterContext();
	
	std::string timeStamp;

	QDateTime time = QDateTime::currentDateTime();
	timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();


	std::cout << "Segment bile network: " << std::endl;
	std::cout << " dir: " << m_pathChannelDPPIV << std::endl;
	std::cout << " file: " << m_filenameChannelDPPIV << std::endl;
	std::cout << " ext: " << m_fileExtensionChannelDPPIV << std::endl;

	std::cout << " dir2: " << m_pathNecroticRegion << std::endl;
	std::cout << " file2: " << m_filenameNecroticRegion << std::endl;
	std::cout << " ext2: " << m_fileExtensionNecroticRegion << std::endl;

	int entry_point_counter = m_entryPoint;

	OrImageFilterType::Pointer orMaskImagesFilter;
	AdaptiveOtsuThresholdImageFilterType::Pointer adapOtsuFilter;
	OtsuThresholdImageFilterType::Pointer otsuFilter;
	InvertIntensityImageFilterType::Pointer invertFilter;
	ThresholdFilterType::Pointer thresDPPIVFilter;
	ScalarVoReaderType::Pointer thresDPPIVReader;
	HoleFillingImageFilterType::Pointer holeFillingDPPIVFilter;
	ScalarVoReaderType::Pointer holeFillingDPPIVReader;
	SubtractImageFilterType::Pointer subtractSinusoidsFromDPPIVFilter;
	ScalarVoReaderType::Pointer subtractSinusoidsFromDPPIVReader;
	HoleFillingImageFilterType::Pointer holeFillingInverseBileFilter;
	ScalarVoReaderType::Pointer holeFillingInverseBileReader;
	HoleFillingImageFilterType::Pointer holeFillingBileFilter;
	ScalarVoReaderType::Pointer holeFillingBileReader;
	SubtractImageFilterType::Pointer subtractNecRegFromBileFilter;
	SubtractImageFilterType::Pointer subtractVeinMaskFromBileFilter;
	BinaryDilateFilterType::Pointer dilateVeinMaskFilter;
	ImageToShapeLabelMapFilterType::Pointer imageToBileShaLabMapFilter;
	ShapeOpeningLabelMapFilterType::Pointer shapeOpeningBileLabMapFilter;
	LabelMapToLabelImageFilterType::Pointer bileLabMapToImageFilter;
	LabelOverlayImageFilterType::Pointer bileOverlayImageFilter;
	ThresholdFilterType::Pointer bileImageFilter;
	ScalarVoReaderType::Pointer bileImageReader;
	Thinning3DImageFilterType::Pointer thinningBileFilter;


	//----------READER---------------------------------------------------------------------------------------------------------
	m_spacing.Fill(1);

	if(GetNumberOfDimensions(m_pathChannelDPPIV + m_filenameChannelDPPIV + m_fileExtensionChannelDPPIV) != 3)
	    throw std::string("Please specify a DPPIV channel file with three dimensions.");

	CScalarVoImageType::Pointer imageDPPIV = CScalarVoImageType::New();
	ReadImage(m_pathChannelDPPIV + m_filenameChannelDPPIV + m_fileExtensionChannelDPPIV, imageDPPIV, m_spacing);

	CScalarVoImageType::Pointer imageSin = CScalarVoImageType::New();
	ReadImage(m_pathChannelDPPIV + "_delete_me_temp" + m_fileExtensionChannelDPPIV, imageSin, m_spacing);

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
		//----------FILTER-THRESHOLD-ON-DPPIV-CHANNEL------------------------------------------------------------------------------
	    switch(m_thresholdingMode)
	    {
	    case 0:
	        {
	        adapOtsuFilter = AdaptiveOtsuThresholdImageFilterType::New();
	        adapOtsuFilter->SetInput(imageDPPIV);
	        adapOtsuFilter->SetInsideValue(255);
	        adapOtsuFilter->SetOutsideValue(0);
	        adapOtsuFilter->SetNumberOfHistogramBins(256);
	        adapOtsuFilter->SetSplineOrder(3);
	        adapOtsuFilter->SetNumberOfControlPoints(10);
	        adapOtsuFilter->SetNumberOfLevels(3);
	        adapOtsuFilter->SetNumberOfSamples(m_adapOtsuSamplePoints);
	        adapOtsuFilter->SetRadius(m_adapOtsuRadius);
	        if(m_withCVMask || m_withPVMask)
	            adapOtsuFilter->SetMaskImage(maskVeinImage);
	        break;
	        }
	    case 1:
	        {
	        otsuFilter = OtsuThresholdImageFilterType::New();
	        otsuFilter->SetInput(imageDPPIV);
	        otsuFilter->Update();
	        m_otsuThreshold = (CScalarPixelType)(otsuFilter->GetThreshold());

	        invertFilter = InvertIntensityImageFilterType::New();
	        invertFilter->SetInput(otsuFilter->GetOutput());
	        invertFilter->SetMaximum(255);
	        break;
	        }
	    default:
	        {
	        thresDPPIVFilter = ThresholdFilterType::New();
	        thresDPPIVFilter->ReleaseDataFlagOn();
	        thresDPPIVFilter->SetOutsideValue(0);
	        thresDPPIVFilter->SetInsideValue(255);
	        thresDPPIVFilter->SetLowerThreshold(m_lowerThreshold);
	        thresDPPIVFilter->SetUpperThreshold(m_upperThreshold);
	        thresDPPIVFilter->SetInput(imageDPPIV);
	        break;
	        }
	    }

		if(m_saveEverything) {
			ScalarVoWriterType::Pointer writer1 = ScalarVoWriterType::New();
			writer1->ReleaseDataFlagOn();
			writer1->SetFileName(m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[0] + m_fileExtensionChannelDPPIV);
			switch(m_thresholdingMode)
			{
			case 0:
			    writer1->SetInput(adapOtsuFilter->GetOutput());
			    break;
			case 1:
			    writer1->SetInput(invertFilter->GetOutput());
			    break;
			default:
			    writer1->SetInput(thresDPPIVFilter->GetOutput());
			    break;
			}
#if (ITK_VERSION_MAJOR >= 4)
			writer1->SetImageIO( itk::TIFFImageIO::New() );
#endif
			writer1->Update();
		}

		entry_point_counter++;
	}
	if(entry_point_counter==1) {
		//----------FILTER-HOLE-FILLING-ON-DPPIV-CHANNEL---------------------------------------------------------------------------
		holeFillingDPPIVFilter = HoleFillingImageFilterType::New();
		holeFillingDPPIVFilter->ReleaseDataFlagOn();
		holeFillingDPPIVFilter->SetRadius(m_holeFilling1NeighborhoodRadius);
		holeFillingDPPIVFilter->SetBackgroundValue(0);
		holeFillingDPPIVFilter->SetForegroundValue(255);
		holeFillingDPPIVFilter->SetMajorityThreshold(m_holeFilling1MajThreshold);					//number of foreground neighbors should be at least (3x3x3-1)/2 + majority

		if(m_entryPoint<1) {
	        switch(m_thresholdingMode)
	        {
	        case 0:
	            holeFillingDPPIVFilter->SetInput(adapOtsuFilter->GetOutput());
	            break;
	        case 1:
	            holeFillingDPPIVFilter->SetInput(invertFilter->GetOutput());
	            break;
	        default:
	            holeFillingDPPIVFilter->SetInput(thresDPPIVFilter->GetOutput());
	            break;
	        }
		}
		else {
			thresDPPIVReader = ScalarVoReaderType::New();
			thresDPPIVReader->SetFileName(m_pathChannelLastStep + m_filenameChannelLastStep + m_fileExtensionChannelLastStep);
			thresDPPIVReader->ReleaseDataBeforeUpdateFlagOn();
			thresDPPIVReader->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
			thresDPPIVReader->SetImageIO( itk::TIFFImageIO::New() );
#endif

			holeFillingDPPIVFilter->SetInput(thresDPPIVReader->GetOutput());
		}
		//-------------------------------------------------------------------------------------------------------------------------

		if(m_saveEverything) {
			ScalarVoWriterType::Pointer writer2 = ScalarVoWriterType::New();
			writer2->ReleaseDataFlagOn();
			writer2->SetFileName(m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[1] + m_fileExtensionChannelDPPIV);
			writer2->SetInput(holeFillingDPPIVFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
			writer2->SetImageIO( itk::TIFFImageIO::New() );
#endif
			writer2->Update();
		}

		entry_point_counter++;
	}
	if(entry_point_counter==2) {
		//----------FILTER-SUBSTRACT-SINUSOIDS-FROM-DPPIV-CHANNEL------------------------------------------------------------------
		subtractSinusoidsFromDPPIVFilter = SubtractImageFilterType::New();
		subtractSinusoidsFromDPPIVFilter->ReleaseDataFlagOn();
		subtractSinusoidsFromDPPIVFilter->SetInput2(imageSin);

		if(m_entryPoint<2)
			subtractSinusoidsFromDPPIVFilter->SetInput1(holeFillingDPPIVFilter->GetOutput());
		else {
			holeFillingDPPIVReader = ScalarVoReaderType::New();
			holeFillingDPPIVReader->SetFileName(m_pathChannelLastStep + m_filenameChannelLastStep + m_fileExtensionChannelLastStep);
			holeFillingDPPIVReader->ReleaseDataBeforeUpdateFlagOn();
			holeFillingDPPIVReader->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
			holeFillingDPPIVReader->SetImageIO( itk::TIFFImageIO::New() );
#endif

			subtractSinusoidsFromDPPIVFilter->SetInput1(holeFillingDPPIVReader->GetOutput());
		}
		//-------------------------------------------------------------------------------------------------------------------------

		if(m_saveEverything) {
			ScalarVoWriterType::Pointer writer3 = ScalarVoWriterType::New();
			writer3->ReleaseDataFlagOn();
			writer3->SetFileName(m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[2] + m_fileExtensionChannelDPPIV);
			writer3->SetInput(subtractSinusoidsFromDPPIVFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
            writer3->SetImageIO( itk::TIFFImageIO::New() );
#endif
			writer3->Update();
		}

		entry_point_counter++;
	}
	if(entry_point_counter==3) {
		//----------FILTER-HOLE-FILLING-ON-INVERSE-TO-REMOVE-NOISE-FROM-SUBSTRACTION-----------------------------------------------
		holeFillingInverseBileFilter = HoleFillingImageFilterType::New();
		holeFillingInverseBileFilter->ReleaseDataFlagOn();
		holeFillingInverseBileFilter->SetRadius(m_inverseHoleFillingNeighborhoodRadius);
		holeFillingInverseBileFilter->SetBackgroundValue(255);
		holeFillingInverseBileFilter->SetForegroundValue(0);
		holeFillingInverseBileFilter->SetMajorityThreshold(m_inverseHoleFillingMajThreshold);					//number of foreground neighbors should be at least (3x3x3-1)/2 + majority

		if(m_entryPoint<3)
			holeFillingInverseBileFilter->SetInput(subtractSinusoidsFromDPPIVFilter->GetOutput());
		else {
			subtractSinusoidsFromDPPIVReader = ScalarVoReaderType::New();
			subtractSinusoidsFromDPPIVReader->SetFileName(m_pathChannelLastStep + m_filenameChannelLastStep + m_fileExtensionChannelLastStep);
			subtractSinusoidsFromDPPIVReader->ReleaseDataBeforeUpdateFlagOn();
			subtractSinusoidsFromDPPIVReader->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
			subtractSinusoidsFromDPPIVReader->SetImageIO( itk::TIFFImageIO::New() );
#endif

			holeFillingInverseBileFilter->SetInput(subtractSinusoidsFromDPPIVReader->GetOutput());
		}
		//-------------------------------------------------------------------------------------------------------------------------

		if(m_saveEverything) {
			ScalarVoWriterType::Pointer writer4 = ScalarVoWriterType::New();
			writer4->ReleaseDataFlagOn();
			writer4->SetFileName(m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[3] + m_fileExtensionChannelDPPIV);
			writer4->SetInput(holeFillingInverseBileFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
            writer4->SetImageIO( itk::TIFFImageIO::New() );
#endif
			writer4->Update();
		}

		entry_point_counter++;
	}
//	if(entry_point_counter==4) {
//		//----------FILTER-HOLE-FILLING-ON-DPPIV-WITHOUT-SINUSOIDS-----------------------------------------------------------------
//		holeFillingBileFilter = HoleFillingImageFilterType::New();
//		holeFillingBileFilter->ReleaseDataFlagOn();
//		holeFillingBileFilter->SetRadius(m_holeFilling2NeighborhoodRadius);
//		holeFillingBileFilter->SetBackgroundValue(0);
//		holeFillingBileFilter->SetForegroundValue(255);
//		holeFillingBileFilter->SetMajorityThreshold(m_holeFilling2MajThreshold);							//number of foreground neighbors should be at least (3x3x3-1)/2 + majority
//
//		if(m_entryPoint<4)
//			holeFillingBileFilter->SetInput(holeFillingInverseBileFilter->GetOutput());
//		else {
//			holeFillingInverseBileReader = ScalarVoReaderType::New();
//			holeFillingInverseBileReader->SetFileName(m_pathChannelLastStep + m_filenameChannelLastStep + m_fileExtensionChannelLastStep);
//			holeFillingInverseBileReader->ReleaseDataBeforeUpdateFlagOn();
//			holeFillingInverseBileReader->UseStreamingOn();
//#if (ITK_VERSION_MAJOR >= 4)
//			holeFillingInverseBileReader->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//
//			holeFillingBileFilter->SetInput(holeFillingInverseBileReader->GetOutput());
//		}
//		//-------------------------------------------------------------------------------------------------------------------------
//
//		if(m_saveEverything) {
//			ScalarVoWriterType::Pointer writer5 = ScalarVoWriterType::New();
//			writer5->ReleaseDataFlagOn();
//			writer5->SetFileName(m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[4] + m_fileExtensionChannelDPPIV);
//			writer5->SetInput(holeFillingBileFilter->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//            writer5->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//			writer5->Update();
//		}
//
//		entry_point_counter++;
//	}
	if(entry_point_counter==4) {

	    if(m_hasNR) {
	        //----------FILTER-SUBSTRACT-SINUSOIDS-FROM-DPPIV-CHANNEL------------------------------------------------------------------
	        subtractNecRegFromBileFilter = SubtractImageFilterType::New();
	        subtractNecRegFromBileFilter->ReleaseDataFlagOn();

	        if(m_entryPoint<4)
	            subtractNecRegFromBileFilter->SetInput1(holeFillingInverseBileFilter->GetOutput());
//	            subtractNecRegFromBileFilter->SetInput1(holeFillingBileFilter->GetOutput());
	        else {
	            holeFillingBileReader = ScalarVoReaderType::New();
	            holeFillingBileReader->SetFileName(m_pathChannelLastStep + m_filenameChannelLastStep + m_fileExtensionChannelLastStep);
	            holeFillingBileReader->ReleaseDataBeforeUpdateFlagOn();
	            holeFillingBileReader->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
	            holeFillingBileReader->SetImageIO( itk::TIFFImageIO::New() );
#endif

	            subtractNecRegFromBileFilter->SetInput1(holeFillingBileReader->GetOutput());
	        }
	        subtractNecRegFromBileFilter->SetInput2(imageNecReg);
	        //--------------------------------------------------------------------------------------------------------------------
	    }
	    if(m_withCVMask || m_withPVMask) {
	        StructuringElementType dilateKernel;
	        dilateKernel.SetRadius(m_maskVeinRadius);
	        dilateKernel.CreateStructuringElement();

	        dilateVeinMaskFilter = BinaryDilateFilterType::New();
	        dilateVeinMaskFilter->SetKernel(dilateKernel);
	        dilateVeinMaskFilter->SetInput(maskVeinImage);

	        subtractVeinMaskFromBileFilter = SubtractImageFilterType::New();
	        subtractVeinMaskFromBileFilter->ReleaseDataFlagOn();
	        if(m_hasNR) subtractVeinMaskFromBileFilter->SetInput1(subtractNecRegFromBileFilter->GetOutput());   //entry point stuff not supported, but not used atm anyway
	        else        subtractVeinMaskFromBileFilter->SetInput1(holeFillingInverseBileFilter->GetOutput());
	        subtractVeinMaskFromBileFilter->SetInput2(dilateVeinMaskFilter->GetOutput());
	    }

		//----------FILTER---CREATE-LABEL-MAP-ON-DPPIV-CHANNEL-DO-SOME-STUFF-AND-CONVERT-IT-BACK-INTO-IMAGE------------------------
		imageToBileShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
		if(m_withCVMask || m_withPVMask)    imageToBileShaLabMapFilter->SetInput(subtractVeinMaskFromBileFilter->GetOutput());
		else if(m_hasNR)                    imageToBileShaLabMapFilter->SetInput(subtractNecRegFromBileFilter->GetOutput());
		else {
		    if(m_entryPoint<4)
		        imageToBileShaLabMapFilter->SetInput(holeFillingInverseBileFilter->GetOutput());
//		        imageToBileShaLabMapFilter->SetInput(holeFillingBileFilter->GetOutput());
		    else {
		        holeFillingBileReader = ScalarVoReaderType::New();
		        holeFillingBileReader->SetFileName(m_pathChannelLastStep + m_filenameChannelLastStep + m_fileExtensionChannelLastStep);
		        holeFillingBileReader->ReleaseDataBeforeUpdateFlagOn();
		        holeFillingBileReader->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
		        holeFillingBileReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
		        imageToBileShaLabMapFilter->SetInput(holeFillingBileReader->GetOutput());
		    }
		}
		imageToBileShaLabMapFilter->Update();
		//--------------------------------------------------------------------------------------------------------------------

		//----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------
		shapeOpeningBileLabMapFilter = ShapeOpeningLabelMapFilterType::New();
		shapeOpeningBileLabMapFilter->SetLambda(m_minimalBileSize);											//attribute value
		shapeOpeningBileLabMapFilter->ReverseOrderingOff();													//removes objects with attribute smaller than lambda
		shapeOpeningBileLabMapFilter->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::NUMBER_OF_PIXELS);
		shapeOpeningBileLabMapFilter->SetInput(imageToBileShaLabMapFilter->GetOutput());
		shapeOpeningBileLabMapFilter->Update();
		//--------------------------------------------------------------------------------------------------------------------

		bileLabMapToImageFilter = LabelMapToLabelImageFilterType::New();
		bileLabMapToImageFilter->SetInput(shapeOpeningBileLabMapFilter->GetOutput());
		//-------------------------------------------------------------------------------------------------------------------------

		//----------FILTER---CREATE-LABEL-MAP-ON-DPPIV-CHANNEL-DO-SOME-STUFF-AND-CONVERT-IT-BACK-INTO-IMAGE------------------------
		bileOverlayImageFilter = LabelOverlayImageFilterType::New();
		bileOverlayImageFilter->SetLabelImage(bileLabMapToImageFilter->GetOutput());
		bileOverlayImageFilter->SetOpacity(m_overlayOpacity);
		bileOverlayImageFilter->SetInput(imageDPPIV);
		bileOverlayImageFilter->Update();
		//-------------------------------------------------------------------------------------------------------------------------

		//----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
		bileImageFilter = ThresholdFilterType::New();
		bileImageFilter->ReleaseDataFlagOn();
		bileImageFilter->SetOutsideValue(0);
		bileImageFilter->SetInsideValue(255);
		bileImageFilter->SetLowerThreshold(1);
		bileImageFilter->SetUpperThreshold(255);
		bileImageFilter->SetInput(bileLabMapToImageFilter->GetOutput());
		//-------------------------------------------------------------------------------------------------------------------------

		//----------FILTER-THINNING-ON-SEGMENTED-NETWORK----------------------------------------------------------------------------
		thinningBileFilter = Thinning3DImageFilterType::New();
		thinningBileFilter->ReleaseDataFlagOn();
		thinningBileFilter->SetInput(bileImageFilter->GetOutput());
		//--------------------------------------------------------------------------------------------------------------------------

		RGBVoWriterType::Pointer writer6 = RGBVoWriterType::New();
		writer6->ReleaseDataFlagOn();
		writer6->SetFileName(m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[4] + m_fileExtensionChannelDPPIV);
		writer6->SetInput(bileOverlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
		writer6->SetImageIO( itk::TIFFImageIO::New() );
#endif
		writer6->Update();

		ScalarVoWriterType::Pointer writer7 = ScalarVoWriterType::New();
		writer7->ReleaseDataFlagOn();
		writer7->SetFileName(m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[5] + m_fileExtensionChannelDPPIV);
		writer7->SetInput(bileImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
		writer7->SetImageIO( itk::TIFFImageIO::New() );
#endif
		writer7->Update();

		ScalarVoWriterType::Pointer writer8 = ScalarVoWriterType::New();
		writer8->ReleaseDataFlagOn();
		writer8->SetFileName(m_pathChannelDPPIV + m_filenameSave + m_saveSuffixes[6] + m_fileExtensionChannelDPPIV);
		writer8->SetInput(thinningBileFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
		writer8->SetImageIO( itk::TIFFImageIO::New() );
#endif
		writer8->Update();

		entry_point_counter++;
	}

	WriteLogFile(timeStamp);
	WriteDataSetSummary();

	QFile file(QString::fromStdString(m_pathChannelDPPIV + "_delete_me_temp" + m_fileExtensionChannelDPPIV));
	file.remove();
}
