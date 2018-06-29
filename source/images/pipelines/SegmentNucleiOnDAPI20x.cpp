///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentNucleiOnDAPI20x.cpp			                                 //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-04-18                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "SegmentNucleiOnDAPI20x.h"


#if (ITK_VERSION_MAJOR < 4)
#define NUMBER_OF_PIXELS SIZE
#endif

#include <QDateTime>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



SegmentNucleiOnDAPI20x::SegmentNucleiOnDAPI20x()
{
	m_overlayOpacity = 0.5;
	m_saveSuffixes[0] = "_step1_bin";
	m_saveSuffixes[1] = "_step2_hole";
	m_saveSuffixes[2] = "_step3_opening";
	m_saveSuffixes[3] = "_step4_prelimNuclei";
	m_saveSuffixes[4] = "_step5_invDistMap";
	m_saveSuffixes[5] = "_step6_bin";
	m_saveSuffixes[6] = "_step6_overlay";
}


SegmentNucleiOnDAPI20x::~SegmentNucleiOnDAPI20x()
{
	// TODO Auto-generated destructor stub
}


void SegmentNucleiOnDAPI20x::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

	file.open((m_pathChannelDAPI + m_logFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

	file << "Begin-segment-nuclei20x-------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
	file << "End-segment-nuclei20x---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


void SegmentNucleiOnDAPI20x::WriteDataSetSummary()
{
    ImageAnalysisSummaryFileIO::AddEntry(HepNucleiSegmentationBin, m_pathChannelDAPI, m_pathChannelDAPI + m_filenameSave + m_saveSuffixes[5] + m_fileExtensionChannelDAPI);
    ImageAnalysisSummaryFileIO::AddEntry(HepNucleiSegmentationOverlay, m_pathChannelDAPI, m_pathChannelDAPI + m_filenameSave + m_saveSuffixes[6] + m_fileExtensionChannelDAPI);
}


void SegmentNucleiOnDAPI20x::ParseParameterContext()
{
    if(m_paramContext->findContext("Segment Nuclei in 20x Datasets",0)==NULL) {
        std::cout << "Error: SegmentNuclei20xContext: Invalid parameter context" << std::endl;
        return;
    }

    m_fullFilenameChannelDAPI = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("DAPI channel", 0)->dataPointer()) );
	m_infoFullFilenameChannelDAPI.setFile(m_fullFilenameChannelDAPI);

	if(!m_infoFullFilenameChannelDAPI.exists())
		throw std::string("Please specify DAPI channel");

	m_pathChannelDAPI = (m_infoFullFilenameChannelDAPI.path() + QString("/")).toStdString();
    m_filenameChannelDAPI = m_infoFullFilenameChannelDAPI.baseName().toStdString();
    m_fileExtensionChannelDAPI = (QString(".") + m_infoFullFilenameChannelDAPI.suffix()).toStdString();

	m_hasNR = *(bool*)(m_paramContext->findParameter("Is there a necrotic region", 0)->dataPointer());

    m_fullFilenameNecroticRegion = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Necrotic Region", 0)->dataPointer()) );
	m_infoFullFilenameNecroticRegion.setFile(m_fullFilenameNecroticRegion);

	if(m_hasNR && !m_infoFullFilenameNecroticRegion.exists())
		throw std::string("Please specify necrotic region binary mask");

	m_pathNecroticRegion = (m_infoFullFilenameNecroticRegion.path() + QString("/")).toStdString();
    m_filenameNecroticRegion = m_infoFullFilenameNecroticRegion.baseName().toStdString();
    m_fileExtensionNecroticRegion = (QString(".") + m_infoFullFilenameNecroticRegion.suffix()).toStdString();
	
	m_fullFilenameSegCV = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Central vein segmentation", 0)->dataPointer()) );
	m_infoFullFilenameSegCV.setFile(m_fullFilenameSegCV);
	m_withCVMask = m_infoFullFilenameSegCV.exists();

	m_fullFilenameSegPV = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Portal vein segmentation", 0)->dataPointer()) );
	m_infoFullFilenameSegPV.setFile(m_fullFilenameSegPV);
	m_withPVMask = m_infoFullFilenameSegPV.exists();

    m_spacing[0] = *(double*)(m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    m_spacing[1] = *(double*)(m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    m_spacing[2] = *(double*)(m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    m_medianRadius = ( *(int*)(m_paramContext->findContext("1.1) Preprocess DAPI Channel",0)->findParameter("Median filter kernel radius", 0)->dataPointer()) );

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
    m_lowerThreshold = *(int*)(m_paramContext->findParameter("DAPI manual threshold", 0)->dataPointer());
    //TODO max;
    m_upperThreshold = 255;

    m_holeFillingNeighborhoodRadius[0] = *(int*)(m_paramContext->findContext("1.3) Hole Filling on 1.2",0)->findParameter("Radius x", 0)->dataPointer());
    m_holeFillingNeighborhoodRadius[1] = *(int*)(m_paramContext->findContext("1.3) Hole Filling on 1.2",0)->findParameter("Radius y", 0)->dataPointer());
    m_holeFillingNeighborhoodRadius[2] = *(int*)(m_paramContext->findContext("1.3) Hole Filling on 1.2",0)->findParameter("Radius z", 0)->dataPointer());
    m_holeFillingMajThreshold = *(int*)(m_paramContext->findContext("1.3) Hole Filling on 1.2",0)->findParameter("Majority threshold", 0)->dataPointer());

    int rad = ( *(int*)(m_paramContext->findContext("1.4) Opening on 1.3",0)->findParameter("Kernel radius", 0)->dataPointer()) );
    m_openingNeighborhoodRadius[0] = rad;
    m_openingNeighborhoodRadius[1] = rad;
    m_openingNeighborhoodRadius[2] = rad;

    double diameter = *(double*)(m_paramContext->findParameter("Smallest hepatocyte diameter", 0)->dataPointer());
    m_minimalNucleiRadius = diameter/2.0;
    m_minimalNucleiVolume = itk::Math::pi * 4/3 * pow(m_minimalNucleiRadius, 3);

    diameter = *(double*)(m_paramContext->findParameter("Biggest hepatocyte diameter", 0)->dataPointer());
    m_maximalNucleiRadius = diameter/2.0;
    m_maximalNucleiVolume = itk::Math::pi * 4/3 * pow(m_maximalNucleiRadius, 3);

    m_floodLevel = *(double*)(m_paramContext->findParameter("Alpha", 0)->dataPointer());

    std::string saveMode = ( (CSParameterChoice*)(m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        m_saveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        m_saveEverything = 0;
    m_logFilenameSave = *(std::string*)(m_paramContext->findParameter("Log file name", 0)->dataPointer());
    m_filenameSave = *(std::string*)(m_paramContext->findParameter("Save prefix", 0)->dataPointer());
}


void SegmentNucleiOnDAPI20x::Update()
{
    ParseParameterContext();

	std::string timeStamp;

	QDateTime time = QDateTime::currentDateTime();
	timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

	std::cout << "Segment nuclei: " << std::endl;
	std::cout << " dir: " << m_pathChannelDAPI << std::endl;
	std::cout << " file: " << m_filenameChannelDAPI << std::endl;
	std::cout << " ext: " << m_fileExtensionChannelDAPI << std::endl;

	MedianImageFilterType::Pointer                  medianFilter;
	OrImageFilterType::Pointer                      orMaskImagesFilter;
	AdaptiveOtsuThresholdImageFilterType::Pointer   adapOtsuFilter;
	OtsuThresholdImageFilterType::Pointer			otsuFilter;	
	InvertIntensityImageFilterType::Pointer			invertFilter;
	ThresholdFilterType::Pointer                    thresDAPIFilter;
	HoleFillingImageFilterType::Pointer             holeFillingFilter;
	OpeningImageFilterType::Pointer                 openingFilter;
	ImageToShapeLabelMapFilterType::Pointer         imageToNucleiShaLabMapFilter;
	ShapeOpeningLabelMapFilterType1::Pointer        shapeOpeningNucleiLabMapFilter;
	LabelMapToLabelImageFilterType1::Pointer        nucleiLabMapToImageFilter;
	LabelOverlayImageFilterType::Pointer            nucleiOverlayImageFilter;
	ThresholdFilterType::Pointer                    nucleiImageFilter;
	InvertIntensityImageFilterType::Pointer         invertIntensity1;
	SubtractImageFilterType::Pointer                subtractNecRegFromNucleiFilter;
	MaurerDistanceMapImageFilterType::Pointer       distanceMap;
	RescaleImageFilterType::Pointer                 rescaler;
	IInvertIntensityImageFilterType::Pointer        invertIntensity2;
	MorphoWatershedImageFilterType::Pointer         morphWatershed;
	MaskImageFilterType::Pointer                    maskImage;
	LabelImageToShapeLabelMapFilterType::Pointer    watershedImageToLabelMap;
	ShapeOpeningLabelMapFilterType2::Pointer        watershedShapeOpeningLabMapFilter;
	ShapeOpeningLabelMapFilterType3::Pointer        watershedShapeOpeningLabMapFilter2;
	LabelMapToLabelImageFilterType2::Pointer        watershedLabelMapToImage;
	ThresholdFilterType::Pointer                    watershedImageFilter;
	LabelOverlayImageFilterType::Pointer            watershedOverlayImageFilter;


	//----------READER---------------------------------------------------------------------------------------------------------
	if(GetNumberOfDimensions(m_pathChannelDAPI + m_filenameChannelDAPI + m_fileExtensionChannelDAPI) != 3)
	    throw std::string("Please specify a DAPI channel file with three dimensions.");

    CScalarVoImageType::Pointer readerImage = CScalarVoImageType::New();
    ReadImage(m_pathChannelDAPI + m_filenameChannelDAPI + m_fileExtensionChannelDAPI, readerImage, m_spacing);

    CScalarVoImageType::Pointer necRegImage = CScalarVoImageType::New();
    if(m_hasNR)
        ReadImage(m_pathNecroticRegion + m_filenameNecroticRegion + m_fileExtensionNecroticRegion, necRegImage, m_spacing);

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

    //----------PREPROCESS-WITH-MEDIAN-FILTERING-------------------------------------------------------------------------------
    MedianImageFilterType::InputSizeType radius;						//create a size object (basically an array)
    radius.Fill(m_medianRadius);														//fill it with the radius for the median filtering

    medianFilter = MedianImageFilterType::New();						//instantiation of median filter object
    medianFilter->SetRadius(radius);									//set the radius for the median filter
    medianFilter->SetInput(readerImage);								//set as input for the filter the imge that is stored in the readerImage object

    ScalarVoWriterType::Pointer writerMedian = ScalarVoWriterType::New();
    writerMedian->ReleaseDataFlagOn();
    writerMedian->SetFileName(m_pathChannelDAPI + "nuclei_step0_median" + m_fileExtensionChannelDAPI);
    writerMedian->SetInput(medianFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writerMedian->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerMedian->Update();
    //-------------------------------------------------------------------------------------------------------------------------

	//----------FILTER-THRESHOLD-ON-DAPI-CHANNEL-------------------------------------------------------------------------------
	//Pipeline allows choice between three different threshold methods
	switch(m_thresholdingMode)											//switch for the threshold mode
	{
	case 0:																//adaptive otsu case
	    {
		//this filter calculates locally (at random sample points) (within a region with specified radius) with otsu method thresholds
		//these thresholds are interpolated over data set with B-Splines so that for every pixel an individual threshold value is available
	    adapOtsuFilter = AdaptiveOtsuThresholdImageFilterType::New();	//instantiation of filter object (filter is templated over the image type passed on from the reader -> see *.h file)
	    adapOtsuFilter->SetInput(medianFilter->GetOutput());			//Input for this filter is the output of the median filter
	    adapOtsuFilter->SetInsideValue(255);							//input pixels within the threshold range will be set to this value in the output image (value for segmented objects)
	    adapOtsuFilter->SetOutsideValue(0);								//input pixels outside of threshold range will be set to this value in the output image (background value for segmented image)
	    adapOtsuFilter->SetNumberOfHistogramBins(256);					//specify number of intensity values
	    adapOtsuFilter->SetSplineOrder(3);								//set spline order for spline interpolation
	    adapOtsuFilter->SetNumberOfControlPoints(5);					//set number of control points for spline interpolation
	    adapOtsuFilter->SetNumberOfLevels(3);							//set number of levels for spline interpolation
	    adapOtsuFilter->SetNumberOfSamples(m_adapOtsuSamplePoints);		//set number of sample points at which otsu thresholds will be calculated
	    adapOtsuFilter->SetRadius(m_adapOtsuRadius);					//set radius for region for which otsu threshold is calculated
	    if(m_withCVMask || m_withPVMask)
	        adapOtsuFilter->SetMaskImage(maskVeinImage);
	    break;
	    }
	case 1:																//(standard) otsu case
	    {
	    otsuFilter = OtsuThresholdImageFilterType::New();				//instantiation of otsu-threshold filter
	    otsuFilter->SetInput(medianFilter->GetOutput());				//set as input output of median filter
	    otsuFilter->Update();											//call update = process pipeline: read -> median -> otsu threshold: an output image is now available
	    m_otsuThreshold = (CScalarPixelType)(otsuFilter->GetThreshold());

	    invertFilter = InvertIntensityImageFilterType::New();			//necessary to invert result, because foreground value is 0, but should be 255
	    invertFilter->SetInput(otsuFilter->GetOutput());
	    invertFilter->SetMaximum(255);
	    break;
	    }
	default:															//(manual) binary threshold case
	    {
		thresDAPIFilter = ThresholdFilterType::New();					//instantiation of binary threshold filter
		thresDAPIFilter->SetInput(medianFilter->GetOutput());			//set as input for this filter the output of the median filter
		thresDAPIFilter->SetOutsideValue(0);							//set background value of thresholded image
		thresDAPIFilter->SetInsideValue(255);							//set foreground value of thresholded image
		thresDAPIFilter->SetLowerThreshold(m_lowerThreshold);			//set lower threshold
		thresDAPIFilter->SetUpperThreshold(m_upperThreshold);			//set upper threshold
		break;
	    }
	}

	if(m_saveEverything) {												//if save flag is set
		ScalarVoWriterType::Pointer writer1 = ScalarVoWriterType::New();		//instantiation of writer
		writer1->ReleaseDataFlagOn();										
		writer1->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixes[0] + m_fileExtensionChannelDAPI);		//set file name
		switch(m_thresholdingMode)												//depending on threshold mode set as input for writer the output of the threshold filter
		{
		case 0:
		    writer1->SetInput(adapOtsuFilter->GetOutput());
		    break;
		case 1:
		    writer1->SetInput(invertFilter->GetOutput());
		    break;
		default:
		    writer1->SetInput(thresDAPIFilter->GetOutput());
		    break;
		}
#if (ITK_VERSION_MAJOR >= 4)
		writer1->SetImageIO( itk::TIFFImageIO::New() );
#endif
		writer1->Update();														//invoke the execution of the pipeline (reader -> median -> threshold -> writer): image is saved
	}
	//-------------------------------------------------------------------------------------------------------------------------

	//----------FILTER-HOLE-FILLING-ON-DAPI-CHANNEL----------------------------------------------------------------------------
	holeFillingFilter = HoleFillingImageFilterType::New();
	holeFillingFilter->SetRadius(m_holeFillingNeighborhoodRadius);
	holeFillingFilter->SetBackgroundValue(0);
	holeFillingFilter->SetForegroundValue(255);
	holeFillingFilter->SetMajorityThreshold(m_holeFillingMajThreshold);					//number of foreground neighbors should be at least (3x3x3-1)/2 + majority
	holeFillingFilter->SetMaximumNumberOfIterations(30);
	switch(m_thresholdingMode)
	{
	case 0:
	    holeFillingFilter->SetInput(adapOtsuFilter->GetOutput());
	    break;
	case 1:
	    holeFillingFilter->SetInput(invertFilter->GetOutput());
	    break;
	default:
	    holeFillingFilter->SetInput(thresDAPIFilter->GetOutput());
	    break;
	}
	//-------------------------------------------------------------------------------------------------------------------------

	if(m_saveEverything) {
		ScalarVoWriterType::Pointer writer2 = ScalarVoWriterType::New();
		writer2->ReleaseDataFlagOn();
		writer2->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixes[1] + m_fileExtensionChannelDAPI);
		writer2->SetInput(holeFillingFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer2->SetImageIO( itk::TIFFImageIO::New() );
#endif
		writer2->Update();
	}
	//-------------------------------------------------------------------------------------------------------------------------

	//----------FILTER-OPENING-ON-DAPI-CHANNEL---------------------------------------------------------------------------------
	openingFilter = OpeningImageFilterType::New();
	openingFilter->SetRadius(m_openingNeighborhoodRadius);
	openingFilter->SetBackgroundValue(0);
	openingFilter->SetForegroundValue(255);
	openingFilter->SetInput(holeFillingFilter->GetOutput());
	openingFilter->Update();
	//-------------------------------------------------------------------------------------------------------------------------

	//----------FILTER-SUBSTRACT-NECREGION-FROM-DAPI-CHANNEL-------------------------------------------------------------------
	itk::SmartPointer<CScalarVoImageType> inputForDistanceMapImage;

	if(m_hasNR) {
	    subtractNecRegFromNucleiFilter = SubtractImageFilterType::New();
	    subtractNecRegFromNucleiFilter->ReleaseDataFlagOn();
	    subtractNecRegFromNucleiFilter->SetInput1(openingFilter->GetOutput());
	    subtractNecRegFromNucleiFilter->SetInput2(necRegImage);
	    subtractNecRegFromNucleiFilter->Update();

	    inputForDistanceMapImage = subtractNecRegFromNucleiFilter->GetOutput();
	}
	else
	    inputForDistanceMapImage = openingFilter->GetOutput();

	inputForDistanceMapImage->DisconnectPipeline();
	inputForDistanceMapImage->SetSpacing(m_spacing);

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer3 = ScalarVoWriterType::New();
        writer3->ReleaseDataFlagOn();
        writer3->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixes[2] + m_fileExtensionChannelDAPI);
        writer3->SetInput(inputForDistanceMapImage);
#if (ITK_VERSION_MAJOR >= 4)
        writer3->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer3->Update();
    }
	//--------------------------------------------------------------------------------------------------------------------

	//----------DISTANCE-MAP---------------------------------------------------------------------------------------------------
	distanceMap = MaurerDistanceMapImageFilterType::New();
	distanceMap->SetUseImageSpacing(true);
	distanceMap->SquaredDistanceOff();
	distanceMap->SetBackgroundValue(255);
	distanceMap->SetInsideIsPositive(true);
	distanceMap->SetInput(inputForDistanceMapImage);
	distanceMap->Update();
	//-------------------------------------------------------------------------------------------------------------------------

	//----------FILTER-INVERT-IMAGE-FILTER-------------------------------------------------------------------------------------
	rescaler = RescaleImageFilterType::New();
	rescaler->ReleaseDataFlagOn();
	rescaler->SetInput(distanceMap->GetOutput());

	if(m_saveEverything) {
	    SScalarVoWriterType::Pointer writer5 = SScalarVoWriterType::New();
	    writer5->ReleaseDataFlagOn();
	    writer5->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixes[4] + m_fileExtensionChannelDAPI);
	    writer5->SetInput(rescaler->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
	    writer5->SetImageIO( itk::TIFFImageIO::New() );
#endif
	    writer5->Update();
	}
	//-------------------------------------------------------------------------------------------------------------------------

	//----------WATERSHED-----------------------------------------------------------------------------------------------------
	morphWatershed = MorphoWatershedImageFilterType::New();
	morphWatershed->SetLevel(m_floodLevel);
	morphWatershed->FullyConnectedOn();
	morphWatershed->ReleaseDataFlagOn();
	morphWatershed->SetInput(distanceMap->GetOutput());
	//-------------------------------------------------------------------------------------------------------------------------

	//----------MASK-THE-WATERSHED-SEGMENTS------------------------------------------------------------------------------------
	maskImage = MaskImageFilterType::New();
	maskImage->ReleaseDataFlagOn();
	maskImage->SetInput1(morphWatershed->GetOutput());
	maskImage->SetInput2(inputForDistanceMapImage);
	maskImage->Update();

	itk::SmartPointer<IScalarVoImageType> maskedImage = maskImage->GetOutput();
	maskedImage->DisconnectPipeline();

	maskedImage->SetSpacing(m_spacing);
	//-------------------------------------------------------------------------------------------------------------------------

	//----------TO-LABEL-MAP-AND-BACK-----------------------------------------------------------------------------------------
	watershedImageToLabelMap = LabelImageToShapeLabelMapFilterType::New();
	watershedImageToLabelMap->ReleaseDataFlagOn();
	watershedImageToLabelMap->ComputeFeretDiameterOn();
	watershedImageToLabelMap->SetNumberOfThreads(8);
	watershedImageToLabelMap->SetInput(maskedImage);
	watershedImageToLabelMap->Update();

	itk::SmartPointer<LabelImageToShapeLabelMapFilterType::OutputImageType> hepNucleiLabelMap = watershedImageToLabelMap->GetOutput();
	hepNucleiLabelMap->DisconnectPipeline();
	std::cout << "hep nuclei label map spacing (after feret diameter calculation) = " << hepNucleiLabelMap->GetSpacing() << std::endl;

	//----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------------
	double minHepDiameter = 2.0*m_minimalNucleiRadius;
	double maxHepDiameter = 2.0*m_maximalNucleiRadius;

	std::cout << "after watersheding, before classification: hep nuclei " << hepNucleiLabelMap->GetNumberOfLabelObjects() << std::endl;

	std::vector<unsigned long> labelsToRemoveFromHepNucleiLabelMap;

	for(unsigned int i=0; i<hepNucleiLabelMap->GetNumberOfLabelObjects(); i++) {
	    float distToBorder = 0;

	    for(unsigned int j=0; j < hepNucleiLabelMap->GetNthLabelObject(i)->GetNumberOfPixels(); j++) {
	        itk::Index<3> pos = hepNucleiLabelMap->GetNthLabelObject(i)->GetIndex(j);
	        float d = fabs( distanceMap->GetOutput()->GetPixel(pos) );

	        if(d > distToBorder)
	            distToBorder = d;
	    }

	    double inscSphereDiam = 2. * distToBorder;
	    double roundness = (2. * hepNucleiLabelMap->GetNthLabelObject(i)->GetEquivalentSphericalRadius()) / hepNucleiLabelMap->GetNthLabelObject(i)->GetFeretDiameter();

	    bool isHep = true;

	    if( inscSphereDiam < minHepDiameter || maxHepDiameter < inscSphereDiam ) {
	        labelsToRemoveFromHepNucleiLabelMap.push_back(hepNucleiLabelMap->GetNthLabelObject(i)->GetLabel());
	        isHep = false;
	    }
	}

	std::cout << "labelsToRemoveFromHepNucleiLabelMap.size(): " << labelsToRemoveFromHepNucleiLabelMap.size() << std::endl;
	for(unsigned int i=0; i<labelsToRemoveFromHepNucleiLabelMap.size(); i++)
	    hepNucleiLabelMap->RemoveLabel(labelsToRemoveFromHepNucleiLabelMap[i]);
	std::cout << "after watersheding: hep nuclei " << hepNucleiLabelMap->GetNumberOfLabelObjects() << std::endl;

	LabelMapToLabelImageFilterType2::Pointer hepNucleiLabelMapToImage = LabelMapToLabelImageFilterType2::New();
	hepNucleiLabelMapToImage->SetInput(hepNucleiLabelMap);
	hepNucleiLabelMapToImage->Update();
	//-------------------------------------------------------------------------------------------------------------------------

	//----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
    watershedImageFilter = ThresholdFilterType::New();
    watershedImageFilter->SetOutsideValue(0);
    watershedImageFilter->SetInsideValue(255);
    watershedImageFilter->SetLowerThreshold(1);
    watershedImageFilter->SetUpperThreshold(255);
    watershedImageFilter->SetInput(hepNucleiLabelMapToImage->GetOutput());

    ScalarVoWriterType::Pointer writer6 = ScalarVoWriterType::New();
    writer6->ReleaseDataFlagOn();
    writer6->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixes[5] + m_fileExtensionChannelDAPI);
    writer6->SetInput(watershedImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer6->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer6->Update();
	//-------------------------------------------------------------------------------------------------------------------------

    //----------VISULALIZATION-ONLY-FILTER---NOT-NEEDED-BY-OTHER-PIPELINE-FILTERS----------------------------------------------
    watershedOverlayImageFilter = LabelOverlayImageFilterType::New();
    watershedOverlayImageFilter->SetOpacity(m_overlayOpacity);
    watershedOverlayImageFilter->ReleaseDataFlagOn();
    watershedOverlayImageFilter->SetInput(readerImage);
    watershedOverlayImageFilter->SetLabelImage(hepNucleiLabelMapToImage->GetOutput());

    RGBVoWriterType::Pointer writer7 = RGBVoWriterType::New();
    writer7->ReleaseDataFlagOn();
    writer7->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixes[6] + m_fileExtensionChannelDAPI);
    writer7->SetInput(watershedOverlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer7->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer7->Update();
    //-------------------------------------------------------------------------------------------------------------------------
	
	WriteLogFile(timeStamp);
	WriteDataSetSummary();
}
