///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentNucleiOnDAPI60x.cpp                                           //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-07                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "SegmentNucleiOnDAPI60x.h"

#if (ITK_VERSION_MAJOR < 4)
#define NUMBER_OF_PIXELS SIZE
#endif

#include <QDateTime>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"


#include "itkPointSet.h"


template< unsigned int VImageDimension > SegmentNucleiOnDAPI60x< VImageDimension >::SegmentNucleiOnDAPI60x()
{
    m_overlayOpacity = 0.5;

    m_saveSuffixesForFinals[0] = "_step6_nonHepNuclei_bin";
    m_saveSuffixesForFinals[1] = "_step6_hepNuclei_bin";
    m_saveSuffixesForFinals[2] = "_step6_nonHepNuclei_overlay";
    m_saveSuffixesForFinals[3] = "_step6_hepNuclei_overlay";
}


template< unsigned int VImageDimension > SegmentNucleiOnDAPI60x< VImageDimension >::~SegmentNucleiOnDAPI60x()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void SegmentNucleiOnDAPI60x< VImageDimension >::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    this->m_paramContext->dump(parameter);

    file.open((m_pathChannelDAPI + m_logFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-segment-nuclei60x-------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-segment-nuclei60x---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


template< unsigned int VImageDimension > void SegmentNucleiOnDAPI60x< VImageDimension >::WriteDataSetSummary()
{
    ImageAnalysisSummaryFileIO::AddEntry(NonHepNucleiSegmentationBin, m_pathChannelDAPI, m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[0] + m_fileExtensionChannelDAPI);
    ImageAnalysisSummaryFileIO::AddEntry(NonHepNucleiSegmentationOverlay, m_pathChannelDAPI, m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[2] + m_fileExtensionChannelDAPI);
    ImageAnalysisSummaryFileIO::AddEntry(HepNucleiSegmentationBin, m_pathChannelDAPI, m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[1] + m_fileExtensionChannelDAPI);
    ImageAnalysisSummaryFileIO::AddEntry(HepNucleiSegmentationOverlay, m_pathChannelDAPI, m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[3] + m_fileExtensionChannelDAPI);
}


template< unsigned int VImageDimension > void SegmentNucleiOnDAPI60x< VImageDimension >::ParseParameterContext()
{
    if(this->m_paramContext->findContext("Segment Nuclei in 60x Datasets",0)==NULL) {
        std::cout << "Error: SegmentNuclei60xContext: Invalid parameter context" << std::endl;
        return;
    }

    m_fullFilenameChannelDAPI = QString::fromStdString( *(std::string*)(this->m_paramContext->findParameter("DAPI channel", 0)->dataPointer()) );
	m_infoFullFilenameChannelDAPI.setFile(m_fullFilenameChannelDAPI);

	if(!m_infoFullFilenameChannelDAPI.exists())
		throw std::string("Please specify DAPI channel");

	m_pathChannelDAPI = (m_infoFullFilenameChannelDAPI.path() + QString("/")).toStdString();
    m_filenameChannelDAPI = m_infoFullFilenameChannelDAPI.completeBaseName().toStdString();
    m_fileExtensionChannelDAPI = (QString(".") + m_infoFullFilenameChannelDAPI.suffix()).toStdString();

    m_fullFilenameSegCV = QString::fromStdString( *(std::string*)(this->m_paramContext->findParameter("Central vein segmentation", 0)->dataPointer()) );
    m_infoFullFilenameSegCV.setFile(m_fullFilenameSegCV);
    m_withCVMask = m_infoFullFilenameSegCV.exists();

    m_fullFilenameSegPV = QString::fromStdString( *(std::string*)(this->m_paramContext->findParameter("Portal vein segmentation", 0)->dataPointer()) );
    m_infoFullFilenameSegPV.setFile(m_fullFilenameSegPV);
    m_withPVMask = m_infoFullFilenameSegPV.exists();

    m_spacing[0] = *(double*)(this->m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    m_spacing[1] = *(double*)(this->m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    if(ImageDimension==3)   m_spacing[2] = *(double*)(this->m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    m_medianRadius = ( *(int*)(this->m_paramContext->findContext("1.1) Preprocess DAPI Channel",0)->findParameter("Median filter kernel radius", 0)->dataPointer()) );
    m_greyscaleOpeningRadius = ( *(int*)(this->m_paramContext->findContext("1.1) Preprocess DAPI Channel",0)->findParameter("Greyscale opening kernel radius", 0)->dataPointer()) );

    std::string thresMode = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Threshold mode", 0)->dataPointer()) )->currentString();
    if( thresMode.compare("Adaptive Otsu-Threshold")==0 )
        m_thresholdingMode = 0;
    else if( thresMode.compare("Otsu-Threshold")==0 )
        m_thresholdingMode = 1;
    else if( thresMode.compare("Manual Threshold")==0 )
        m_thresholdingMode = 2;

    m_adapOtsuRadius[0] = *(int*)(this->m_paramContext->findParameter("Sample region size x", 0)->dataPointer());
    m_adapOtsuRadius[1] = *(int*)(this->m_paramContext->findParameter("Sample region size y", 0)->dataPointer());
    if(ImageDimension==3)   m_adapOtsuRadius[2] = *(int*)(this->m_paramContext->findParameter("Sample region size z", 0)->dataPointer());
    m_adapOtsuSamplePoints = *(int*)(this->m_paramContext->findParameter("Number samples", 0)->dataPointer());
    m_lowerThreshold = *(int*)(this->m_paramContext->findParameter("DAPI manual min threshold", 0)->dataPointer());
    m_upperThreshold = *(int*)(this->m_paramContext->findParameter("DAPI manual max threshold", 0)->dataPointer());

    m_invHoleFillingNeighborhoodRadius[0] = *(int*)(this->m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2",0)->findParameter("Radius x", 0)->dataPointer());
    m_invHoleFillingNeighborhoodRadius[1] = *(int*)(this->m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2",0)->findParameter("Radius y", 0)->dataPointer());
    if(ImageDimension==3)   m_invHoleFillingNeighborhoodRadius[2] = *(int*)(this->m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2",0)->findParameter("Radius z", 0)->dataPointer());
    m_invHoleFillingMajThreshold = *(int*)(this->m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2",0)->findParameter("Majority threshold", 0)->dataPointer());

    m_cavWithRescaling = *(bool*)(this->m_paramContext->findContext("1.4) Cavity Filling on 1.3", 0)->findParameter("Fast (less accurate)", 0)->dataPointer());
    m_cavitiyFillingRadius = *(int*)(this->m_paramContext->findContext("1.4) Cavity Filling on 1.3",0)->findParameter("Radius", 0)->dataPointer());
    m_cavitiyFillingMinThreshold = *(double*)(this->m_paramContext->findContext("1.4) Cavity Filling on 1.3",0)->findParameter("Minimal fraction of surrounding foreground", 0)->dataPointer());
    m_cavitiyFillingMaxThreshold = *(double*)(this->m_paramContext->findContext("1.4) Cavity Filling on 1.3",0)->findParameter("Maximal fraction of surrounding foreground", 0)->dataPointer());

    typename CScalarImageType::SizeType rad;
    rad[0] = ( *(int*)(this->m_paramContext->findContext("1.5) Closing on 1.4", 0)->findParameter("Kernel radius x", 0)->dataPointer()) );
    rad[1] = ( *(int*)(this->m_paramContext->findContext("1.5) Closing on 1.4", 0)->findParameter("Kernel radius y", 0)->dataPointer()) );
    if(ImageDimension==3)   rad[2] = ( *(int*)(this->m_paramContext->findContext("1.5) Closing on 1.4", 0)->findParameter("Kernel radius z", 0)->dataPointer()) );
    m_closingStructuringElement.SetRadius(rad);
    m_closingStructuringElement.CreateStructuringElement();

    rad[0] = ( *(int*)(this->m_paramContext->findContext("1.6) Opening on 1.5", 0)->findParameter("Kernel radius x", 0)->dataPointer()) );
    rad[1] = ( *(int*)(this->m_paramContext->findContext("1.6) Opening on 1.5", 0)->findParameter("Kernel radius y", 0)->dataPointer()) );
    if(ImageDimension==3)   rad[2] = ( *(int*)(this->m_paramContext->findContext("1.6) Opening on 1.5", 0)->findParameter("Kernel radius z", 0)->dataPointer()) );
    m_openingStructuringElement.SetRadius(rad);
    m_openingStructuringElement.CreateStructuringElement();

    double diameter = *(double*)(this->m_paramContext->findParameter("Smallest non-hepatocyte diameter", 0)->dataPointer());
    m_minNonHepRadius = diameter/2.0;
    if(ImageDimension==2)       m_minNonHepVolume = itk::Math::pi * pow(m_minNonHepRadius, 2);              //Area in this case
    else if(ImageDimension==3)  m_minNonHepVolume = itk::Math::pi * 4.0/3.0 * pow(m_minNonHepRadius, 3);

    diameter = *(double*)(this->m_paramContext->findParameter("Biggest non-hepatocyte diameter", 0)->dataPointer());
    m_maxNonHepRadius = diameter/2.0;
    if(ImageDimension==2)       m_maxNonHepVolume = itk::Math::pi * pow(m_maxNonHepRadius, 2);              //Area in this case
    else if(ImageDimension==3)  m_maxNonHepVolume = itk::Math::pi * 4.0/3.0 * pow(m_maxNonHepRadius, 3);

    diameter = *(double*)(this->m_paramContext->findParameter("Smallest hepatocyte diameter", 0)->dataPointer());
    m_minHepRadius = diameter/2.0;
    if(ImageDimension==2)       m_minHepVolume = itk::Math::pi * pow(m_minHepRadius, 2);              //Area in this case
    else if(ImageDimension==3)  m_minHepVolume = itk::Math::pi * 4.0/3.0 * pow(m_minHepRadius, 3);

    diameter = *(double*)(this->m_paramContext->findParameter("Biggest hepatocyte diameter", 0)->dataPointer());
    m_maxHepRadius = diameter/2.0;
    if(ImageDimension==2)   m_maxHepVolume = itk::Math::pi * pow(m_maxHepRadius, 2);              //Area in this case
    if(ImageDimension==3)   m_maxHepVolume = itk::Math::pi * 4.0/3.0 * pow(m_maxHepRadius, 3);

    m_roundnessThreshold = *(double*)(this->m_paramContext->findParameter("Roundness of hepatocyte nuclei [0,1]", 0)->dataPointer());

    m_floodLevel = *(double*)(this->m_paramContext->findParameter("Alpha", 0)->dataPointer());

    std::string saveMode = ( (CSParameterChoice*)(this->m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        m_saveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        m_saveEverything = 0;
    m_logFilenameSave = *(std::string*)(this->m_paramContext->findParameter("Log file name", 0)->dataPointer());
    m_filenameSave = *(std::string*)(this->m_paramContext->findParameter("Save prefix", 0)->dataPointer());
}


template< unsigned int VImageDimension > void SegmentNucleiOnDAPI60x< VImageDimension >::ClassifyNuclei2D(itk::SmartPointer<ShapeLabelMapType> nucleiLabelMap,
        std::vector<unsigned long> &labelsToRemoveFromHepNucleiLabelMap, std::vector<unsigned long> &labelsToRemoveFromNonHepNucleiLabelMap, bool withNucleiAnalysisOutput)
{
    for(unsigned int i=0; i<nucleiLabelMap->GetNumberOfLabelObjects(); i++) {
        double area = nucleiLabelMap->GetNthLabelObject(i)->GetPhysicalSize();
        double elongation = nucleiLabelMap->GetNthLabelObject(i)->GetElongation();

        if( elongation > m_roundnessThreshold || area < m_minHepVolume || m_maxHepVolume < area)
            labelsToRemoveFromHepNucleiLabelMap.push_back(nucleiLabelMap->GetNthLabelObject(i)->GetLabel());
        else
            labelsToRemoveFromNonHepNucleiLabelMap.push_back(nucleiLabelMap->GetNthLabelObject(i)->GetLabel());
    }
    if(withNucleiAnalysisOutput)
        ;   //do something
}


template< unsigned int VImageDimension > void SegmentNucleiOnDAPI60x< VImageDimension >::ClassifyNuclei3D(itk::SmartPointer<ShapeLabelMapType> nucleiLabelMap, itk::SmartPointer<FScalarImageType> distanceMap,
        std::vector<unsigned long> &labelsToRemoveFromHepNucleiLabelMap, std::vector<unsigned long> &labelsToRemoveFromNonHepNucleiLabelMap, bool withNucleiAnalysisOutput)
{
    double minNonHepDiameter = 2.0*m_minNonHepRadius;
    double maxNonHepDiameter = 2.0*m_maxNonHepRadius;
    double minHepDiameter = 2.0*m_minHepRadius;
    double maxHepDiameter = 2.0*m_maxHepRadius;

    std::fstream file1;
    if(withNucleiAnalysisOutput) {
        file1.open((m_pathChannelDAPI + "nuclei_roundness.txt").c_str(), std::fstream::out);
        file1.width(20);
        file1 << "roundness";
        file1.width(20);
        file1 << "elongation";
        file1.width(20);
        file1 << "feretDia";
        file1.width(20);
        file1 << "equSphDia";
        file1.width(20);
        file1 << "inscrSphDia";
        file1.width(20);
        file1 << "class" << std::endl;
    }

    for(unsigned int i=0; i<nucleiLabelMap->GetNumberOfLabelObjects(); i++) {
        float distToBorder = 0;

        for(unsigned int j=0; j < nucleiLabelMap->GetNthLabelObject(i)->GetNumberOfPixels(); j++) {
            itk::Index<ImageDimension> pos = nucleiLabelMap->GetNthLabelObject(i)->GetIndex(j);
            float d = fabs( distanceMap->GetPixel(pos) );

            if(d > distToBorder)
                distToBorder = d;
        }

        double inscSphereDiam = 2. * distToBorder;
        double roundness = (2. * nucleiLabelMap->GetNthLabelObject(i)->GetEquivalentSphericalRadius()) / nucleiLabelMap->GetNthLabelObject(i)->GetFeretDiameter();

        if(withNucleiAnalysisOutput) {
            file1.width(20);
            file1 << roundness;
            file1.width(20);
            file1 << nucleiLabelMap->GetNthLabelObject(i)->GetElongation();
            file1.width(20);
            file1 << nucleiLabelMap->GetNthLabelObject(i)->GetFeretDiameter();
            file1.width(20);
            file1 << (2. * nucleiLabelMap->GetNthLabelObject(i)->GetEquivalentSphericalRadius());
            file1.width(20);
            file1 << inscSphereDiam;
        }

        bool isHep = true, isNonHep = true;

        if( inscSphereDiam < minNonHepDiameter || maxNonHepDiameter < inscSphereDiam ) {
            labelsToRemoveFromNonHepNucleiLabelMap.push_back(nucleiLabelMap->GetNthLabelObject(i)->GetLabel());
            isNonHep = false;
        }

        if( inscSphereDiam < minHepDiameter || maxHepDiameter < inscSphereDiam ) {
            labelsToRemoveFromHepNucleiLabelMap.push_back(nucleiLabelMap->GetNthLabelObject(i)->GetLabel());
            isHep = false;
        }

        if( minHepDiameter <= inscSphereDiam && inscSphereDiam <= maxNonHepDiameter ) {
            if( roundness > m_roundnessThreshold ) {
                labelsToRemoveFromNonHepNucleiLabelMap.push_back(nucleiLabelMap->GetNthLabelObject(i)->GetLabel());
                isNonHep = false;
            }
            else {
                labelsToRemoveFromHepNucleiLabelMap.push_back(nucleiLabelMap->GetNthLabelObject(i)->GetLabel());
                isHep = false;
            }
        }

        if(withNucleiAnalysisOutput) {
            file1.width(20);
            if(!isHep && !isNonHep)
                file1 << "1" << std::endl;
            else if(!isHep && isNonHep)
                file1 << "2" << std::endl;
            else if(isHep && !isNonHep)
                file1 << "3" << std::endl;
            else if(isHep && isNonHep)
                file1 << "4" << std::endl;
        }
    }

    if(withNucleiAnalysisOutput)
        file1.close();
}


template< unsigned int VImageDimension > void SegmentNucleiOnDAPI60x< VImageDimension >::Update()
{
    ParseParameterContext();

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    std::cout << "Segment nuclei: " << std::endl;
    std::cout << " dir: " << m_pathChannelDAPI << std::endl;
    std::cout << " file: " << m_filenameChannelDAPI << std::endl;
    std::cout << " ext: " << m_fileExtensionChannelDAPI << std::endl;


    typename OrImageFilterType::Pointer                      orMaskImagesFilter;
    typename AdaptiveOtsuThresholdImageFilterType::Pointer   adapOtsuFilter;
    typename OtsuThresholdImageFilterType::Pointer           otsuFilter;
    typename InvertIntensityImageFilterType::Pointer         invertFilter;
    typename NearestNeighborExtrapolatorImageFunctionType::Pointer   extrapolator;
    typename NearestNeighborInterpolatorImageFunctionType::Pointer   interpolator;
    typename ResampleImageFilterType::Pointer                downscalingFilter;
    typename CavityFillingImageFilter::Pointer               cavityFillingFilter;
    typename ResampleImageFilterType::Pointer                upscalingFilter;
    typename ThresholdFilterType::Pointer                    thresDAPIFilter;
    typename HoleFillingImageFilterType::Pointer             holeFillingFilter;
    typename ClosingImageFilterType::Pointer                 closingFilter;
    typename OpeningImageFilterType::Pointer                 openingFilter;
    typename ErodeImageFilterType::Pointer                   erodeFilter;
    typename ImageToShapeLabelMapFilterType::Pointer         imageToNucleiShaLabMapFilter;
    typename ShapeOpeningLabelMapFilterType::Pointer         shapeOpeningNucleiLabMapFilter;
    typename LabelMapToLabelImageFilterType::Pointer         nucleiLabMapToImageFilter;
    typename LabelOverlayImageFilterType::Pointer            nucleiOverlayImageFilter;
    typename ThresholdFilterType::Pointer                    nucleiImageFilter;
    typename InvertIntensityImageFilterType::Pointer         invertIntensity1;
    typename MaurerDistanceMapImageFilterType::Pointer       distanceMap;
    typename RescaleImageFilterType::Pointer                 rescaler1, rescaler2;
    typename FInvertIntensityImageFilterType::Pointer        invertIntensity2;
    typename MorphoWatershedImageFilterType::Pointer         morphWatershed;
    typename MaskImageFilterType::Pointer                    maskImageFilter;
    typename LabelImageToShapeLabelMapFilterType::Pointer    watershedImageToLabelMap;
    typename ShapeOpeningLabelMapFilterType::Pointer         watershedShapeOpeningLabMapFilter;
    typename LabelMapToLabelImageFilterType::Pointer         watershedLabelMapToImage;
    typename ThresholdFilterType::Pointer                    watershedImageFilter;
    typename LabelOverlayImageFilterType::Pointer            watershedOverlayImageFilter;


    //----------READER---------------------------------------------------------------------------------------------------------
    if(this->GetNumberOfDimensions(m_pathChannelDAPI + m_filenameChannelDAPI + m_fileExtensionChannelDAPI) != ImageDimension)
        throw std::string("Data " + m_pathChannelDAPI + m_filenameChannelDAPI + m_fileExtensionChannelDAPI + " has different dimensionality than expected by pipeline.");

    typename CScalarImageType::Pointer readerDAPIImage = CScalarImageType::New();
    this->ReadImage(m_pathChannelDAPI + m_filenameChannelDAPI + m_fileExtensionChannelDAPI, readerDAPIImage, m_spacing);

    typename CScalarImageType::Pointer segCVImage = CScalarImageType::New();
    if(m_withCVMask)
        this->ReadImage(m_fullFilenameSegCV.toStdString(), segCVImage, m_spacing);

    typename CScalarImageType::Pointer segPVImage = CScalarImageType::New();
    if(m_withPVMask)
        this->ReadImage(m_fullFilenameSegPV.toStdString(), segPVImage, m_spacing);

    typename CScalarImageType::Pointer maskVeinImage = CScalarImageType::New();
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
    typename MedianImageFilterType::InputSizeType medianRadius;
    medianRadius.Fill(m_medianRadius);

    typename MedianImageFilterType::Pointer medianFilter = MedianImageFilterType::New();
    medianFilter->SetRadius(medianRadius);
    medianFilter->SetInput(readerDAPIImage);

    if(m_saveEverything) {
        typename CScalarImageWriterType::Pointer writerMedian1 = CScalarImageWriterType::New();
        writerMedian1->ReleaseDataFlagOn();
        writerMedian1->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step0a_median" + m_fileExtensionChannelDAPI);
        writerMedian1->SetInput(medianFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerMedian1->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerMedian1->Update();
    }

    StructuringElementType greyOpeningKernel;
    greyOpeningKernel.SetRadius(m_greyscaleOpeningRadius);
    greyOpeningKernel.CreateStructuringElement();

    typename GrayscaleErodeImageFilterType::Pointer greyErodeFilter = GrayscaleErodeImageFilterType::New();
    greyErodeFilter->SetInput(medianFilter->GetOutput());
    greyErodeFilter->SetKernel(greyOpeningKernel);

    if(m_saveEverything) {
        typename CScalarImageWriterType::Pointer writerErode = CScalarImageWriterType::New();
        writerErode->ReleaseDataFlagOn();
        writerErode->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step0b_erode" + m_fileExtensionChannelDAPI);
        writerErode->SetInput(greyErodeFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerErode->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerErode->Update();
    }

    typename GrayscaleDilateImageFilterType::Pointer dilateFilter = GrayscaleDilateImageFilterType::New();
    dilateFilter->SetInput(greyErodeFilter->GetOutput());
    dilateFilter->SetKernel(greyOpeningKernel);

    if(m_saveEverything) {
        typename CScalarImageWriterType::Pointer writerDilate = CScalarImageWriterType::New();
        writerDilate->ReleaseDataFlagOn();
        writerDilate->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step0c_dilate" + m_fileExtensionChannelDAPI);
        writerDilate->SetInput(dilateFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerDilate->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerDilate->Update();
    }

    //----------FILTER-THRESHOLD-ON-DAPI-CHANNEL-------------------------------------------------------------------------------
    switch(m_thresholdingMode)
    {
    case 0:
        {
        adapOtsuFilter = AdaptiveOtsuThresholdImageFilterType::New();
        adapOtsuFilter->SetInput(dilateFilter->GetOutput());
        adapOtsuFilter->SetInsideValue(255);
        adapOtsuFilter->SetOutsideValue(0);
        adapOtsuFilter->SetNumberOfHistogramBins(256);
        adapOtsuFilter->SetSplineOrder(3);
        adapOtsuFilter->SetNumberOfControlPoints(5);
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
        otsuFilter->SetInput(dilateFilter->GetOutput());
        otsuFilter->Update();
        m_otsuThreshold = (CScalarPixelType)(otsuFilter->GetThreshold());

        invertFilter = InvertIntensityImageFilterType::New();
        invertFilter->SetInput(otsuFilter->GetOutput());
        invertFilter->SetMaximum(255);
        break;
        }
    default:
        {
        thresDAPIFilter = ThresholdFilterType::New();
        thresDAPIFilter->SetOutsideValue(0);
        thresDAPIFilter->SetInsideValue(255);
        thresDAPIFilter->SetLowerThreshold(m_lowerThreshold);
        thresDAPIFilter->SetUpperThreshold(m_upperThreshold);
        thresDAPIFilter->SetInput(dilateFilter->GetOutput());
        break;
        }
    }

    if(m_saveEverything) {
        typename CScalarImageWriterType::Pointer writer1 = CScalarImageWriterType::New();
        writer1->ReleaseDataFlagOn();
        writer1->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step1_bin" + m_fileExtensionChannelDAPI);
        switch(m_thresholdingMode)
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
        writer1->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    typename HoleFillingImageFilterType::Pointer invHoleFillingFilter = HoleFillingImageFilterType::New();
    invHoleFillingFilter->SetRadius(m_invHoleFillingNeighborhoodRadius);
    invHoleFillingFilter->SetBackgroundValue(255);
    invHoleFillingFilter->SetForegroundValue(0);
    invHoleFillingFilter->SetMajorityThreshold(m_invHoleFillingMajThreshold);
    invHoleFillingFilter->SetMaximumNumberOfIterations(20);
    switch(m_thresholdingMode)
    {
    case 0:
        invHoleFillingFilter->SetInput(adapOtsuFilter->GetOutput());
        break;
    case 1:
        invHoleFillingFilter->SetInput(invertFilter->GetOutput());
        break;
    default:
        invHoleFillingFilter->SetInput(thresDAPIFilter->GetOutput());
        break;
    }

    if(m_saveEverything) {
        typename CScalarImageWriterType::Pointer writerIHF = CScalarImageWriterType::New();
        writerIHF->ReleaseDataFlagOn();
        writerIHF->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step2_invHole" + m_fileExtensionChannelDAPI);
        writerIHF->SetInput(invHoleFillingFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerIHF->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerIHF->Update();
    }

    if(m_cavitiyFillingRadius!=0) {
        //----------FILTER-HOLE-FILLING-ON-DAPI-CHANNEL----------------------------------------------------------------------------
        if(m_cavWithRescaling) {
            typename CScalarImageType::SizeType inputSize = readerDAPIImage->GetLargestPossibleRegion().GetSize();
            typename CScalarImageType::SpacingType inputSpacing = readerDAPIImage->GetSpacing();

            std::cout << "Input size: " << inputSize << std::endl;
            std::cout << "Input spacing: " << inputSpacing << std::endl;

            typename CScalarImageType::SizeType outputSize;
            outputSize[0] = std::ceil((double)inputSize[0]/2.);
            outputSize[1] = std::ceil((double)inputSize[1]/2.);
            if(ImageDimension==3)   outputSize[2] = std::ceil((double)inputSize[2]/2.);

            typename CScalarImageType::SpacingType outputSpacing;
            outputSpacing[0] = inputSpacing[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
            outputSpacing[1] = inputSpacing[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
            if(ImageDimension==3)   outputSpacing[2] = inputSpacing[2] * (static_cast<double>(inputSize[2]) / static_cast<double>(outputSize[2]));

            typename CScalarImageType::SpacingType scaledSpacing;
            scaledSpacing[0] = m_spacing[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
            scaledSpacing[1] = m_spacing[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
            if(ImageDimension==3)   scaledSpacing[2] = m_spacing[2] * (static_cast<double>(inputSize[2]) / static_cast<double>(outputSize[2]));

            std::cout << "Output size: " << outputSize << std::endl;
            std::cout << "Output spacing: " << outputSpacing << std::endl;
            std::cout << "Scaled spacing: [" << scaledSpacing << "]" << std::endl;

            extrapolator = NearestNeighborExtrapolatorImageFunctionType::New();
            interpolator = NearestNeighborInterpolatorImageFunctionType::New();

            downscalingFilter = ResampleImageFilterType::New();
            downscalingFilter->SetSize(outputSize);
            downscalingFilter->SetOutputSpacing(outputSpacing);
            downscalingFilter->SetExtrapolator(extrapolator);
            downscalingFilter->SetInterpolator(interpolator);
            downscalingFilter->ReleaseDataFlagOn();
            downscalingFilter->SetInput(invHoleFillingFilter->GetOutput());
            downscalingFilter->Update();
            std::cout << "spacing of downscaling filter " << downscalingFilter->GetOutput()->GetSpacing() << std::endl;

            cavityFillingFilter = CavityFillingImageFilter::New();
            cavityFillingFilter->SetRadius((int)std::ceil((double)m_cavitiyFillingRadius/2.));
            cavityFillingFilter->SetLowerThreshold(m_cavitiyFillingMinThreshold);
            cavityFillingFilter->SetUpperThreshold(m_cavitiyFillingMaxThreshold);
            cavityFillingFilter->SetSpacing(scaledSpacing);
            cavityFillingFilter->SetMaximumNumberOfIterations(6);
            cavityFillingFilter->SetInput(downscalingFilter->GetOutput());

            upscalingFilter = ResampleImageFilterType::New();
            upscalingFilter->SetSize(inputSize);
            upscalingFilter->SetOutputSpacing(inputSpacing);
            upscalingFilter->SetExtrapolator(extrapolator);
            upscalingFilter->SetInterpolator(interpolator);
            upscalingFilter->ReleaseDataFlagOn();
            upscalingFilter->SetInput(cavityFillingFilter->GetOutput());
            upscalingFilter->Update();
            std::cout << "spacing of downscaling filter " << upscalingFilter->GetOutput()->GetSpacing() << std::endl;
        }
        else {
            double *spacing = new double[ImageDimension];
            spacing[0] = m_spacing[0];
            spacing[1] = m_spacing[1];
            if(ImageDimension==3)   spacing[2] = m_spacing[2];

            cavityFillingFilter = CavityFillingImageFilter::New();
            cavityFillingFilter->SetRadius(m_cavitiyFillingRadius);
            cavityFillingFilter->SetLowerThreshold(m_cavitiyFillingMinThreshold);
            cavityFillingFilter->SetUpperThreshold(m_cavitiyFillingMaxThreshold);
            cavityFillingFilter->SetSpacing(spacing);
            cavityFillingFilter->SetMaximumNumberOfIterations(6);
            cavityFillingFilter->SetInput(invHoleFillingFilter->GetOutput());
        }

        if(m_saveEverything) {
            typename CScalarImageWriterType::Pointer writer2 = CScalarImageWriterType::New();
            writer2->ReleaseDataFlagOn();
            writer2->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step3_cavFill" + m_fileExtensionChannelDAPI);
            if(m_cavWithRescaling)  writer2->SetInput(upscalingFilter->GetOutput());
            else                    writer2->SetInput(cavityFillingFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
            writer2->SetImageIO( itk::TIFFImageIO::New() );
#endif
            writer2->Update();
        }
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-CLOSING-ON-DAPI-CHANNEL---------------------------------------------------------------------------------
    closingFilter = ClosingImageFilterType::New();
    closingFilter->ReleaseDataFlagOn();
    closingFilter->SetKernel(m_closingStructuringElement);
    closingFilter->SetForegroundValue(255);
    if(m_cavitiyFillingRadius!=0) {
        if(m_cavWithRescaling)  closingFilter->SetInput(upscalingFilter->GetOutput());
        else                    closingFilter->SetInput(cavityFillingFilter->GetOutput());
    }
    else
        closingFilter->SetInput(invHoleFillingFilter->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-OPENING-ON-DAPI-CHANNEL---------------------------------------------------------------------------------
    openingFilter = OpeningImageFilterType::New();
    openingFilter->ReleaseDataFlagOn();
    openingFilter->SetKernel(m_openingStructuringElement);
    openingFilter->SetBackgroundValue(0);
    openingFilter->SetForegroundValue(255);
    openingFilter->SetInput(closingFilter->GetOutput());
    openingFilter->Update();

    itk::SmartPointer<CScalarImageType> openingImage = openingFilter->GetOutput();
    openingImage->DisconnectPipeline();

    openingImage->SetSpacing(m_spacing);

    if(m_saveEverything) {
        typename CScalarImageWriterType::Pointer writer3 = CScalarImageWriterType::New();
        writer3->ReleaseDataFlagOn();
        writer3->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step4_clopening" + m_fileExtensionChannelDAPI);
        writer3->SetInput(openingImage);
#if (ITK_VERSION_MAJOR >= 4)
        writer3->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer3->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------DISTANCE-MAP---------------------------------------------------------------------------------------------------
    distanceMap = MaurerDistanceMapImageFilterType::New();
    distanceMap->SetUseImageSpacing(true);
    distanceMap->SquaredDistanceOff();
    distanceMap->SetBackgroundValue(255);
    distanceMap->SetInsideIsPositive(true);
    distanceMap->SetInput(openingImage);
    distanceMap->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-INVERT-IMAGE-FILTER-------------------------------------------------------------------------------------
    rescaler2 = RescaleImageFilterType::New();
    rescaler2->ReleaseDataFlagOn();
    rescaler2->SetInput(distanceMap->GetOutput());

    if(m_saveEverything) {
        typename SScalarWriterType::Pointer writer5 = SScalarWriterType::New();
        writer5->ReleaseDataFlagOn();
        writer5->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step5_distMap" + m_fileExtensionChannelDAPI);
        writer5->SetInput(rescaler2->GetOutput());
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
    maskImageFilter = MaskImageFilterType::New();
    maskImageFilter->ReleaseDataFlagOn();
    maskImageFilter->SetInput1(morphWatershed->GetOutput());
    maskImageFilter->SetInput2(openingImage);
    maskImageFilter->Update();

    itk::SmartPointer<IScalarImageType> maskImage = maskImageFilter->GetOutput();
    maskImage->DisconnectPipeline();

    maskImage->SetSpacing(m_spacing);

    openingImage->ReleaseData();
    openingImage = NULL;
    std::cout << "mask spacing (before feret diameter calculation) = " << maskImage->GetSpacing() << std::endl;
    //-------------------------------------------------------------------------------------------------------------------------

    //----------TO-LABEL-MAP-AND-BACK-----------------------------------------------------------------------------------------
    watershedImageToLabelMap = LabelImageToShapeLabelMapFilterType::New();
    if(ImageDimension == 3) watershedImageToLabelMap->ComputeFeretDiameterOn();
    watershedImageToLabelMap->SetInput(maskImage);
    watershedImageToLabelMap->Update();

    itk::SmartPointer<ShapeLabelMapType> hepNucleiLabelMap = watershedImageToLabelMap->GetOutput();
    hepNucleiLabelMap->DisconnectPipeline();
    std::cout << "hep nuclei label map spacing (after feret diameter calculation) = " << hepNucleiLabelMap->GetSpacing() << std::endl;
    std::cout << "after watersheding, before classification: hep nuclei " << hepNucleiLabelMap->GetNumberOfLabelObjects() << std::endl;

    watershedImageToLabelMap->Update();

    itk::SmartPointer<ShapeLabelMapType> nonHepNucleiLabelMap = watershedImageToLabelMap->GetOutput();
    nonHepNucleiLabelMap->DisconnectPipeline();
    std::cout << "non hep nuclei label map spacing (after feret diameter calculation) = " << nonHepNucleiLabelMap->GetSpacing() << std::endl;
    std::cout << "after watersheding, before classification: non-hep nuclei " << nonHepNucleiLabelMap->GetNumberOfLabelObjects() << std::endl;

    maskImage->ReleaseData();
    maskImage = NULL;

    std::vector<unsigned long> labelsToRemoveFromHepNucleiLabelMap, labelsToRemoveFromNonHepNucleiLabelMap;

    if(ImageDimension==2)
        ClassifyNuclei2D(hepNucleiLabelMap, labelsToRemoveFromHepNucleiLabelMap, labelsToRemoveFromNonHepNucleiLabelMap, false);
    if(ImageDimension==3)
        ClassifyNuclei3D(hepNucleiLabelMap, distanceMap->GetOutput(), labelsToRemoveFromHepNucleiLabelMap, labelsToRemoveFromNonHepNucleiLabelMap, false);

    std::cout << "labelsToRemoveFromHepNucleiLabelMap.size(): " << labelsToRemoveFromHepNucleiLabelMap.size() << std::endl;
    for(unsigned int i=0; i<labelsToRemoveFromHepNucleiLabelMap.size(); i++)
        hepNucleiLabelMap->RemoveLabel(labelsToRemoveFromHepNucleiLabelMap[i]);
    std::cout << "after watersheding: hep nuclei " << hepNucleiLabelMap->GetNumberOfLabelObjects() << std::endl;

    std::cout << "labelsToRemoveFromNonHepNucleiLabelMap.size(): " << labelsToRemoveFromNonHepNucleiLabelMap.size() << std::endl;
    for(unsigned int i=0; i<labelsToRemoveFromNonHepNucleiLabelMap.size(); i++)
        nonHepNucleiLabelMap->RemoveLabel(labelsToRemoveFromNonHepNucleiLabelMap[i]);
    std::cout << "after watersheding: non-hep nuclei " << nonHepNucleiLabelMap->GetNumberOfLabelObjects() << std::endl;

    typename LabelMapToLabelImageFilterType::Pointer nonHepNucleiLabelMapToImage = LabelMapToLabelImageFilterType::New();
    nonHepNucleiLabelMapToImage->SetInput(nonHepNucleiLabelMap);
//    nonHepNucleiLabelMapToImage->Update();

    typename LabelMapToLabelImageFilterType::Pointer hepNucleiLabelMapToImage = LabelMapToLabelImageFilterType::New();
    hepNucleiLabelMapToImage->SetInput(hepNucleiLabelMap);
//    hepNucleiLabelMapToImage->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
    typename ThresholdFilterType::Pointer thresholdSmallNucleiImageFilter = ThresholdFilterType::New();
    thresholdSmallNucleiImageFilter->ReleaseDataFlagOn();
    thresholdSmallNucleiImageFilter->SetOutsideValue(0);
    thresholdSmallNucleiImageFilter->SetInsideValue(255);
    thresholdSmallNucleiImageFilter->SetLowerThreshold(1);
    thresholdSmallNucleiImageFilter->SetUpperThreshold(255);
    thresholdSmallNucleiImageFilter->SetInput(nonHepNucleiLabelMapToImage->GetOutput());

    typename CScalarImageWriterType::Pointer writer6a = CScalarImageWriterType::New();
    writer6a->ReleaseDataFlagOn();
    writer6a->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[0] + m_fileExtensionChannelDAPI);
    writer6a->SetInput(thresholdSmallNucleiImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer6a->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer6a->Update();


    typename ThresholdFilterType::Pointer thresholdBigNucleiImageFilter = ThresholdFilterType::New();
    thresholdBigNucleiImageFilter->ReleaseDataFlagOn();
    thresholdBigNucleiImageFilter->SetOutsideValue(0);
    thresholdBigNucleiImageFilter->SetInsideValue(255);
    thresholdBigNucleiImageFilter->SetLowerThreshold(1);
    thresholdBigNucleiImageFilter->SetUpperThreshold(255);
    thresholdBigNucleiImageFilter->SetInput(hepNucleiLabelMapToImage->GetOutput());

    typename CScalarImageWriterType::Pointer writer6b = CScalarImageWriterType::New();
    writer6b->ReleaseDataFlagOn();
    writer6b->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[1] + m_fileExtensionChannelDAPI);
    writer6b->SetInput(thresholdBigNucleiImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer6b->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer6b->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------VISULALIZATION-ONLY-FILTER---NOT-NEEDED-BY-OTHER-PIPELINE-FILTERS----------------------------------------------
    typename LabelOverlayImageFilterType::Pointer smallNucleiOverlayImageFilter = LabelOverlayImageFilterType::New();
    smallNucleiOverlayImageFilter->SetOpacity(m_overlayOpacity);
    smallNucleiOverlayImageFilter->ReleaseDataFlagOn();
    smallNucleiOverlayImageFilter->SetInput(readerDAPIImage);
    smallNucleiOverlayImageFilter->SetLabelImage(nonHepNucleiLabelMapToImage->GetOutput());

    typename RGBWriterType::Pointer writer7a = RGBWriterType::New();
    writer7a->ReleaseDataFlagOn();
    writer7a->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[2] + m_fileExtensionChannelDAPI);
    writer7a->SetInput(smallNucleiOverlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer7a->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer7a->Update();

    typename LabelOverlayImageFilterType::Pointer bigNucleiOverlayImageFilter = LabelOverlayImageFilterType::New();
    bigNucleiOverlayImageFilter->SetOpacity(m_overlayOpacity);
    bigNucleiOverlayImageFilter->ReleaseDataFlagOn();
    bigNucleiOverlayImageFilter->SetInput(readerDAPIImage);
    bigNucleiOverlayImageFilter->SetLabelImage(hepNucleiLabelMapToImage->GetOutput());

    typename RGBWriterType::Pointer writer7b = RGBWriterType::New();
    writer7b->ReleaseDataFlagOn();
    writer7b->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[3] + m_fileExtensionChannelDAPI);
    writer7b->SetInput(bigNucleiOverlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer7b->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer7b->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    hepNucleiLabelMap->ReleaseData();
    hepNucleiLabelMap = NULL;

    nonHepNucleiLabelMap->ReleaseData();
    nonHepNucleiLabelMap = NULL;

    readerDAPIImage->ReleaseData();
    readerDAPIImage = NULL;

    maskVeinImage->ReleaseData();
    maskVeinImage = NULL;

    segCVImage->ReleaseData();
    segCVImage = NULL;

    segPVImage->ReleaseData();
    segPVImage = NULL;

    WriteLogFile(timeStamp);
    WriteDataSetSummary();
}
