///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentNucleiOnDAPIWithHough.cpp                                     //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-23-04                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "SegmentNucleiOnDAPIWithHough.h"

#if (ITK_VERSION_MAJOR < 4)
#define NUMBER_OF_PIXELS SIZE
#endif

#include <time.h>
#include <QDateTime>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"


SegmentNucleiOnDAPIWithHough::SegmentNucleiOnDAPIWithHough()
{
    m_overlayOpacity = 0.5;

    m_saveSuffixesForFinals[0] = "_step1_hough_bin";
    m_saveSuffixesForFinals[1] = "_step1_hough_overlay";
}


SegmentNucleiOnDAPIWithHough::~SegmentNucleiOnDAPIWithHough()
{
    // TODO Auto-generated destructor stub
}


void SegmentNucleiOnDAPIWithHough::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

    file.open((m_pathChannelDAPI + m_logFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-segment-nuclei-with-hough--------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-segment-nuclei-with-hough-------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


void SegmentNucleiOnDAPIWithHough::WriteDataSetSummary()
{
    ImageAnalysisSummaryFileIO::AddEntry(HepNucleiSegmentationBin, m_pathChannelDAPI, m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[0] + m_fileExtensionChannelDAPI);
    ImageAnalysisSummaryFileIO::AddEntry(HepNucleiSegmentationOverlay, m_pathChannelDAPI, m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[1] + m_fileExtensionChannelDAPI);
}


void SegmentNucleiOnDAPIWithHough::ParseParameterContext()
{
    if(m_paramContext->findContext("Segment Nuclei using Hough Transformation",0)==NULL) {
        std::cout << "Error: SegmentNucleiWithHoughContext: Invalid parameter context" << std::endl;
        return;
    }

    m_fullFilenameChannelDAPI = *(std::string*)(m_paramContext->findParameter("DAPI channel", 0)->dataPointer());

    m_spacing[0] = *(double*)(m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    m_spacing[1] = *(double*)(m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    m_spacing[2] = *(double*)(m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    m_medianRadius = ( *(int*)(m_paramContext->findContext("1.1) Preprocess DAPI Channel",0)->findParameter("Median filter radius", 0)->dataPointer()) );
    m_greyscaleOpeningRadius = ( *(int*)(m_paramContext->findContext("1.1) Preprocess DAPI Channel",0)->findParameter("Greyscale opening radius", 0)->dataPointer()) );

    m_houghNumSpheres = ( *(int*)(m_paramContext->findParameter("Maximal number of spheres", 0)->dataPointer()) );
    m_houghThreshold = ( *(int*)(m_paramContext->findParameter("Threshold on input image", 0)->dataPointer()) );
    m_houghGradientThreshold = ( *(int*)(m_paramContext->findParameter("Gradient Threshold", 0)->dataPointer()) );
    m_houghOutputThreshold = ( *(double*)(m_paramContext->findParameter("Output Threshold", 0)->dataPointer()) );
    m_houghMinRadius = ( *(double*)(m_paramContext->findParameter("Minimal Radius", 0)->dataPointer()) );
    m_houghMaxRadius = ( *(double*)(m_paramContext->findParameter("Maximal Radius", 0)->dataPointer()) );
    m_houghSigmaGradient = ( *(double*)(m_paramContext->findParameter("Sigma Gradient (in microns)", 0)->dataPointer()) );
    m_houghVariance = ( *(double*)(m_paramContext->findParameter("Variance (in microns)", 0)->dataPointer()) );
    m_houghSphereRadiusRatio = ( *(double*)(m_paramContext->findParameter("Sphere Radius Ratio [0,1]", 0)->dataPointer()) );
    m_houghVotingRadiusRatio = ( *(double*)(m_paramContext->findParameter("Voting Radius Ratio [0,1]", 0)->dataPointer()) );
    m_houghSamplingRatio = ( *(double*)(m_paramContext->findParameter("Sampling Ratio [0,1]", 0)->dataPointer()) );

    m_sphereToImageResampleFactor = ( *(int*)(m_paramContext->findParameter("Resample Factor", 0)->dataPointer()) );
    m_sphereToImageDilationRadius = ( *(int*)(m_paramContext->findParameter("Dilation Radius", 0)->dataPointer()) );

//    double diameter = *(double*)(m_paramContext->findParameter("Smallest non-hepatocyte diameter", 0)->dataPointer());
//    m_minNonHepRadius = diameter/2.0;
//
//    diameter = *(double*)(m_paramContext->findParameter("Smallest hepatocyte diameter", 0)->dataPointer());
//    m_minHepRadius = diameter/2.0;
//
//    diameter = *(double*)(m_paramContext->findParameter("Biggest hepatocyte diameter", 0)->dataPointer());
//    m_maxHepRadius = diameter/2.0;

    std::string saveMode = ( (CSParameterChoice*)(m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        m_saveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        m_saveEverything = 0;
    m_logFilenameSave = *(std::string*)(m_paramContext->findParameter("Log file name", 0)->dataPointer());
    m_filenameSave = *(std::string*)(m_paramContext->findParameter("Save prefix", 0)->dataPointer());
}


void SegmentNucleiOnDAPIWithHough::Update()
{
    ParseParameterContext();

    bool f = FilenameParser::ParseFilename(m_fullFilenameChannelDAPI, m_pathChannelDAPI, m_filenameChannelDAPI, m_fileExtensionChannelDAPI);

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    if(!f) {
        std::cout << "Error: SegmentNucleiOnDAPIWithHough: Could not execute pipeline, because input is invalid: " << m_fullFilenameChannelDAPI << std::endl;
        return;
    }

    std::cout << "Segment nuclei: " << std::endl;
    std::cout << " dir: " << m_pathChannelDAPI << std::endl;
    std::cout << " file: " << m_filenameChannelDAPI << std::endl;
    std::cout << " ext: " << m_fileExtensionChannelDAPI << std::endl;


    //----------READER---------------------------------------------------------------------------------------------------------
    ScalarVoReaderType::Pointer                 reader;
    HoughTransformFilterType::Pointer           houghFilter;
    HoughTransformFilterType::SpheresListType   circles;
    SpheresListType::const_iterator             itSpheres;
    SpatialObjectToImageFilterType::Pointer     drawSphereFilter;
    GroupType::Pointer                          group;
    LabelOverlayImageFilterType::Pointer        labelOverlayFilter;

    LabelImageToShapeLabelMapFilterType::Pointer    sphereImageToLabelMap;

    RGBVoWriterType::Pointer                    rgbWriter;
    ScalarVoWriterType::Pointer                 writer1;


    clock_t beginT, endT;
    beginT = clock();

    reader = ScalarVoReaderType::New();
    reader->SetFileName(m_pathChannelDAPI + m_filenameChannelDAPI + m_fileExtensionChannelDAPI);
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    reader->Update();

    CScalarVoImageType::Pointer localImage = reader->GetOutput();
    localImage->DisconnectPipeline();
    localImage->SetSpacing(m_spacing);
    //-------------------------------------------------------------------------------------------------------------------------

    //----------PREPROCESSING--------------------------------------------------------------------------------------------------
    MedianImageFilterType::InputSizeType medianRadius;
    medianRadius.Fill(m_medianRadius);

    MedianImageFilterType::Pointer medianFilter = MedianImageFilterType::New();
    medianFilter->SetRadius(medianRadius);
    medianFilter->SetInput(localImage);

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writerMedian1 = ScalarVoWriterType::New();
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

    GrayscaleErodeImageFilterType::Pointer erodeFilter = GrayscaleErodeImageFilterType::New();
    erodeFilter->SetInput(medianFilter->GetOutput());
    erodeFilter->SetKernel(greyOpeningKernel);

    GrayscaleDilateImageFilterType::Pointer dilateFilter = GrayscaleDilateImageFilterType::New();
    dilateFilter->SetInput(erodeFilter->GetOutput());
    dilateFilter->SetKernel(greyOpeningKernel);
    dilateFilter->Update();

    CScalarVoImageType::Pointer tempImage = dilateFilter->GetOutput();
    tempImage->DisconnectPipeline();
    tempImage->SetSpacing(m_spacing);

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writerDilate = ScalarVoWriterType::New();
        writerDilate->ReleaseDataFlagOn();
        writerDilate->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step0b_opening" + m_fileExtensionChannelDAPI);
        writerDilate->SetInput(tempImage);
#if (ITK_VERSION_MAJOR >= 4)
        writerDilate->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerDilate->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------HOUGH-TRANSFORMATION-------------------------------------------------------------------------------------------
    std::cout << "Computing Hough Map" << std::endl;

    houghFilter = HoughTransformFilterType::New();
    houghFilter->SetInput( tempImage );
    houghFilter->SetNumberOfSpheres( m_houghNumSpheres );
    houghFilter->SetMinimumRadius( m_houghMinRadius );
    houghFilter->SetMaximumRadius( m_houghMaxRadius );
    houghFilter->SetSigmaGradient( m_houghSigmaGradient );
    houghFilter->SetVariance( m_houghVariance );
    houghFilter->SetSphereRadiusRatio( m_houghSphereRadiusRatio );
    houghFilter->SetVotingRadiusRatio( m_houghVotingRadiusRatio );
    houghFilter->SetThreshold( m_houghThreshold );
    houghFilter->SetOutputThreshold( m_houghOutputThreshold );
    houghFilter->SetGradientThreshold( m_houghGradientThreshold );
    houghFilter->SetSamplingRatio( m_houghSamplingRatio );
    houghFilter->SetNbOfThreads( 8 );
    houghFilter->Update();

    circles = houghFilter->GetSpheres( );
    itSpheres = circles.begin();

    endT = clock();
    std::cout << "hough time " << ( ( endT - beginT )/CLOCKS_PER_SEC ) << std::endl;
    std::cout << "Found " << circles.size() << " circle(s)." << std::endl;
    //-------------------------------------------------------------------------------------------------------------------------

    //----------CONVERT-HOUGH-SPHERES-TO-IMAGE---------------------------------------------------------------------------------
    beginT = clock();

    // Computing the circles output
    CScalarVoImageType::RegionType region;
    region.SetSize( localImage->GetLargestPossibleRegion().GetSize() );
    region.SetIndex( localImage->GetLargestPossibleRegion().GetIndex() );

    CScalarVoImageType::Pointer  localOutputImage = CScalarVoImageType::New();
    localOutputImage->SetRegions( region );
    localOutputImage->SetOrigin(localImage->GetOrigin());
    localOutputImage->SetSpacing(localImage->GetSpacing());
    localOutputImage->Allocate();
    localOutputImage->FillBuffer(0);

    itk::FixedArray<double, Dimension>  rad;

    group = GroupType::New();
    drawSphereFilter = SpatialObjectToImageFilterType::New();

    unsigned int count = 1;
    while( itSpheres != circles.end() )
    {
        std::cout << "Center: ";
        std::cout << (*itSpheres)->GetObjectToParentTransform()->GetOffset() << std::endl;
        std::cout << "Radius: " << (*itSpheres)->GetRadius() << std::endl;

        rad[0] = (*itSpheres)->GetRadius()[0] / m_spacing[0];
        rad[1] = (*itSpheres)->GetRadius()[1] / m_spacing[1];
        rad[2] = (*itSpheres)->GetRadius()[2] / m_spacing[2];

        SphereType::Pointer sphere = SphereType::New();
        TransformType::Pointer transform = TransformType::New();
        TransformType::OutputVectorType translation;

        sphere->SetRadius(rad);
        sphere->SetId(count-1);

        transform->SetIdentity();

        translation[0] = (*itSpheres)->GetObjectToParentTransform()->GetOffset()[0];
        translation[1] = (*itSpheres)->GetObjectToParentTransform()->GetOffset()[1];
        translation[2] = (*itSpheres)->GetObjectToParentTransform()->GetOffset()[2];
        transform->Translate(translation, false);

        sphere->SetObjectToParentTransform(transform);

        group->AddSpatialObject(sphere);

        sphere->SetDefaultInsideValue(count);

        itSpheres++;
        count++;
    }

    drawSphereFilter->SetSize(localOutputImage->GetLargestPossibleRegion().GetSize());
    drawSphereFilter->SetUseObjectValue(true);
    drawSphereFilter->SetOutsideValue(0);
    drawSphereFilter->SetNumberOfThreads(8);
    drawSphereFilter->SetMaskResampleFactor(m_sphereToImageResampleFactor);
    drawSphereFilter->SetMaskDilationSize(m_sphereToImageDilationRadius);
    drawSphereFilter->SetInput(group);
    drawSphereFilter->Update();

    CScalarVoImageType::Pointer outputImage = drawSphereFilter->GetOutput();
    outputImage->DisconnectPipeline();
    outputImage->SetSpacing(m_spacing);

    endT = clock();
    std::cout << "draw sphere time " << ( ( endT - beginT )/CLOCKS_PER_SEC ) << std::endl;
    //-------------------------------------------------------------------------------------------------------------------------

    //----------LABEL-FUN------------------------------------------------------------------------------------------------------
    labelOverlayFilter = LabelOverlayImageFilterType::New();
    labelOverlayFilter->SetInput( localImage );
    labelOverlayFilter->SetLabelImage( outputImage );

    sphereImageToLabelMap = LabelImageToShapeLabelMapFilterType::New();
    sphereImageToLabelMap->ComputeFeretDiameterOn();
    sphereImageToLabelMap->ReleaseDataFlagOn();
    sphereImageToLabelMap->SetInput(outputImage);
    sphereImageToLabelMap->Update();

    itk::SmartPointer<LabelImageToShapeLabelMapFilterType::OutputImageType> labelMapImage2 = sphereImageToLabelMap->GetOutput();
    labelMapImage2->DisconnectPipeline();

    std::cout << "spacing = " << labelMapImage2->GetSpacing() << std::endl;
    std::cout << "after labeling: " << labelMapImage2->GetNumberOfLabelObjects() << " objects" << std::endl;

    for(unsigned int i=0; i<labelMapImage2->GetNumberOfLabelObjects(); i++) {
        std::cout << "Label object " << i << ":" << std::endl;

        LabelImageToShapeLabelMapFilterType::OutputImageType::IndexType index;

        PointSetPointType point = labelMapImage2->GetNthLabelObject(i)->GetCentroid();
        labelMapImage2->TransformPhysicalPointToIndex(point, index);

        std::cout << " - index of centroid : " << index[0] << ", " << index[1] << ", " << index[2] << std::endl;
        std::cout << " - ~diameter : " << 2.0 * pow((labelMapImage2->GetNthLabelObject(i)->GetPhysicalSize() * 3.0/4.0 / itk::Math::pi), 1.0/3.0) << std::endl;
        std::cout << " - feret diameter : " << labelMapImage2->GetNthLabelObject(i)->GetFeretDiameter() << std::endl;
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------WRITER---------------------------------------------------------------------------------------------------------
    rgbWriter = RGBVoWriterType::New();
    rgbWriter->SetInput( labelOverlayFilter->GetOutput() );
    rgbWriter->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[1] + m_fileExtensionChannelDAPI);
#if (ITK_VERSION_MAJOR >= 4)
    rgbWriter->SetImageIO( itk::TIFFImageIO::New() );
#endif
    rgbWriter->Update();

    writer1 = ScalarVoWriterType::New();
    writer1->SetFileName(m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[0] + m_fileExtensionChannelDAPI);
    writer1->SetInput( outputImage );
#if (ITK_VERSION_MAJOR >= 4)
    writer1->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer1->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    WriteLogFile(timeStamp);
    WriteDataSetSummary();
}
