///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentAndClassifyStellateCells.h                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-01-14                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "SegmentAndClassifyStellateCells.h"

#if (ITK_VERSION_MAJOR < 4)
#define NUMBER_OF_PIXELS SIZE
#endif

#include <fstream>

#include <QDateTime>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



SegmentAndClassifyStellateCells::SegmentAndClassifyStellateCells()
{
    m_overlayOpacity = 0.4;
}


SegmentAndClassifyStellateCells::~SegmentAndClassifyStellateCells()
{
}


void SegmentAndClassifyStellateCells::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

    file.open((m_pathDesminChannel + m_logFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-segment-and-classify-stellate-cells-----------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-segment--and-classify-stellate-cells------------------------------------------------------------------------------------------------\n";

    file.close();
}


void SegmentAndClassifyStellateCells::WriteDataSetSummary()
{
    //todo ??
//    ImageAnalysisSummaryFileIO::AddEntry(SinusoidSegmentationBin, m_pathChannel1, m_pathChannel1 + m_filenameSave + m_finalSaveSuffixes[0] + m_fileExtensionChannel1);
}


void SegmentAndClassifyStellateCells::ParseParameterContext()
{
    if(m_paramContext->findContext("Segment Stellate Cell Cytoskeleton",0)==NULL) {
        std::cout << "Error: SegmentAndClassifyStellateCells: Invalid parameter context" << std::endl;
        return;
    }

    m_fullFilenameDesminChannel = *(std::string*)(m_paramContext->findParameter("Desmin channel", 0)->dataPointer());
    m_fullFilenameNonHepNucleiSeg = *(std::string*)(m_paramContext->findParameter("Non-Hepatic nuclei segmentation", 0)->dataPointer());
    m_fullFilenameHepNucleiSeg = *(std::string*)(m_paramContext->findParameter("Hepatic nuclei segmentation", 0)->dataPointer());

    m_spacing[0] = *(double*)(m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    m_spacing[1] = *(double*)(m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    m_spacing[2] = *(double*)(m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    m_medianRadius = ( *(int*)(m_paramContext->findContext("1.1) Preprocess Desmin Channel",0)->findParameter("Median filter radius", 0)->dataPointer()) );
    m_greyscaleOpeningRadius = *(int*)(m_paramContext->findContext("1.1) Preprocess Desmin Channel", 0)->findParameter("Greyscale opening radius", 0)->dataPointer());
    m_withConvexFilterPreprocessing = *(bool*)(m_paramContext->findContext("1.1) Preprocess Desmin Channel", 0)->findParameter("With convex image filter", 0)->dataPointer());
    m_convexFilterHeight = *(int*)(m_paramContext->findContext("1.1) Preprocess Desmin Channel", 0)->findParameter("Convex image filter height", 0)->dataPointer());

    std::string thresMode = ( (CSParameterChoice*)(m_paramContext->findParameter("Threshold mode", 0)->dataPointer()) )->currentString();
    if( thresMode.compare("Adaptive Otsu-Threshold")==0 )
        m_thresholdingMode = 0;
    else if( thresMode.compare("Otsu-Threshold")==0 )
        m_thresholdingMode = 1;
    else if( thresMode.compare("Manual Threshold")==0 )
        m_thresholdingMode = 2;

    m_adapOtsuRadius[0] = *(int*)(m_paramContext->findContext("1.2) Binary Threshold on 1.1", 0)->findParameter("Sample region size x", 0)->dataPointer());
    m_adapOtsuRadius[1] = *(int*)(m_paramContext->findContext("1.2) Binary Threshold on 1.1", 0)->findParameter("Sample region size y", 0)->dataPointer());
    m_adapOtsuRadius[2] = *(int*)(m_paramContext->findContext("1.2) Binary Threshold on 1.1", 0)->findParameter("Sample region size z", 0)->dataPointer());
    m_adapOtsuSamplePoints = *(int*)(m_paramContext->findContext("1.2) Binary Threshold on 1.1", 0)->findParameter("Number samples", 0)->dataPointer());
    m_lowerThreshold = *(int*)(m_paramContext->findContext("1.2) Binary Threshold on 1.1", 0)->findParameter("Desmin manual threshold", 0)->dataPointer());
    //TODO max;
    m_upperThreshold = 255;

    m_inverseHoleFillingNeighborhoodRadius[0] = *(int*)(m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2", 0)->findParameter("Radius x", 0)->dataPointer());
    m_inverseHoleFillingNeighborhoodRadius[1] = *(int*)(m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2", 0)->findParameter("Radius y", 0)->dataPointer());
    m_inverseHoleFillingNeighborhoodRadius[2] = *(int*)(m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2", 0)->findParameter("Radius z", 0)->dataPointer());
    m_inverseHoleFillingMajThreshold = *(int*)(m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2", 0)->findParameter("Majority threshold", 0)->dataPointer());

    m_minimalStellateCellSize = *(int*)(m_paramContext->findContext("1.4) Remove objects on 1.3", 0)->findParameter("Minimal structure size", 0)->dataPointer());

    m_nucleiBodyOverlap = *(double*)(m_paramContext->findContext("1.5) Classify stellate nuclei on 1.4", 0)->findParameter("Non-hepatic nuclei - stellate body superposition", 0)->dataPointer());

    std::string saveMode = ( (CSParameterChoice*)(m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        m_saveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        m_saveEverything = 0;
    m_logFilenameSave = *(std::string*)(m_paramContext->findParameter("Log file name", 0)->dataPointer());
    m_filenameSave = *(std::string*)(m_paramContext->findParameter("Save prefix", 0)->dataPointer());
}


void SegmentAndClassifyStellateCells::Update()
{
    ParseParameterContext();

    bool f1 = FilenameParser::ParseFilename(m_fullFilenameDesminChannel, m_pathDesminChannel, m_filenameDesminChannel, m_fileExtensionDesminChannel);
    bool f2 = FilenameParser::ParseFilename(m_fullFilenameNonHepNucleiSeg, m_pathNonHepNucleiSeg, m_filenameNonHepNucleiSeg, m_fileExtensionNonHepNucleiSeg);
    bool f3 = FilenameParser::ParseFilename(m_fullFilenameHepNucleiSeg, m_pathHepNucleiSeg, m_filenameHepNucleiSeg, m_fileExtensionHepNucleiSeg);

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    if(!f1 || !f2 || !f3) {
        std::cout << "f1 = " << f1 << std::endl;
        std::cout << "Error: SegmentAndClassifyStellateCells: Could not execute pipeline, because input is invalid" << std::endl;
        return;
    }

    std::cout << "Segment and classify stellate cells: " << std::endl;
    std::cout << " dir1: " << m_pathDesminChannel << std::endl;
    std::cout << " file1: " << m_filenameDesminChannel << std::endl;
    std::cout << " ext1: " << m_fileExtensionDesminChannel << std::endl;


    ScalarVoReaderType::Pointer readerDesmin;
    ScalarVoReaderType::Pointer readerNonHepNuclei;
    ScalarVoReaderType::Pointer readerHepNuclei;
    MedianImageFilterType::Pointer medianFilter;
    GreyscaleErodeImageFilterType::Pointer erodeFilter;
    GreyscaleDilateImageFilterType::Pointer dilateFilter;
    ConvexImageFilterType::Pointer convexFilter;
    AdaptiveOtsuThresholdImageFilterType::Pointer adapOtsuDesminFilter;
    OtsuThresholdImageFilterType::Pointer otsuFilter;
    InvertIntensityImageFilterType::Pointer invertFilter;
    ThresholdFilterType::Pointer thresDesminFilter;
    HoleFillingImageFilterType::Pointer holeFillingInverseDesminFilter;
    ImageToShapeLabelMapFilterType::Pointer imageToSCShaLabMapFilter;
    ShapeOpeningLabelMapFilterType::Pointer scShapeOpeningLabMapFilter;
    LabelMapToLabelImageFilterType::Pointer scLabMapToImageFilter;
    LabelOverlayImageFilterType::Pointer scOverlayImageFilter;
    ThresholdFilterType::Pointer scImageFilter;
    ImageToShapeLabelMapFilterType::Pointer nucleiImageToShaLabMapFilter;
    MaskImageFilterType::Pointer maskFilter;
    ImageToShapeLabelMapFilterType::Pointer maskedNucleiImageToShaLabMapFilter;
    LabelMapToLabelImageFilterType::Pointer registeredNucleiLabMapToImageFilter;
    LabelOverlayImageFilterType::Pointer registeredNucleiOverlayImageFilter;
    ThresholdFilterType::Pointer registeredNucleiImageFilter;

    //----------READER--------------------------------------------------------------------------------------------------------
    readerDesmin = ScalarVoReaderType::New();
    readerDesmin->SetFileName(m_pathDesminChannel + m_filenameDesminChannel + m_fileExtensionDesminChannel);
    readerDesmin->ReleaseDataBeforeUpdateFlagOn();
    readerDesmin->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerDesmin->SetImageIO( itk::TIFFImageIO::New() );
#endif

    readerNonHepNuclei = ScalarVoReaderType::New();
    readerNonHepNuclei->SetFileName(m_pathNonHepNucleiSeg + m_filenameNonHepNucleiSeg + m_fileExtensionNonHepNucleiSeg);
    readerNonHepNuclei->ReleaseDataBeforeUpdateFlagOn();
    readerNonHepNuclei->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerNonHepNuclei->SetImageIO( itk::TIFFImageIO::New() );
#endif

    readerHepNuclei = ScalarVoReaderType::New();
    readerHepNuclei->SetFileName(m_pathHepNucleiSeg + m_filenameHepNucleiSeg + m_fileExtensionHepNucleiSeg);
    readerHepNuclei->ReleaseDataBeforeUpdateFlagOn();
    readerHepNuclei->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerHepNuclei->SetImageIO( itk::TIFFImageIO::New() );
#endif
    //-------------------------------------------------------------------------------------------------------------------------

    //----------PREPROCESSING-OF-DESMIN-CHANNEL--------------------------------------------------------------------------------
    MedianImageFilterType::InputSizeType medianRadius;
    medianRadius.Fill(m_medianRadius);

    medianFilter = MedianImageFilterType::New();
    medianFilter->SetRadius(medianRadius);
    medianFilter->SetInput(readerDesmin->GetOutput());

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writerMedian = ScalarVoWriterType::New();
        writerMedian->ReleaseDataFlagOn();
        writerMedian->SetFileName(m_pathDesminChannel + m_filenameSave + "_step0_Desmin_median" + m_fileExtensionDesminChannel);
        writerMedian->SetInput(medianFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerMedian->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerMedian->Update();
    }

    StructuringElementType greyscaleOpeningKernel;
    greyscaleOpeningKernel.SetRadius(m_greyscaleOpeningRadius);
    greyscaleOpeningKernel.CreateStructuringElement();

    erodeFilter = GreyscaleErodeImageFilterType::New();
    erodeFilter->SetInput(medianFilter->GetOutput());
    erodeFilter->SetKernel(greyscaleOpeningKernel);

    dilateFilter = GreyscaleDilateImageFilterType::New();
    dilateFilter->SetInput(erodeFilter->GetOutput());
    dilateFilter->SetKernel(greyscaleOpeningKernel);

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writerOpening = ScalarVoWriterType::New();
        writerOpening->ReleaseDataFlagOn();
        writerOpening->SetFileName(m_pathDesminChannel + m_filenameSave + "_step0_Desmin_opening" + m_fileExtensionDesminChannel);
        writerOpening->SetInput(dilateFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerOpening->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerOpening->Update();
    }

    //----------FILTER-THRESHOLD-ON-DAPI-CHANNEL-------------------------------------------------------------------------------
    if(m_withConvexFilterPreprocessing) {
        convexFilter = ConvexImageFilterType::New();
        convexFilter->SetInput(dilateFilter->GetOutput());
        convexFilter->SetHeight(m_convexFilterHeight);
        convexFilter->SetFullyConnected(true);

        if(m_saveEverything) {
            ScalarVoWriterType::Pointer writerConvex = ScalarVoWriterType::New();
            writerConvex->ReleaseDataFlagOn();
            writerConvex->SetFileName(m_pathDesminChannel + m_filenameSave + "_step0_convex" + m_fileExtensionDesminChannel);
            writerConvex->SetInput(convexFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
            writerConvex->SetImageIO( itk::TIFFImageIO::New() );
#endif
            writerConvex->Update();
        }
    }

    switch(m_thresholdingMode)
    {
    case 0:
    {
        adapOtsuDesminFilter = AdaptiveOtsuThresholdImageFilterType::New();
        if(m_withConvexFilterPreprocessing) adapOtsuDesminFilter->SetInput(convexFilter->GetOutput());
        else                                adapOtsuDesminFilter->SetInput(dilateFilter->GetOutput());
        adapOtsuDesminFilter->SetInsideValue(255);
        adapOtsuDesminFilter->SetOutsideValue(0);
        adapOtsuDesminFilter->SetNumberOfHistogramBins(256);
        adapOtsuDesminFilter->SetSplineOrder(3);
        adapOtsuDesminFilter->SetNumberOfControlPoints(5);
        adapOtsuDesminFilter->SetNumberOfLevels(3);
        adapOtsuDesminFilter->SetNumberOfSamples(m_adapOtsuSamplePoints);
        adapOtsuDesminFilter->SetRadius(m_adapOtsuRadius);
        break;
    }
    case 1:
    {
        otsuFilter = OtsuThresholdImageFilterType::New();
        if(m_withConvexFilterPreprocessing) otsuFilter->SetInput(convexFilter->GetOutput());
        else                                otsuFilter->SetInput(dilateFilter->GetOutput());
        otsuFilter->Update();
        m_otsuThreshold = (CScalarPixelType)(otsuFilter->GetThreshold());

        invertFilter = InvertIntensityImageFilterType::New();
        invertFilter->SetInput(otsuFilter->GetOutput());
        invertFilter->SetMaximum(255);
        break;
    }
    default:
    {
        thresDesminFilter = ThresholdFilterType::New();
        if(m_withConvexFilterPreprocessing) thresDesminFilter->SetInput(convexFilter->GetOutput());
        else                                thresDesminFilter->SetInput(dilateFilter->GetOutput());
        thresDesminFilter->SetOutsideValue(0);
        thresDesminFilter->SetInsideValue(255);
        thresDesminFilter->SetLowerThreshold(m_lowerThreshold);
        thresDesminFilter->SetUpperThreshold(m_upperThreshold);
        break;
    }
    }

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer1 = ScalarVoWriterType::New();
        writer1->ReleaseDataFlagOn();
        writer1->SetFileName(m_pathDesminChannel + m_filenameSave + "_step1_bin" + m_fileExtensionDesminChannel);
        switch(m_thresholdingMode)
        {
        case 0:
            writer1->SetInput(adapOtsuDesminFilter->GetOutput());
            break;
        case 1:
            writer1->SetInput(invertFilter->GetOutput());
            break;
        default:
            writer1->SetInput(thresDesminFilter->GetOutput());
            break;
        }
#if (ITK_VERSION_MAJOR >= 4)
    writer1->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer1->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-HOLE-FILLING-ON-INVERSE-TO-REMOVE-NOISE-FROM-SUBSTRACTION-----------------------------------------------
    holeFillingInverseDesminFilter = HoleFillingImageFilterType::New();
    holeFillingInverseDesminFilter->ReleaseDataFlagOn();
    holeFillingInverseDesminFilter->SetRadius(m_inverseHoleFillingNeighborhoodRadius);
    holeFillingInverseDesminFilter->SetBackgroundValue(255);
    holeFillingInverseDesminFilter->SetForegroundValue(0);
    holeFillingInverseDesminFilter->SetMajorityThreshold(m_inverseHoleFillingMajThreshold);                   //number of foreground neighbors should be at least (3x3x3-1)/2 + majority
    switch(m_thresholdingMode)
    {
    case 0:
        holeFillingInverseDesminFilter->SetInput(adapOtsuDesminFilter->GetOutput());
        break;
    case 1:
        holeFillingInverseDesminFilter->SetInput(invertFilter->GetOutput());
        break;
    default:
        holeFillingInverseDesminFilter->SetInput(thresDesminFilter->GetOutput());
        break;
    }

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer4 = ScalarVoWriterType::New();
        writer4->ReleaseDataFlagOn();
        writer4->SetFileName(m_pathDesminChannel + m_filenameSave + "_step2_invHole" + m_fileExtensionDesminChannel);
        writer4->SetInput(holeFillingInverseDesminFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer4->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer4->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER---CREATE-LABEL-MAP-ON-Desmin-CHANNEL-DO-SOME-STUFF-AND-CONVERT-IT-BACK-INTO-IMAGE-----------------------
    imageToSCShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    imageToSCShaLabMapFilter->SetInput(holeFillingInverseDesminFilter->GetOutput());
    imageToSCShaLabMapFilter->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------------
    std::cout << "remove all objects with less than " << m_minimalStellateCellSize << " pixels " << std::endl;

    std::cout << "before removal " << imageToSCShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects" << std::endl;

    scShapeOpeningLabMapFilter = ShapeOpeningLabelMapFilterType::New();
    scShapeOpeningLabMapFilter->SetLambda(m_minimalStellateCellSize);                     //attribute value
    scShapeOpeningLabMapFilter->ReverseOrderingOff();                                     //removes objects with attribute smaller than lambda
    scShapeOpeningLabMapFilter->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::NUMBER_OF_PIXELS);
    scShapeOpeningLabMapFilter->SetInput(imageToSCShaLabMapFilter->GetOutput());
    scShapeOpeningLabMapFilter->Update();

    std::cout << "after removal " << scShapeOpeningLabMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects" << std::endl;

    scLabMapToImageFilter = LabelMapToLabelImageFilterType::New();
    scLabMapToImageFilter->SetInput(scShapeOpeningLabMapFilter->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------VISULALIZATION-ONLY-FILTER---NOT-NEEDED-BY-OTHER-PIPELINE-FILTERS----------------------------------------------
    scOverlayImageFilter = LabelOverlayImageFilterType::New();
    scOverlayImageFilter->SetLabelImage(scLabMapToImageFilter->GetOutput());
    scOverlayImageFilter->SetOpacity(m_overlayOpacity);
    scOverlayImageFilter->ReleaseDataFlagOn();
    scOverlayImageFilter->SetInput(readerDesmin->GetOutput());
    scOverlayImageFilter->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
    scImageFilter = ThresholdFilterType::New();
    scImageFilter->SetOutsideValue(0);
    scImageFilter->SetInsideValue(255);
    scImageFilter->SetLowerThreshold(1);
    scImageFilter->SetUpperThreshold(255);
    scImageFilter->SetInput(scLabMapToImageFilter->GetOutput());
    //--------------------------------------------------------------------------------------------------------------------------

    //----------WRITER----------------------------------------------------------------------------------------------------------
    RGBVoWriterType::Pointer writer6 = RGBVoWriterType::New();
    writer6->ReleaseDataFlagOn();
    writer6->SetFileName(m_pathDesminChannel + m_filenameSave + "_step3_overlay" + m_fileExtensionDesminChannel);
    writer6->SetInput(scOverlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer6->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer6->Update();

    ScalarVoWriterType::Pointer writer7 = ScalarVoWriterType::New();
    writer7->ReleaseDataFlagOn();
    writer7->SetFileName(m_pathDesminChannel + m_filenameSave + "_step3_bin" + m_fileExtensionDesminChannel);
    writer7->SetInput(scImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer7->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer7->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER---CREATE-LABEL-MAP-ON-Desmin-CHANNEL-DO-SOME-STUFF-AND-CONVERT-IT-BACK-INTO-IMAGE-----------------------
    nucleiImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    nucleiImageToShaLabMapFilter->SetInput(readerNonHepNuclei->GetOutput());
    nucleiImageToShaLabMapFilter->Update();

    maskFilter = MaskImageFilterType::New();
    maskFilter->SetInput(readerNonHepNuclei->GetOutput());
    maskFilter->SetMaskImage(scImageFilter->GetOutput());

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer5 = ScalarVoWriterType::New();
        writer5->ReleaseDataFlagOn();
        writer5->SetFileName(m_pathDesminChannel + "stellateNuclei" + "_step1_masked_nuclei" + m_fileExtensionDesminChannel);
        writer5->SetInput(maskFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer5->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer5->Update();
    }

    maskedNucleiImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    maskedNucleiImageToShaLabMapFilter->SetInput(maskFilter->GetOutput());
    maskedNucleiImageToShaLabMapFilter->Update();
    //-------------------------------------------------------------------------------------------------------------------------
    std::map<int, double> nucleusToCoveredFraction;
    std::vector<unsigned long> nucleiToRemove;

    for(unsigned int i=0; i<maskedNucleiImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
        for(unsigned int j=0; j<nucleiImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); j++) {
            if( nucleiImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(j)->HasIndex( maskedNucleiImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetIndex(0) ) ) {
                int nucleiPixel = nucleiImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(j)->GetNumberOfPixels();
                int maskedPixel = maskedNucleiImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixels();

                if(nucleusToCoveredFraction.count(j)==0)
                    nucleusToCoveredFraction[j] = (double)maskedPixel/(double)nucleiPixel;
                else
                    nucleusToCoveredFraction[j] += (double)maskedPixel/(double)nucleiPixel;
            }
        }
    }
    for(unsigned int i=0; i<nucleusToCoveredFraction.size(); i++) {
        std::cout << "nucleus " << i << " has coincidence fraction of " << nucleusToCoveredFraction[i] << std::endl;
        if(nucleusToCoveredFraction[i] < m_nucleiBodyOverlap)
            nucleiToRemove.push_back(nucleiImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetLabel());
    }

    for(unsigned int i=0; i<nucleiToRemove.size(); i++)
        nucleiImageToShaLabMapFilter->GetOutput()->RemoveLabel(nucleiToRemove[i]);

    registeredNucleiLabMapToImageFilter = LabelMapToLabelImageFilterType::New();
    registeredNucleiLabMapToImageFilter->SetInput(nucleiImageToShaLabMapFilter->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------VISULALIZATION-ONLY-FILTER---NOT-NEEDED-BY-OTHER-PIPELINE-FILTERS----------------------------------------------
    registeredNucleiOverlayImageFilter = LabelOverlayImageFilterType::New();
    registeredNucleiOverlayImageFilter->SetLabelImage(registeredNucleiLabMapToImageFilter->GetOutput());
    registeredNucleiOverlayImageFilter->SetOpacity(m_overlayOpacity);
    registeredNucleiOverlayImageFilter->ReleaseDataFlagOn();
    registeredNucleiOverlayImageFilter->SetInput(readerDesmin->GetOutput());
    registeredNucleiOverlayImageFilter->Update();

    //----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
    registeredNucleiImageFilter = ThresholdFilterType::New();
    registeredNucleiImageFilter->SetOutsideValue(0);
    registeredNucleiImageFilter->SetInsideValue(255);
    registeredNucleiImageFilter->SetLowerThreshold(1);
    registeredNucleiImageFilter->SetUpperThreshold(255);
    registeredNucleiImageFilter->SetInput(registeredNucleiLabMapToImageFilter->GetOutput());
    //--------------------------------------------------------------------------------------------------------------------------

    //----------WRITER----------------------------------------------------------------------------------------------------------
    RGBVoWriterType::Pointer writer8 = RGBVoWriterType::New();
    writer8->ReleaseDataFlagOn();
    writer8->SetFileName(m_pathDesminChannel + "stellateNuclei" + "_step2_overlay" + m_fileExtensionDesminChannel);
    writer8->SetInput(registeredNucleiOverlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer8->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer8->Update();

    ScalarVoWriterType::Pointer writer9 = ScalarVoWriterType::New();
    writer9->ReleaseDataFlagOn();
    writer9->SetFileName(m_pathDesminChannel + "stellateNuclei" + "_step2_bin" + m_fileExtensionDesminChannel);
    writer9->SetInput(registeredNucleiImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer9->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer9->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    WriteLogFile(timeStamp);
    WriteDataSetSummary();
}
