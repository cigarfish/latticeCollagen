///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Skeletonization.cpp                                                  //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-12-14                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "Skeletonization.h"

#include <QDateTime>

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"


SkeletonizationPipeline::SkeletonizationPipeline()
{
    // TODO Auto-generated constructor stub

}


SkeletonizationPipeline::~SkeletonizationPipeline()
{
    // TODO Auto-generated destructor stub
}


void SkeletonizationPipeline::WriteLogFile(std::string timeStamp)
{
    //no logging in here so far
    //all pipelines that need skeletonization implement it directly
}


void SkeletonizationPipeline::ParseParameterContext()
{
    if(m_paramContext->findContext("Skeletonization",0)==NULL) {
        std::cout << "Error: Skeletonization: Invalid parameter context: " << std::endl;
        return;
    }
    m_fullFilename = *(std::string*)(m_paramContext->findParameter("Image filename", 0)->dataPointer());
}


void SkeletonizationPipeline::Update()
{
    ParseParameterContext();

    bool f = FilenameParser::ParseFilename(m_fullFilename, m_path, m_filename, m_fileExtension);

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    if(!f) {
        std::cout << "Error: Skeletonization: Could not execute pipeline, because input is invalid: " << m_fullFilename << std::endl;
        return;
    }

    std::cout << "Skeletonization: " << std::endl;
    std::cout << " dir: " << m_path << std::endl;
    std::cout << " file: " << m_filename << std::endl;
    std::cout << " ext: " << m_fileExtension << std::endl;


    ScalarVoReaderType::Pointer         reader;
    ThinningImageFilterType::Pointer    thinningFilter;
    ScalarVoWriterType::Pointer         writer;


    reader = ScalarVoReaderType::New();
    reader->SetFileName(m_path + m_filename + m_fileExtension);
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif

    thinningFilter = ThinningImageFilterType::New();
    thinningFilter->ReleaseDataFlagOn();
    thinningFilter->SetInput(reader->GetOutput());

    writer = ScalarVoWriterType::New();
    writer->ReleaseDataFlagOn();
    writer->SetFileName(m_path + m_filename + "_skeleton" + m_fileExtension);
    writer->SetInput(thinningFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();


    //Code for outward flux skeletonization (from Insight Journal)
    //class could be extended with different skeltonization procedures

//    typedef float                                                               FPixelType;
//    typedef itk::Image<FPixelType, 2>                                           FScalarImageType;
//
//    //Test------------------------------------------------------------------------------------------------------------------------
//    double sigma = 0.01, threshold = 200000.0;
//
//    typedef unsigned char                                                           CPixelType;
//    typedef float                                                                   FPixelType;
//    typedef itk::Image<CPixelType, 2>                                               CScalarImageType;
//    typedef itk::Image<FPixelType, 2>                                               FScalarImageType;
//
//    typedef itk::ImageFileReader<FScalarImageType> ScalarReaderType;
//    ScalarReaderType::Pointer reader = ScalarReaderType::New();
//    reader->SetFileName(path + filename + fileExtension);
//#if (ITK_VERSION_MAJOR >= 4)
//    reader->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//    reader->Update();
//
////                typedef itk::DanielssonDistanceMapImageFilter<FScalarImageType, FScalarImageType> DanielssonDistanceMapImageFilterType;
////                DanielssonDistanceMapImageFilterType::Pointer distanceMap = DanielssonDistanceMapImageFilterType::New();
////                distanceMap->SetInput(reader->GetOutput());
//
//    itk::Vector<double, 2> spacing, spacing20x;
//    spacing20x[0] = 10; spacing20x[1] = 10;// spacing20x[2] = 10;
//
//    itk::SmartPointer<FScalarImageType> segImage = reader->GetOutput();
//    segImage->DisconnectPipeline();
//    segImage->SetSpacing(spacing20x);
//
//    typedef itk::SignedMaurerDistanceMapImageFilter<FScalarImageType, FScalarImageType> SignedMaurerDistanceMapImageFilterType;
//    SignedMaurerDistanceMapImageFilterType::Pointer distanceMap = SignedMaurerDistanceMapImageFilterType::New();
//    distanceMap->SquaredDistanceOff();
//    distanceMap->UseImageSpacingOn();
//    distanceMap->SetInput(segImage);
//
//    typedef itk::DiscreteGaussianImageFilter<FScalarImageType, FScalarImageType> DiscreteGaussianImageFilterType;
//    DiscreteGaussianImageFilterType::Pointer gaussianFilter = DiscreteGaussianImageFilterType::New();
//    gaussianFilter->SetVariance(sigma);
//    gaussianFilter->SetUseImageSpacingOn();
//    gaussianFilter->SetInput(distanceMap->GetOutput());
//
//    typedef itk::GradientImageFilter<FScalarImageType, FPixelType, FPixelType> GradientFilterType;
//    GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
//    gradientFilter->SetInput(gaussianFilter->GetOutput());
//
//    // Compute the average outward flux.
//    typedef itk::AverageOutwardFluxImageFilter<FScalarImageType, FPixelType, GradientFilterType::OutputImageType::PixelType> AOFFilterType;
//    AOFFilterType::Pointer aofFilter = AOFFilterType::New();
//    aofFilter->SetInput(gaussianFilter->GetOutput());
//    aofFilter->SetGradientImage(gradientFilter->GetOutput());
//
//    typedef itk::MedialCurveImageFilter<FScalarImageType> MedialCurveFilter;
//    MedialCurveFilter::Pointer medialFilter = MedialCurveFilter::New();
//    medialFilter->SetInput(gaussianFilter->GetOutput());
//    medialFilter->SetAverageOutwardFluxImage( aofFilter->GetOutput() );
//    medialFilter->SetThreshold( threshold );
//    //Test------------------------------------------------------------------------------------------------------------------------

    WriteLogFile(timeStamp);
}

