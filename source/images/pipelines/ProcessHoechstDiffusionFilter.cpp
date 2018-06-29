///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ProcessHoechstDiffusionFilter.cpp                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-12-06                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "ProcessHoechstDiffusionFilter.h"

#include <fstream>
#include <string>
#include <sstream>

#include <QDateTime>

#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIOFactory.h>
#endif
#include <itkNrrdImageIO.h>

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"


ProcessHoechstDiffusionFilter::ProcessHoechstDiffusionFilter()
{
    mFrame = 0;

#if (ITK_VERSION_MAJOR >= 4)
    itk::TIFFImageIOFactory::RegisterOneFactory();
#endif
}


ProcessHoechstDiffusionFilter::~ProcessHoechstDiffusionFilter()
{
    // TODO Auto-generated destructor stub
}


void ProcessHoechstDiffusionFilter::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

    file.open((mPath + "log.txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-Hoechst Spheroid Processing----------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-processing-----------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


void ProcessHoechstDiffusionFilter::ParseParameterContext()
{
    if(m_paramContext->findContext("Analyze Hoechst Diffusion",0)==NULL) {
        std::cout << "Error: Analyze Hoechst Diffusion: Invalid parameter context" << std::endl;
        return;
    }

    mDatasetID = *(std::string*)(m_paramContext->findParameter("Dataset ID", 0)->dataPointer());

    mFullFilename = *(std::string*)(m_paramContext->findParameter("First dataset slide", 0)->dataPointer());

    mVoxSpacing[0] = *(double*)(m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    mVoxSpacing[1] = *(double*)(m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());

    mTimeStep = *(double*)(m_paramContext->findParameter("Time step", 0)->dataPointer());
}


void ProcessHoechstDiffusionFilter::BuildFilename()
{
    std::stringstream filename;

    std::cout << "mFilename 0 = " << mFilename << std::endl;
    mFilename = mFilename.substr(0, mFilename.length()-3);
    std::cout << "mFilename 1 = " << mFilename << std::endl;

    filename << mFilename;

    if(mFrame/10 == 0) {
        filename << 0 << 0 << mFrame;
    }
    else if(mFrame/10 < 10) {
        filename << 0 << mFrame;
    }
    else
        filename << mFrame;

    mFilename = filename.str();
    std::cout << "mFilename 2 = " << mFilename << std::endl;
}


void ProcessHoechstDiffusionFilter::WriteDistBlueValueTableFile()
{
    std::fstream file;
    std::stringstream filename;
    filename << mPath << "/analysis/" << "distBlue_Frame" << mFrame << ".txt";

    file.open(filename.str().c_str(), std::fstream::out);

    file.flags(std::fstream::left | std::fstream::scientific);
    file.width(15);
    file << "nr";
    file.width(15);
    file << "dist";
    file.width(15);
    file << "blue" << std::endl;

    for(unsigned int i=0; i<mDists.size(); i++) {
        file.width(15);
        file << i;
        file.width(15);
        file << mDists[i];
        file.width(15);
        file << (int)mBlueValues[i] << std::endl;
    }

    file.close();
}


void ProcessHoechstDiffusionFilter::Update()
{
    ParseParameterContext();

    bool f = FilenameParser::ParseFilename(mFullFilename, mPath, mFilename, mFileExtension);

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    if(!f) {
        std::cout << "Error: ProcessHoechstDiffusionFilter: Could not execute pipeline, because input is invalid: " << mFullFilename << std::endl;
        return;
    }

    std::cout << "Process Hoechst Diffusion Filter: " << std::endl;
    std::cout << " dir: " << mPath << std::endl;
    std::cout << " file: " << mFilename << std::endl;
    std::cout << " ext: " << mFileExtension << std::endl;

    CRGBReaderType::Pointer                         reader;
    CScalarWriterType::Pointer                      writer1, writer2, writer3, writer4;
    RGBVoWriterType::Pointer                        writer5;
    FScalarWriterType::Pointer                      writer6;
    BlueAdaptorType::Pointer                        blueAdaptor;
    RedAdaptorType::Pointer                         redAdaptor;
    BlueRescaleIntensityImageFilterType::Pointer    rescaleBlue;
    RedRescaleIntensityImageFilterType::Pointer     rescaleRed;
    ThresholdFilterType::Pointer                    thresholdFilter;
    ClosingImageFilterType::Pointer                 closing;
    ImageToShapeLabelMapFilterType::Pointer         imageToShapeLabMapFilter;
    ShapeOpeningLabelMapFilterType::Pointer         shapeOpeningLabMapFilter;
    LabelMapToLabelImageFilterType::Pointer         labMapToImageFilter;
    LabelOverlayImageFilterType::Pointer            overlayImageFilter;
    DistanceMapImageFilterType::Pointer             distanceMapFilter;


    mThreshold = 20;
    itk::Size<2> rad;
    rad[0] = 12;
    rad[1] = 12;
    mClosingStructuringElement.SetRadius(rad);
    mClosingStructuringElement.CreateStructuringElement();
    mMinimalSpheroidSize = 50000;

    while(true)
    {
        BuildFilename();

        std::cout << "complete path: " << mPath + mFilename + mFileExtension << std::endl;

        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO((mPath + mFilename + mFileExtension).c_str(), itk::ImageIOFactory::ReadMode);
        std::cout << "1";
        if(imageIO->CanReadFile((mPath + mFilename + mFileExtension).c_str()))
        {
            std::cout << "2";
            reader = CRGBReaderType::New();
            reader->SetFileName(mPath + mFilename + mFileExtension);
#if (ITK_VERSION_MAJOR >= 4)
            reader->SetImageIO( itk::TIFFImageIO::New() );
#endif
            reader->Update();

            CRGBImageType::Pointer originalImage = reader->GetOutput();
            originalImage->DisconnectPipeline();
            originalImage->SetSpacing(mVoxSpacing);

            blueAdaptor = BlueAdaptorType::New();                                   //Adaptor to present only blue channel of RGB input image
            blueAdaptor->SetImage(originalImage);

            redAdaptor = RedAdaptorType::New();                                   //Adaptor to present only blue channel of RGB input image
            redAdaptor->SetImage(originalImage);

            rescaleBlue = BlueRescaleIntensityImageFilterType::New();
            rescaleBlue->SetInput(blueAdaptor);
            rescaleBlue->SetOutputMinimum(0);
            rescaleBlue->SetOutputMaximum(255);

            rescaleRed = RedRescaleIntensityImageFilterType::New();
            rescaleRed->SetInput(redAdaptor);
            rescaleRed->SetOutputMinimum(0);
            rescaleRed->SetOutputMaximum(255);

            thresholdFilter = ThresholdFilterType::New();
            thresholdFilter->SetInput(rescaleRed->GetOutput());
            thresholdFilter->SetOutsideValue(0);
            thresholdFilter->SetInsideValue(255);
            thresholdFilter->SetLowerThreshold(mThreshold);
            thresholdFilter->SetUpperThreshold(255);

            closing = ClosingImageFilterType::New();
            closing->SetKernel(mClosingStructuringElement);
            closing->SetInput(thresholdFilter->GetOutput());

            imageToShapeLabMapFilter = ImageToShapeLabelMapFilterType::New();
            imageToShapeLabMapFilter->SetInput(closing->GetOutput());
            imageToShapeLabMapFilter->Update();

            //----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------
            shapeOpeningLabMapFilter = ShapeOpeningLabelMapFilterType::New();
            shapeOpeningLabMapFilter->SetLambda(mMinimalSpheroidSize);                         //attribute value
            shapeOpeningLabMapFilter->ReverseOrderingOff();                                     //removes objects with attribute smaller than lambda
            shapeOpeningLabMapFilter->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::NUMBER_OF_PIXELS);
            shapeOpeningLabMapFilter->SetInput(imageToShapeLabMapFilter->GetOutput());
            shapeOpeningLabMapFilter->Update();

            labMapToImageFilter = LabelMapToLabelImageFilterType::New();
            labMapToImageFilter->SetInput(shapeOpeningLabMapFilter->GetOutput());
            labMapToImageFilter->Update();
            //-------------------------------------------------------------------------------------------------------------------------

            overlayImageFilter = LabelOverlayImageFilterType::New();
            overlayImageFilter->ReleaseDataFlagOn();
            overlayImageFilter->SetLabelImage(labMapToImageFilter->GetOutput());
            overlayImageFilter->SetOpacity(0.5);
            overlayImageFilter->SetInput(rescaleRed->GetOutput());
            overlayImageFilter->Update();

            distanceMapFilter = DistanceMapImageFilterType::New();
            distanceMapFilter->ReleaseDataFlagOn();
            distanceMapFilter->UseImageSpacingOn();
            distanceMapFilter->SquaredDistanceOff();
            distanceMapFilter->InsideIsPositiveOn();
            distanceMapFilter->SetBackgroundValue(0);
            distanceMapFilter->SetInput(labMapToImageFilter->GetOutput());
            distanceMapFilter->Update();

            FScalarImageType::Pointer distMap = distanceMapFilter->GetOutput();
            distMap->DisconnectPipeline();

            mDists.clear();
            mBlueValues.clear();

            unsigned int numPxl = shapeOpeningLabMapFilter->GetOutput()->GetNthLabelObject(0)->GetNumberOfPixels();
            std::cout << "numPxl = " << numPxl << std::endl;
            for(unsigned int i=0; i<numPxl; i++) {
                mDists.push_back( distMap->GetPixel( shapeOpeningLabMapFilter->GetOutput()->GetNthLabelObject(0)->GetIndex(i) ) );
//                std::cout << "mDists[" << i << "] = " << mDists[mDists.size()-1] << std::endl;
                mBlueValues.push_back( blueAdaptor->GetPixel( shapeOpeningLabMapFilter->GetOutput()->GetNthLabelObject(0)->GetIndex(i) ) );
//                std::cout << "mBlueValues[" << i << "] = " << (int)mBlueValues[mBlueValues.size()-1] << std::endl;
            }

            WriteDistBlueValueTableFile();

//            writer1 = CScalarWriterType::New();
//            writer1->SetFileName(mPath + mFilename + "_blue_test" + mFileExtension);
//            writer1->SetInput(rescaleBlue->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//            writer1->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//            writer1->Update();
//
//            writer2 = CScalarWriterType::New();
//            writer2->SetFileName(mPath + mFilename + "_red_step0_test" + mFileExtension);
//            writer2->SetInput(rescaleRed->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//            writer2->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//            writer2->Update();
//
//            writer3 = CScalarWriterType::New();
//            writer3->SetFileName(mPath + mFilename + "_red_step1_thres" + mFileExtension);
//            writer3->SetInput(thresholdFilter->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//            writer3->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//            writer3->Update();
//
//            writer4 = CScalarWriterType::New();
//            writer4->SetFileName(mPath + mFilename + "_red_step2_closing" + mFileExtension);
//            writer4->SetInput(closing->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//            writer4->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//            writer4->Update();

            writer5 = RGBVoWriterType::New();
            writer5->ReleaseDataFlagOn();
            writer5->SetFileName(mPath + "/overlay/" + mFilename + "_red_step3_overlay" + mFileExtension);
            writer5->SetInput(overlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
            writer5->SetImageIO( itk::TIFFImageIO::New() );
#endif
            writer5->Update();

            mFrame++;

//            writer6 = FScalarWriterType::New();
//            writer6->SetFileName(mPath + mFilename + "_red_step4_distMap" + ".nrrd");
//            writer6->SetInput(distanceMapFilter->GetOutput());
//            writer6->SetImageIO( itk::NrrdImageIO::New() );
//            writer6->Update();
        }
        else
            break;
    }

    WriteLogFile(timeStamp);
}
