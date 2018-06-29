///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  BasePipeline.tpp                                                     //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-06-10                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "BasePipeline.h"

#include <fstream>
#include <iostream>

#include "itkTIFFImageIO.h"


template< unsigned int VImageDimension > void BasePipeline< VImageDimension >::ReadImage(std::string filename, CScalarImagePointerType &image, CScalarSpacingType spacing)
{
    if(BasePipeline::GetNumberOfDimensions(filename) != ImageDimension)
        throw std::string("Data " + filename + " has different dimensionality than expected by pipeline.");

    typename ScalarVoReaderType::Pointer reader = ScalarVoReaderType::New();
    reader->SetFileName(filename);
    reader->ReleaseDataBeforeUpdateFlagOn();
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    reader->Update();

    image = reader->GetOutput();
    image->DisconnectPipeline();
    image->SetSpacing(spacing);
}


template< unsigned int VImageDimension > void BasePipeline< VImageDimension >::ReadImage(std::string filename, RGBCImagePointerType &image, CScalarSpacingType spacing)
{
    if(BasePipeline::GetNumberOfDimensions(filename) != ImageDimension)
        throw std::string("Data " + filename + " has different dimensionality than expected by pipeline.");

    typename RGBVoReaderType::Pointer reader = RGBVoReaderType::New();
    reader->SetFileName(filename);
    reader->ReleaseDataBeforeUpdateFlagOn();
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    reader->Update();

    image = reader->GetOutput();
    image->DisconnectPipeline();
    image->SetSpacing(spacing);
}


template< unsigned int VImageDimension > void BasePipeline< VImageDimension >::ReadImageLayer(std::string filename, CScalarImagePointerType &image, CScalarSpacingType spacing, CScalarPixelType foreground)
{
    ReadImage(filename, image, spacing);

    typename ThresholdFilterType1::Pointer thresholdFilter = ThresholdFilterType1::New();
    thresholdFilter->ReleaseDataFlagOn();
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->SetInsideValue(255);
    thresholdFilter->SetLowerThreshold(foreground);
    thresholdFilter->SetUpperThreshold(foreground);
    thresholdFilter->SetInput(image);
    thresholdFilter->Update();

    image = thresholdFilter->GetOutput();
    image->DisconnectPipeline();
    image->SetSpacing(spacing);
}


template< unsigned int VImageDimension > void BasePipeline< VImageDimension >::CreateImage(CScalarImagePointerType &image, CScalarRegionType region, CScalarPixelType constant, CScalarSpacingType spacing)
{
    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(constant);
    image->SetSpacing(spacing);
}


template< unsigned int VImageDimension > void BasePipeline< VImageDimension >::CreateImage(FScalarImagePointerType &image, CScalarRegionType region, FScalarPixelType constant, CScalarSpacingType spacing)
{
    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(constant);
    image->SetSpacing(spacing);
}


template< unsigned int VImageDimension > void BasePipeline< VImageDimension >::BuildDistanceMap(CScalarImagePointerType &inImage, FScalarImagePointerType &distMap, CScalarSpacingType spacing, bool distanceToSegmentation)
{
    inImage->SetSpacing(spacing);

    typename SignedMaurerDistanceMapImageFilterType::Pointer distanceMapFilter = SignedMaurerDistanceMapImageFilterType::New();
    distanceMapFilter->ReleaseDataFlagOn();
    distanceMapFilter->UseImageSpacingOn();
    distanceMapFilter->SquaredDistanceOff();
    if(distanceToSegmentation)
        distanceMapFilter->SetBackgroundValue(0);
    else
        distanceMapFilter->SetBackgroundValue(255);
    distanceMapFilter->SetInput(inImage);
    distanceMapFilter->Update();

    distMap = distanceMapFilter->GetOutput();
    distMap->DisconnectPipeline();
}


template< unsigned int VImageDimension > void BasePipeline< VImageDimension >::BuildNonZeroDistanceMap(CScalarImagePointerType &inImage, FScalarImagePointerType &distMap, FScalarPixelType constant, CScalarSpacingType spacing, bool distanceToSegmentation)
{
    typename CMinMaxCalculatorType::Pointer minMaxCalc = CMinMaxCalculatorType::New();
    minMaxCalc->SetImage(inImage);
    minMaxCalc->Compute();

    if(minMaxCalc->GetMaximum() != 0) {
        BuildDistanceMap(inImage, distMap, spacing, distanceToSegmentation);
    }
    else {
        distMap = FScalarVoImageType::New();
        distMap->SetRegions(inImage->GetLargestPossibleRegion());
        distMap->Allocate();
        distMap->FillBuffer(constant);
    }
}


template< unsigned int VImageDimension > void BasePipeline< VImageDimension >::BuildLabelImage(CScalarImagePointerType &inImage, LScalarImagePointerType &labelImage, CScalarSpacingType spacing, bool fullyConneceted)
{
    typename ImageToShapeLabelMapFilterType::Pointer imageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    imageToShaLabMapFilter->SetInput(inImage);
    imageToShaLabMapFilter->SetFullyConnected(fullyConneceted);

    typename LabelMapToLabelImageFilterType::Pointer labelMapToImage = LabelMapToLabelImageFilterType::New();
    labelMapToImage->SetInput(imageToShaLabMapFilter->GetOutput());
    labelMapToImage->Update();

    labelImage = labelMapToImage->GetOutput();
    labelImage->DisconnectPipeline();
    labelImage->SetSpacing(spacing);
}


template< unsigned int VImageDimension > void BasePipeline< VImageDimension >::LabelToBinImage(CScalarImagePointerType inLabelImage, unsigned char labelToExtract, CScalarImagePointerType &outBinImage)
{
    typename ThresholdFilterType1::Pointer thresholdFilter = ThresholdFilterType1::New();
    thresholdFilter->ReleaseDataFlagOn();
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->SetInsideValue(255);
    thresholdFilter->SetLowerThreshold(labelToExtract);
    thresholdFilter->SetUpperThreshold(labelToExtract);
    thresholdFilter->SetInput(inLabelImage);
    thresholdFilter->Update();

    outBinImage = thresholdFilter->GetOutput();
    outBinImage->DisconnectPipeline();
}


template< unsigned int VImageDimension > void BasePipeline< VImageDimension >::LabelToBinImage(LScalarImagePointerType inLabelImage, unsigned long labelToExtract, CScalarImagePointerType &outBinImage)
{
    typename ThresholdFilterType2::Pointer thresholdFilter = ThresholdFilterType2::New();
    thresholdFilter->ReleaseDataFlagOn();
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->SetInsideValue(255);
    thresholdFilter->SetLowerThreshold(labelToExtract);
    thresholdFilter->SetUpperThreshold(labelToExtract);
    thresholdFilter->SetInput(inLabelImage);
    thresholdFilter->Update();

    outBinImage = thresholdFilter->GetOutput();
    outBinImage->DisconnectPipeline();
}


template< unsigned int VImageDimension > unsigned int BasePipeline< VImageDimension >::GetNumberOfDimensions(std::string file)
{
    unsigned int dim;

    itk::TIFFImageIO::Pointer imageIO = itk::TIFFImageIO::New();
    imageIO->SetFileName(file);

    if(imageIO->CanReadFile(file.c_str())) {
        imageIO->ReadImageInformation();
        dim = imageIO->GetNumberOfDimensions();
    }
    else
        dim = 0;

    return dim;
}


template< unsigned int VImageDimension > long double BasePipeline< VImageDimension >::ComputeVolumeOfRegions(std::string filename, CScalarSpacingType spacing, long double voxelVolume, CScalarPixelType foreground)
{
    long double volume = 0;

    std::ifstream file(filename.c_str());
    if(file.good())
    {
        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::ImageIOFactory::ReadMode);
        imageIO->ReadImageInformation();

        typename CScalarVoImageType::Pointer image;
        ReadImage(filename, image, spacing);

        volume = ComputeVolumeOfRegions(image, voxelVolume, foreground);
    }
    return volume;
}


template< unsigned int VImageDimension > long double BasePipeline< VImageDimension >::ComputeVolumeOfRegions(CScalarImagePointerType &image, long double voxelVolume, CScalarPixelType foreground)
{
    long double volume = 0;

    typename ThresholdFilterType1::Pointer thresholdFilter = ThresholdFilterType1::New();
    thresholdFilter->ReleaseDataFlagOn();
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->SetInsideValue(255);
    thresholdFilter->SetLowerThreshold(foreground);
    thresholdFilter->SetUpperThreshold(foreground);
    thresholdFilter->SetInput(image);
    thresholdFilter->Update();

    typename ImageToShapeLabelMapFilterType::Pointer imageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    imageToShaLabMapFilter->SetInput(thresholdFilter->GetOutput());
    imageToShaLabMapFilter->Update();

    int numPixel = 0;

    int numLabObjects = imageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
    for(int i=0; i<numLabObjects; ++i) {
        imageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->Optimize();
        numPixel += imageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixels();
    }

    volume = numPixel * voxelVolume;

    return volume;
}
