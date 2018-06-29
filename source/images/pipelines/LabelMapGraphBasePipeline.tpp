///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelMapGraphBasePipeline.tpp                                        //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-04-01                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "LabelMapGraphBasePipeline.h"

#include <itkImage.h>
#include <itkNrrdImageIO.h>

#include <vtkGraphReader.h>
#include <vtkGraphWriter.h>



template< unsigned int VImageDimension > LabelMapGraphBasePipeline< VImageDimension >::LabelMapGraphBasePipeline()
{
    mObjectGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
}


template< unsigned int VImageDimension > LabelMapGraphBasePipeline< VImageDimension >::~LabelMapGraphBasePipeline()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void LabelMapGraphBasePipeline< VImageDimension >::ReadOriginalImage(std::string filename)
{
    typename CRGBReaderType::Pointer reader = CRGBReaderType::New();
    reader->SetFileName(filename.c_str());
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    reader->Update();

    mOriginalImage = reader->GetOutput();
    mOriginalImage->DisconnectPipeline();
    mOriginalImage->SetSpacing(mVoxelSpacing);
}


template< unsigned int VImageDimension > void LabelMapGraphBasePipeline< VImageDimension >::ReadObjectImage(std::string filename)
{
    typename IScalarReaderType::Pointer objectReader = IScalarReaderType::New();
    objectReader->SetFileName(filename.c_str());
#if (ITK_VERSION_MAJOR >= 4)
    objectReader->SetImageIO( itk::NrrdImageIO::New() );
#endif
    objectReader->Update();

    typename IMinMaxCalculatorType::Pointer minMaxCalc1 = IMinMaxCalculatorType::New();
    minMaxCalc1->SetImage(objectReader->GetOutput());
    minMaxCalc1->Compute();

    typename CastILImageFilterType::Pointer castFilter = CastILImageFilterType::New();
    castFilter->SetInput(objectReader->GetOutput());
    castFilter->Update();

    mObjectImage = castFilter->GetOutput();
    mObjectImage->DisconnectPipeline();
    mObjectImage->SetSpacing(mVoxelSpacing);

    //TODO: get rid of this duplication
    castFilter->Update();

    mObjectImage2 = castFilter->GetOutput();
    mObjectImage2->DisconnectPipeline();
    mObjectImage2->SetSpacing(mVoxelSpacing);

    typename LMinMaxCalculatorType::Pointer minMaxCalc2 = LMinMaxCalculatorType::New();
    minMaxCalc2->SetImage(mObjectImage);
    minMaxCalc2->Compute();

    if(minMaxCalc1->GetMaximum() != minMaxCalc2->GetMaximum()) {
        std::cout << "Error: lost superpixels during casting operation" << std::endl;
        return;
    }
}


template< unsigned int VImageDimension > void LabelMapGraphBasePipeline< VImageDimension >::ReadObjectGraph(std::string filename)
{
    vtkSmartPointer<vtkGraphReader> graphReader = vtkSmartPointer<vtkGraphReader>::New();
    graphReader->SetFileName(filename.c_str());
    graphReader->Update();
    graphReader->GetOutput()->ToUndirectedGraph(mObjectGraph);
}


template< unsigned int VImageDimension > void LabelMapGraphBasePipeline< VImageDimension >::SaveObjectGraph(std::string filename)
{
    vtkSmartPointer<vtkGraphWriter> graphWriter = vtkSmartPointer<vtkGraphWriter>::New();
    graphWriter->SetFileName(filename.c_str());
    graphWriter->SetInput(mObjectGraph);
    graphWriter->Update();
}


template< unsigned int VImageDimension > void LabelMapGraphBasePipeline< VImageDimension >::SaveObjectImage(std::string filename)
{
    typename LMinMaxCalculatorType::Pointer minMaxCalc = LMinMaxCalculatorType::New();
    minMaxCalc->SetImage(mObjectImage);
    minMaxCalc->Compute();

    if(std::numeric_limits<unsigned int>::max() < minMaxCalc->GetMaximum()) {
        //TODO some neat error handling here
        std::cout << "SuperpixelPreparation: Can't save superpixels. LabelIDs exceed maximal unsigned int maximum." << std::endl;
        return;
    }

    typename CastLIImageFilterType::Pointer castFilter = CastLIImageFilterType::New();
    castFilter->SetInput(mObjectImage);

    typename IScalarWriterType::Pointer labelWriter = IScalarWriterType::New();
    labelWriter->SetFileName(filename);
    labelWriter->SetInput(castFilter->GetOutput());
    labelWriter->SetImageIO( itk::NrrdImageIO::New() );
    labelWriter->Update();
}



template< unsigned int VImageDimension > void LabelMapGraphBasePipeline< VImageDimension >::SaveObjectContourOverlayImage(std::string filename)
{
    typename LabelContourImageFilterType::Pointer superpixelContour = LabelContourImageFilterType::New();
    superpixelContour->SetFullyConnected(false);
    superpixelContour->SetInput(mObjectImage);

    CRGBPixelType pxl;
    pxl.Fill(100);

    typename MaskNegatedFilterType::Pointer maskFilter = MaskNegatedFilterType::New();
    maskFilter->SetInput(mOriginalImage);
    maskFilter->SetMaskImage(superpixelContour->GetOutput());
    maskFilter->SetOutsideValue(pxl);

    typename CRGBWriterType::Pointer outlineWriter = CRGBWriterType::New();
    outlineWriter->SetFileName(filename);
    outlineWriter->SetInput(maskFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    outlineWriter->SetImageIO( itk::TIFFImageIO::New() );
#endif
    outlineWriter->Update();
}


template< unsigned int VImageDimension > void LabelMapGraphBasePipeline< VImageDimension >::SaveObjectOverlayImage(std::string filename)
{
    typename BlueAdaptorType::Pointer BlueAdaptor = BlueAdaptorType::New();
    BlueAdaptor->SetImage(mOriginalImage);

    typename LabelOverlayImageFilterType::Pointer overlayImage = LabelOverlayImageFilterType::New();
    overlayImage->ReleaseDataFlagOn();
    overlayImage->SetOpacity(0.35);
    overlayImage->SetInput(BlueAdaptor);
    overlayImage->SetLabelImage(mObjectImage);

    typename CRGBWriterType::Pointer overlayWriter = CRGBWriterType::New();
    overlayWriter->SetFileName(filename);
    overlayWriter->SetInput(overlayImage->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    overlayWriter->SetImageIO( itk::TIFFImageIO::New() );
#endif
    overlayWriter->Update();
}


template< unsigned int VImageDimension > void LabelMapGraphBasePipeline< VImageDimension >::InitLabelMap()
{
    std::cout << "Init LabelMapGraph" << std::endl;

    typename LabelImageToLabelMapGraphFilterType::Pointer labelImageToLabelMapGraph = LabelImageToLabelMapGraphFilterType::New();           //Build label map from slic label image
    labelImageToLabelMapGraph->SetOriginalImage(mOriginalImage);
    labelImageToLabelMapGraph->SetInput(mObjectImage);
    labelImageToLabelMapGraph->SetLabelImage(mObjectImage2);
    labelImageToLabelMapGraph->UsePrecalculatedGraph(mObjectGraph);
    labelImageToLabelMapGraph->SetBackgroundValue(0);
    labelImageToLabelMapGraph->Update();

    mLabelMap = labelImageToLabelMapGraph->GetOutput();
    mLabelMap->DisconnectPipeline();

    std::cout << "label objects in label map graph = " << mLabelMap->GetNumberOfLabelObjects() << std::endl;
    std::cout << "-----number label objects/vertices = " << mLabelMap->GetLabelMapGraph()->GetNumberOfVertices() << ", number label object pairs/edges = " << mLabelMap->GetLabelMapGraph()->GetNumberOfEdges() << std::endl;

}
