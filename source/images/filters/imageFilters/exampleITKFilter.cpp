///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  exampleITKFilter.cpp                                                 //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2011-10-27 18:45:16                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  main code by Adrian Friebel <adrian.friebel@uni-leipzig.de>                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "exampleITKFilter.h"


ExampleITKFilter::ExampleITKFilter()
{
	mp_rgbReader = RGBReaderType::New();

	mp_redChannel	= RGBChannelImageAdaptorType::New();
	mp_greenChannel	= RGBChannelImageAdaptorType::New();
	mp_blueChannel 	= RGBChannelImageAdaptorType::New();

	mp_rescaleFilterRed 	= RescalIntensityFilterType::New();
	mp_rescaleFilterGreen	= RescalIntensityFilterType::New();
	mp_rescaleFilterBlue 	= RescalIntensityFilterType::New();

	mp_thresFilterRed 	= ThresholdFilterType::New();
	mp_thresFilterGreen = ThresholdFilterType::New();
	mp_thresFilterBlue 	= ThresholdFilterType::New();

	mp_scalarWriter = ScalarWriterType::New();

	filter();
}


ExampleITKFilter::~ExampleITKFilter()
{}


void
ExampleITKFilter::Update()
{
	std::cout << "Filter wird jetzt ausgeführt" << std::endl;
	mp_scalarWriter->Update();
}

void
ExampleITKFilter::filter()
{
	std::string input_path = "/home/ls1/friebel/Workspace/2011.09.13__cell_core_segmenter/Data/WholeSlideScans/";
	std::string input_filename = "test";
	std::string input_file_suffix = ".tiff";
	//----------READER--------------------------------------------------------------------------------------------------------
	mp_rgbReader->SetFileName(input_path+input_filename+input_file_suffix);
	//-------------------------------------------------------------------------------------------------------------------------


	//----------SET-UP-ADAPTORS-TO-PRESENT-CHANNELS-SEPARATELY-----------------------------------------------------------------
	mp_redChannel->SelectNthElement(0);
	mp_redChannel->SetImage(mp_rgbReader->GetOutput());

	mp_greenChannel->SelectNthElement(1);
	mp_greenChannel->SetImage(mp_rgbReader->GetOutput());

	mp_blueChannel->SelectNthElement(2);
	mp_blueChannel->SetImage(mp_rgbReader->GetOutput());
	//-------------------------------------------------------------------------------------------------------------------------


	//----------FILTER---RESCALE-INTENSITY-FILTER-SEPARATELY-ON-CHANNELS-------------------------------------------------------
	mp_rescaleFilterRed->SetOutputMinimum(0);
	mp_rescaleFilterRed->SetOutputMaximum(255);
	mp_rescaleFilterRed->SetInput(mp_redChannel);

	mp_rescaleFilterGreen->SetOutputMinimum(0);
	mp_rescaleFilterGreen->SetOutputMaximum(255);
	mp_rescaleFilterGreen->SetInput(mp_greenChannel);

	mp_rescaleFilterBlue->SetOutputMinimum(0);
	mp_rescaleFilterBlue->SetOutputMaximum(255);
	mp_rescaleFilterBlue->SetInput(mp_blueChannel);
	//-------------------------------------------------------------------------------------------------------------------------


	//----------FILTER---THRESHOLD-CHANNELS-SEPARATELY-------------------------------------------------------------------------
	const PixelComponentType outsideValue = 0;
	const PixelComponentType insideValue  = 255;
	const PixelComponentType lowerThresholdRed 		= 60;
	const PixelComponentType upperThresholdRed 		= 155;
	const PixelComponentType lowerThresholdGreen 	= 0;
	const PixelComponentType upperThresholdGreen 	= 255;
	const PixelComponentType lowerThresholdBlue 	= 110;
	const PixelComponentType upperThresholdBlue 	= 175;


	mp_thresFilterRed->SetOutsideValue(outsideValue);
	mp_thresFilterRed->SetInsideValue(insideValue);
	mp_thresFilterRed->SetLowerThreshold(lowerThresholdRed);
	mp_thresFilterRed->SetUpperThreshold(upperThresholdRed);
	mp_thresFilterRed->SetNumberOfThreads(2);
	mp_thresFilterRed->SetInput(mp_rescaleFilterRed->GetOutput());

	mp_thresFilterGreen->SetOutsideValue(0);											//no green at all
	mp_thresFilterGreen->SetInsideValue(0);
	mp_thresFilterGreen->SetLowerThreshold(lowerThresholdGreen);
	mp_thresFilterGreen->SetUpperThreshold(upperThresholdGreen);
	mp_thresFilterGreen->SetNumberOfThreads(2);
	mp_thresFilterGreen->SetInput(mp_rescaleFilterGreen->GetOutput());

	mp_thresFilterBlue->SetOutsideValue(outsideValue);
	mp_thresFilterBlue->SetInsideValue(insideValue);
	mp_thresFilterBlue->SetLowerThreshold(lowerThresholdBlue);
	mp_thresFilterBlue->SetUpperThreshold(upperThresholdBlue);
	mp_thresFilterBlue->SetNumberOfThreads(2);
	mp_thresFilterBlue->SetInput(mp_rescaleFilterBlue->GetOutput());
	//-------------------------------------------------------------------------------------------------------------------------

	//----------WRITER---------------------------------------------------------------------------------------------------------
	std::string output_path 			= input_path;
	std::string output_filename_prefix 	= "CS";
	std::string output_filename 		= input_filename;
	std::string output_filename_suffix 	= "thresBlue";
	std::string output_file_suffix 		= input_file_suffix;

	mp_scalarWriter->SetFileName(output_path+output_filename_prefix+output_filename+output_filename_suffix+output_file_suffix);
	mp_scalarWriter->SetInput(mp_thresFilterBlue->GetOutput());
	//-------------------------------------------------------------------------------------------------------------------------

	//----------VTK-INTERFACE--------------------------------------------------------------------------------------------------
//	ITKVTKConnectorType::Pointer itkvtkConnectorOriginalImage = ITKVTKConnectorType::New();
//	itkvtkConnectorOriginalImage->SetInput(sliceFilter->GetOutput());
//	itkvtkConnectorOriginalImage->ReleaseDataFlagOn();
//
//	vtkSmartPointer<vtkImageFlip> flipYFilter = vtkSmartPointer<vtkImageFlip>::New();
//	flipYFilter->SetFilteredAxis(1); 											// flip y axis
//	flipYFilter->SetInput(itkvtkConnectorOriginalImage->GetOutput());
//	flipYFilter->ReleaseDataFlagOn();
//
//	vtkSmartPointer<vtkImageViewer> image_view = vtkSmartPointer<vtkImageViewer>::New();
//	image_view->SetInput(flipYFilter->GetOutput());
//
//	widget.SetRenderWindow(image_view->GetRenderWindow());
//	image_view->SetupInteractor(widget.GetRenderWindow()->GetInteractor());
//	image_view->SetSize(300,100);
//
//	image_view->GetRenderer()->GetRenderWindow()->Render();
//
//	widget.show();
//	app.exec();
}
