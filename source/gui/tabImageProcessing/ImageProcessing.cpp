///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ImageProcessing.cpp                                                  //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-11-09                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "ImageProcessing.h"

#include <QFileDialog>
#include <QLabel>
#include <QLCDNumber>
#include <QMessageBox>
#include <QModelIndex>
#include <QPushButton>
#include <QtGui/QWidget>

#include <itkImageToVTKImageFilter.h>
#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIO.h>
#endif
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4 || ITK_VERSION_MAJOR > 4)
#include <itkSCIFIOImageIO.h>
#endif

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkImageFlip.h>
#include <vtkPOVExporter.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include "../../Core.h"
#include "GraphViewer.h"
#include "ImageProcessingWorker.h"
#include "JobManagerModel.h"
#include "../SeedPicker.h"
#include "../tools/ImageAnalysisSummaryFileIO.h"
#include "../QCSParameterModel.h"
#include "../QCSParameterDelegate.h"
#include "../../images/filters/graphFilters/AnalyzeBileNetworkFilter.h"
#include "../../images/pipelines/ExtractGraph.h"
#include "../../images/pipelines/ObjectBasedSegmentation.h"
#include "../../images/tools/FileFormatConverter.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



const std::string ImageProcessing::PipelineName[] =
{
    "Analyze Cells",
    "Approximate Cell Shape",
    "CLAHE",
    "Crop",
    "Extract and Analyze Graph",
    "Segment Necrotic Region",
    "Segment Nuclei in 20x Datasets",
    "Segment Nuclei in 60x Datasets",
    "Segment Sinusoids + Bile Canaliculi in 20x Datasets",
    "Segment Sinusoids + Bile Canaliculi in 60x Datasets",
    "Segment Veins",
#ifndef CS_TI_QUANT_ONLY
    "Add Images",
    "Analyze Lobule",
    "Analyze Nuclei",
    "Approximate Lobule Shape",
    "Background Elimination",
    "HConvex Image Filter",
    "Analyze Hoechst Diffusion",
    "Segment Nuclei using Hough Transformation",
    "Visualize 3D",
    "Superpixel Preparation",
    "Object-based Segmentation",
    "Train Classifier",
    "Segment Cell Membrane",
    "Objects to Matrix",
    "Compare Segmentations",
    "Segment Stellate Cell Cytoskeleton",
    "Analyze Stellate Cells",
    "Skeletonization",
    "File Format Converter",
    "Test Pipeline",
#endif
    "" // be sure this empty string is the last entry
};

// Constructor
ImageProcessing::ImageProcessing(QWidget *parent)
	: QWidget(parent),
	  m_isBusy(false)
{
#ifdef CS_BUILD_PYTHON
	Core *core;
	core->initPythonInterpreter();
#endif

	// Setting up the Qt Designer code
	setupUi(this);

	m_testOutput = true;

	jobManager = new JobManager();
    jobTableView->setModel( new JobManagerModel( jobManager ) );
    jobTableView->setSelectionBehavior(QAbstractItemView::SelectRows);

	connect( selectPipelineComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(selectionPipelineChanged(int)) );
	connect( startPipelineButton, SIGNAL(clicked(bool)), this, SLOT(startPipelineButtonClicked()) );
	connect( loadJobsButton, SIGNAL(clicked(bool)), this, SLOT(loadJobsButtonClicked()) );
	connect( addJobButton, SIGNAL(clicked(bool)), this, SLOT(addJobButtonClicked()) );
	connect( startQueuedJobsButton, SIGNAL(clicked(bool)), this, SLOT(startQueuedJobsButtonClicked()) );
    connect( jobTableView, SIGNAL(clicked(const QModelIndex &)), this, SLOT(toggleDeleteJobsButton()) );
    connect( deleteJobsButton, SIGNAL(clicked()), this, SLOT(deleteSelectedJobs()));

    std::vector<std::string> imageDataDim;
    imageDataDim.push_back("2D");
    imageDataDim.push_back("3D");

	std::vector<std::string> thresholdModes;
	thresholdModes.push_back("Adaptive Otsu-Threshold");
	thresholdModes.push_back("Otsu-Threshold");
	thresholdModes.push_back("Manual Threshold");

	std::vector<std::string> analysisMode;
	analysisMode.push_back("Without Analysis");
	analysisMode.push_back("With Analysis");

	std::vector<std::string> saveModes;
	saveModes.push_back("Save everything");
	saveModes.push_back("Save only essentials");

    mpParameterRoot = new CSParameterContext("Root");

	//Convert + Analyze Network Graph from Skeleton Image*************************************************************************
	//TODO: make choice dialog for magnification selection -> change predefined voxel size values

	std::vector<std::string> networkType;
	networkType.push_back("Bile Canaliculi Network");
	networkType.push_back("Sinusoidal Network");
	m_extractGraph_NetworkType = new CSParameterChoice(networkType, 0);
	std::vector<std::string> inputType;
	inputType.push_back("Skeleton Image");
	inputType.push_back("Final Graph");
	m_extractGraph_InputType = new CSParameterChoice(inputType, 0);
	m_extractGraph_InputType->setControlSubContexts(true);

	m_extractGraph_Folder = "/path/to/dataset_containing_folder";
	m_extractGraph_DataSetFileList = "/path/to/file.ias";
	m_extractGraph_NetworkSegmentationFile = "/path/to/network_bin_file.tif";
	m_extractGraph_SkeletonInputFile = "/path/to/skeleton_image.tif";
	m_extractGraph_GraphInputFile = "/path/to/graph_0_file.txt";

	m_extractGraph_ResamplingFactor = 1;
	m_extractGraph_ResamplingMaxDist = 0.00;
	m_extractGraph_RemoveDeadEndsThreshold = 2.00;
	m_extractGraph_RemoveAllDeadEnds = false;
	m_extractGraph_CollapseIsecNodes = 1.50;
	m_extractGraph_GeomPruning = 10.00;

	m_extractGraph_AnalysisMode = new CSParameterChoice(analysisMode, 0);
	m_extractGraph_AnalysisMode->setControlSubContexts(true);
	m_extractGraph_AnalysisFilePrefix = "Bile_";
	m_extractGraph_DataSetID = "YourDataSetID";
	m_extractGraph_CellSegmentationFile = "/path/to/cell_bin_file.tif";
	m_extractGraph_CVSegmentationFile = "/path/to/vein_central_bin_file.tif";
	m_extractGraph_PVSegmentationFile = "/path/to/vein_portal_bin_file.tif";
	m_extractGraph_MaskDeadEnds = false;
	m_extractGraph_CVMaskFile = "/path/to/vein_central_mask_bin_file.tif";
	m_extractGraph_PVMaskFile = "/path/to/vein_portal_mask_bin_file.tif";
	m_extractGraph_NecRegSegFile = "/path/to/necReg_bin_file.tif";
	m_extractGraph_VoxSpacing[0] = 0.207;
	m_extractGraph_VoxSpacing[1] = 0.207;
	m_extractGraph_VoxSpacing[2] = 0.54;
	m_extractGraph_StaticRadius = 0.50;
	m_extractGraph_WriteToAnalysisFile = true;

	m_extractGraph_Display_VanillaGraph = false;
	m_extractGraph_Display_ResampledGraph = false;
	m_extractGraph_Display_CollapsedGraph = false;
	m_extractGraph_Display_FinalGraph = false;

	m_extractGraph_Filename = "graph";
	m_extractGraph_Save_Mode = new CSParameterChoice(saveModes, 1);

	m_extractGraphContext = mpParameterRoot->addContext( PipelineName[PipelineExtractAndAnalyzeGraph] );
	m_extractGraphContext->setDebug(true);

	m_extractGraphContext->addParameter("Network type", CSParameter::Choice, m_extractGraph_NetworkType, "");

	CSParameterContext *extractGraph_imageReader_subContext = m_extractGraphContext->addContext("Reader");
	extractGraph_imageReader_subContext->addParameter("Dataset containing folder", CSParameter::DirName, &m_extractGraph_Folder, "");
	extractGraph_imageReader_subContext->addParameter("Dataset file list", CSParameter::FileName, &m_extractGraph_DataSetFileList, "")->setAttribute("FileTypes","Dataset file list (*.ias)");
	extractGraph_imageReader_subContext->addParameter("Input type", CSParameter::Choice, m_extractGraph_InputType, "");
	CSParameterContext *extractGraph_skeletonReader_subContext = extractGraph_imageReader_subContext->addContext("Skeleton Reader");
	extractGraph_skeletonReader_subContext->addParameter("Skeleton image", CSParameter::FileName, &m_extractGraph_SkeletonInputFile, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	CSParameterContext *extractGraph_graphReader_subContext = extractGraph_imageReader_subContext->addContext("Graph Reader");
	extractGraph_graphReader_subContext->addParameter("Graph file", CSParameter::FileName, &m_extractGraph_GraphInputFile, "")->setAttribute("FileTypes","Graph file (*.txt)");
	extractGraph_graphReader_subContext->setVisible(false);
	extractGraph_imageReader_subContext->addParameter("Voxel spacing x", CSParameter::Double, &m_extractGraph_VoxSpacing[0], "micron");
	extractGraph_imageReader_subContext->addParameter("Voxel spacing y", CSParameter::Double, &m_extractGraph_VoxSpacing[1], "micron");
	extractGraph_imageReader_subContext->addParameter("Voxel spacing z", CSParameter::Double, &m_extractGraph_VoxSpacing[2], "micron");
	extractGraph_imageReader_subContext->addParameter("Specify network segmentation file", CSParameter::FileName, &m_extractGraph_NetworkSegmentationFile, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");

	CSParameterContext *extractGraph_parameter_subContext = m_extractGraphContext->addContext("Graph Extraction Parameter");
	CSParameterContext *extractGraph_resampling_subContext = extractGraph_parameter_subContext->addContext("Resampling Filter");
	extractGraph_resampling_subContext->addParameter("Resampling factor", CSParameter::Int, &m_extractGraph_ResamplingFactor, "");
	extractGraph_resampling_subContext->addParameter("Maximal resampling distance", CSParameter::Double, &m_extractGraph_ResamplingMaxDist, "micron");
	CSParameterContext *extractGraph_pruning_subContext = extractGraph_parameter_subContext->addContext("Remove Dead End Filter");
	extractGraph_pruning_subContext->addParameter("Length threshold", CSParameter::Double, &m_extractGraph_RemoveDeadEndsThreshold, "micron");
#ifndef CS_TI_QUANT_ONLY
	extractGraph_pruning_subContext->addParameter("Remove all dead ends (regardless of length)", CSParameter::Bool, &m_extractGraph_RemoveAllDeadEnds, "");
#endif
	CSParameterContext *extractGraph_collapseIsecs_subContext = extractGraph_parameter_subContext->addContext("Collapse Intersection Nodes Filter");
	extractGraph_collapseIsecs_subContext->addParameter("Distance threshold", CSParameter::Double, &m_extractGraph_CollapseIsecNodes, "micron");
	CSParameterContext *extractGraph_geomPruning_subContext = extractGraph_parameter_subContext->addContext("Geometric Pruning Filter");
	extractGraph_geomPruning_subContext->addParameter("Drop angle", CSParameter::Double, &m_extractGraph_GeomPruning, "degree");

	CSParameterContext *extractGraph_analysis_subContext = m_extractGraphContext->addContext("Network Analysis");
	extractGraph_analysis_subContext->addParameter("Analysis mode", CSParameter::Choice, m_extractGraph_AnalysisMode, "");
	CSParameterContext *extractGraph_withoutAnalysis_subContext = extractGraph_analysis_subContext->addContext("Without Analysis");
	CSParameterContext *extractGraph_withAnalysis_subContext = extractGraph_analysis_subContext->addContext("With Analysis");
	extractGraph_withAnalysis_subContext->addParameter("Output file prefix", CSParameter::String, &m_extractGraph_AnalysisFilePrefix, "");
	extractGraph_withAnalysis_subContext->addParameter("Specify data set name", CSParameter::String, &m_extractGraph_DataSetID, "");
	extractGraph_withAnalysis_subContext->addParameter("Specify cell segmentation file", CSParameter::FileName, &m_extractGraph_CellSegmentationFile, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	extractGraph_withAnalysis_subContext->addParameter("Specify central vein segmentation file", CSParameter::FileName, &m_extractGraph_CVSegmentationFile, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	extractGraph_withAnalysis_subContext->addParameter("Specify portal vein segmentation file", CSParameter::FileName, &m_extractGraph_PVSegmentationFile, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	extractGraph_withAnalysis_subContext->addParameter("Specify necrotic region segmentation file", CSParameter::FileName, &m_extractGraph_NecRegSegFile, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	extractGraph_withAnalysis_subContext->addParameter("Assumed average radius of branch", CSParameter::Double, &m_extractGraph_StaticRadius, "micron");
	extractGraph_withAnalysis_subContext->addParameter("Write results to analysis file", CSParameter::Bool, &m_extractGraph_WriteToAnalysisFile, "");
	extractGraph_withAnalysis_subContext->setVisible(false);

	CSParameterContext *extractGraph_display_subContext = m_extractGraphContext->addContext("Display Options");
	extractGraph_display_subContext->addParameter("Display vanilla graph", CSParameter::Bool, &m_extractGraph_Display_VanillaGraph, "");
	extractGraph_display_subContext->addParameter("Display resampled graph", CSParameter::Bool, &m_extractGraph_Display_ResampledGraph, "");
	extractGraph_display_subContext->addParameter("Display pruned graph", CSParameter::Bool, &m_extractGraph_Display_PrunedGraph, "");
	extractGraph_display_subContext->addParameter("Display collapsed graph", CSParameter::Bool, &m_extractGraph_Display_CollapsedGraph, "");
	extractGraph_display_subContext->addParameter("Display final graph", CSParameter::Bool, &m_extractGraph_Display_FinalGraph, "");

	CSParameterContext *extractGraph_save_subContext = m_extractGraphContext->addContext("Save Selector");
	extractGraph_save_subContext->addParameter("Save name", CSParameter::String, &m_extractGraph_Filename, "");
	extractGraph_save_subContext->addParameter("Save mode", CSParameter::Choice, m_extractGraph_Save_Mode, "");
	//****************************************************************************************************************************


	//Segment sinusoids and bile canaliculi in 20x datasets***********************************************************************
	//TODO: entry point stuff atm not useable

	m_segmentSinBile20x_DPPIV_Filename = "/path/to/DPPIV_channel_file.tif";
	m_segmentSinBile20x_DMs_Filename = "/path/to/DMs_channel_file.tif";
	m_segmentSinBile20x_UseNecReg = false;
	m_segmentSinBile20x_NecReg_Filename = "/path/to/necrotic_segmentation_file.tif";
	m_segmentSinBile20x_OtherInput_Filename = "/path/to/intermediate_result_file.tif";
	m_segmentSinBile20x_CVSeg_Filename = "/path/to/CV_bin_file.tif --optional";
	m_segmentSinBile20x_PVSeg_Filename = "/path/to/PV_bin_file.tif --optional";
	m_segmentSin20x_DPPIV_Threshold_Mode = new CSParameterChoice(thresholdModes, 0);
	m_segmentSin20x_DPPIV_Threshold_Mode->setControlSubContexts(true);
	m_segmentSin20x_DPPIV_AdapOtsuThreshold_RegionSize[0] = 100;
	m_segmentSin20x_DPPIV_AdapOtsuThreshold_RegionSize[1] = 100;
	m_segmentSin20x_DPPIV_AdapOtsuThreshold_RegionSize[2] = 10;
	m_segmentSin20x_DPPIV_AdapOtsuThreshold_NumSamples = 10000;
	m_segmentSin20x_DPPIV_ManualThreshold = 60;
	m_segmentSin20x_DMs_Threshold_Mode = new CSParameterChoice(thresholdModes, 0);
	m_segmentSin20x_DMs_Threshold_Mode->setControlSubContexts(true);
	m_segmentSin20x_DMs_AdapOtsuThreshold_RegionSize[0] = 100;
	m_segmentSin20x_DMs_AdapOtsuThreshold_RegionSize[1] = 100;
	m_segmentSin20x_DMs_AdapOtsuThreshold_RegionSize[2] = 10;
	m_segmentSin20x_DMs_AdapOtsuThreshold_NumSamples = 10000;
	m_segmentSin20x_DMs_ManualThreshold = 60;
	m_segmentSin20x_HoleFilling_RegionSize[0] = 1;
	m_segmentSin20x_HoleFilling_RegionSize[1] = 1;
	m_segmentSin20x_HoleFilling_RegionSize[2] = 1;
	m_segmentSin20x_HoleFilling_MajThreshold = 4;
	m_segmentSin20x_Closing_RadStrucElem = 2;
	m_segmentSin20x_Opening_RadStrucElem = 2;
	m_segmentSin20x_MaskVein_RegionSize[0] = 2;
	m_segmentSin20x_MaskVein_RegionSize[1] = 2;
	m_segmentSin20x_MaskVein_RegionSize[2] = 2;
	m_segmentSin20x_MinObjSize = 1000;

	m_segmentBile20x_DPPIV_Threshold_Mode = new CSParameterChoice(thresholdModes, 0);
	m_segmentBile20x_DPPIV_Threshold_Mode->setControlSubContexts(true);
	m_segmentBile20x_DPPIV_AdapOtsuThreshold_RegionSize[0] = 100;
	m_segmentBile20x_DPPIV_AdapOtsuThreshold_RegionSize[1] = 100;
	m_segmentBile20x_DPPIV_AdapOtsuThreshold_RegionSize[2] = 10;
	m_segmentBile20x_DPPIV_AdapOtsuThreshold_NumSamples = 10000;
	m_segmentBile20x_DPPIV_ManualThreshold = 60;
	m_segmentBile20x_HoleFilling_RegionSize[0] = 1;
	m_segmentBile20x_HoleFilling_RegionSize[1] = 1;
	m_segmentBile20x_HoleFilling_RegionSize[2] = 1;
	m_segmentBile20x_HoleFilling_MajThreshold = 8;
	m_segmentBile20x_InvHoleFilling_RegionSize[0] = 1;
	m_segmentBile20x_InvHoleFilling_RegionSize[1] = 1;
	m_segmentBile20x_InvHoleFilling_RegionSize[2] = 1;
	m_segmentBile20x_InvHoleFilling_MajThreshold = 5;
	m_segmentBile20x_MaskVein_RegionSize[0] = 2;
	m_segmentBile20x_MaskVein_RegionSize[1] = 2;
	m_segmentBile20x_MaskVein_RegionSize[2] = 2;
	m_segmentBile20x_MinObjSize = 300;

	m_segmentSinBile20x_Log_Filename = "log";
	m_segmentSin20x_Filename_Prefix = "sinus";
	m_segmentSin20x_Save_Mode = new CSParameterChoice(saveModes, 0);
	m_segmentBile20x_Filename_Prefix = "bile";
	m_segmentBile20x_Save_Mode = new CSParameterChoice(saveModes, 0);

	m_segmentSinBile20xContext = mpParameterRoot->addContext( PipelineName[PipelineSegmentSinusoidsAndBile20x] );
	m_segmentSinBile20xContext->setDebug(true);

	CSParameterContext *segmentSinBile20x_imageReader_subContext = m_segmentSinBile20xContext->addContext("Image Reader");
	segmentSinBile20x_imageReader_subContext->addParameter("DPPIV channel", CSParameter::FileName, &m_segmentSinBile20x_DPPIV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentSinBile20x_imageReader_subContext->addParameter("DMs channel", CSParameter::FileName, &m_segmentSinBile20x_DMs_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentSinBile20x_imageReader_subContext->addParameter("Central vein segmentation", CSParameter::FileName, &m_segmentSinBile20x_CVSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentSinBile20x_imageReader_subContext->addParameter("Portal vein segmentation", CSParameter::FileName, &m_segmentSinBile20x_PVSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentSinBile20x_imageReader_subContext->addParameter("Is there a necrotic region", CSParameter::Bool, &m_segmentSinBile20x_UseNecReg, "");
	segmentSinBile20x_imageReader_subContext->addParameter("Necrotic Region", CSParameter::FileName, &m_segmentSinBile20x_NecReg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");

	CSParameterContext *segmentSin20x_subContext = m_segmentSinBile20xContext->addContext("Segment Sinusoids on DPPIV + DMs Channel");

	CSParameterContext *segmentSin20x_binThreshold_subContext = segmentSin20x_subContext->addContext("1.1) Multi-Channel Binary Threshold Filter on DPPIV + DMs Channel");
	CSParameterContext *segmentSin20x_binThresholdDPPIV_subContext = segmentSin20x_binThreshold_subContext->addContext("1.1a) Binary Threshold on DPPIV Channel");
	segmentSin20x_binThresholdDPPIV_subContext->addParameter("DPPIV channel threshold mode", CSParameter::Choice, m_segmentSin20x_DPPIV_Threshold_Mode, "");
	CSParameterContext *segmentSin20x_adapOtsuThresholdDPPIV_subContext = segmentSin20x_binThresholdDPPIV_subContext->addContext("Adaptive Otsu-Threshold");
	segmentSin20x_adapOtsuThresholdDPPIV_subContext->addParameter("Sample region size x", CSParameter::Int, &m_segmentSin20x_DPPIV_AdapOtsuThreshold_RegionSize[0], "pixel");
	segmentSin20x_adapOtsuThresholdDPPIV_subContext->addParameter("Sample region size y", CSParameter::Int, &m_segmentSin20x_DPPIV_AdapOtsuThreshold_RegionSize[1], "pixel");
	segmentSin20x_adapOtsuThresholdDPPIV_subContext->addParameter("Sample region size z", CSParameter::Int, &m_segmentSin20x_DPPIV_AdapOtsuThreshold_RegionSize[2], "pixel");
	segmentSin20x_adapOtsuThresholdDPPIV_subContext->addParameter("Number samples", CSParameter::Int, &m_segmentSin20x_DPPIV_AdapOtsuThreshold_NumSamples, "");
	CSParameterContext *segmentSin20x_otsuThresholdDPPIV_subContext = segmentSin20x_binThresholdDPPIV_subContext->addContext("Otsu-Threshold");
	segmentSin20x_otsuThresholdDPPIV_subContext->setVisible(false);
	CSParameterContext *segmentSin20x_manualThresholdDPPIV_subContext = segmentSin20x_binThresholdDPPIV_subContext->addContext("Manual Threshold");
	segmentSin20x_manualThresholdDPPIV_subContext->addParameter("DPPIV manual threshold", CSParameter::Int, &m_segmentSin20x_DPPIV_ManualThreshold, "");
	segmentSin20x_manualThresholdDPPIV_subContext->setVisible(false);
	CSParameterContext *segmentSin20x_binThresholdDMs_subContext = segmentSin20x_binThreshold_subContext->addContext("1.1b) Binary Threshold on DMs Channel");
	segmentSin20x_binThresholdDMs_subContext->addParameter("DMs channel threshold mode", CSParameter::Choice, m_segmentSin20x_DMs_Threshold_Mode, "");
	CSParameterContext *segmentSin20x_adapOtsuThresholdDMs_subContext = segmentSin20x_binThresholdDMs_subContext->addContext("Adaptive Otsu-Threshold");
	segmentSin20x_adapOtsuThresholdDMs_subContext->addParameter("Sample region size x", CSParameter::Int, &m_segmentSin20x_DMs_AdapOtsuThreshold_RegionSize[0], "pixel");
	segmentSin20x_adapOtsuThresholdDMs_subContext->addParameter("Sample region size y", CSParameter::Int, &m_segmentSin20x_DMs_AdapOtsuThreshold_RegionSize[1], "pixel");
	segmentSin20x_adapOtsuThresholdDMs_subContext->addParameter("Sample region size z", CSParameter::Int, &m_segmentSin20x_DMs_AdapOtsuThreshold_RegionSize[2], "pixel");
	segmentSin20x_adapOtsuThresholdDMs_subContext->addParameter("Number samples", CSParameter::Int, &m_segmentSin20x_DMs_AdapOtsuThreshold_NumSamples, "");
	CSParameterContext *segmentSin20x_otsuThresholdDMs_subContext = segmentSin20x_binThresholdDMs_subContext->addContext("Otsu-Threshold");
	segmentSin20x_otsuThresholdDMs_subContext->setVisible(false);
	CSParameterContext *segmentSin20x_manualThresholdDMs_subContext = segmentSin20x_binThresholdDMs_subContext->addContext("Manual Threshold");
	segmentSin20x_manualThresholdDMs_subContext->addParameter("DMs manual threshold", CSParameter::Int, &m_segmentSin20x_DMs_ManualThreshold, "");
	segmentSin20x_manualThresholdDMs_subContext->setVisible(false);
	CSParameterContext *segmentSin20x_holeFilling_subContext = segmentSin20x_subContext->addContext("1.2) Hole Filling on 1.1");
	segmentSin20x_holeFilling_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentSin20x_HoleFilling_RegionSize[0], "pixel");
	segmentSin20x_holeFilling_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentSin20x_HoleFilling_RegionSize[1], "pixel");
	segmentSin20x_holeFilling_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentSin20x_HoleFilling_RegionSize[2], "pixel");
	segmentSin20x_holeFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &m_segmentSin20x_HoleFilling_MajThreshold, "");
	CSParameterContext *segmentSin20x_closing_subContext = segmentSin20x_subContext->addContext("1.3) Closing on 1.2");
	segmentSin20x_closing_subContext->addParameter("Kernel radius", CSParameter::Int, &m_segmentSin20x_Closing_RadStrucElem, "pixel");
	CSParameterContext *segmentSin20x_opening_subContext = segmentSin20x_subContext->addContext("1.4) Opening on 1.3");
	segmentSin20x_opening_subContext->addParameter("Kernel radius", CSParameter::Int, &m_segmentSin20x_Opening_RadStrucElem, "pixel");
	CSParameterContext *segmentSin20x_maskVein_subContext = segmentSin20x_subContext->addContext("1.5) Optional vein masking on 1.4");
	segmentSin20x_maskVein_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentSin20x_MaskVein_RegionSize[0], "pixel");
	segmentSin20x_maskVein_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentSin20x_MaskVein_RegionSize[1], "pixel");
	segmentSin20x_maskVein_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentSin20x_MaskVein_RegionSize[2], "pixel");
	CSParameterContext *segmentSin20x_removeObjects_subContext = segmentSin20x_subContext->addContext("1.6) Remove objects on 1.5");
	segmentSin20x_removeObjects_subContext->addParameter("Smallest allowed object size", CSParameter::Int, &m_segmentSin20x_MinObjSize, "voxel");
	CSParameterContext *segmentSin20x_saveSelector_subContext = segmentSin20x_subContext->addContext("Save Selector");
	segmentSin20x_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_segmentSinBile20x_Log_Filename, "");
	segmentSin20x_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_segmentSin20x_Filename_Prefix, "");
	segmentSin20x_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_segmentSin20x_Save_Mode, "");

	CSParameterContext *segmentBile20x_subContext = m_segmentSinBile20xContext->addContext("Segment Bile Canaliculi on DPPIV Channel");

	CSParameterContext *segmentBile20x_binThreshold_subContext = segmentBile20x_subContext->addContext("2.1) Binary Threshold Filter on DPPIV Channel");
	segmentBile20x_binThreshold_subContext->addParameter("DPPIV channel threshold mode", CSParameter::Choice, m_segmentBile20x_DPPIV_Threshold_Mode, "");
	CSParameterContext *segmentBile20x_adapOtsuThreshold_subContext = segmentBile20x_binThreshold_subContext->addContext("Adaptive Otsu-Threshold");
	segmentBile20x_adapOtsuThreshold_subContext->addParameter("Sample region size x", CSParameter::Int, &m_segmentBile20x_DPPIV_AdapOtsuThreshold_RegionSize[0], "pixel");
	segmentBile20x_adapOtsuThreshold_subContext->addParameter("Sample region size y", CSParameter::Int, &m_segmentBile20x_DPPIV_AdapOtsuThreshold_RegionSize[1], "pixel");
	segmentBile20x_adapOtsuThreshold_subContext->addParameter("Sample region size z", CSParameter::Int, &m_segmentBile20x_DPPIV_AdapOtsuThreshold_RegionSize[2], "pixel");
	segmentBile20x_adapOtsuThreshold_subContext->addParameter("Number samples", CSParameter::Int, &m_segmentBile20x_DPPIV_AdapOtsuThreshold_NumSamples, "");
	CSParameterContext *segmentBile20x_otsuThreshold_subContext = segmentBile20x_binThreshold_subContext->addContext("Otsu-Threshold");
	segmentBile20x_otsuThreshold_subContext->setVisible(false);
	CSParameterContext *segmentBile20x_manualThreshold_subContext = segmentBile20x_binThreshold_subContext->addContext("Manual Threshold");
	segmentBile20x_manualThreshold_subContext->addParameter("DPPIV manual threshold", CSParameter::Int, &m_segmentBile20x_DPPIV_ManualThreshold, "");
	segmentBile20x_manualThreshold_subContext->setVisible(false);
	CSParameterContext *segmentBile20x_holeFilling_subContext = segmentBile20x_subContext->addContext("2.2) Hole Filling on 2.1");
	segmentBile20x_holeFilling_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentBile20x_HoleFilling_RegionSize[0], "pixel");
	segmentBile20x_holeFilling_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentBile20x_HoleFilling_RegionSize[1], "pixel");
	segmentBile20x_holeFilling_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentBile20x_HoleFilling_RegionSize[2], "pixel");
	segmentBile20x_holeFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &m_segmentBile20x_HoleFilling_MajThreshold, "");
	CSParameterContext *segmentBile20x_invHoleFilling_subContext = segmentBile20x_subContext->addContext("2.3) Inverse Hole Filling on 2.2");
	segmentBile20x_invHoleFilling_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentBile20x_InvHoleFilling_RegionSize[0], "pixel");
	segmentBile20x_invHoleFilling_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentBile20x_InvHoleFilling_RegionSize[1], "pixel");
	segmentBile20x_invHoleFilling_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentBile20x_InvHoleFilling_RegionSize[2], "pixel");
	segmentBile20x_invHoleFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &m_segmentBile20x_InvHoleFilling_MajThreshold, "");
	CSParameterContext *segmentBile20x_maskVein_subContext = segmentBile20x_subContext->addContext("2.4) Optional vein masking on 2.3");
	segmentBile20x_maskVein_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentBile20x_MaskVein_RegionSize[0], "pixel");
	segmentBile20x_maskVein_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentBile20x_MaskVein_RegionSize[1], "pixel");
	segmentBile20x_maskVein_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentBile20x_MaskVein_RegionSize[2], "pixel");
	CSParameterContext *segmentBile20x_removeObjects_subContext = segmentBile20x_subContext->addContext("2.5) Remove objects on 2.4");
	segmentBile20x_removeObjects_subContext->addParameter("Smallest allowed object size", CSParameter::Int, &m_segmentBile20x_MinObjSize, "voxel");
	CSParameterContext *segmentBile20x_saveSelector_subContext = segmentBile20x_subContext->addContext("Save Selector");
	segmentBile20x_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_segmentSinBile20x_Log_Filename, "");
	segmentBile20x_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_segmentBile20x_Filename_Prefix, "");
	segmentBile20x_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_segmentBile20x_Save_Mode, "");
	//****************************************************************************************************************************

	//Segment sinusoids and bile canaliculi in 60x datasets***********************************************************************
	m_segmentSinBile60x_DPPIV_Filename = "/path/to/DPPIV_channel_file.tif";
	m_segmentSinBile60x_DMs_Filename = "/path/to/DMs_channel_file.tif";
	m_segmentSinBile60x_CVSeg_Filename = "/path/to/CV_bin_file.tif --optional";
	m_segmentSinBile60x_PVSeg_Filename = "/path/to/PV_bin_file.tif --optional";
	m_segmentSinBile60x_VoxelSpacing[0] = 0.207;
	m_segmentSinBile60x_VoxelSpacing[1] = 0.207;
	m_segmentSinBile60x_VoxelSpacing[2] = 0.54;
	m_segmentSinBile60x_UseNecReg = false;
	m_segmentSinBile60x_NecReg_Filename = "/path/to/necrotic_segmentation_file.tif";
	m_segmentSinBile60x_OtherInput_Filename = "/path/to/intermediate_result_file.tif";
	m_segmentSin60x_Closing_RadStrucElem = 0;
	m_segmentSin60x_DPPIV_Threshold_Mode = new CSParameterChoice(thresholdModes, 0);
	m_segmentSin60x_DPPIV_Threshold_Mode->setControlSubContexts(true);
	m_segmentSin60x_DPPIV_AdapOtsuThreshold_RegionSize[0] = 100;
	m_segmentSin60x_DPPIV_AdapOtsuThreshold_RegionSize[1] = 100;
	m_segmentSin60x_DPPIV_AdapOtsuThreshold_RegionSize[2] = 10;
	m_segmentSin60x_DPPIV_AdapOtsuThreshold_NumSamples = 10000;
	m_segmentSin60x_DPPIV_ManualThreshold = 60;
	m_segmentSin60x_DMs_Threshold_Mode = new CSParameterChoice(thresholdModes, 0);
	m_segmentSin60x_DMs_Threshold_Mode->setControlSubContexts(true);
	m_segmentSin60x_DMs_AdapOtsuThreshold_RegionSize[0] = 100;
	m_segmentSin60x_DMs_AdapOtsuThreshold_RegionSize[1] = 100;
	m_segmentSin60x_DMs_AdapOtsuThreshold_RegionSize[2] = 10;
	m_segmentSin60x_DMs_AdapOtsuThreshold_NumSamples = 10000;
	m_segmentSin60x_DMs_ManualThreshold = 60;
	m_segmentSin60x_InvHoleFilling_RegionSize[0] = 1;
	m_segmentSin60x_InvHoleFilling_RegionSize[1] = 1;
	m_segmentSin60x_InvHoleFilling_RegionSize[2] = 1;
	m_segmentSin60x_InvHoleFilling_MajThreshold = 4;
	m_segmentSin60x_CavityFilling_WithRescaling = true;
	m_segmentSin60x_CavityFilling_Radius = 20;
	m_segmentSin60x_CavityFilling_MinThreshold = 0.7;
	m_segmentSin60x_CavityFilling_MaxThreshold = 1.0;
	m_segmentSin60x_Closing_xRadStrucElem = 3;
	m_segmentSin60x_Closing_yRadStrucElem = 3;
	m_segmentSin60x_Closing_zRadStrucElem = 1;
	m_segmentSin60x_Opening_xRadStrucElem = 3;
	m_segmentSin60x_Opening_yRadStrucElem = 3;
	m_segmentSin60x_Opening_zRadStrucElem = 1;
	m_segmentSin60x_MaskVein_RegionSize[0] = 10;
	m_segmentSin60x_MaskVein_RegionSize[1] = 10;
	m_segmentSin60x_MaskVein_RegionSize[2] = 5;
	m_segmentSin60x_MinObjSize = 10000;

	m_segmentBile60x_Median_RadStrucElem = 0;
	m_segmentBile60x_GreyOpening_RadStrucElem = 0;
	m_segmentBile60x_DPPIV_Threshold_Mode = new CSParameterChoice(thresholdModes, 0);
	m_segmentBile60x_DPPIV_Threshold_Mode->setControlSubContexts(true);
	m_segmentBile60x_DPPIV_AdapOtsuThreshold_RegionSize[0] = 100;
	m_segmentBile60x_DPPIV_AdapOtsuThreshold_RegionSize[1] = 100;
	m_segmentBile60x_DPPIV_AdapOtsuThreshold_RegionSize[2] = 10;
	m_segmentBile60x_DPPIV_AdapOtsuThreshold_NumSamples = 10000;
	m_segmentBile60x_DPPIV_ManualThreshold = 60;
	m_segmentBile60x_HoleFilling_RegionSize[0] = 1;
	m_segmentBile60x_HoleFilling_RegionSize[1] = 1;
	m_segmentBile60x_HoleFilling_RegionSize[2] = 1;
	m_segmentBile60x_HoleFilling_MajThreshold = 1;
	m_segmentBile60x_BinOpening_RadStrucElem = 0;
	m_segmentBile60x_InvHoleFilling_RegionSize[0] = 1;
	m_segmentBile60x_InvHoleFilling_RegionSize[1] = 1;
	m_segmentBile60x_InvHoleFilling_RegionSize[2] = 1;
	m_segmentBile60x_InvHoleFilling_MajThreshold = 6;
	m_segmentBile60x_MaskVein_RegionSize[0] = 10;
	m_segmentBile60x_MaskVein_RegionSize[1] = 10;
	m_segmentBile60x_MaskVein_RegionSize[2] = 5;
	m_segmentBile60x_MinObjSize = 3000;

	m_segmentSinBile60x_Log_Filename = "log";
	m_segmentSin60x_Filename_Prefix = "sinus";
	m_segmentSin60x_Save_Mode = new CSParameterChoice(saveModes, 0);
	m_segmentBile60x_Filename_Prefix = "bile";
	m_segmentBile60x_Save_Mode = new CSParameterChoice(saveModes, 0);

	m_segmentSinBile60xContext = mpParameterRoot->addContext( PipelineName[PipelineSegmentSinusoidsAndBile60x] );
	m_segmentSinBile60xContext->setDebug(true);

	CSParameterContext *segmentSinBile60x_imageReader_subContext = m_segmentSinBile60xContext->addContext("Image Reader");
	segmentSinBile60x_imageReader_subContext->addParameter("DPPIV channel", CSParameter::FileName, &m_segmentSinBile60x_DPPIV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentSinBile60x_imageReader_subContext->addParameter("DMs channel", CSParameter::FileName, &m_segmentSinBile60x_DMs_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentSinBile60x_imageReader_subContext->addParameter("Voxel spacing x", CSParameter::Double, &m_segmentSinBile60x_VoxelSpacing[0], "micron");
	segmentSinBile60x_imageReader_subContext->addParameter("Voxel spacing y", CSParameter::Double, &m_segmentSinBile60x_VoxelSpacing[1], "micron");
	segmentSinBile60x_imageReader_subContext->addParameter("Voxel spacing z", CSParameter::Double, &m_segmentSinBile60x_VoxelSpacing[2], "micron");
	segmentSinBile60x_imageReader_subContext->addParameter("Central vein segmentation", CSParameter::FileName, &m_segmentSinBile60x_CVSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentSinBile60x_imageReader_subContext->addParameter("Portal vein segmentation", CSParameter::FileName, &m_segmentSinBile60x_PVSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentSinBile60x_imageReader_subContext->addParameter("Is there a necrotic region", CSParameter::Bool, &m_segmentSinBile60x_UseNecReg, "");
	segmentSinBile60x_imageReader_subContext->addParameter("Necrotic Region", CSParameter::FileName, &m_segmentSinBile60x_NecReg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");

	CSParameterContext *segmentSin60x_subContext = m_segmentSinBile60xContext->addContext("Segment Sinusoids on DPPIV + DMs Channel");

	CSParameterContext *segmentSin60x_preprocessing_subContext = segmentSin60x_subContext->addContext("1.1) Preprocess DPPIV + DMs Channel");
	segmentSin60x_preprocessing_subContext->addParameter("Greyscale closing kernel radius", CSParameter::Int, &m_segmentSin60x_Closing_RadStrucElem, "pixel");
	CSParameterContext *segmentSin60x_binThreshold_subContext = segmentSin60x_subContext->addContext("1.2) Multi-Channel Binary Threshold Filter on 1.1");
	CSParameterContext *segmentSin60x_binThresholdDPPIV_subContext = segmentSin60x_binThreshold_subContext->addContext("1.2a) Binary Threshold on DPPIV Channel");
	segmentSin60x_binThresholdDPPIV_subContext->addParameter("DPPIV channel threshold mode", CSParameter::Choice, m_segmentSin60x_DPPIV_Threshold_Mode, "");
	CSParameterContext *segmentSin60x_adapOtsuThresholdDPPIV_subContext = segmentSin60x_binThresholdDPPIV_subContext->addContext("Adaptive Otsu-Threshold");
	segmentSin60x_adapOtsuThresholdDPPIV_subContext->addParameter("Sample region size x", CSParameter::Int, &m_segmentSin60x_DPPIV_AdapOtsuThreshold_RegionSize[0], "pixel");
	segmentSin60x_adapOtsuThresholdDPPIV_subContext->addParameter("Sample region size y", CSParameter::Int, &m_segmentSin60x_DPPIV_AdapOtsuThreshold_RegionSize[1], "pixel");
	segmentSin60x_adapOtsuThresholdDPPIV_subContext->addParameter("Sample region size z", CSParameter::Int, &m_segmentSin60x_DPPIV_AdapOtsuThreshold_RegionSize[2], "pixel");
	segmentSin60x_adapOtsuThresholdDPPIV_subContext->addParameter("Number samples", CSParameter::Int, &m_segmentSin60x_DPPIV_AdapOtsuThreshold_NumSamples, "");
	CSParameterContext *segmentSin60x_otsuThresholdDPPIV_subContext = segmentSin60x_binThresholdDPPIV_subContext->addContext("Otsu-Threshold");
	segmentSin60x_otsuThresholdDPPIV_subContext->setVisible(false);
	CSParameterContext *segmentSin60x_manualThresholdDPPIV_subContext = segmentSin60x_binThresholdDPPIV_subContext->addContext("Manual Threshold");
	segmentSin60x_manualThresholdDPPIV_subContext->addParameter("DPPIV manual threshold", CSParameter::Int, &m_segmentSin60x_DPPIV_ManualThreshold, "");
	segmentSin60x_manualThresholdDPPIV_subContext->setVisible(false);
	CSParameterContext *segmentSin60x_binThresholdDMs_subContext = segmentSin60x_binThreshold_subContext->addContext("1.2b) Binary Threshold on DMs Channel");
	segmentSin60x_binThresholdDMs_subContext->addParameter("DMs channel threshold mode", CSParameter::Choice, m_segmentSin60x_DMs_Threshold_Mode, "");
	CSParameterContext *segmentSin60x_adapOtsuThresholdDMs_subContext = segmentSin60x_binThresholdDMs_subContext->addContext("Adaptive Otsu-Threshold");
	segmentSin60x_adapOtsuThresholdDMs_subContext->addParameter("Sample region size x", CSParameter::Int, &m_segmentSin60x_DMs_AdapOtsuThreshold_RegionSize[0], "pixel");
	segmentSin60x_adapOtsuThresholdDMs_subContext->addParameter("Sample region size y", CSParameter::Int, &m_segmentSin60x_DMs_AdapOtsuThreshold_RegionSize[1], "pixel");
	segmentSin60x_adapOtsuThresholdDMs_subContext->addParameter("Sample region size z", CSParameter::Int, &m_segmentSin60x_DMs_AdapOtsuThreshold_RegionSize[2], "pixel");
	segmentSin60x_adapOtsuThresholdDMs_subContext->addParameter("Number samples", CSParameter::Int, &m_segmentSin60x_DMs_AdapOtsuThreshold_NumSamples, "");
	CSParameterContext *segmentSin60x_otsuThresholdDMs_subContext = segmentSin60x_binThresholdDMs_subContext->addContext("Otsu-Threshold");
	segmentSin60x_otsuThresholdDMs_subContext->setVisible(false);
	CSParameterContext *segmentSin60x_manualThresholdDMs_subContext = segmentSin60x_binThresholdDMs_subContext->addContext("Manual Threshold");
	segmentSin60x_manualThresholdDMs_subContext->addParameter("DMs manual threshold", CSParameter::Int, &m_segmentSin60x_DMs_ManualThreshold, "");
	segmentSin60x_manualThresholdDMs_subContext->setVisible(false);
	CSParameterContext *segmentSin60x_invHoleFilling_subContext = segmentSin60x_subContext->addContext("1.3) Inverse Hole Filling on 1.2");
	segmentSin60x_invHoleFilling_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentSin60x_InvHoleFilling_RegionSize[0], "pixel");
	segmentSin60x_invHoleFilling_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentSin60x_InvHoleFilling_RegionSize[1], "pixel");
	segmentSin60x_invHoleFilling_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentSin60x_InvHoleFilling_RegionSize[2], "pixel");
	segmentSin60x_invHoleFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &m_segmentSin60x_InvHoleFilling_MajThreshold, "");
	CSParameterContext *segmentSin60x_cavityFilling_subContext = segmentSin60x_subContext->addContext("1.4) Cavity Filling on 1.3");
	segmentSin60x_cavityFilling_subContext->addParameter("Fast (less accurate)", CSParameter::Bool, &m_segmentSin60x_CavityFilling_WithRescaling, "");
	segmentSin60x_cavityFilling_subContext->addParameter("Radius", CSParameter::Int, &m_segmentSin60x_CavityFilling_Radius, "pixel");
	segmentSin60x_cavityFilling_subContext->addParameter("Minimal fraction of surrounding foreground", CSParameter::Double, &m_segmentSin60x_CavityFilling_MinThreshold, "");
	segmentSin60x_cavityFilling_subContext->addParameter("Maximal fraction of surrounding foreground", CSParameter::Double, &m_segmentSin60x_CavityFilling_MaxThreshold, "");
	CSParameterContext *segmentSin60x_opening_subContext = segmentSin60x_subContext->addContext("1.5) Closing & Opening on 1.4");
	segmentSin60x_opening_subContext->addParameter("Closing kernel radius x", CSParameter::Int, &m_segmentSin60x_Closing_xRadStrucElem, "pixel");
	segmentSin60x_opening_subContext->addParameter("Closing kernel radius y", CSParameter::Int, &m_segmentSin60x_Closing_yRadStrucElem, "pixel");
	segmentSin60x_opening_subContext->addParameter("Closing kernel radius z", CSParameter::Int, &m_segmentSin60x_Closing_zRadStrucElem, "pixel");
	segmentSin60x_opening_subContext->addParameter("Opening kernel radius x", CSParameter::Int, &m_segmentSin60x_Opening_xRadStrucElem, "pixel");
	segmentSin60x_opening_subContext->addParameter("Opening kernel radius y", CSParameter::Int, &m_segmentSin60x_Opening_yRadStrucElem, "pixel");
	segmentSin60x_opening_subContext->addParameter("Opening kernel radius z", CSParameter::Int, &m_segmentSin60x_Opening_zRadStrucElem, "pixel");
	CSParameterContext *segmentSin60x_maskVein_subContext = segmentSin60x_subContext->addContext("1.6) Optional vein masking on 1.5");
	segmentSin60x_maskVein_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentSin60x_MaskVein_RegionSize[0], "pixel");
	segmentSin60x_maskVein_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentSin60x_MaskVein_RegionSize[1], "pixel");
	segmentSin60x_maskVein_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentSin60x_MaskVein_RegionSize[2], "pixel");
	CSParameterContext *segmentSin60x_removeObjects_subContext = segmentSin60x_subContext->addContext("1.7) Remove objects on 1.6");
	segmentSin60x_removeObjects_subContext->addParameter("Smallest allowed object size", CSParameter::Int, &m_segmentSin60x_MinObjSize, "voxel");
	CSParameterContext *segmentSin60x_saveSelector_subContext = segmentSin60x_subContext->addContext("Save Selector");
	segmentSin60x_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_segmentSinBile60x_Log_Filename, "");
	segmentSin60x_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_segmentSin60x_Filename_Prefix, "");
	segmentSin60x_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_segmentSin60x_Save_Mode, "");

	CSParameterContext *segmentBile60x_subContext = m_segmentSinBile60xContext->addContext("Segment Bile Canaliculi on DPPIV Channel");

	CSParameterContext *segmentBile60x_preprocessing_subContext = segmentBile60x_subContext->addContext("2.1) Preprocess DPPIV Channel");
	segmentBile60x_preprocessing_subContext->addParameter("Median filter kernel radius", CSParameter::Int, &m_segmentBile60x_Median_RadStrucElem, "pixel");
	segmentBile60x_preprocessing_subContext->addParameter("Greyscale opening kernel radius", CSParameter::Int, &m_segmentBile60x_GreyOpening_RadStrucElem, "pixel");
	CSParameterContext *segmentBile60x_binThreshold_subContext = segmentBile60x_subContext->addContext("2.2) Binary Threshold Filter on 2.1");
	segmentBile60x_binThreshold_subContext->addParameter("DPPIV channel threshold mode", CSParameter::Choice, m_segmentBile60x_DPPIV_Threshold_Mode, "");
	CSParameterContext *segmentBile60x_adapOtsuThreshold_subContext = segmentBile60x_binThreshold_subContext->addContext("Adaptive Otsu-Threshold");
	segmentBile60x_adapOtsuThreshold_subContext->addParameter("Sample region size x", CSParameter::Int, &m_segmentBile60x_DPPIV_AdapOtsuThreshold_RegionSize[0], "pixel");
	segmentBile60x_adapOtsuThreshold_subContext->addParameter("Sample region size y", CSParameter::Int, &m_segmentBile60x_DPPIV_AdapOtsuThreshold_RegionSize[1], "pixel");
	segmentBile60x_adapOtsuThreshold_subContext->addParameter("Sample region size z", CSParameter::Int, &m_segmentBile60x_DPPIV_AdapOtsuThreshold_RegionSize[2], "pixel");
	segmentBile60x_adapOtsuThreshold_subContext->addParameter("Number samples", CSParameter::Int, &m_segmentBile60x_DPPIV_AdapOtsuThreshold_NumSamples, "");
	CSParameterContext *segmentBile60x_otsuThreshold_subContext = segmentBile60x_binThreshold_subContext->addContext("Otsu-Threshold");
	segmentBile60x_otsuThreshold_subContext->setVisible(false);
	CSParameterContext *segmentBile60x_manualThreshold_subContext = segmentBile60x_binThreshold_subContext->addContext("Manual Threshold");
	segmentBile60x_manualThreshold_subContext->addParameter("DPPIV manual threshold", CSParameter::Int, &m_segmentBile60x_DPPIV_ManualThreshold, "");
	segmentBile60x_manualThreshold_subContext->setVisible(false);
	CSParameterContext *segmentBile60x_holeFilling_subContext = segmentBile60x_subContext->addContext("2.3) Hole Filling on 2.2");
	segmentBile60x_holeFilling_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentBile60x_HoleFilling_RegionSize[0], "pixel");
	segmentBile60x_holeFilling_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentBile60x_HoleFilling_RegionSize[1], "pixel");
	segmentBile60x_holeFilling_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentBile60x_HoleFilling_RegionSize[2], "pixel");
	segmentBile60x_holeFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &m_segmentBile60x_HoleFilling_MajThreshold, "");
	CSParameterContext *segmentBile60x_opening_subContext = segmentBile60x_subContext->addContext("2.4) Opening on 2.3");
	segmentBile60x_opening_subContext->addParameter("Kernel radius", CSParameter::Int, &m_segmentBile60x_BinOpening_RadStrucElem, "pixel");
	CSParameterContext *segmentBile60x_invHoleFilling_subContext = segmentBile60x_subContext->addContext("2.5) Inverse Hole Filling on 2.4");
	segmentBile60x_invHoleFilling_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentBile60x_InvHoleFilling_RegionSize[0], "pixel");
	segmentBile60x_invHoleFilling_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentBile60x_InvHoleFilling_RegionSize[1], "pixel");
	segmentBile60x_invHoleFilling_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentBile60x_InvHoleFilling_RegionSize[2], "pixel");
	segmentBile60x_invHoleFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &m_segmentBile60x_InvHoleFilling_MajThreshold, "");
	CSParameterContext *segmentBile60x_maskVein_subContext = segmentBile60x_subContext->addContext("2.6) Optional vein masking on 2.5");
	segmentBile60x_maskVein_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentBile60x_MaskVein_RegionSize[0], "pixel");
	segmentBile60x_maskVein_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentBile60x_MaskVein_RegionSize[1], "pixel");
	segmentBile60x_maskVein_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentBile60x_MaskVein_RegionSize[2], "pixel");
	CSParameterContext *segmentBile60x_removeObjects_subContext = segmentBile60x_subContext->addContext("2.7) Remove objects on 2.6");
	segmentBile60x_removeObjects_subContext->addParameter("Smallest allowed object size", CSParameter::Int, &m_segmentBile60x_MinObjSize, "voxel");
	CSParameterContext *segmentBile60x_saveSelector_subContext = segmentBile60x_subContext->addContext("Save Selector");
	segmentBile60x_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_segmentSinBile60x_Log_Filename, "");
	segmentBile60x_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_segmentBile60x_Filename_Prefix, "");
	segmentBile60x_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_segmentBile60x_Save_Mode, "");
    //****************************************************************************************************************************


	//Segment nuclei in 20x datasets**********************************************************************************************
	m_segmentNuclei20x_DAPI_Filename = "/path/to/DAPI_channel_file.tif";
	m_segmentNuclei20x_CVSeg_Filename = "/path/to/CV_bin_file.tif --optional";
	m_segmentNuclei20x_PVSeg_Filename = "/path/to/PV_bin_file.tif --optional";
	m_segmentNuclei20x_UseNecReg = false;
	m_segmentNuclei20x_NecReg_Filename = "/path/to/necrotic_segmentation_file.tif";
	m_segmentNuclei20x_VoxelSpacing[0] = 0.620996;
	m_segmentNuclei20x_VoxelSpacing[1] = 0.620996;
	m_segmentNuclei20x_VoxelSpacing[2] = 0.54;
	m_segmentNuclei20x_Median_Radius = 0;
	m_segmentNuclei20x_DAPI_Threshold_Mode = new CSParameterChoice(thresholdModes, 0);
	m_segmentNuclei20x_DAPI_Threshold_Mode->setControlSubContexts(true);
	m_segmentNuclei20x_DAPI_AdapOtsuThreshold_RegionSize[0] = 100;
	m_segmentNuclei20x_DAPI_AdapOtsuThreshold_RegionSize[1] = 100;
	m_segmentNuclei20x_DAPI_AdapOtsuThreshold_RegionSize[2] = 10;
	m_segmentNuclei20x_DAPI_AdapOtsuThreshold_NumSamples = 10000;
	m_segmentNuclei20x_DAPI_ManualThreshold = 10;
	m_segmentNuclei20x_HoleFilling_RegionSize[0] = 1;
	m_segmentNuclei20x_HoleFilling_RegionSize[1] = 1;
	m_segmentNuclei20x_HoleFilling_RegionSize[2] = 1;
	m_segmentNuclei20x_HoleFilling_MajThreshold = 1;
	m_segmentNuclei20x_Opening_RadStrucElem = 1;
	m_segmentNuclei20x_MinDiameterThreshold = 6.0;
	m_segmentNuclei20x_MaxDiameterThreshold = 20.0;
	m_segmentNuclei20x_WatershedFloodLevel = 0.18;

	m_segmentNuclei20x_Log_Filename = "log";
	m_segmentNuclei20x_Filename_Prefix = "nuclei";
	m_segmentNuclei20x_Save_Mode = new CSParameterChoice(saveModes, 0);

	m_segmentNuclei20xContext = mpParameterRoot->addContext( PipelineName[PipelineSegmentNuclei20x] );
	m_segmentNuclei20xContext->setDebug(true);

	CSParameterContext *segmentNuclei20x_imageReader_subContext = m_segmentNuclei20xContext->addContext("Image Reader");
	segmentNuclei20x_imageReader_subContext->addParameter("DAPI channel", CSParameter::FileName, &m_segmentNuclei20x_DAPI_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentNuclei20x_imageReader_subContext->addParameter("Voxel spacing x", CSParameter::Double, &m_segmentNuclei20x_VoxelSpacing[0], "micron");
	segmentNuclei20x_imageReader_subContext->addParameter("Voxel spacing y", CSParameter::Double, &m_segmentNuclei20x_VoxelSpacing[1], "micron");
	segmentNuclei20x_imageReader_subContext->addParameter("Voxel spacing z", CSParameter::Double, &m_segmentNuclei20x_VoxelSpacing[2], "micron");
	segmentNuclei20x_imageReader_subContext->addParameter("Central vein segmentation", CSParameter::FileName, &m_segmentNuclei20x_CVSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentNuclei20x_imageReader_subContext->addParameter("Portal vein segmentation", CSParameter::FileName, &m_segmentNuclei20x_PVSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentNuclei20x_imageReader_subContext->addParameter("Is there a necrotic region", CSParameter::Bool, &m_segmentNuclei20x_UseNecReg, "");
	segmentNuclei20x_imageReader_subContext->addParameter("Necrotic Region", CSParameter::FileName, &m_segmentNuclei20x_NecReg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	CSParameterContext *segmentNuclei20x_preprocessing_subContext = m_segmentNuclei20xContext->addContext("1.1) Preprocess DAPI Channel");
	segmentNuclei20x_preprocessing_subContext->addParameter("Median filter kernel radius", CSParameter::Int, &m_segmentNuclei20x_Median_Radius, "pixel");
	CSParameterContext *segmentNuclei20x_binThreshold_subContext = m_segmentNuclei20xContext->addContext("1.2) Binary Threshold Filter on DAPI Channel on 1.1");
	segmentNuclei20x_binThreshold_subContext->addParameter("Threshold mode", CSParameter::Choice, m_segmentNuclei20x_DAPI_Threshold_Mode, "");
	CSParameterContext *segmentNuclei20x_adapOtsuThreshold_subContext = segmentNuclei20x_binThreshold_subContext->addContext("Adaptive Otsu-Threshold");
	segmentNuclei20x_adapOtsuThreshold_subContext->addParameter("Sample region size x", CSParameter::Int, &m_segmentNuclei20x_DAPI_AdapOtsuThreshold_RegionSize[0], "pixel");
	segmentNuclei20x_adapOtsuThreshold_subContext->addParameter("Sample region size y", CSParameter::Int, &m_segmentNuclei20x_DAPI_AdapOtsuThreshold_RegionSize[1], "pixel");
	segmentNuclei20x_adapOtsuThreshold_subContext->addParameter("Sample region size z", CSParameter::Int, &m_segmentNuclei20x_DAPI_AdapOtsuThreshold_RegionSize[2], "pixel");
	segmentNuclei20x_adapOtsuThreshold_subContext->addParameter("Number samples", CSParameter::Int, &m_segmentNuclei20x_DAPI_AdapOtsuThreshold_NumSamples, "");
	CSParameterContext *segmentNuclei20x_otsuThreshold_subContext = segmentNuclei20x_binThreshold_subContext->addContext("Otsu-Threshold");
	segmentNuclei20x_otsuThreshold_subContext->setVisible(false);
	CSParameterContext *segmentNuclei20x_manualThreshold_subContext = segmentNuclei20x_binThreshold_subContext->addContext("Manual Threshold");
	segmentNuclei20x_manualThreshold_subContext->addParameter("DAPI manual threshold", CSParameter::Int, &m_segmentNuclei20x_DAPI_ManualThreshold, "");
	segmentNuclei20x_manualThreshold_subContext->setVisible(false);
	CSParameterContext *segmentNuclei20x_holeFilling_subContext = m_segmentNuclei20xContext->addContext("1.3) Hole Filling on 1.2");
	segmentNuclei20x_holeFilling_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentNuclei20x_HoleFilling_RegionSize[0], "pixel");
	segmentNuclei20x_holeFilling_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentNuclei20x_HoleFilling_RegionSize[1], "pixel");
	segmentNuclei20x_holeFilling_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentNuclei20x_HoleFilling_RegionSize[2], "pixel");
	segmentNuclei20x_holeFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &m_segmentNuclei20x_HoleFilling_MajThreshold, "");
	CSParameterContext *segmentNuclei20x_opening_subContext = m_segmentNuclei20xContext->addContext("1.4) Opening on 1.3");
	segmentNuclei20x_opening_subContext->addParameter("Kernel radius", CSParameter::Int, &m_segmentNuclei20x_Opening_RadStrucElem, "pixel");
	CSParameterContext *segmentNuclei20x_watershed_subContext = m_segmentNuclei20xContext->addContext("1.5) Nuclei Separation on 1.4");
	segmentNuclei20x_watershed_subContext->addParameter("Alpha", CSParameter::Double, &m_segmentNuclei20x_WatershedFloodLevel, "");
	CSParameterContext *segmentNuclei20x_removeObjects_subContext = m_segmentNuclei20xContext->addContext("1.6) Remove objects on 1.5");
	segmentNuclei20x_removeObjects_subContext->addParameter("Smallest hepatocyte diameter", CSParameter::Double, &m_segmentNuclei20x_MinDiameterThreshold, "micron");
	segmentNuclei20x_removeObjects_subContext->addParameter("Biggest hepatocyte diameter", CSParameter::Double, &m_segmentNuclei20x_MaxDiameterThreshold, "micron");

	CSParameterContext *segmentNuclei20x_saveSelector_subContext = m_segmentNuclei20xContext->addContext("Save Selector");
	segmentNuclei20x_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_segmentNuclei20x_Log_Filename, "");
	segmentNuclei20x_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_segmentNuclei20x_Filename_Prefix, "");
	segmentNuclei20x_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_segmentNuclei20x_Save_Mode, "");
	//****************************************************************************************************************************


	//Segment Nuclei in 60x Datasets**********************************************************************************************
	m_segmentNuclei60x_DataDimMode = new CSParameterChoice(imageDataDim, 1);
	m_segmentNuclei60x_DataDimMode->setControlSubContexts(false);

	m_segmentNuclei60x_DAPI_Filename = "/path/to/DAPI_channel_file.tif";
	m_segmentNuclei60x_CVSeg_Filename = "/path/to/CV_bin_file.tif --optional";
	m_segmentNuclei60x_PVSeg_Filename = "/path/to/PV_bin_file.tif --optional";
	m_segmentNuclei60x_VoxelSpacing[0] = 0.207;     //IfADo: x,y: 0.207; z: 0.54
	m_segmentNuclei60x_VoxelSpacing[1] = 0.207;     //Iryna: x,y: 0.438; z: 0.2
	m_segmentNuclei60x_VoxelSpacing[2] = 0.54;
	m_segmentNuclei60x_Median_Radius = 0;
	m_segmentNuclei60x_Opening1_RadStrucElem = 0;
	m_segmentNuclei60x_DAPI_Threshold_Mode = new CSParameterChoice(thresholdModes, 0);
	m_segmentNuclei60x_DAPI_Threshold_Mode->setControlSubContexts(true);
	m_segmentNuclei60x_DAPI_AdapOtsuThreshold_RegionSize[0] = 100;
	m_segmentNuclei60x_DAPI_AdapOtsuThreshold_RegionSize[1] = 100;
	m_segmentNuclei60x_DAPI_AdapOtsuThreshold_RegionSize[2] = 10;
	m_segmentNuclei60x_DAPI_AdapOtsuThreshold_NumSamples = 10000;
	m_segmentNuclei60x_DAPI_ManualMinThreshold = 10;
	m_segmentNuclei60x_DAPI_ManualMaxThreshold = 255;
	m_segmentNuclei60x_InvHoleFilling_RegionSize[0] = 1;
	m_segmentNuclei60x_InvHoleFilling_RegionSize[1] = 1;
	m_segmentNuclei60x_InvHoleFilling_RegionSize[2] = 1;
	m_segmentNuclei60x_InvHoleFilling_MajThreshold = 4;
	m_segmentNuclei60x_CavityFilling_WithResampling = true;
	m_segmentNuclei60x_CavityFilling_Radius = 8;
	m_segmentNuclei60x_CavityFilling_MinThreshold = 0.75;
	m_segmentNuclei60x_CavityFilling_MaxThreshold = 1.0;
	m_segmentNuclei60x_Closing_xRadStrucElem = 2;
	m_segmentNuclei60x_Closing_yRadStrucElem = 2;
	m_segmentNuclei60x_Closing_zRadStrucElem = 1;
	m_segmentNuclei60x_Opening_xRadStrucElem = 2;
	m_segmentNuclei60x_Opening_yRadStrucElem = 2;
	m_segmentNuclei60x_Opening_zRadStrucElem = 1;
	m_segmentNuclei60x_MinNonHepDiameterThreshold = 2.5;
	m_segmentNuclei60x_MaxNonHepDiameterThreshold = 7.0;
	m_segmentNuclei60x_MinHepDiameterThreshold = 6.0;
	m_segmentNuclei60x_MaxHepDiameterThreshold = 20.0;
	m_segmentNuclei60x_roundness = 0.85;
	m_segmentNuclei60x_WatershedFloodLevel = 0.18;

	m_segmentNuclei60x_Log_Filename = "log";
	m_segmentNuclei60x_Filename_Prefix = "nuclei";
	m_segmentNuclei60x_Save_Mode = new CSParameterChoice(saveModes, 0);

	m_segmentNuclei60xContext = mpParameterRoot->addContext( PipelineName[PipelineSegmentNuclei60x] );
	m_segmentNuclei60xContext->setDebug(true);

	CSParameterContext *segmentNuclei60x_imageReader_subContext = m_segmentNuclei60xContext->addContext("Image Reader");
	segmentNuclei60x_imageReader_subContext->addParameter("Data dimensionality", CSParameter::Choice, m_segmentNuclei60x_DataDimMode, "");
	segmentNuclei60x_imageReader_subContext->addParameter("DAPI channel", CSParameter::FileName, &m_segmentNuclei60x_DAPI_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentNuclei60x_imageReader_subContext->addParameter("Voxel spacing x", CSParameter::Double, &m_segmentNuclei60x_VoxelSpacing[0], "micron");
	segmentNuclei60x_imageReader_subContext->addParameter("Voxel spacing y", CSParameter::Double, &m_segmentNuclei60x_VoxelSpacing[1], "micron");
	segmentNuclei60x_imageReader_subContext->addParameter("Voxel spacing z", CSParameter::Double, &m_segmentNuclei60x_VoxelSpacing[2], "micron");
	segmentNuclei60x_imageReader_subContext->addParameter("Central vein segmentation", CSParameter::FileName, &m_segmentNuclei60x_CVSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentNuclei60x_imageReader_subContext->addParameter("Portal vein segmentation", CSParameter::FileName, &m_segmentNuclei60x_PVSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	CSParameterContext *segmentNuclei60x_preprocessing_subContext = m_segmentNuclei60xContext->addContext("1.1) Preprocess DAPI Channel");
	segmentNuclei60x_preprocessing_subContext->addParameter("Median filter kernel radius", CSParameter::Int, &m_segmentNuclei60x_Median_Radius, "pixel");
	segmentNuclei60x_preprocessing_subContext->addParameter("Greyscale opening kernel radius", CSParameter::Int, &m_segmentNuclei60x_Opening1_RadStrucElem, "pixel");
	CSParameterContext *segmentNuclei60x_binThreshold_subContext = m_segmentNuclei60xContext->addContext("1.2) Binary Threshold on 1.1");
	segmentNuclei60x_binThreshold_subContext->addParameter("Threshold mode", CSParameter::Choice, m_segmentNuclei60x_DAPI_Threshold_Mode, "");
	CSParameterContext *segmentNuclei60x_adapOtsuThreshold_subContext = segmentNuclei60x_binThreshold_subContext->addContext("Adaptive Otsu-Threshold");
	segmentNuclei60x_adapOtsuThreshold_subContext->addParameter("Sample region size x", CSParameter::Int, &m_segmentNuclei60x_DAPI_AdapOtsuThreshold_RegionSize[0], "pixel");
	segmentNuclei60x_adapOtsuThreshold_subContext->addParameter("Sample region size y", CSParameter::Int, &m_segmentNuclei60x_DAPI_AdapOtsuThreshold_RegionSize[1], "pixel");
	segmentNuclei60x_adapOtsuThreshold_subContext->addParameter("Sample region size z", CSParameter::Int, &m_segmentNuclei60x_DAPI_AdapOtsuThreshold_RegionSize[2], "pixel");
	segmentNuclei60x_adapOtsuThreshold_subContext->addParameter("Number samples", CSParameter::Int, &m_segmentNuclei60x_DAPI_AdapOtsuThreshold_NumSamples, "");
	CSParameterContext *segmentNuclei60x_otsuThreshold_subContext = segmentNuclei60x_binThreshold_subContext->addContext("Otsu-Threshold");
	segmentNuclei60x_otsuThreshold_subContext->setVisible(false);
	CSParameterContext *segmentNuclei60x_manualThreshold_subContext = segmentNuclei60x_binThreshold_subContext->addContext("Manual Threshold");
	segmentNuclei60x_manualThreshold_subContext->addParameter("DAPI manual min threshold", CSParameter::Int, &m_segmentNuclei60x_DAPI_ManualMinThreshold, "");
	segmentNuclei60x_manualThreshold_subContext->addParameter("DAPI manual max threshold", CSParameter::Int, &m_segmentNuclei60x_DAPI_ManualMaxThreshold, "");
	segmentNuclei60x_manualThreshold_subContext->setVisible(false);
	CSParameterContext *segmentNuclei60x_invHoleFilling_subContext = m_segmentNuclei60xContext->addContext("1.3) Inverse Hole Filling on 1.2");
	segmentNuclei60x_invHoleFilling_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentNuclei60x_InvHoleFilling_RegionSize[0], "pixel");
	segmentNuclei60x_invHoleFilling_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentNuclei60x_InvHoleFilling_RegionSize[1], "pixel");
	segmentNuclei60x_invHoleFilling_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentNuclei60x_InvHoleFilling_RegionSize[2], "pixel");
	segmentNuclei60x_invHoleFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &m_segmentNuclei60x_InvHoleFilling_MajThreshold, "");
	CSParameterContext *segmentNuclei60x_cavityFilling_subContext = m_segmentNuclei60xContext->addContext("1.4) Cavity Filling on 1.3");
	segmentNuclei60x_cavityFilling_subContext->addParameter("Fast (less accurate)", CSParameter::Bool, &m_segmentNuclei60x_CavityFilling_WithResampling, "");
	segmentNuclei60x_cavityFilling_subContext->addParameter("Radius", CSParameter::Int, &m_segmentNuclei60x_CavityFilling_Radius, "pixel");
	segmentNuclei60x_cavityFilling_subContext->addParameter("Minimal fraction of surrounding foreground", CSParameter::Double, &m_segmentNuclei60x_CavityFilling_MinThreshold, "");
	segmentNuclei60x_cavityFilling_subContext->addParameter("Maximal fraction of surrounding foreground", CSParameter::Double, &m_segmentNuclei60x_CavityFilling_MaxThreshold, "");
	CSParameterContext *segmentNuclei60x_closing_subContext = m_segmentNuclei60xContext->addContext("1.5) Closing on 1.4");
	segmentNuclei60x_closing_subContext->addParameter("Kernel radius x", CSParameter::Int, &m_segmentNuclei60x_Closing_xRadStrucElem, "pixel");
	segmentNuclei60x_closing_subContext->addParameter("Kernel radius y", CSParameter::Int, &m_segmentNuclei60x_Closing_yRadStrucElem, "pixel");
	segmentNuclei60x_closing_subContext->addParameter("Kernel radius z", CSParameter::Int, &m_segmentNuclei60x_Closing_zRadStrucElem, "pixel");
	CSParameterContext *segmentNuclei60x_opening_subContext = m_segmentNuclei60xContext->addContext("1.6) Opening on 1.5");
	segmentNuclei60x_opening_subContext->addParameter("Kernel radius x", CSParameter::Int, &m_segmentNuclei60x_Opening_xRadStrucElem, "pixel");
	segmentNuclei60x_opening_subContext->addParameter("Kernel radius y", CSParameter::Int, &m_segmentNuclei60x_Opening_yRadStrucElem, "pixel");
	segmentNuclei60x_opening_subContext->addParameter("Kernel radius z", CSParameter::Int, &m_segmentNuclei60x_Opening_zRadStrucElem, "pixel");
	CSParameterContext *segmentNuclei60x_watershed_subContext = m_segmentNuclei60xContext->addContext("1.7) Nuclei Separation on 1.6");
	segmentNuclei60x_watershed_subContext->addParameter("Alpha", CSParameter::Double, &m_segmentNuclei60x_WatershedFloodLevel, "");
	CSParameterContext *segmentNuclei60x_removeObjects_subContext = m_segmentNuclei60xContext->addContext("1.8) Remove and group objects on 1.7");
	segmentNuclei60x_removeObjects_subContext->addParameter("Smallest non-hepatocyte diameter", CSParameter::Double, &m_segmentNuclei60x_MinNonHepDiameterThreshold, "micron");
	segmentNuclei60x_removeObjects_subContext->addParameter("Biggest non-hepatocyte diameter", CSParameter::Double, &m_segmentNuclei60x_MaxNonHepDiameterThreshold, "micron");
	segmentNuclei60x_removeObjects_subContext->addParameter("Smallest hepatocyte diameter", CSParameter::Double, &m_segmentNuclei60x_MinHepDiameterThreshold, "micron");
	segmentNuclei60x_removeObjects_subContext->addParameter("Biggest hepatocyte diameter", CSParameter::Double, &m_segmentNuclei60x_MaxHepDiameterThreshold, "micron");
	segmentNuclei60x_removeObjects_subContext->addParameter("Roundness of hepatocyte nuclei [0,1]", CSParameter::Double, &m_segmentNuclei60x_roundness, "");

	CSParameterContext *segmentNuclei60x_saveSelector_subContext = m_segmentNuclei60xContext->addContext("Save Selector");
	segmentNuclei60x_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_segmentNuclei60x_Log_Filename, "");
	segmentNuclei60x_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_segmentNuclei60x_Filename_Prefix, "");
	segmentNuclei60x_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_segmentNuclei60x_Save_Mode, "");
	//****************************************************************************************************************************


#ifndef CS_TI_QUANT_ONLY
	//Nuclei Segmentation with Hough**********************************************************************************************
	//TODO: no explicit size dependent filtering so far. Only by Hough.

	m_segmentNucleiWithHough_DAPI_Filename = "/path/to/DAPI_channel_file.tif";
	m_segmentNucleiWithHough_VoxelSpacing[0] = 0.438;
	m_segmentNucleiWithHough_VoxelSpacing[1] = 0.438;
	m_segmentNucleiWithHough_VoxelSpacing[2] = 0.2;
	m_segmentNucleiWithHough_Median_Radius = 3;
	m_segmentNucleiWithHough_Opening1_RadStrucElem = 2;
	m_segmentNucleiWithHough_Hough_NumSpheres = 400;
	m_segmentNucleiWithHough_Hough_Threshold = 90;
	m_segmentNucleiWithHough_Hough_GradientThreshold = 1;
	m_segmentNucleiWithHough_Hough_OutputThreshold = 0.1;
	m_segmentNucleiWithHough_Hough_MinRadius = 2.0;
	m_segmentNucleiWithHough_Hough_MaxRadius = 15.0;
	m_segmentNucleiWithHough_Hough_SigmaGradient = 3;
	m_segmentNucleiWithHough_Hough_Variance = 3;
	m_segmentNucleiWithHough_Hough_SphereRadiusRatio = 0.7;
	m_segmentNucleiWithHough_Hough_VotingRadiusRatio = 0.3;
	m_segmentNucleiWithHough_Hough_SamplingRatio = 1;
	m_segmentNucleiWithHough_SphereToImage_ResampleFactor = 4;
	m_segmentNucleiWithHough_SphereToImage_DilationRadius = 1;

	m_segmentNucleiWithHough_MinNonHepDiameterThreshold = 5.0;
	m_segmentNucleiWithHough_MinHepDiameterThreshold = 9.0;
	m_segmentNucleiWithHough_MaxHepDiameterThreshold = 20.0;

	m_segmentNucleiWithHough_Log_Filename = "log";
	m_segmentNucleiWithHough_Filename_Prefix = "nuclei";
	m_segmentNucleiWithHough_Save_Mode = new CSParameterChoice(saveModes, 0);

	m_segmentNucleiWithHoughContext = mpParameterRoot->addContext( PipelineName[PipelineSegmentNucleiHoughTransformation] );
	m_segmentNucleiWithHoughContext->setDebug(true);


	CSParameterContext *segmentNucleiWithHough_imageReader_subContext = m_segmentNucleiWithHoughContext->addContext("Image Reader");
	segmentNucleiWithHough_imageReader_subContext->addParameter("DAPI channel", CSParameter::FileName, &m_segmentNucleiWithHough_DAPI_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentNucleiWithHough_imageReader_subContext->addParameter("Voxel spacing x", CSParameter::Double, &m_segmentNucleiWithHough_VoxelSpacing[0], "micron");
	segmentNucleiWithHough_imageReader_subContext->addParameter("Voxel spacing y", CSParameter::Double, &m_segmentNucleiWithHough_VoxelSpacing[1], "micron");
	segmentNucleiWithHough_imageReader_subContext->addParameter("Voxel spacing z", CSParameter::Double, &m_segmentNucleiWithHough_VoxelSpacing[2], "micron");
	CSParameterContext *segmentNucleiWithHough_preprocessing_subContext = m_segmentNucleiWithHoughContext->addContext("1.1) Preprocess DAPI Channel");
	segmentNucleiWithHough_preprocessing_subContext->addParameter("Median filter radius", CSParameter::Int, &m_segmentNucleiWithHough_Median_Radius, "");
	segmentNucleiWithHough_preprocessing_subContext->addParameter("Greyscale opening radius", CSParameter::Int, &m_segmentNucleiWithHough_Opening1_RadStrucElem, "");
	CSParameterContext *segmentNucleiWithHough_hough_subContext = m_segmentNucleiWithHoughContext->addContext("1.2) Hough Transformation on 1.1");
	segmentNucleiWithHough_hough_subContext->addParameter("Maximal number of spheres", CSParameter::Int, &m_segmentNucleiWithHough_Hough_NumSpheres, "");
	segmentNucleiWithHough_hough_subContext->addParameter("Threshold on input image", CSParameter::Int, &m_segmentNucleiWithHough_Hough_Threshold, "");
	segmentNucleiWithHough_hough_subContext->addParameter("Gradient Threshold", CSParameter::Int, &m_segmentNucleiWithHough_Hough_GradientThreshold, "");
	segmentNucleiWithHough_hough_subContext->addParameter("Output Threshold", CSParameter::Double, &m_segmentNucleiWithHough_Hough_OutputThreshold, "");
	segmentNucleiWithHough_hough_subContext->addParameter("Minimal Radius", CSParameter::Double, &m_segmentNucleiWithHough_Hough_MinRadius, "");
	segmentNucleiWithHough_hough_subContext->addParameter("Maximal Radius", CSParameter::Double, &m_segmentNucleiWithHough_Hough_MaxRadius, "");
	segmentNucleiWithHough_hough_subContext->addParameter("Sigma Gradient", CSParameter::Double, &m_segmentNucleiWithHough_Hough_SigmaGradient, "micron");
	segmentNucleiWithHough_hough_subContext->addParameter("Variance", CSParameter::Double, &m_segmentNucleiWithHough_Hough_Variance, "micron");
	segmentNucleiWithHough_hough_subContext->addParameter("Sphere Radius Ratio [0,1]", CSParameter::Double, &m_segmentNucleiWithHough_Hough_SphereRadiusRatio, "");
	segmentNucleiWithHough_hough_subContext->addParameter("Voting Radius Ratio [0,1]", CSParameter::Double, &m_segmentNucleiWithHough_Hough_VotingRadiusRatio, "");
	segmentNucleiWithHough_hough_subContext->addParameter("Sampling Ratio [0,1]", CSParameter::Double, &m_segmentNucleiWithHough_Hough_SamplingRatio, "");
	CSParameterContext *segmentNucleiWithHough_houghToImage_subContext = m_segmentNucleiWithHoughContext->addContext("1.3) Projecting spheres to image 1.2");
	segmentNucleiWithHough_houghToImage_subContext->addParameter("Resample Factor", CSParameter::Int, &m_segmentNucleiWithHough_SphereToImage_ResampleFactor, "");
	segmentNucleiWithHough_houghToImage_subContext->addParameter("Dilation Radius", CSParameter::Int, &m_segmentNucleiWithHough_SphereToImage_DilationRadius, "");
//	CSParameterContext *segmentNucleiWithHough_removeObjects_subContext = m_segmentNucleiWithHoughContext->addContext("1.4) Remove and group objects on 1.3");
//	segmentNucleiWithHough_removeObjects_subContext->addParameter("Smallest non-hepatocyte diameter", CSParameter::Double, &m_segmentNucleiWithHough_MinNonHepDiameterThreshold, "micron");
//	segmentNucleiWithHough_removeObjects_subContext->addParameter("Smallest hepatocyte diameter", CSParameter::Double, &m_segmentNucleiWithHough_MinHepDiameterThreshold, "micron");
//	segmentNucleiWithHough_removeObjects_subContext->addParameter("Biggest hepatocyte diameter", CSParameter::Double, &m_segmentNucleiWithHough_MaxHepDiameterThreshold, "micron");

	CSParameterContext *segmentNucleiWithHough_saveSelector_subContext = m_segmentNucleiWithHoughContext->addContext("Save Selector");
	segmentNucleiWithHough_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_segmentNucleiWithHough_Log_Filename, "");
	segmentNucleiWithHough_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_segmentNucleiWithHough_Filename_Prefix, "");
	segmentNucleiWithHough_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_segmentNucleiWithHough_Save_Mode, "");
	//****************************************************************************************************************************

#endif // CS_TI_QUANT_ONLY

	//Estimate Cell Shape*********************************************************************************************************
	m_estimateCellShape_DataSetFileList = "/path/to/file.ias";
	m_estimateCellShape_DPPIV_Filename = "/path/to/dppiv_channel_file.tif";
	m_estimateCellShape_Nuclei_Filename = "/path/to/hep_nuclei_bin_file.tif";
	m_estimateCellShape_Bile_Filename = "/path/to/bile_bin_file.tif";
	m_estimateCellShape_Sinusoid_Filename = "/path/to/sin_bin_file.tif";
	m_estimateCellShape_CV_Filename = "/path/to/CV_bin_file.tif --optional";
	m_estimateCellShape_PV_Filename = "/path/to/PV_bin_file.tif --optional";
	m_estimateCellShape_SinusoidSkeleton_Filename = "/path/to/sin_skeleton_file.tif";
	m_estimateCellShape_UseNecReg = false;
	m_estimateCellShape_NecReg_Filename = "/path/to/necrotic_segmentation_file.tif";
	m_estimateCellShape_VoxelSpacing[0] = 0.207;
	m_estimateCellShape_VoxelSpacing[1] = 0.207; //20x: 0.620996; 60x: 0.207
	m_estimateCellShape_VoxelSpacing[2] = 0.54;
	m_estimateCellShape_Bile_Weight = 0.9;
	m_estimateCellShape_MinCellDiameterThreshold = 10.0;
	m_estimateCellShape_WatershedFloodLevel = 5;

	m_estimateCellShape_WithPositionFilter = false;
	m_estimateCellShape_NumberCells = 5;
	m_estimateCellShape_Pos[0] = 0;
	m_estimateCellShape_Pos[1] = 0;
	m_estimateCellShape_Pos[2] = 0;

	m_estimateCellShape_Log_Filename = "log";
	m_estimateCellShape_Filename_Prefix = "cellShape";
	m_estimateCellShape_Save_Mode = new CSParameterChoice(saveModes, 0);
	m_estimateCellShape_WriteCellOutlineToFile = false;
	m_estimateCellShape_WriteSinusoidOutlineToFile = false;
	m_estimateCellShape_WriteSinusoidGraphToFile = false;

	m_estimateCellShapeContext = mpParameterRoot->addContext( PipelineName[PipelineApproximateCellShape] );
	m_estimateCellShapeContext->setDebug(true);

	CSParameterContext *estimateCellShape_imageReader_subContext = m_estimateCellShapeContext->addContext("Image Reader");
	estimateCellShape_imageReader_subContext->addParameter("Dataset file list", CSParameter::FileName, &m_estimateCellShape_DataSetFileList, "")->setAttribute("FileTypes","Dataset file list (*.ias)");
	estimateCellShape_imageReader_subContext->addParameter("DPPIV channel", CSParameter::FileName, &m_estimateCellShape_DPPIV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateCellShape_imageReader_subContext->addParameter("Hepatic nuclei segmentation", CSParameter::FileName, &m_estimateCellShape_Nuclei_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateCellShape_imageReader_subContext->addParameter("Bile segmentation", CSParameter::FileName, &m_estimateCellShape_Bile_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateCellShape_imageReader_subContext->addParameter("Sinusoid segmentation", CSParameter::FileName, &m_estimateCellShape_Sinusoid_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateCellShape_imageReader_subContext->addParameter("Central vein segmentation", CSParameter::FileName, &m_estimateCellShape_CV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateCellShape_imageReader_subContext->addParameter("Portal vein segmentation", CSParameter::FileName, &m_estimateCellShape_PV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateCellShape_imageReader_subContext->addParameter("Is there a necrotic region", CSParameter::Bool, &m_estimateCellShape_UseNecReg, "");
	estimateCellShape_imageReader_subContext->addParameter("Necrotic region", CSParameter::FileName, &m_estimateCellShape_NecReg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateCellShape_imageReader_subContext->addParameter("Voxel spacing x", CSParameter::Double, &m_estimateCellShape_VoxelSpacing[0], "micron");
	estimateCellShape_imageReader_subContext->addParameter("Voxel spacing y", CSParameter::Double, &m_estimateCellShape_VoxelSpacing[1], "micron");
	estimateCellShape_imageReader_subContext->addParameter("Voxel spacing z", CSParameter::Double, &m_estimateCellShape_VoxelSpacing[2], "micron");
	CSParameterContext *estimateCellShape_parameter_subContext = m_estimateCellShapeContext->addContext("Parameter");
	estimateCellShape_parameter_subContext->addParameter("Bile weight [0,1]", CSParameter::Double, &m_estimateCellShape_Bile_Weight, "");
	estimateCellShape_parameter_subContext->addParameter("Alpha", CSParameter::Double, &m_estimateCellShape_WatershedFloodLevel, "");
	estimateCellShape_parameter_subContext->addParameter("Minimal cell diameter", CSParameter::Double, &m_estimateCellShape_MinCellDiameterThreshold, "micron");
	CSParameterContext *estimateCellShape_WithPositionFilter_subContext = m_estimateCellShapeContext->addContext("Position Filter");
	estimateCellShape_WithPositionFilter_subContext->addParameter("With position filter", CSParameter::Bool, &m_estimateCellShape_WithPositionFilter, "");
	estimateCellShape_WithPositionFilter_subContext->addParameter("Number of cells to keep", CSParameter::Int, &m_estimateCellShape_NumberCells, "");
	estimateCellShape_WithPositionFilter_subContext->addParameter("Coordinate x", CSParameter::Double, &m_estimateCellShape_Pos[0], "");
	estimateCellShape_WithPositionFilter_subContext->addParameter("Coordinate y", CSParameter::Double, &m_estimateCellShape_Pos[1], "");
	estimateCellShape_WithPositionFilter_subContext->addParameter("Coordinate z", CSParameter::Double, &m_estimateCellShape_Pos[2], "");

	CSParameterContext *estimateCellShape_saveSelector_subContext = m_estimateCellShapeContext->addContext("Save Selector");
	estimateCellShape_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_estimateCellShape_Log_Filename, "");
	estimateCellShape_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_estimateCellShape_Filename_Prefix, "");
	estimateCellShape_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_estimateCellShape_Save_Mode, "");
	estimateCellShape_saveSelector_subContext->addParameter("Write cell outline to file", CSParameter::Bool, &m_estimateCellShape_WriteCellOutlineToFile, "");
	estimateCellShape_saveSelector_subContext->addParameter("Write sinusoid outline to file", CSParameter::Bool, &m_estimateCellShape_WriteSinusoidOutlineToFile, "");
#ifndef CS_TI_QUANT_ONLY
	estimateCellShape_saveSelector_subContext->addParameter("Write sinusoid graph to file", CSParameter::Bool, &m_estimateCellShape_WriteSinusoidGraphToFile, "");
	estimateCellShape_saveSelector_subContext->addParameter("Sinusoid skeleton image", CSParameter::FileName, &m_estimateCellShape_SinusoidSkeleton_Filename, "");
#endif
	//****************************************************************************************************************************


#ifndef CS_TI_QUANT_ONLY

	//Estimate Lobule Shape*********************************************************************************************************
	m_estimateLobuleShape_DataDimMode = new CSParameterChoice(imageDataDim, 1);
	m_estimateLobuleShape_DataDimMode->setControlSubContexts(false);

	m_estimateLobuleShape_RawData_Filename = "/path/to/raw_data_file.tif";
	m_estimateLobuleShape_CV_Filename = "/path/to/CV_bin_file.tif";
	m_estimateLobuleShape_PV_Filename = "/path/to/PV_bin_file.tif";
	m_estimateLobuleShape_Tissue_Filename = "/path/to/tissue_bin_file.tif --optional";
	m_estimateLobuleShape_CV_Intensity = 50;
	m_estimateLobuleShape_PV_Intensity = 100;
	m_estimateLobuleShape_Tissue_Intensity = 150;
	m_estimateLobuleShape_Void_Intensity = 0;
	m_estimateLobuleShape_VoxelSpacing[0] = 0.46;
	m_estimateLobuleShape_VoxelSpacing[1] = 0.46; //20x: 0.620996; 60x: 0.207
	m_estimateLobuleShape_VoxelSpacing[2] = 0.54;
	m_estimateLobuleShape_PV_Weight = 0.5;
	m_estimateLobuleShape_MinLobuleDiameterThreshold = 200.0;
	m_estimateLobuleShape_MaxLobuleDiameterThreshold = 1100.0;
	m_estimateLobuleShape_WatershedFloodLevel = 5;

	m_estimateLobuleShape_AnalysisMode = new CSParameterChoice(analysisMode, 0);
	m_estimateLobuleShape_AnalysisMode->setControlSubContexts(true);
	m_estimateLobuleShape_DataSetID = "YourDataSetID";

	m_estimateLobuleShape_Log_Filename = "log";
	m_estimateLobuleShape_Filename_Prefix = "lobuleShape";
	m_estimateLobuleShape_Save_Mode = new CSParameterChoice(saveModes, 0);

	m_estimateLobuleShapeContext = mpParameterRoot->addContext( PipelineName[PipelineApproximateLobuleShape] );
	m_estimateLobuleShapeContext->setDebug(true);

	CSParameterContext *estimateLobuleShape_imageReader_subContext = m_estimateLobuleShapeContext->addContext("Image Reader");
	estimateLobuleShape_imageReader_subContext->addParameter("Data dimensionality", CSParameter::Choice, m_estimateLobuleShape_DataDimMode, "");
	estimateLobuleShape_imageReader_subContext->addParameter("Raw data", CSParameter::FileName, &m_estimateLobuleShape_RawData_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateLobuleShape_imageReader_subContext->addParameter("Central vein segmentation", CSParameter::FileName, &m_estimateLobuleShape_CV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateLobuleShape_imageReader_subContext->addParameter("Central vein intensity", CSParameter::Int, &m_estimateLobuleShape_CV_Intensity, "");
	estimateLobuleShape_imageReader_subContext->addParameter("Portal vein segmentation", CSParameter::FileName, &m_estimateLobuleShape_PV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateLobuleShape_imageReader_subContext->addParameter("Portal vein intensity", CSParameter::Int, &m_estimateLobuleShape_PV_Intensity, "");
	estimateLobuleShape_imageReader_subContext->addParameter("Tissue segmentation", CSParameter::FileName, &m_estimateLobuleShape_Tissue_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	estimateLobuleShape_imageReader_subContext->addParameter("Tissue intensity", CSParameter::Int, &m_estimateLobuleShape_Tissue_Intensity, "");
	estimateLobuleShape_imageReader_subContext->addParameter("Void intensity", CSParameter::Int, &m_estimateLobuleShape_Void_Intensity, "");
	estimateLobuleShape_imageReader_subContext->addParameter("Voxel spacing x", CSParameter::Double, &m_estimateLobuleShape_VoxelSpacing[0], "micron");
	estimateLobuleShape_imageReader_subContext->addParameter("Voxel spacing y", CSParameter::Double, &m_estimateLobuleShape_VoxelSpacing[1], "micron");
	estimateLobuleShape_imageReader_subContext->addParameter("Voxel spacing z", CSParameter::Double, &m_estimateLobuleShape_VoxelSpacing[2], "micron");
	CSParameterContext *estimateLobuleShape_parameter_subContext = m_estimateLobuleShapeContext->addContext("Parameter");
	estimateLobuleShape_parameter_subContext->addParameter("Portal vein weight [0,1]", CSParameter::Double, &m_estimateLobuleShape_PV_Weight, "");
	estimateLobuleShape_parameter_subContext->addParameter("Alpha", CSParameter::Double, &m_estimateLobuleShape_WatershedFloodLevel, "");
	estimateLobuleShape_parameter_subContext->addParameter("Minimal lobule diameter", CSParameter::Double, &m_estimateLobuleShape_MinLobuleDiameterThreshold, "micron");
	estimateLobuleShape_parameter_subContext->addParameter("Maximal lobule diameter", CSParameter::Double, &m_estimateLobuleShape_MaxLobuleDiameterThreshold, "micron");

	CSParameterContext *estimateLobuleShape_analysis_subContext = m_estimateLobuleShapeContext->addContext("Lobule analysis");
	estimateLobuleShape_analysis_subContext->addParameter("Analysis mode", CSParameter::Choice, m_estimateLobuleShape_AnalysisMode, "");
	CSParameterContext *estimateLobuleShape_WithoutAnalysis_subContext = estimateLobuleShape_analysis_subContext->addContext("Without Analysis");
	CSParameterContext *estimateLobuleShape_WithAnalysis_subContext = estimateLobuleShape_analysis_subContext->addContext("With Analysis");
	estimateLobuleShape_WithAnalysis_subContext->addParameter("Specify data set name", CSParameter::String, &m_estimateLobuleShape_DataSetID, "");
	estimateLobuleShape_WithAnalysis_subContext->setVisible(false);

	CSParameterContext *estimateLobuleShape_saveSelector_subContext = m_estimateLobuleShapeContext->addContext("Save Selector");
	estimateLobuleShape_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_estimateLobuleShape_Log_Filename, "");
	estimateLobuleShape_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_estimateLobuleShape_Filename_Prefix, "");
	estimateLobuleShape_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_estimateLobuleShape_Save_Mode, "");
	//****************************************************************************************************************************

	//Segment Stellate Cells******************************************************************************************************
	m_segmentStellateCells_Desmin_Filename = "/path/to/Desmin_channel_file.tif";
	m_segmentStellateCells_DAPI_Filename = "/path/to/DAPI_channel_file.tif";
	m_segmentStellateCells_NonHepNuclei_Filename = "/path/to/non_hep_nuclei_bin_file.tif";
	m_segmentStellateCells_HepNuclei_Filename = "/path/to/hep_nuclei_bin_file.tif";
	m_segmentStellateCells_VoxelSpacing[0] = 0.438;
	m_segmentStellateCells_VoxelSpacing[1] = 0.438;
	m_segmentStellateCells_VoxelSpacing[2] = 0.2;
	m_segmentStellateCells_Median_Radius = 2;
	m_segmentStellateCells_Opening_RadStrucElem = 2;
	m_segmentStellateCells_WithConvexFilter = false;
	m_segmentStellateCells_ConvexFilter_Height = 170;
	m_segmentStellateCells_Threshold_Mode = new CSParameterChoice(thresholdModes, 0);
	m_segmentStellateCells_Threshold_Mode->setControlSubContexts(true);
	m_segmentStellateCells_AdapOtsuThreshold_RegionSize[0] = 100;
	m_segmentStellateCells_AdapOtsuThreshold_RegionSize[1] = 100;
	m_segmentStellateCells_AdapOtsuThreshold_RegionSize[2] = 10;
	m_segmentStellateCells_AdapOtsuThreshold_NumSamples = 10000;
	m_segmentStellateCells_Desmin_ManualThreshold = 5;
	m_segmentStellateCells_InvHoleFilling_RegionSize[0] = 1;
	m_segmentStellateCells_InvHoleFilling_RegionSize[1] = 1;
	m_segmentStellateCells_InvHoleFilling_RegionSize[2] = 1;
	m_segmentStellateCells_InvHoleFilling_MajThreshold = 4;
	m_segmentStellateCells_MinSizeThreshold = 1000;
	m_segmentStellateCells_NucleiBodySuperposition = 0.2;

	m_segmentStellateCells_Log_Filename = "log";
	m_segmentStellateCells_Filename_Prefix = "stellate";
	m_segmentStellateCells_Save_Mode = new CSParameterChoice(saveModes, 0);

	m_segmentStellateCellsContext = mpParameterRoot->addContext( PipelineName[PipelineSegmentStellateCellCytoskeleton] );
	m_segmentStellateCellsContext->setDebug(true);

	CSParameterContext *segmentStellateCells_imageReader_subContext = m_segmentStellateCellsContext->addContext("Image Reader");
	segmentStellateCells_imageReader_subContext->addParameter("Desmin channel", CSParameter::FileName, &m_segmentStellateCells_Desmin_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentStellateCells_imageReader_subContext->addParameter("DAPI channel", CSParameter::FileName, &m_segmentStellateCells_DAPI_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentStellateCells_imageReader_subContext->addParameter("Non-Hepatic nuclei segmentation", CSParameter::FileName, &m_segmentStellateCells_NonHepNuclei_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentStellateCells_imageReader_subContext->addParameter("Hepatic nuclei segmentation", CSParameter::FileName, &m_segmentStellateCells_HepNuclei_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentStellateCells_imageReader_subContext->addParameter("Voxel spacing x", CSParameter::Double, &m_segmentStellateCells_VoxelSpacing[0], "micron");
	segmentStellateCells_imageReader_subContext->addParameter("Voxel spacing y", CSParameter::Double, &m_segmentStellateCells_VoxelSpacing[1], "micron");
	segmentStellateCells_imageReader_subContext->addParameter("Voxel spacing z", CSParameter::Double, &m_segmentStellateCells_VoxelSpacing[2], "micron");
	CSParameterContext *segmentStellateCells_preprocessing_subContext = m_segmentStellateCellsContext->addContext("1.1) Preprocess Desmin Channel");
	segmentStellateCells_preprocessing_subContext->addParameter("Median filter radius", CSParameter::Int, &m_segmentStellateCells_Median_Radius, "");
	segmentStellateCells_preprocessing_subContext->addParameter("Greyscale opening radius", CSParameter::Int, &m_segmentStellateCells_Opening_RadStrucElem, "");
	segmentStellateCells_preprocessing_subContext->addParameter("With convex image filter", CSParameter::Bool, &m_segmentStellateCells_WithConvexFilter, "");
	segmentStellateCells_preprocessing_subContext->addParameter("Convex image filter height", CSParameter::Int, &m_segmentStellateCells_ConvexFilter_Height, "");
	CSParameterContext *segmentStellateCells_binThreshold_subContext = m_segmentStellateCellsContext->addContext("1.2) Binary Threshold on 1.1");
	segmentStellateCells_binThreshold_subContext->addParameter("Threshold mode", CSParameter::Choice, m_segmentStellateCells_Threshold_Mode, "");
	CSParameterContext *segmentStellateCells_adapOtsuThreshold_subContext = segmentStellateCells_binThreshold_subContext->addContext("Adaptive Otsu-Threshold");
	segmentStellateCells_adapOtsuThreshold_subContext->addParameter("Sample region size x", CSParameter::Int, &m_segmentStellateCells_AdapOtsuThreshold_RegionSize[0], "");
	segmentStellateCells_adapOtsuThreshold_subContext->addParameter("Sample region size y", CSParameter::Int, &m_segmentStellateCells_AdapOtsuThreshold_RegionSize[1], "");
	segmentStellateCells_adapOtsuThreshold_subContext->addParameter("Sample region size z", CSParameter::Int, &m_segmentStellateCells_AdapOtsuThreshold_RegionSize[2], "");
	segmentStellateCells_adapOtsuThreshold_subContext->addParameter("Number samples", CSParameter::Int, &m_segmentStellateCells_AdapOtsuThreshold_NumSamples, "");
	CSParameterContext *segmentStellateCells_otsuThreshold_subContext = segmentStellateCells_binThreshold_subContext->addContext("Otsu-Threshold");
	segmentStellateCells_otsuThreshold_subContext->setVisible(false);
	CSParameterContext *segmentStellateCells_manualThreshold_subContext = segmentStellateCells_binThreshold_subContext->addContext("Manual Threshold");
	segmentStellateCells_manualThreshold_subContext->addParameter("Desmin manual threshold", CSParameter::Int, &m_segmentStellateCells_Desmin_ManualThreshold, "");
	segmentStellateCells_manualThreshold_subContext->setVisible(false);
	CSParameterContext *segmentStellateCells_invHoleFilling_subContext = m_segmentStellateCellsContext->addContext("1.3) Inverse Hole Filling on 1.2");
	segmentStellateCells_invHoleFilling_subContext->addParameter("Radius x", CSParameter::Int, &m_segmentStellateCells_InvHoleFilling_RegionSize[0], "");
	segmentStellateCells_invHoleFilling_subContext->addParameter("Radius y", CSParameter::Int, &m_segmentStellateCells_InvHoleFilling_RegionSize[1], "");
	segmentStellateCells_invHoleFilling_subContext->addParameter("Radius z", CSParameter::Int, &m_segmentStellateCells_InvHoleFilling_RegionSize[2], "");
	segmentStellateCells_invHoleFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &m_segmentStellateCells_InvHoleFilling_MajThreshold, "");
	CSParameterContext *segmentStellateCells_removeObjects_subContext = m_segmentStellateCellsContext->addContext("1.4) Remove objects on 1.3");
	segmentStellateCells_removeObjects_subContext->addParameter("Minimal structure size", CSParameter::Int, &m_segmentStellateCells_MinSizeThreshold, "pixel");
	CSParameterContext *segmentStellateCells_classifyStellateNuclei_subContext = m_segmentStellateCellsContext->addContext("1.5) Classify stellate nuclei on 1.4");
	segmentStellateCells_classifyStellateNuclei_subContext->addParameter("Non-hepatic nuclei - stellate body superposition", CSParameter::Double, &m_segmentStellateCells_NucleiBodySuperposition, "percent");

	CSParameterContext *segmentStellateCells_saveSelector_subContext = m_segmentStellateCellsContext->addContext("Save Selector");
	segmentStellateCells_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_segmentStellateCells_Log_Filename, "");
	segmentStellateCells_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_segmentStellateCells_Filename_Prefix, "");
	segmentStellateCells_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_segmentStellateCells_Save_Mode, "");
	//****************************************************************************************************************************

#endif // CS_TI_QUANT_ONLY

	//Analyze Cells***************************************************************************************************************
	mAnalyzeCellsDataSetID = "Dataset ID";
	mAnalyzeCellsDataSetFileList = "/path/to/file.ias";
	mAnalyzeCellsCellShapeFilename = "/path/to/cellShape_bin_file.tif";
	mAnalyzeCellsNucleiFilename = "/path/to/nuclei_bin_file.tif";
	mAnalyzeCellsSinusoidFilename = "/path/to/sinusoid_bin_file.tif";
	mAnalyzeCellsBileFilename = "/path/to/bile_bin_file.tif";
	mAnalyzeCellsCVFilename = "/path/to/CV_bin_file.tif --optional";
	mAnalyzeCellsPVFilename = "/path/to/PV_bin_file.tif --optional";

	mAnalyzeCellsVoxSpacing[0] = 0.46;
	mAnalyzeCellsVoxSpacing[1] = 0.46;
	mAnalyzeCellsVoxSpacing[2] = 0.54;

	mAnalyzeCellsSaveAsGraph = false;

	mAnalyzeCellsContext = mpParameterRoot->addContext( PipelineName[PipelineAnalyzeCells] );
	mAnalyzeCellsContext->setDebug(true);

	mAnalyzeCellsContext->addParameter("Dataset ID", CSParameter::String, &mAnalyzeCellsDataSetID, "");
	mAnalyzeCellsContext->addParameter("Dataset file list", CSParameter::FileName, &mAnalyzeCellsDataSetFileList, "")->setAttribute("FileTypes","Dataset file list (*.ias)");
	mAnalyzeCellsContext->addParameter("Cell shape binary image", CSParameter::FileName, &mAnalyzeCellsCellShapeFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeCellsContext->addParameter("Nuclei segmentation", CSParameter::FileName, &mAnalyzeCellsNucleiFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeCellsContext->addParameter("Sinusoid segmentation", CSParameter::FileName, &mAnalyzeCellsSinusoidFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeCellsContext->addParameter("Bile segmentation", CSParameter::FileName, &mAnalyzeCellsBileFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeCellsContext->addParameter("CV segmentation", CSParameter::FileName, &mAnalyzeCellsCVFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeCellsContext->addParameter("PV segmentation", CSParameter::FileName, &mAnalyzeCellsPVFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeCellsContext->addParameter("Voxel spacing x", CSParameter::Double, &mAnalyzeCellsVoxSpacing[0], "micron");
	mAnalyzeCellsContext->addParameter("Voxel spacing y", CSParameter::Double, &mAnalyzeCellsVoxSpacing[1], "micron");
	mAnalyzeCellsContext->addParameter("Voxel spacing z", CSParameter::Double, &mAnalyzeCellsVoxSpacing[2], "micron");
	mAnalyzeCellsContext->addParameter("Save analysis result as graph", CSParameter::Bool, &mAnalyzeCellsSaveAsGraph, "");
	//****************************************************************************************************************************


#ifndef CS_TI_QUANT_ONLY

	//Analyze Nuclei**************************************************************************************************************
	mAnalyzeNucleiDataDimMode = new CSParameterChoice(imageDataDim, 1);
	mAnalyzeNucleiDataDimMode->setControlSubContexts(false);
	mAnalyzeNucleiDataSetID = "Dataset ID";
	mAnalyzeNucleiHepNucleiFilename = "/path/to/hepNuclei_bin_file.tif";
	mAnalyzeNucleiNonHepNucleiFilename = "/path/to/nonHepNuclei_bin_file.tif";
	mAnalyzeNucleiNonProlifNucleiIntensity = 200;
	mAnalyzeNucleiProlifNucleiIntensity = 100;
	mAnalyzeNucleiCVFilename = "/path/to/CV_bin_file.tif --optional";
	mAnalyzeNucleiCVIntensity = 50;
	mAnalyzeNucleiPVFilename = "/path/to/PV_bin_file.tif --optional";
	mAnalyzeNucleiPVIntensity = 100;
	mAnalyzeNucleiTissueFilename = "/path/to/tissue_bin_file.tif --optional";
	mAnalyzeNucleiTissueIntensity = 150;
	mAnalyzeNucleiNecRegionFilename = "/path/to/necroticRegion_bin_file.tif --optional";
	mAnalyzeNucleiNecRegionIntensity = 200;
	mAnalyzeNucleiVoidIntensity = 0;
	mAnalyzeNucleiLobuleFilename = "/path/to/lobule_bin_file.tif --optional";
	mAnalyzeNucleiLobulesIntensity = 150;
	mAnalyzeNucleiSaveAsGraph = false;

	mAnalyzeNucleiVoxSpacing[0] = 0.46;
	mAnalyzeNucleiVoxSpacing[1] = 0.46;
	mAnalyzeNucleiVoxSpacing[2] = 0.54;

	mAnalyzeNucleiContext = mpParameterRoot->addContext( PipelineName[PipelineAnalyzeNuclei] );
	mAnalyzeNucleiContext->setDebug(true);

	mAnalyzeNucleiContext->addParameter("Data dimensionality", CSParameter::Choice, mAnalyzeNucleiDataDimMode, "");
	mAnalyzeNucleiContext->addParameter("Dataset ID", CSParameter::String, &mAnalyzeNucleiDataSetID, "");
	mAnalyzeNucleiContext->addParameter("Hepatic nuclei segmentation", CSParameter::FileName, &mAnalyzeNucleiHepNucleiFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeNucleiContext->addParameter("Non-hepatic nuclei segmentation", CSParameter::FileName, &mAnalyzeNucleiNonHepNucleiFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeNucleiContext->addParameter("Non-proliferating nuclei intensity", CSParameter::Int, &mAnalyzeNucleiNonProlifNucleiIntensity, "");
	mAnalyzeNucleiContext->addParameter("Proliferating nuclei intensity", CSParameter::Int, &mAnalyzeNucleiProlifNucleiIntensity, "");
	mAnalyzeNucleiContext->addParameter("CV segmentation", CSParameter::FileName, &mAnalyzeNucleiCVFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeNucleiContext->addParameter("CV intensity", CSParameter::Int, &mAnalyzeNucleiCVIntensity, "");
	mAnalyzeNucleiContext->addParameter("PV segmentation", CSParameter::FileName, &mAnalyzeNucleiPVFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeNucleiContext->addParameter("PV intensity", CSParameter::Int, &mAnalyzeNucleiPVIntensity, "");
	mAnalyzeNucleiContext->addParameter("Tissue segmentation", CSParameter::FileName, &mAnalyzeNucleiTissueFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeNucleiContext->addParameter("Tissue intensity", CSParameter::Int, &mAnalyzeNucleiTissueIntensity, "");
	mAnalyzeNucleiContext->addParameter("Necrotic Region segmentation", CSParameter::FileName, &mAnalyzeNucleiNecRegionFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeNucleiContext->addParameter("Necrotic Region intensity", CSParameter::Int, &mAnalyzeNucleiNecRegionIntensity, "");
	mAnalyzeNucleiContext->addParameter("Void intensity", CSParameter::Int, &mAnalyzeNucleiVoidIntensity, "");
	mAnalyzeNucleiContext->addParameter("Lobule segmentation", CSParameter::FileName, &mAnalyzeNucleiLobuleFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeNucleiContext->addParameter("Lobule intensity", CSParameter::Int, &mAnalyzeNucleiLobulesIntensity, "");
	mAnalyzeNucleiContext->addParameter("Voxel spacing x", CSParameter::Double, &mAnalyzeNucleiVoxSpacing[0], "micron");
	mAnalyzeNucleiContext->addParameter("Voxel spacing y", CSParameter::Double, &mAnalyzeNucleiVoxSpacing[1], "micron");
	mAnalyzeNucleiContext->addParameter("Voxel spacing z", CSParameter::Double, &mAnalyzeNucleiVoxSpacing[2], "micron");
	mAnalyzeNucleiContext->addParameter("Save analysis result as graph", CSParameter::Bool, &mAnalyzeNucleiSaveAsGraph, "");
	//****************************************************************************************************************************

	//Analyze Lobule**************************************************************************************************************
	mAnalyzeLobuleCVFilename = "/path/to/CV_bin_file.tif";
	mAnalyzeLobulePVFilename = "/path/to/PV_bin_file.tif";
	mAnalyzeLobuleVoxSpacing[0] = 0.621;
	mAnalyzeLobuleVoxSpacing[1] = 0.621;
	mAnalyzeLobuleVoxSpacing[2] = 0.54;
	mAnalyzeLobuleSaveMode = new CSParameterChoice(saveModes, 1);

	mAnalyzeLobuleContext = mpParameterRoot->addContext( PipelineName[PipelineAnalyzeLobule] );
	mAnalyzeLobuleContext->setDebug(true);

	mAnalyzeLobuleContext->addParameter("Central vein segmentation", CSParameter::FileName, &mAnalyzeLobuleCVFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeLobuleContext->addParameter("Portal vein segmentation", CSParameter::FileName, &mAnalyzeLobulePVFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeLobuleContext->addParameter("Voxel spacing x", CSParameter::Double, &mAnalyzeLobuleVoxSpacing[0], "micron");
	mAnalyzeLobuleContext->addParameter("Voxel spacing y", CSParameter::Double, &mAnalyzeLobuleVoxSpacing[1], "micron");
	mAnalyzeLobuleContext->addParameter("Voxel spacing z", CSParameter::Double, &mAnalyzeLobuleVoxSpacing[2], "micron");
	mAnalyzeLobuleContext->addParameter("Save mode", CSParameter::Choice, mAnalyzeLobuleSaveMode, "");
	//****************************************************************************************************************************

	//Analyze Hoechst Diffusion***************************************************************************************************
	mAnalyzeHoechstDiffDataSetID = "Mov_1";
	mAnalyzeHoechstDiffFilename = "/home/ls1/friebel/Workspace/Data/Confocal/Bile_Canaliculi_Analysis/HoechstSpheroid/1000_CELLS_40UM INTO THE TISSUE_Subset_t015_t000.tif";

	mAnalyzeHoechstDiffVoxSpacing[0] = 0.692;
	mAnalyzeHoechstDiffVoxSpacing[1] = 0.692;
	mAnalyzeHoechstDiffTimeStep = 31.8;

	mAnalyzeHoechstDiffContext = mpParameterRoot->addContext( PipelineName[PipelineAnalyzeHoechstDiffusion] );
	mAnalyzeHoechstDiffContext->setDebug(true);

	mAnalyzeHoechstDiffContext->addParameter("Dataset ID", CSParameter::String, &mAnalyzeHoechstDiffDataSetID, "");
	mAnalyzeHoechstDiffContext->addParameter("First dataset slide", CSParameter::FileName, &mAnalyzeHoechstDiffFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	mAnalyzeHoechstDiffContext->addParameter("Voxel spacing x", CSParameter::Double, &mAnalyzeHoechstDiffVoxSpacing[0], "micron");
	mAnalyzeHoechstDiffContext->addParameter("Voxel spacing y", CSParameter::Double, &mAnalyzeHoechstDiffVoxSpacing[1], "micron");
	mAnalyzeHoechstDiffContext->addParameter("Time step", CSParameter::Double, &mAnalyzeHoechstDiffTimeStep, "seconds");
	//****************************************************************************************************************************

	//Analyze Stellate Cells******************************************************************************************************
	m_analyzeStellateCells_dataSetID = "YourDataSetID";
    m_analyzeStellateCells_HepSeg_Filename = "/path/to/hep_nuclei_bin_file.tif";
    m_analyzeStellateCells_NonHepSeg_Filename = "/path/to/nonHep_nuclei_bin_file.tif";
    m_analyzeStellateCells_StCSeg_Filename = "/path/to/stc_nuclei_bin_file.tif";
    m_analyzeStellateCells_StCBodies_Filename = "/path/to/stc_bin_file.tif";
    m_analyzeStellateCells_Hepatocyte_Filename = "/path/to/hep_bin_file.tif";
    m_analyzeStellateCells_CV_Filename = "/path/to/cv_bin_file.tif";
    m_analyzeStellateCells_PV_Filename = "/path/to/pv_bin_file.tif";

    m_analyzeStellateCells_VoxSpacing[0] = 0.438;
    m_analyzeStellateCells_VoxSpacing[1] = 0.438;
    m_analyzeStellateCells_VoxSpacing[2] = 0.2;

    m_analyzeStellateCellsContext = mpParameterRoot->addContext( PipelineName[PipelineAnalyzeStellateCells] );
    m_analyzeStellateCellsContext->setDebug(true);

    m_analyzeStellateCellsContext->addParameter("Dataset ID", CSParameter::String, &m_analyzeStellateCells_dataSetID, "");
    m_analyzeStellateCellsContext->addParameter("Hepatic nuclei segmentation", CSParameter::FileName, &m_analyzeStellateCells_HepSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    m_analyzeStellateCellsContext->addParameter("Non-Hepatic nuclei segmentation", CSParameter::FileName, &m_analyzeStellateCells_NonHepSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    m_analyzeStellateCellsContext->addParameter("Stellate nuclei segmentation", CSParameter::FileName, &m_analyzeStellateCells_StCSeg_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    m_analyzeStellateCellsContext->addParameter("Stellate cell body segmentation", CSParameter::FileName, &m_analyzeStellateCells_StCBodies_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    m_analyzeStellateCellsContext->addParameter("Hepatocyte segmentation", CSParameter::FileName, &m_analyzeStellateCells_Hepatocyte_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    m_analyzeStellateCellsContext->addParameter("Central Vein segmentation", CSParameter::FileName, &m_analyzeStellateCells_CV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    m_analyzeStellateCellsContext->addParameter("Portal Vein segmentation", CSParameter::FileName, &m_analyzeStellateCells_PV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    m_analyzeStellateCellsContext->addParameter("Voxel spacing x", CSParameter::Double, &m_analyzeStellateCells_VoxSpacing[0], "micron");
    m_analyzeStellateCellsContext->addParameter("Voxel spacing y", CSParameter::Double, &m_analyzeStellateCells_VoxSpacing[1], "micron");
    m_analyzeStellateCellsContext->addParameter("Voxel spacing z", CSParameter::Double, &m_analyzeStellateCells_VoxSpacing[2], "micron");
    //****************************************************************************************************************************

#endif // CS_TI_QUANT_ONLY


	//Segment Necrotic Region*****************************************************************************************************
	m_segmentNecroticRegion_DPPIV_Filename = "/path/to/DPPIV_channel_file.tif";
	m_segmentNecroticRegion_DMs_Filename = "/path/to/DM_channel_file.tif";
	m_segmentNecroticRegion_DMs_Threshold_Mode = new CSParameterChoice(thresholdModes, 0);
	m_segmentNecroticRegion_DMs_Threshold_Mode->setControlSubContexts(true);
	m_segmentNecroticRegion_DMs_AdapOtsuThreshold_RegionSize[0] = 100;
	m_segmentNecroticRegion_DMs_AdapOtsuThreshold_RegionSize[1] = 100;
	m_segmentNecroticRegion_DMs_AdapOtsuThreshold_RegionSize[2] = 10;
	m_segmentNecroticRegion_DMs_AdapOtsuThreshold_NumSamples = 10000;
	m_segmentNecroticRegion_DMs_ManualThreshold = 10;
	m_segmentNecroticRegion_Erode1_RadStrucElem = 2;
	m_segmentNecroticRegion_Dilate1_RadStrucElem = 5;
	m_segmentNecroticRegion_Erode2_RadStrucElem = 8;
	m_segmentNecroticRegion_Dilate2_RadStrucElem = 10;
	m_segmentNecroticRegion_RemvObjSizeThreshold = 500000;

	m_segmentNecroticRegion_Log_Filename = "log";
	m_segmentNecroticRegion_Filename_Prefix = "necroticRegion";
	m_segmentNecroticRegion_Save_Mode = new CSParameterChoice(saveModes, 0);

	m_segmentNecroticRegionContext = mpParameterRoot->addContext( PipelineName[PipelineSegmentNecroticRegion] );
	m_segmentNecroticRegionContext->setDebug(true);

	CSParameterContext *segmentNecroticRegion_imageReader_subContext = m_segmentNecroticRegionContext->addContext("Image Reader");
	segmentNecroticRegion_imageReader_subContext->addParameter("DPPIV channel", CSParameter::FileName, &m_segmentNecroticRegion_DPPIV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentNecroticRegion_imageReader_subContext->addParameter("DMs channel", CSParameter::FileName, &m_segmentNecroticRegion_DMs_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	CSParameterContext *segmentNecroticRegion_binThreshold_subContext = m_segmentNecroticRegionContext->addContext("1.1) Binary Threshold Filter on DMs Channel");
	segmentNecroticRegion_binThreshold_subContext->addParameter("Threshold mode", CSParameter::Choice, m_segmentNecroticRegion_DMs_Threshold_Mode, "");
	CSParameterContext *segmentNecroticRegion_adapOtsuThreshold_subContext = segmentNecroticRegion_binThreshold_subContext->addContext("Adaptive Otsu-Threshold");
	segmentNecroticRegion_adapOtsuThreshold_subContext->addParameter("Sample region size x", CSParameter::Int, &m_segmentNecroticRegion_DMs_AdapOtsuThreshold_RegionSize[0], "pixel");
	segmentNecroticRegion_adapOtsuThreshold_subContext->addParameter("Sample region size y", CSParameter::Int, &m_segmentNecroticRegion_DMs_AdapOtsuThreshold_RegionSize[1], "pixel");
	segmentNecroticRegion_adapOtsuThreshold_subContext->addParameter("Sample region size z", CSParameter::Int, &m_segmentNecroticRegion_DMs_AdapOtsuThreshold_RegionSize[2], "pixel");
	segmentNecroticRegion_adapOtsuThreshold_subContext->addParameter("Number samples", CSParameter::Int, &m_segmentNecroticRegion_DMs_AdapOtsuThreshold_NumSamples, "");
	CSParameterContext *segmentNecroticRegion_otsuThreshold_subContext = segmentNecroticRegion_binThreshold_subContext->addContext("Otsu-Threshold");
	segmentNecroticRegion_otsuThreshold_subContext->setVisible(false);
	CSParameterContext *segmentNecroticRegion_manualThreshold_subContext = segmentNecroticRegion_binThreshold_subContext->addContext("Manual Threshold");
	segmentNecroticRegion_manualThreshold_subContext->addParameter("DMs manual threshold", CSParameter::Int, &m_segmentNecroticRegion_DMs_ManualThreshold, "");
	segmentNecroticRegion_manualThreshold_subContext->setVisible(false);
	CSParameterContext *segmentNecroticRegion_erode1_subContext = m_segmentNecroticRegionContext->addContext("1.2) Erode A on 1.1");
	segmentNecroticRegion_erode1_subContext->addParameter("Kernel radius", CSParameter::Int, &m_segmentNecroticRegion_Erode1_RadStrucElem, "pixel");
	CSParameterContext *segmentNecroticRegion_dilate1_subContext = m_segmentNecroticRegionContext->addContext("1.3) Dilate A on 1.2");
	segmentNecroticRegion_dilate1_subContext->addParameter("Kernel radius", CSParameter::Int, &m_segmentNecroticRegion_Dilate1_RadStrucElem, "pixel");
	CSParameterContext *segmentNecroticRegion_erode2_subContext = m_segmentNecroticRegionContext->addContext("1.4) Erode B on 1.3");
	segmentNecroticRegion_erode2_subContext->addParameter("Kernel radius", CSParameter::Int, &m_segmentNecroticRegion_Erode2_RadStrucElem, "pixel");
	CSParameterContext *segmentNecroticRegion_dilate2_subContext = m_segmentNecroticRegionContext->addContext("1.5) Dilate B on 1.4");
	segmentNecroticRegion_dilate2_subContext->addParameter("Kernel radius", CSParameter::Int, &m_segmentNecroticRegion_Dilate2_RadStrucElem, "pixel");
	CSParameterContext *segmentNecroticRegion_remSmallObj_subContext = m_segmentNecroticRegionContext->addContext("1.6) Remove Small Objects on 1.5");
	segmentNecroticRegion_remSmallObj_subContext->addParameter("Remove objects with volume smaller than", CSParameter::Int, &m_segmentNecroticRegion_RemvObjSizeThreshold, "voxel");

	CSParameterContext *segmentNecroticRegion_saveSelector_subContext = m_segmentNecroticRegionContext->addContext("Save Selector");
	segmentNecroticRegion_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &m_segmentNecroticRegion_Log_Filename, "");
	segmentNecroticRegion_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &m_segmentNecroticRegion_Filename_Prefix, "");
	segmentNecroticRegion_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, m_segmentNecroticRegion_Save_Mode, "");
	//****************************************************************************************************************************


#ifndef CS_TI_QUANT_ONLY

	//Skeletonization Pipeline****************************************************************************************************
	m_skeletonization_Filename = "/path/to/dataSet_file.tif";

	m_skeletonizationContext = mpParameterRoot->addContext( PipelineName[PipelineSkeletonization] );
	m_skeletonizationContext->setDebug(true);

	m_skeletonizationContext->addParameter("Image filename", CSParameter::FileName, &m_skeletonization_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	//****************************************************************************************************************************

	//3D Visualization parameters*************************************************************************************************
	m_visualize3DBileFilename = "/path/to/bile_bin_file.tif";
	m_visualize3DSinusoidFilename = "/path/to/sinusoid_bin_file.tif";
	m_visualize3DNucleiFilename = "/path/to/nuclei_bin_file.tif";
	m_visualize3DGraphFilename = "/path/to/graph0_file.txt";
	m_visualize3DCellShapeFilename = "/path/to/cellShape_bin_file.tif";

	m_visualize3DContext = mpParameterRoot->addContext( PipelineName[PipelineVisualize3D] );
	m_visualize3DContext->setDebug(true);

	m_visualize3DContext->addParameter("Input Bile Image:", CSParameter::FileName, &m_visualize3DBileFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	m_visualize3DContext->addParameter("Input Sinusoid Image:", CSParameter::FileName, &m_visualize3DSinusoidFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	m_visualize3DContext->addParameter("Input Nuclei Image:", CSParameter::FileName, &m_visualize3DNucleiFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	m_visualize3DContext->addParameter("Input Graph Input:", CSParameter::FileName, &m_visualize3DGraphFilename, "")->setAttribute("FileTypes","Graph file (*.txt)");
	m_visualize3DContext->addParameter("Input Cell Shape Image:", CSParameter::FileName, &m_visualize3DCellShapeFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	//****************************************************************************************************************************

#endif // CS_TI_QUANT_ONLY


	//Pipeline: Vein Segmentation*************************************************************************************************
	//TODO: add table view element to show seed point number / and or coords
	//TODO: make mask dilation radius either adjustable or dependent on vox size
	seedPointSelectionPushButton->hide();

	connect( seedPointSelectionPushButton, SIGNAL(clicked(bool)), this, SLOT(seedPointSelectionButtonClicked()) );

	m_segmentVeins_DAPI_Filename = "/path/to/DAPI_channel_file.tif  --  optional";
	m_segmentVeins_DPPIV_Filename = "/path/to/DPPIV_channel_file.tif";
	m_segmentVeins_DMs_Filename = "/path/to/DMs_channel_file.tif";
	m_segmentVeins_DAPI_Threshold = 20;
	m_segmentVeins_DPPIV_Threshold = 20;
	m_segmentVeins_DMs_Threshold = 20;
	m_segmentVeins_OpeningRadius = 6;
	m_segmentVeins_ClosingRadius = 6;
	m_segmentVeins_With_Resampling = false;
	m_segmentVeins_With_Segmentation = true;
	m_segmentVeins_With_LobulePrediction = false;
	m_numberCVSeedPoints = 0;
	m_numberPVSeedPoints = 0;

	m_segmentVeinsContext = mpParameterRoot->addContext( PipelineName[PipelineSegmentVeins] );
	m_segmentVeinsContext->setDebug(true);

	CSParameterContext *segmentVeins_imageReader_subContext = m_segmentVeinsContext->addContext("Image Reader");
	segmentVeins_imageReader_subContext->addParameter("DAPI channel", CSParameter::FileName, &m_segmentVeins_DAPI_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentVeins_imageReader_subContext->addParameter("DPPIV channel", CSParameter::FileName, &m_segmentVeins_DPPIV_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	segmentVeins_imageReader_subContext->addParameter("DMs channel", CSParameter::FileName, &m_segmentVeins_DMs_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	CSParameterContext *segmentVeins_threshold_subContext = m_segmentVeinsContext->addContext("Thresholds");
	segmentVeins_threshold_subContext->addParameter("DAPI threshold", CSParameter::Int, &m_segmentVeins_DAPI_Threshold, "");
	segmentVeins_threshold_subContext->addParameter("DPPIV threshold", CSParameter::Int, &m_segmentVeins_DPPIV_Threshold, "");
	segmentVeins_threshold_subContext->addParameter("DMs threshold", CSParameter::Int, &m_segmentVeins_DMs_Threshold, "");
	CSParameterContext *segmentVeins_postProcessing_subContext = m_segmentVeinsContext->addContext("Postprocessing");
	segmentVeins_postProcessing_subContext->addParameter("Opening kernel radius", CSParameter::Int, &m_segmentVeins_OpeningRadius, "pixel");
	segmentVeins_postProcessing_subContext->addParameter("Closing kernel radius", CSParameter::Int, &m_segmentVeins_ClosingRadius, "pixel");
	CSParameterContext *segmentVeins_cvSeeds_subContext = m_segmentVeinsContext->addContext("Central Vein Seed Points");
	segmentVeins_cvSeeds_subContext->addParameter("Number seed points", CSParameter::Int, &m_numberCVSeedPoints, "");
	CSParameterContext *segmentVeins_pvSeeds_subContext = m_segmentVeinsContext->addContext("Portal Vein Seed Points");
	segmentVeins_pvSeeds_subContext->addParameter("Number seed points", CSParameter::Int, &m_numberPVSeedPoints, "");
	CSParameterContext *segmentVeins_options_subContext = m_segmentVeinsContext->addContext("Options");
	segmentVeins_options_subContext->addParameter("With resampling", CSParameter::Bool, &m_segmentVeins_With_Resampling, "");
	segmentVeins_options_subContext->addParameter("With segmentation", CSParameter::Bool, &m_segmentVeins_With_Segmentation, "");
	segmentVeins_options_subContext->addParameter("With lobule prediction", CSParameter::Bool, &m_segmentVeins_With_LobulePrediction, "");
	//****************************************************************************************************************************


	//Pipeline: Crop image********************************************************************************************************
	//TODO: connect file name specified event to a method that loads the image and sets the end values to the size of the image
	m_cropImageContext_DataDimMode = new CSParameterChoice(imageDataDim, 1);
	m_cropImageContext_DataDimMode->setControlSubContexts(false);
	m_cropImage_Filename = "/path/to/dataSet_file.tif";
	m_cropImage[0][0] = 0;
	m_cropImage[0][1] = 1024;
	m_cropImage[1][0] = 0;
	m_cropImage[1][1] = 1024;
	m_cropImage[2][0] = 0;
	m_cropImage[2][1] = 100;

	m_cropImageContext = mpParameterRoot->addContext( PipelineName[PipelineCrop] );
	m_cropImageContext->setDebug(true);

	m_cropImageContext->addParameter("Data dimensionality", CSParameter::Choice, m_cropImageContext_DataDimMode, "");
	m_cropImageContext->addParameter("Image filename", CSParameter::FileName, &m_cropImage_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	CSParameterContext *cropImage_subContext = m_cropImageContext->addContext("RegionOfInterest");
	cropImage_subContext->addParameter("x start", CSParameter::Int, &m_cropImage[0][0], "pixel");
	cropImage_subContext->addParameter("x end", CSParameter::Int, &m_cropImage[0][1], "pixel");
	cropImage_subContext->addParameter("y start", CSParameter::Int, &m_cropImage[1][0], "pixel");
	cropImage_subContext->addParameter("y end", CSParameter::Int, &m_cropImage[1][1], "pixel");
	cropImage_subContext->addParameter("z start", CSParameter::Int, &m_cropImage[2][0], "pixel");
	cropImage_subContext->addParameter("z end", CSParameter::Int, &m_cropImage[2][1], "pixel");
	//****************************************************************************************************************************


	//Pipeline: CLAHE*************************************************************************************************************
	m_CLAHE_DataDimMode = new CSParameterChoice(imageDataDim, 1);
	m_CLAHE_DataDimMode->setControlSubContexts(false);
	m_CLAHE_Filename = "/path/to/dataSet_file.tif";
	m_CLAHE_HistWindRad = 20;
	m_CLAHE_StepSize = 2;
	m_CLAHE_ClipLevel = 0.012;

	m_CLAHEContext = mpParameterRoot->addContext( PipelineName[PipelineCLAHE] );
	m_CLAHEContext->setDebug(true);

	m_CLAHEContext->addParameter("Data dimensionality", CSParameter::Choice, m_CLAHE_DataDimMode, "");
	m_CLAHEContext->addParameter("Image filename", CSParameter::FileName, &m_CLAHE_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	CSParameterContext *CLAHE_subContext = m_CLAHEContext->addContext("CLAHE Parameters");
	CLAHE_subContext->addParameter("Histogram Window Size", CSParameter::Int, &m_CLAHE_HistWindRad, "pixel");
	CLAHE_subContext->addParameter("Step Size", CSParameter::Int, &m_CLAHE_StepSize, "pixel");
	CLAHE_subContext->addParameter("Clip Level", CSParameter::Double, &m_CLAHE_ClipLevel, "");
	//****************************************************************************************************************************

#ifndef CS_TI_QUANT_ONLY

	//Pipeline: Add Images********************************************************************************************************
	m_AddImages_DataDimMode = new CSParameterChoice(imageDataDim, 1);
	m_AddImages_DataDimMode->setControlSubContexts(false);
	m_AddImages_Filename1 = "/path/to/dataSet_file1.tif";
	m_AddImages_Filename2 = "/path/to/dataSet_file2.tif";
	m_AddImages_Intensity1 = 200;
	m_AddImages_Intensity2 = 100;

	m_AddImagesContext = mpParameterRoot->addContext( PipelineName[PipelineAddImages] );
	m_AddImagesContext->setDebug(true);

	m_AddImagesContext->addParameter("Data dimensionality", CSParameter::Choice, m_AddImages_DataDimMode, "");
	m_AddImagesContext->addParameter("Image filename 1", CSParameter::FileName, &m_AddImages_Filename1, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	m_AddImagesContext->addParameter("Image filename 2", CSParameter::FileName, &m_AddImages_Filename2, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	m_AddImagesContext->addParameter("Intensity 1", CSParameter::Int, &m_AddImages_Intensity1, "");
	m_AddImagesContext->addParameter("Intensity 2", CSParameter::Int, &m_AddImages_Intensity2, "");
	//****************************************************************************************************************************

	//Pipeline: Background Eliminiation*******************************************************************************************
	m_backgroundElimination_DataDimMode = new CSParameterChoice(imageDataDim, 1);
	m_backgroundElimination_DataDimMode->setControlSubContexts(false);
	m_backgroundElimination_Filename = "/path/to/dataSet_file.tif";
	m_backgroundElimination_MedianKernelRadius = 2;
	m_backgroundElimination_KernelRadius = 10;

	m_backgroundEliminationContext = mpParameterRoot->addContext( PipelineName[PipelineBackgroundElimination] );
	m_backgroundEliminationContext->setDebug(true);

	m_backgroundEliminationContext->addParameter("Data dimensionality", CSParameter::Choice, m_backgroundElimination_DataDimMode, "");
	m_backgroundEliminationContext->addParameter("Image filename", CSParameter::FileName, &m_backgroundElimination_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
	CSParameterContext *backgroundElimination_subContext = m_backgroundEliminationContext->addContext("Background Elimination Parameters");
	backgroundElimination_subContext->addParameter("Median kernel radius", CSParameter::Int, &m_backgroundElimination_MedianKernelRadius, "");
	backgroundElimination_subContext->addParameter("Background Elimination kernel radius", CSParameter::Int, &m_backgroundElimination_KernelRadius, "");
	//****************************************************************************************************************************

    //Pipeline: Background Eliminiation*******************************************************************************************
	m_convexFilter_DataDimMode = new CSParameterChoice(imageDataDim, 1);
	m_convexFilter_DataDimMode->setControlSubContexts(false);
	m_convexFilter_Filename = "/path/to/dataSet_file.tif";
	m_convexFilter_Height = 150;
	m_convexFilter_FullyConnected = true;

	m_convexFilterContext = mpParameterRoot->addContext( PipelineName[PipelineHConvexImageFilter] );
	m_convexFilterContext->setDebug(true);

	m_convexFilterContext->addParameter("Data dimensionality", CSParameter::Choice, m_convexFilter_DataDimMode, "");
	m_convexFilterContext->addParameter("Image filename", CSParameter::FileName, &m_convexFilter_Filename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    CSParameterContext *convexFilter_subContext = m_convexFilterContext->addContext("Convex Filter Parameters");
    convexFilter_subContext->addParameter("Height Level for identification of local maxima", CSParameter::Int, &m_convexFilter_Height, "");
    convexFilter_subContext->addParameter("Fully Connected", CSParameter::Bool, &m_convexFilter_FullyConnected, "");
    //****************************************************************************************************************************

    //Pipeline: Superpixel Preparation********************************************************************************************
    std::vector<std::string> startingPoint;
    startingPoint.push_back("With Superpixel-Partitioning");
    startingPoint.push_back("Add Features only");

    std::vector<std::string> superpixelModes;
    superpixelModes.push_back("Specify size");
    superpixelModes.push_back("Specify number");

    mSuperpixelPreparationDataDimMode = new CSParameterChoice(imageDataDim, 1);
    mSuperpixelPreparationDataDimMode->setControlSubContexts(false);
    mSuperpixelPreparationStartingPoint = new CSParameterChoice(startingPoint, 0);
    mSuperpixelPreparationStartingPoint->setControlSubContexts(true);
    mSuperpixelPreparationRawImageFilename = "/path/to/rawImage_file.tif";
    mSuperpixelPreparationObjectImageFilename = "/path/to/objectImage_file.nrrd";
    mSuperpixelPreparationObjectGraphFilename = "/path/to/objectGraph_file.txt";
    mSuperpixelPreparationSuperpixelMode = new CSParameterChoice(superpixelModes, 0);
    mSuperpixelPreparationSuperpixelMode->setControlSubContexts(true);
    mSuperpixelPreparationVoxSpacing[0] = 0.620996;
    mSuperpixelPreparationVoxSpacing[1] = 0.620996;
    mSuperpixelPreparationVoxSpacing[2] = 0.54;
    mSuperpixelPreparationSuperpixelSpacing[0] = 5.0;
    mSuperpixelPreparationSuperpixelSpacing[1] = 5.0;
    mSuperpixelPreparationSuperpixelSpacing[2] = 5.0;
    mSuperpixelPreparationSuperpixelNumber = 5000;
    mSuperpixelPreparationWithFeature = new bool[LabelMapType::NumberUnaryFeatures+LabelMapType::NumberBinaryFeatures];
    for(unsigned int i=0; i<LabelMapType::NumberUnaryFeatures+LabelMapType::NumberBinaryFeatures; i++)
        mSuperpixelPreparationWithFeature[i] = true;

    mSuperpixelPreparationLogFilename = "log";
    mSuperpixelPreparationFilenamePrefix = "superPrep";

    mpSuperpixelPreparationContext = mpParameterRoot->addContext( PipelineName[PipelineSuperpixelPreparation] );
    mpSuperpixelPreparationContext->setDebug(true);

    CSParameterContext *superpixelPreparationReaderContext = mpSuperpixelPreparationContext->addContext("Reader");
    superpixelPreparationReaderContext->addParameter("Data dimensionality", CSParameter::Choice, mSuperpixelPreparationDataDimMode, "");
    superpixelPreparationReaderContext->addParameter("Image filename", CSParameter::FileName, &mSuperpixelPreparationRawImageFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    superpixelPreparationReaderContext->addParameter("Image voxel spacing x", CSParameter::Double, &mSuperpixelPreparationVoxSpacing[0], "micron");
    superpixelPreparationReaderContext->addParameter("Image voxel spacing y", CSParameter::Double, &mSuperpixelPreparationVoxSpacing[1], "micron");
    superpixelPreparationReaderContext->addParameter("Image voxel spacing z", CSParameter::Double, &mSuperpixelPreparationVoxSpacing[2], "micron");

    CSParameterContext *superpixelPreparationEntrypointContext = mpSuperpixelPreparationContext->addContext("Start configuration");
    superpixelPreparationEntrypointContext->addParameter("Start configuration", CSParameter::Choice, mSuperpixelPreparationStartingPoint, "");
    CSParameterContext *superpixelPreparationWithSuperpixelContext = superpixelPreparationEntrypointContext->addContext("With Superpixel Partitioning");
    CSParameterContext *superpixelPreparationSuperpixelSubContext = superpixelPreparationWithSuperpixelContext->addContext("Superpixel Parameters");
    superpixelPreparationSuperpixelSubContext->addParameter("Superpixel mode", CSParameter::Choice, mSuperpixelPreparationSuperpixelMode, "");
    CSParameterContext *superpixelPreparationBySizeSubContext = superpixelPreparationSuperpixelSubContext->addContext("Specify superpixel size");
    superpixelPreparationBySizeSubContext->addParameter("Superpixel spacing x", CSParameter::Double, &mSuperpixelPreparationSuperpixelSpacing[0], "micron");
    superpixelPreparationBySizeSubContext->addParameter("Superpixel spacing y", CSParameter::Double, &mSuperpixelPreparationSuperpixelSpacing[1], "micron");
    superpixelPreparationBySizeSubContext->addParameter("Superpixel spacing z", CSParameter::Double, &mSuperpixelPreparationSuperpixelSpacing[2], "micron");
    CSParameterContext *superpixelPreparationByNumberSubContext = superpixelPreparationSuperpixelSubContext->addContext("Specify superpixel number");
    superpixelPreparationByNumberSubContext->addParameter("Superpixel number", CSParameter::Int, &mSuperpixelPreparationSuperpixelNumber, "");
    superpixelPreparationByNumberSubContext->setVisible(false);
    CSParameterContext *superpixelPreparationWithoutSuperpixelContext = superpixelPreparationEntrypointContext->addContext("Without Superpixel Partitioning");
    CSParameterContext *superpixelPreparationReaderSubContext = superpixelPreparationWithoutSuperpixelContext->addContext("Superpixel Partitioning Reader");
    superpixelPreparationReaderSubContext->addParameter("Object dataset filename", CSParameter::FileName, &mSuperpixelPreparationObjectImageFilename, "")->setAttribute("FileTypes","NRRD Images (*.nrrd)");
    superpixelPreparationReaderSubContext->addParameter("Object graph filename", CSParameter::FileName, &mSuperpixelPreparationObjectGraphFilename, "")->setAttribute("FileTypes","Graph file (*.txt)");
    superpixelPreparationWithoutSuperpixelContext->setVisible(false);

    CSParameterContext *superpixelPreparationFeatureSelectionSubContext = superpixelPreparationWithoutSuperpixelContext->addContext("Feature Selection");
    CSParameterContext *superpixelPreparationUnaryFeatureSelectionSubContext = superpixelPreparationFeatureSelectionSubContext->addContext("Unary Features");
    unsigned int superObjI=0;
    for(superObjI=0; superObjI<LabelMapType::NumberUnaryFeatures; superObjI++)
        superpixelPreparationUnaryFeatureSelectionSubContext->addParameter("With " + LabelMapType::UnaryFeatureName[superObjI], CSParameter::Bool, &mSuperpixelPreparationWithFeature[superObjI], "");
    CSParameterContext *superpixelPreparationBinaryFeatureSelectionSubContext = superpixelPreparationFeatureSelectionSubContext->addContext("Binary Features");
    for(unsigned int j=0; j<LabelMapType::NumberBinaryFeatures; j++)
        superpixelPreparationBinaryFeatureSelectionSubContext->addParameter("With " + LabelMapType::BinaryFeatureName[j], CSParameter::Bool, &mSuperpixelPreparationWithFeature[superObjI+j], "");

    CSParameterContext *superpixelPreparationSaveSelectorSubContext = mpSuperpixelPreparationContext->addContext("Save Selector");
    superpixelPreparationSaveSelectorSubContext->addParameter("Log file name", CSParameter::String, &mSuperpixelPreparationLogFilename, "");
    superpixelPreparationSaveSelectorSubContext->addParameter("Save prefix", CSParameter::String, &mSuperpixelPreparationFilenamePrefix, "");
    //****************************************************************************************************************************

    //Pipeline: Object-based Segmentation*****************************************************************************************
    std::vector<std::string> segmentationTargetModes;
    for(int i=0; i<ObjectBasedSegmentation<3>::NumberClasses; i++)
        segmentationTargetModes.push_back(ObjectBasedSegmentation<3>::ClassName[i]);

    std::vector<std::string> segmentationModes;
    for(int i=0; i<ObjectBasedSegmentation<3>::NumberSegmentationModes; i++)
        segmentationModes.push_back(ObjectBasedSegmentation<3>::SegmentationModeName[i]);

    std::vector<std::string> classifierChoice;
    for(int i=0; i<ObjectBasedSegmentation<3>::NumberClassifiers; i++)
        classifierChoice.push_back(ObjectBasedSegmentation<3>::ClassifierName[i]);

    mObjectBasedSegmentationDataDimMode = new CSParameterChoice(imageDataDim, 1);
    mObjectBasedSegmentationDataDimMode->setControlSubContexts(false);
    mObjectBasedSegmentationRawImageFilename = "/path/to/rawImage_file.tif";
    mObjectBasedSegmentationObjectImageFilename = "/path/to/objectImage_file.nrrd";
    mObjectBasedSegmentationObjectGraphFilename = "/path/to/objectGraph_file.txt";
    mObjectBasedSegmentationTrainingDatabaseFilename = "/path/to/trainingDatabase_file.csv";
    mObjectBasedSegmentationVoxSpacing[0] = 0.620996;
    mObjectBasedSegmentationVoxSpacing[1] = 0.620996;
    mObjectBasedSegmentationVoxSpacing[2] = 0.54;
    mObjectBasedSegmentationSegmentationMode = new CSParameterChoice(segmentationModes, 0);
    mObjectBasedSegmentationSegmentationMode->setControlSubContexts(false);
    mObjectBasedSegmentationClassifierChoice = new CSParameterChoice(classifierChoice, 0);
    mObjectBasedSegmentationClassifierChoice->setControlSubContexts(true);
    mObjectBasedSegmentationSegmentationTargetMode = new CSParameterChoice(segmentationTargetModes, 0);
    mObjectBasedSegmentationSegmentationTargetMode->setControlSubContexts(false);
    mObjectBasedSegmentationKeepOnlyTargetClassObjects = false;
    mObjectBasedSegmentationWithFeature = new bool[LabelMapType::NumberUnaryFeatures+LabelMapType::NumberBinaryFeatures];
    for(unsigned int i=0; i<LabelMapType::NumberUnaryFeatures+LabelMapType::NumberBinaryFeatures; i++)
        mObjectBasedSegmentationWithFeature[i] = true;

    mObjectBasedSegmentationLogFilename = "log";
    mObjectBasedSegmentationFilenamePrefix = "objBasedSeg";
    mpObjectBasedSegmentationSaveMode = new CSParameterChoice(saveModes, 1);

    mpObjectBasedSegmentationContext = mpParameterRoot->addContext( PipelineName[PipelineObjectBasedSegmentation] );
    mpObjectBasedSegmentationContext->setDebug(true);

    CSParameterContext *objectBasedSegReaderSubContext = mpObjectBasedSegmentationContext->addContext("Superpixel Partitioning Reader");
    objectBasedSegReaderSubContext->addParameter("Data dimensionality", CSParameter::Choice, mObjectBasedSegmentationDataDimMode, "");
    objectBasedSegReaderSubContext->addParameter("Image voxel spacing x", CSParameter::Double, &mObjectBasedSegmentationVoxSpacing[0], "micron");
    objectBasedSegReaderSubContext->addParameter("Image voxel spacing y", CSParameter::Double, &mObjectBasedSegmentationVoxSpacing[1], "micron");
    objectBasedSegReaderSubContext->addParameter("Image voxel spacing z", CSParameter::Double, &mObjectBasedSegmentationVoxSpacing[2], "micron");
    objectBasedSegReaderSubContext->addParameter("Image filename", CSParameter::FileName, &mObjectBasedSegmentationRawImageFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    objectBasedSegReaderSubContext->addParameter("Object dataset filename", CSParameter::FileName, &mObjectBasedSegmentationObjectImageFilename, "")->setAttribute("FileTypes","NRRD Images (*.nrrd)");
    objectBasedSegReaderSubContext->addParameter("Object graph filename", CSParameter::FileName, &mObjectBasedSegmentationObjectGraphFilename, "")->setAttribute("FileTypes","Graph file (*.txt)");

    CSParameterContext *objectBasedSegSegmentationSubContext = mpObjectBasedSegmentationContext->addContext("Segmentation");
    objectBasedSegSegmentationSubContext->addParameter("Training database filename", CSParameter::FileName, &mObjectBasedSegmentationTrainingDatabaseFilename, "")->setAttribute("FileTypes","CSV Files (*.csv)");
    objectBasedSegSegmentationSubContext->addParameter("Segmentation mode", CSParameter::Choice, mObjectBasedSegmentationSegmentationMode, "");
    objectBasedSegSegmentationSubContext->addParameter("Classifier", CSParameter::Choice, mObjectBasedSegmentationClassifierChoice, "");
    CSParameterContext *objectBasedSegWithoutSegContext = objectBasedSegSegmentationSubContext->addContext("Without Segmentation");
    CSParameterContext *objectBasedSegKNNBasedMergingContext = objectBasedSegSegmentationSubContext->addContext("KNN-Merging Segmentation");
    objectBasedSegKNNBasedMergingContext->setVisible(false);
    CSParameterContext *objectBasedSegPythonSVMContext = objectBasedSegSegmentationSubContext->addContext("Python-based SVM-Merging Segmentation");
    objectBasedSegPythonSVMContext->setVisible(false);
    CSParameterContext *objectBasedSegFilterLargeObjectsContext = objectBasedSegSegmentationSubContext->addContext("Filter large objects Segmentation");
    objectBasedSegFilterLargeObjectsContext->setVisible(false);
    objectBasedSegSegmentationSubContext->addParameter("Segmentation target class", CSParameter::Choice, mObjectBasedSegmentationSegmentationTargetMode, "");
    objectBasedSegSegmentationSubContext->addParameter("Keep only target class label objects", CSParameter::Bool, &mObjectBasedSegmentationKeepOnlyTargetClassObjects, "");

    CSParameterContext *objectBasedSegFeatureSelectionSubContext = mpObjectBasedSegmentationContext->addContext("Feature Selection");
    CSParameterContext *objectBasedSegUnaryFeatureSelectionSubContext = objectBasedSegFeatureSelectionSubContext->addContext("Unary Features");
    unsigned int objectObjI=0;
    for(objectObjI=0; objectObjI<LabelMapType::NumberUnaryFeatures; objectObjI++)
        objectBasedSegUnaryFeatureSelectionSubContext->addParameter("With " + LabelMapType::UnaryFeatureName[objectObjI], CSParameter::Bool, &mObjectBasedSegmentationWithFeature[objectObjI], "");
    CSParameterContext *objectBasedSegBinaryFeatureSelectionSubContext = objectBasedSegFeatureSelectionSubContext->addContext("Binary Features");
    for(unsigned int j=0; j<LabelMapType::NumberBinaryFeatures; j++)
        objectBasedSegBinaryFeatureSelectionSubContext->addParameter("With " + LabelMapType::BinaryFeatureName[j], CSParameter::Bool, &mObjectBasedSegmentationWithFeature[objectObjI+j], "");

    CSParameterContext *objectBasedSegSaveSelectorSubContext = mpObjectBasedSegmentationContext->addContext("Save Selector");
    objectBasedSegSaveSelectorSubContext->addParameter("Log file name", CSParameter::String, &mObjectBasedSegmentationLogFilename, "");
    objectBasedSegSaveSelectorSubContext->addParameter("Save prefix", CSParameter::String, &mObjectBasedSegmentationFilenamePrefix, "");
    objectBasedSegSaveSelectorSubContext->addParameter("Save mode", CSParameter::Choice, mpObjectBasedSegmentationSaveMode, "");
    //****************************************************************************************************************************

    //Pipeline: Train Classifier Pipeline*****************************************************************************************
    std::vector<std::string> trainingActionSelection;
    trainingActionSelection.push_back("Add training dataset to training database");
    trainingActionSelection.push_back("Train Classifier");
    std::vector<std::string> classifierSelection;
    classifierSelection.push_back("KNN");

    std::vector<std::string> classSelection;

    for(int i=0; i<ObjectBasedSegmentation<3>::NumberClasses; i++)
        classSelection.push_back(ObjectBasedSegmentation<3>::ClassName[i]);

    mTrainDataDimMode = new CSParameterChoice(imageDataDim, 1);
    mTrainDataDimMode->setControlSubContexts(false);
    mTrainActionChoice = new CSParameterChoice(trainingActionSelection, 0);
    mTrainActionChoice->setControlSubContexts(true);
    mTrainClassifierChoice = new CSParameterChoice(classifierSelection, 0);
    mTrainClassifierChoice->setControlSubContexts(false);
    mTrainClassLabelChoice = new CSParameterChoice(classSelection, 0);
    mTrainClassLabelChoice->setControlSubContexts(false);

    mTrainClassifierOrigImageFilename = "/path/to/rawImage_file.tif";
    mTrainClassifierObjectImageFilename = "/path/to/objectImage_file.nrrd";
    mTrainClassifierObjectGraphFilename = "/path/to/objectGraph_file.txt";
    mTrainClassifierTrainingDatasetFilename = "/path/to/training_dataset.tif";
    mTrainClassifierTrainingDatabaseFilename = "/path/to/training_database.csv";
    mTrainClassifierWithIncremental = false;
    mTrainClassifierWithFeature = new bool[LabelMapType::NumberUnaryFeatures+LabelMapType::NumberBinaryFeatures];
    for(unsigned int i=0; i<LabelMapType::NumberUnaryFeatures+LabelMapType::NumberBinaryFeatures; i++)
        mTrainClassifierWithFeature[i] = true;

    mpTrainClassifierContext = mpParameterRoot->addContext( PipelineName[PipelineTrainClassifier] );
    mpTrainClassifierContext->setDebug(true);

    CSParameterContext *trainAtionSelectionSubContext = mpTrainClassifierContext->addContext("Training Action Selection");
    trainAtionSelectionSubContext->addParameter("Training action", CSParameter::Choice, mTrainActionChoice, "");
    CSParameterContext *trainToDatabaseSelectionSubContext = trainAtionSelectionSubContext->addContext("Add training dataset to training database");
    CSParameterContext *trainClassifierSelectionSubContext = trainAtionSelectionSubContext->addContext("Classifier Selection");
    trainClassifierSelectionSubContext->addParameter("Classifier", CSParameter::Choice, mTrainClassifierChoice, "");
    trainClassifierSelectionSubContext->setVisible(false);

    CSParameterContext *trainClassifierDataSubContext = mpTrainClassifierContext->addContext("Data");
    trainClassifierDataSubContext->addParameter("Data dimensionality", CSParameter::Choice, mTrainDataDimMode, "");
    CSParameterContext *trainClassifierObjectBasedDataSubContext = trainClassifierDataSubContext->addContext("Object-based data");
    trainClassifierObjectBasedDataSubContext->addParameter("Image filename", CSParameter::FileName, &mTrainClassifierOrigImageFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    trainClassifierObjectBasedDataSubContext->addParameter("Object dataset filename", CSParameter::FileName, &mTrainClassifierObjectImageFilename, "")->setAttribute("FileTypes","NRRD Images (*.nrrd)");
    trainClassifierObjectBasedDataSubContext->addParameter("Object graph filename", CSParameter::FileName, &mTrainClassifierObjectGraphFilename, "")->setAttribute("FileTypes","Graph file (*.txt)");
    CSParameterContext *trainClassifierTrainingDataSubContext = trainClassifierDataSubContext->addContext("Training data");
    trainClassifierTrainingDataSubContext->addParameter("Training dataset filename", CSParameter::FileName, &mTrainClassifierTrainingDatasetFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    trainClassifierTrainingDataSubContext->addParameter("Training database filename", CSParameter::FileName, &mTrainClassifierTrainingDatabaseFilename, "")->setAttribute("FileTypes","CSV Files (*.csv)");
    trainClassifierTrainingDataSubContext->addParameter("Class selection", CSParameter::Choice, mTrainClassLabelChoice, "");

    CSParameterContext *trainClassifierOptionsSubContext = mpTrainClassifierContext->addContext("Options");
    trainClassifierOptionsSubContext->addParameter("Incremental training", CSParameter::Bool, &mTrainClassifierWithIncremental, "");

    CSParameterContext *trainClassifierFeatureSelectionSubContext = mpTrainClassifierContext->addContext("Feature Selection");
    unsigned int traI=0;
    for(traI=0; traI<LabelMapType::NumberUnaryFeatures; traI++)
        trainClassifierFeatureSelectionSubContext->addParameter("With " + LabelMapType::UnaryFeatureName[traI], CSParameter::Bool, &mTrainClassifierWithFeature[traI], "");
    for(unsigned int j=0; j<LabelMapType::NumberBinaryFeatures; j++)
        trainClassifierFeatureSelectionSubContext->addParameter("With " + LabelMapType::BinaryFeatureName[j], CSParameter::Bool, &mTrainClassifierWithFeature[traI+j], "");
    //****************************************************************************************************************************

    //Pipeline: Segment cell membrane from bcat staining**************************************************************************
    mSegmentCellMembraneBCatFilename = "/path/to/BCat_channel_file.tif";
    mSegmentCellMembraneSinFilename = "/path/to/sinusoid_bin_file.tif";
    mSegmentCellMembraneNucFilename = "/path/to/hepNuclei_bin_file.tif";
    mSegmentCellMembraneCVFilename = "/path/to/CV_bin_file.tif";
    mSegmentCellMembranePVFilename = "/path/to/PV_bin_file.tif";

    mSegmentCellMembraneVoxelSpacing[0] = 0.207;
    mSegmentCellMembraneVoxelSpacing[1] = 0.207;
    mSegmentCellMembraneVoxelSpacing[2] = 0.54;

    mSegmentCellMembraneBCatThresholdMode = new CSParameterChoice(thresholdModes, 0);
    mSegmentCellMembraneBCatThresholdMode->setControlSubContexts(true);
    mSegmentCellMembraneBCatAdapOtsuThresholdRegionSize[0] = 50;
    mSegmentCellMembraneBCatAdapOtsuThresholdRegionSize[1] = 50;
    mSegmentCellMembraneBCatAdapOtsuThresholdRegionSize[2] = 5;
    mSegmentCellMembraneBCatAdapOtsuThresholdNumSamples = 50000;
    mSegmentCellMembraneBCatManualMinThreshold = 35;
    mSegmentCellMembraneBCatManualMaxThreshold = 255;
    mSegmentCellMembraneHoleFilling_RegionSize[0] = 1;
    mSegmentCellMembraneHoleFilling_RegionSize[1] = 1;
    mSegmentCellMembraneHoleFilling_RegionSize[2] = 1;
    mSegmentCellMembraneHoleFilling_MajThreshold = 2;
    mSegmentCellMembraneInvHoleFilling_RegionSize[0] = 1;
    mSegmentCellMembraneInvHoleFilling_RegionSize[1] = 1;
    mSegmentCellMembraneInvHoleFilling_RegionSize[2] = 1;
    mSegmentCellMembraneInvHoleFilling_MajThreshold = 6;
    mSegmentCellMembraneClosingRadStrucElem[0] = 0;
    mSegmentCellMembraneClosingRadStrucElem[1] = 0;
    mSegmentCellMembraneClosingRadStrucElem[2] = 0;
    mSegmentCellMembraneOpeningRadStrucElem[0] = 8;
    mSegmentCellMembraneOpeningRadStrucElem[1] = 8;
    mSegmentCellMembraneOpeningRadStrucElem[2] = 4;
    mSegmentCellMembraneMinCellDiameterThreshold = 10.0;
    mSegmentCellMembraneWatershedFloodLevel = 1.5;

    mSegmentCellMembraneLogFilename = "log";
    mSegmentCellMembraneFilenamePrefix = "cellShapeSeg";
    mSegmentCellMembraneSaveMode = new CSParameterChoice(saveModes, 0);

    mpSegmentCellMembraneContext = mpParameterRoot->addContext( PipelineName[PipelineSegmentCellMembrane] );
    mpSegmentCellMembraneContext->setDebug(true);

    CSParameterContext *segmentCellMembrane_imageReader_subContext = mpSegmentCellMembraneContext->addContext("Image Reader");
    segmentCellMembrane_imageReader_subContext->addParameter("BCat channel", CSParameter::FileName, &mSegmentCellMembraneBCatFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    segmentCellMembrane_imageReader_subContext->addParameter("Sinusoid segmentation", CSParameter::FileName, &mSegmentCellMembraneSinFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    segmentCellMembrane_imageReader_subContext->addParameter("Hepatic nuclei segmentation", CSParameter::FileName, &mSegmentCellMembraneNucFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    segmentCellMembrane_imageReader_subContext->addParameter("Central vein segmentation", CSParameter::FileName, &mSegmentCellMembraneCVFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    segmentCellMembrane_imageReader_subContext->addParameter("Portal vein segmentation", CSParameter::FileName, &mSegmentCellMembranePVFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    segmentCellMembrane_imageReader_subContext->addParameter("Voxel spacing x", CSParameter::Double, &mSegmentCellMembraneVoxelSpacing[0], "micron");
    segmentCellMembrane_imageReader_subContext->addParameter("Voxel spacing y", CSParameter::Double, &mSegmentCellMembraneVoxelSpacing[1], "micron");
    segmentCellMembrane_imageReader_subContext->addParameter("Voxel spacing z", CSParameter::Double, &mSegmentCellMembraneVoxelSpacing[2], "micron");
    CSParameterContext *segmentCellMembrane_binThreshold_subContext = mpSegmentCellMembraneContext->addContext("1.1) Binary Threshold");
    segmentCellMembrane_binThreshold_subContext->addParameter("Threshold mode", CSParameter::Choice, mSegmentCellMembraneBCatThresholdMode, "");
    CSParameterContext *segmentCellMembrane_adapOtsuThreshold_subContext = segmentCellMembrane_binThreshold_subContext->addContext("Adaptive Otsu-Threshold");
    segmentCellMembrane_adapOtsuThreshold_subContext->addParameter("Sample region size x", CSParameter::Int, &mSegmentCellMembraneBCatAdapOtsuThresholdRegionSize[0], "pixel");
    segmentCellMembrane_adapOtsuThreshold_subContext->addParameter("Sample region size y", CSParameter::Int, &mSegmentCellMembraneBCatAdapOtsuThresholdRegionSize[1], "pixel");
    segmentCellMembrane_adapOtsuThreshold_subContext->addParameter("Sample region size z", CSParameter::Int, &mSegmentCellMembraneBCatAdapOtsuThresholdRegionSize[2], "pixel");
    segmentCellMembrane_adapOtsuThreshold_subContext->addParameter("Number samples", CSParameter::Int, &mSegmentCellMembraneBCatAdapOtsuThresholdNumSamples, "");
    CSParameterContext *segmentCellMembrane_otsuThreshold_subContext = segmentCellMembrane_binThreshold_subContext->addContext("Otsu-Threshold");
    segmentCellMembrane_otsuThreshold_subContext->setVisible(false);
    CSParameterContext *segmentCellMembrane_manualThreshold_subContext = segmentCellMembrane_binThreshold_subContext->addContext("Manual Threshold");
    segmentCellMembrane_manualThreshold_subContext->addParameter("BCat manual min threshold", CSParameter::Int, &mSegmentCellMembraneBCatManualMinThreshold, "");
    segmentCellMembrane_manualThreshold_subContext->addParameter("BCat manual max threshold", CSParameter::Int, &mSegmentCellMembraneBCatManualMaxThreshold, "");
    segmentCellMembrane_manualThreshold_subContext->setVisible(false);
    CSParameterContext *segmentCellMembrane_holeFilling_subContext = mpSegmentCellMembraneContext->addContext("1.2) Hole Filling on 1.1");
    segmentCellMembrane_holeFilling_subContext->addParameter("Radius x", CSParameter::Int, &mSegmentCellMembraneHoleFilling_RegionSize[0], "pixel");
    segmentCellMembrane_holeFilling_subContext->addParameter("Radius y", CSParameter::Int, &mSegmentCellMembraneHoleFilling_RegionSize[1], "pixel");
    segmentCellMembrane_holeFilling_subContext->addParameter("Radius z", CSParameter::Int, &mSegmentCellMembraneHoleFilling_RegionSize[2], "pixel");
    segmentCellMembrane_holeFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &mSegmentCellMembraneHoleFilling_MajThreshold, "");
    CSParameterContext *segmentCellMembrane_invHoleFilling_subContext = mpSegmentCellMembraneContext->addContext("1.3) Inverse Hole Filling on 1.2");
    segmentCellMembrane_invHoleFilling_subContext->addParameter("Radius x", CSParameter::Int, &mSegmentCellMembraneInvHoleFilling_RegionSize[0], "pixel");
    segmentCellMembrane_invHoleFilling_subContext->addParameter("Radius y", CSParameter::Int, &mSegmentCellMembraneInvHoleFilling_RegionSize[1], "pixel");
    segmentCellMembrane_invHoleFilling_subContext->addParameter("Radius z", CSParameter::Int, &mSegmentCellMembraneInvHoleFilling_RegionSize[2], "pixel");
    segmentCellMembrane_invHoleFilling_subContext->addParameter("Majority threshold", CSParameter::Int, &mSegmentCellMembraneInvHoleFilling_MajThreshold, "");
    CSParameterContext *segmentCellMembrane_closing_subContext = mpSegmentCellMembraneContext->addContext("1.4) Closing on 1.3");
    segmentCellMembrane_closing_subContext->addParameter("Kernel radius x", CSParameter::Int, &mSegmentCellMembraneClosingRadStrucElem[0], "pixel");
    segmentCellMembrane_closing_subContext->addParameter("Kernel radius y", CSParameter::Int, &mSegmentCellMembraneClosingRadStrucElem[1], "pixel");
    segmentCellMembrane_closing_subContext->addParameter("Kernel radius z", CSParameter::Int, &mSegmentCellMembraneClosingRadStrucElem[2], "pixel");
    CSParameterContext *segmentCellMembrane_opening_subContext = mpSegmentCellMembraneContext->addContext("1.5) Opening on 1.4");
    segmentCellMembrane_opening_subContext->addParameter("Kernel radius x", CSParameter::Int, &mSegmentCellMembraneOpeningRadStrucElem[0], "pixel");
    segmentCellMembrane_opening_subContext->addParameter("Kernel radius y", CSParameter::Int, &mSegmentCellMembraneOpeningRadStrucElem[1], "pixel");
    segmentCellMembrane_opening_subContext->addParameter("Kernel radius z", CSParameter::Int, &mSegmentCellMembraneOpeningRadStrucElem[2], "pixel");
    CSParameterContext *segmentCellMembrane_watershed_subContext = mpSegmentCellMembraneContext->addContext("1.6) Cell Segmentation on 1.5");
    segmentCellMembrane_watershed_subContext->addParameter("Watershed flood level", CSParameter::Double, &mSegmentCellMembraneWatershedFloodLevel, "");
    segmentCellMembrane_watershed_subContext->addParameter("Smallest cell diameter", CSParameter::Double, &mSegmentCellMembraneMinCellDiameterThreshold, "micron");

    CSParameterContext *segmentCellMembrane_saveSelector_subContext = mpSegmentCellMembraneContext->addContext("Save Selector");
    segmentCellMembrane_saveSelector_subContext->addParameter("Log file name", CSParameter::String, &mSegmentCellMembraneLogFilename, "");
    segmentCellMembrane_saveSelector_subContext->addParameter("Save prefix", CSParameter::String, &mSegmentCellMembraneFilenamePrefix, "");
    segmentCellMembrane_saveSelector_subContext->addParameter("Save mode", CSParameter::Choice, mSegmentCellMembraneSaveMode, "");
    //****************************************************************************************************************************

    //Pipeline: Objects to Matrix Filter******************************************************************************************
    mObjectsToMatrixImageFilename = "/path/to/segmentation_bin_file.tif";
    mObjectsToMatrixColoredImageFilename = "/path/to/segmentation_overlay_file.tif";
    mObjectsToMatrixRows = 5;
    mObjectsToMatrixCols = 5;
    mObjectsToMatrixSpacing[0] = 0.207;
    mObjectsToMatrixSpacing[1] = 0.207;
    mObjectsToMatrixSpacing[2] = 0.54;
    mObjectsToMatrixDim[0] = 1000;
    mObjectsToMatrixDim[1] = 1000;
    mObjectsToMatrixDim[2] = 100;
    mObjectsToMatrixFullyConnected = false;
    mObjectsToMatrixWithFilter = false;
    mObjectsToMatrixFilterBySlide = 0;
    mObjectsToMatrixWriteObjectIdentificationFile = false;

    mpObjectsToMatrixContext = mpParameterRoot->addContext( PipelineName[PipelineObjectsToMatrix] );
    mpObjectsToMatrixContext->setDebug(true);

    mpObjectsToMatrixContext->addParameter("Image filename", CSParameter::FileName, &mObjectsToMatrixImageFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    mpObjectsToMatrixContext->addParameter("Image overlay filename", CSParameter::FileName, &mObjectsToMatrixColoredImageFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    mpObjectsToMatrixContext->addParameter("Voxel spacing x", CSParameter::Double, &mObjectsToMatrixSpacing[0], "micron");
    mpObjectsToMatrixContext->addParameter("Voxel spacing y", CSParameter::Double, &mObjectsToMatrixSpacing[1], "micron");
    mpObjectsToMatrixContext->addParameter("Voxel spacing z", CSParameter::Double, &mObjectsToMatrixSpacing[2], "micron");
    mpObjectsToMatrixContext->addParameter("Matrix rows", CSParameter::Int, &mObjectsToMatrixRows, "");
    mpObjectsToMatrixContext->addParameter("Matrix columns", CSParameter::Int, &mObjectsToMatrixCols, "");
    mpObjectsToMatrixContext->addParameter("Matrix dimension x", CSParameter::Int, &mObjectsToMatrixDim[0], "");
    mpObjectsToMatrixContext->addParameter("Matrix dimension y", CSParameter::Int, &mObjectsToMatrixDim[1], "");
    mpObjectsToMatrixContext->addParameter("Matrix dimension z", CSParameter::Int, &mObjectsToMatrixDim[2], "");
    mpObjectsToMatrixContext->addParameter("Objects fully connected?", CSParameter::Bool, &mObjectsToMatrixFullyConnected, "");
    mpObjectsToMatrixContext->addParameter("Apply filter to objects", CSParameter::Bool, &mObjectsToMatrixWithFilter, "");
    mpObjectsToMatrixContext->addParameter("Filter by slide", CSParameter::Int, &mObjectsToMatrixFilterBySlide, "");
    mpObjectsToMatrixContext->addParameter("Write object identification file", CSParameter::Bool, &mObjectsToMatrixWriteObjectIdentificationFile, "");
    //****************************************************************************************************************************

    //Pipeline: Compare Segmentations Filter**************************************************************************************
    mCompareSegmentationsDatasetID = 0;
    mCompareSegmentationsGoldImageFilename = "/path/to/segmentation_bin_file.tif";
    mCompareSegmentationsEvalImageFilename = "/path/to/segmentation_bin_file.tif";

    mpCompareSegmentationsContext = mpParameterRoot->addContext( PipelineName[PipelineCompareSegmentations] );
    mpCompareSegmentationsContext->setDebug(true);

    mpCompareSegmentationsContext->addParameter("Dataset ID", CSParameter::Int, &mCompareSegmentationsDatasetID, "");
    mpCompareSegmentationsContext->addParameter("Goldstandard segmentation filename", CSParameter::FileName, &mCompareSegmentationsGoldImageFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    mpCompareSegmentationsContext->addParameter("Evaluation segmentation filename", CSParameter::FileName, &mCompareSegmentationsEvalImageFilename, "")->setAttribute("FileTypes","TIFF Images (*.tif *.tiff)");
    //****************************************************************************************************************************

    //Pipeline: File Format Converter Filter**************************************************************************************
    mFileFormatConverterFilename = "/path/to/input_file.tif";

    std::vector<std::string> inputDataFormatNames, outputDataFormatNames;
    std::vector<FileFormatConverter::OutputData> outputDataFormats;

    for(unsigned int i=0; i<FileFormatConverter::NumberInputFileFormats; i++)
        inputDataFormatNames.push_back(FileFormatConverter::InputDataName[i]);

    outputDataFormats = FileFormatConverter::AvailableConversions((FileFormatConverter::InputData)0);

    for(unsigned int i=0; i<outputDataFormats.size(); i++)
        outputDataFormatNames.push_back(FileFormatConverter::OutputDataName[outputDataFormats[i]]);

    mFileFormatConverterInputChoice = new CSParameterChoice(inputDataFormatNames, 0);
    mFileFormatConverterInputChoice->setControlSubContexts(false);
    mFileFormatConverterOutputChoice = new CSParameterChoice(outputDataFormatNames, 0);
    mFileFormatConverterOutputChoice->setControlSubContexts(false);

    mpFileFormatConverterContext = mpParameterRoot->addContext( PipelineName[PipelineFileFormatConverter] );
    mpFileFormatConverterContext->setDebug(true);

    mpFileFormatConverterContext->addParameter("Input filename", CSParameter::FileName, &mFileFormatConverterFilename, "");
    mpFileFormatConverterContext->addParameter("Input file format", CSParameter::Choice, mFileFormatConverterInputChoice, "");
    mpFileFormatConverterContext->addParameter("Output file format", CSParameter::Choice, mFileFormatConverterOutputChoice, "");
    //****************************************************************************************************************************

    //Pipeline: Compare Segmentations Filter**************************************************************************************
    mpTestPipelineContext = mpParameterRoot->addContext( PipelineName[PipelineTest] );
    mpTestPipelineContext->setDebug(true);
    //****************************************************************************************************************************

#endif // CS_TI_QUANT_ONLY


    //Init Tree view with first element
    mAnalyzeCellsContext->setupGUI( pipelineEditorTreeView );

	pipelineEditorTreeView->resizeColumnToContents(0);
	pipelineEditorTreeView->resizeColumnToContents(1);
	pipelineEditorTreeView->setAlternatingRowColors(true);
	pipelineEditorTreeView->expandAll();

	stopPipelineButton->setVisible(false);
	pipelineProgressBar->setVisible(false);

	selectPipelineComboBox->clear();

    QStringList menuItems = QStringList();
    for ( unsigned int i=0; PipelineName[i].size() !=0; ++i )
        menuItems << QApplication::translate("ImageProcessing",
                                             PipelineName[i].c_str(),
                                             0,
                                             QApplication::UnicodeUTF8);
	selectPipelineComboBox->insertItems(0, menuItems);

	connect( ((QAbstractItemModel*)(pipelineEditorTreeView->model())), SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)), this, SLOT(pipelineEditorTreeViewChanged(const QModelIndex&, const QModelIndex&)) );

    mpWorker = new ImageProcessingWorker( mpParameterRoot, this );
    mpWorkerThread = new QThread();
    mpWorker->moveToThread( mpWorkerThread );

    connect( mpWorkerThread, SIGNAL(started()), mpWorker, SLOT(execute()) );
    connect( mpWorker, SIGNAL(finished()), mpWorkerThread, SLOT(quit()) );
    connect( mpWorker, SIGNAL(finished()), this, SLOT(pipelineFinished()) );
    connect( mpWorker, SIGNAL(aborted(QString)), this, SLOT(processPipelineError(QString)) );

    connect( mpWorker, SIGNAL(displayResult(void *, int, bool)), this, SLOT(displayResult(void *, int, bool)) );

    jobQueueProgressBar->hide();
    mpJobManagerThread = new QThread();
    jobManager->moveToThread(mpJobManagerThread);
    connect( mpJobManagerThread, SIGNAL(started()), jobManager, SLOT(Start()) );
    connect( jobManager, SIGNAL(finished()), mpJobManagerThread, SLOT(quit()) );
    connect( jobManager, SIGNAL(finished()), jobQueueProgressBar, SLOT(hide()) );
    connect( jobManager, SIGNAL(aborted(QString)), this, SLOT(processPipelineError(QString)) );
}


ImageProcessing::~ImageProcessing()
{
}


void ImageProcessing::startPipelineButtonClicked()
{
    pipelineProgressBar->show();
    mpWorker->setPipeline( (Pipeline)selectPipelineComboBox->currentIndex() );
    mpWorkerThread->start();
    m_isBusy = true;
    startPipelineButton->setEnabled(false);
}


void ImageProcessing::pipelineFinished()
{
    std::cout << "\n\nPipeline finished!\n\n";
    m_isBusy = false;
    startPipelineButton->setEnabled(true);
    pipelineProgressBar->hide();
}


void ImageProcessing::displayResult( void *pipelineInput, int dm, bool displayGraphAlreadyDisplayed )
{
    std::cout << "In displayResult...\n";
    CSParameterContext * pipelineContext = (CSParameterContext *) mpWorker->mpPipelineContextCopy;

    switch (mpWorker->mPipelineID)
    {
    case PipelineExtractAndAnalyzeGraph:
    {
        ExtractGraph * extractPipeline;
        AnalyzeBileNetworkFilter *analyzePipeline;

        switch (dm)
        {
        case DisplayGraphWithoutAnalysis:
            extractPipeline = static_cast<ExtractGraph * >(pipelineInput);
            break;
        case DisplayGraphWithAnalysis:
            analyzePipeline = static_cast<AnalyzeBileNetworkFilter *>(pipelineInput);
            break;
        default:
            std::cout << "ImageProcessing::displayResult():  pipeline display mode unknown.";
            break;
        }

        std::string readMode = ( (CSParameterChoice*)(pipelineContext->findContext("Reader", 0)->findParameter("Input type", 0)->dataPointer()) )->currentString();
        int entryPoint = 0;
        if(readMode.compare("Final Graph")==0) entryPoint = 1;

        bool displayVanillaGraph = ( *(bool*)(pipelineContext->findParameter("Display vanilla graph", 0)->dataPointer()) );
        bool displayResampledGraph = ( *(bool*)(pipelineContext->findParameter("Display resampled graph", 0)->dataPointer()) );
        bool displayPrunedGraph = ( *(bool*)(pipelineContext->findParameter("Display pruned graph", 0)->dataPointer()) );
        bool displayCollapsedGraph = ( *(bool*)(pipelineContext->findParameter("Display collapsed graph", 0)->dataPointer()) );
        bool displayFinalGraph = ( *(bool*)(pipelineContext->findParameter("Display final graph", 0)->dataPointer()) );
        double voxelSpacing[3];
        voxelSpacing[0] = ( *(double*)(pipelineContext->findParameter("Voxel spacing x", 0)->dataPointer()) );
        voxelSpacing[1] = ( *(double*)(pipelineContext->findParameter("Voxel spacing y", 0)->dataPointer()) );
        voxelSpacing[2] = ( *(double*)(pipelineContext->findParameter("Voxel spacing z", 0)->dataPointer()) );

        std::string analysisMode = ( (CSParameterChoice*)(pipelineContext->findParameter("Analysis mode", 0)->dataPointer()) )->currentString();
        bool withAnalysis = (analysisMode.compare("With Analysis")==0);

        if(dm==DisplayGraphWithoutAnalysis && entryPoint==0 && displayVanillaGraph) {
            GraphViewer *vanillaGraphViewer = new GraphViewer();
            vanillaGraphViewer->SetName("Vanilla graph");
            vanillaGraphViewer->SetGraphs(extractPipeline->GetVanillaGraphs());
            vanillaGraphViewer->SetScaling(voxelSpacing[0], voxelSpacing[1], voxelSpacing[2]);
            vanillaGraphViewer->Show(false);

            displays["Vanilla graph"] = vanillaGraphViewer;
        }
        if(dm==DisplayGraphWithoutAnalysis && entryPoint==0 && displayResampledGraph) {
            GraphViewer *resampledGraphViewer = new GraphViewer();
            resampledGraphViewer->SetName("Resampled graph");
            resampledGraphViewer->SetGraphs(extractPipeline->GetResampledGraphs());
            resampledGraphViewer->SetScaling(voxelSpacing[0], voxelSpacing[1], voxelSpacing[2]);
            resampledGraphViewer->Show(false);

            displays["Resampled graph"] = resampledGraphViewer;
        }
        if(dm==DisplayGraphWithoutAnalysis && entryPoint==0 && displayPrunedGraph) {
            GraphViewer *prunedGraphViewer = new GraphViewer();
            prunedGraphViewer->SetName("Pruned graph");
            prunedGraphViewer->SetGraphs(extractPipeline->GetPrunedGraphs());
            prunedGraphViewer->SetScaling(voxelSpacing[0], voxelSpacing[1], voxelSpacing[2]);
            prunedGraphViewer->Show(false);

            displays["Pruned graph"] = prunedGraphViewer;
        }
        if(dm==DisplayGraphWithoutAnalysis && entryPoint==0 && displayCollapsedGraph) {
            GraphViewer *collapsedGraphViewer = new GraphViewer();
            collapsedGraphViewer->SetName("Collapsed graph");
            collapsedGraphViewer->SetGraphs(extractPipeline->GetCollapsedGraphs());
            collapsedGraphViewer->SetScaling(voxelSpacing[0], voxelSpacing[1], voxelSpacing[2]);
            collapsedGraphViewer->Show(false);

            displays["Collapsed graph"] = collapsedGraphViewer;
        }
        if(displayFinalGraph && ((dm==DisplayGraphWithoutAnalysis && !displayGraphAlreadyDisplayed) || dm==DisplayGraphWithAnalysis)) {
            GraphViewer *finalGraphViewer = new GraphViewer();
            finalGraphViewer->SetName("Final graph");
            if(dm==DisplayGraphWithoutAnalysis && !displayGraphAlreadyDisplayed)
                finalGraphViewer->SetGraphs(extractPipeline->GetFinalGraphs());
            else if(dm==DisplayGraphWithAnalysis)
                finalGraphViewer->SetGraphs(analyzePipeline->GetGraphs());
            finalGraphViewer->SetScaling(voxelSpacing[0], voxelSpacing[1], voxelSpacing[2]);
            finalGraphViewer->Show(false);

            displays["Final graph"] = finalGraphViewer;
        }
    }
    break;

    case PipelineAnalyzeCells:
    case PipelineApproximateCellShape:
    case PipelineCLAHE:
    case PipelineCrop:
    case PipelineSegmentNecroticRegion:
    case PipelineSegmentNuclei20x:
    case PipelineSegmentNuclei60x:
    case PipelineSegmentSinusoidsAndBile20x:
    case PipelineSegmentSinusoidsAndBile60x:
    case PipelineSegmentVeins:
#ifndef CS_TI_QUANT_ONLY
    case PipelineBackgroundElimination:
    case PipelineHConvexImageFilter:
    case PipelineAnalyzeHoechstDiffusion:
    case PipelineSegmentNucleiHoughTransformation:
    case PipelineVisualize3D:
    case PipelineObjectBasedSegmentation:
    {
        unsigned int ImageDimension = 3;

        CSParameter *dataDimParam = pipelineContext->findParameter("Data dimensionality", 0);
        if(dataDimParam != NULL) {
            std::string dataDimMode = ( (CSParameterChoice*)(dataDimParam->dataPointer()) )->currentString();
            if( dataDimMode.compare("2D")==0 )
                ImageDimension = 2;
            else if( dataDimMode.compare("3D")==0 )
                ImageDimension = 3;
            else
                throw std::string("Unknown data dimensionality.");
        }

        double voxelSpacing[3];
        voxelSpacing[0] = ( *(double*)(pipelineContext->findParameter("Image voxel spacing x", 0)->dataPointer()) );
        voxelSpacing[1] = ( *(double*)(pipelineContext->findParameter("Image voxel spacing y", 0)->dataPointer()) );
        if(ImageDimension==3)
            voxelSpacing[2] = ( *(double*)(pipelineContext->findParameter("Image voxel spacing z", 0)->dataPointer()) );
        else
            voxelSpacing[2] = 1.;

        std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs;

        if(ImageDimension==2) {
            ObjectBasedSegmentation<2> * objSegmentationPipeline = static_cast<ObjectBasedSegmentation<2> * >(pipelineInput);
            graphs.push_back(objSegmentationPipeline->GetObjectGraph());
        }
        else if(ImageDimension==3) {
            ObjectBasedSegmentation<3> * objSegmentationPipeline = static_cast<ObjectBasedSegmentation<3> * >(pipelineInput);
            graphs.push_back(objSegmentationPipeline->GetObjectGraph());
        }

        GraphViewer *graphViewer = new GraphViewer();
        graphViewer->SetName("Object graph");
        graphViewer->SetGraphs(graphs);
        graphViewer->SetScaling(voxelSpacing[0], voxelSpacing[1], voxelSpacing[2]);
        graphViewer->Show(false);

        displays["Object Graph"] = graphViewer;
    }
    break;
    case PipelineTrainClassifier:
    case PipelineSegmentCellMembrane:
    case PipelineObjectsToMatrix:
    case PipelineCompareSegmentations:
    case PipelineSegmentStellateCellCytoskeleton:
    case PipelineAnalyzeStellateCells:
    case PipelineSkeletonization:
    case PipelineTest:
#endif
    default:
        std::cout << "ImageProcessing::displayResult():  pipeline results not handled for \""
                  << PipelineName[mpWorker->mPipelineID] <<"\".\n";
        break;
    }
}


void ImageProcessing::processPipelineError( QString errorMsg )
{
    // for testing:
    std::cout << "An error occured:  " << errorMsg.toStdString() << std::endl;

    QMessageBox mb;
    mb.setText(errorMsg);
    mb.exec();

    pipelineFinished();
    mpWorkerThread->quit();
}


void ImageProcessing::selectionPipelineChanged(int index)
{
    seedPointSelectionPushButton->hide();

    switch ( index )
    {
    case PipelineAnalyzeCells:
        mAnalyzeCellsContext->setupGUI( pipelineEditorTreeView );
        break;

    case PipelineApproximateCellShape:
        m_estimateCellShapeContext->setupGUI( pipelineEditorTreeView );
        break;

    case PipelineCLAHE:
        m_CLAHEContext->setupGUI( pipelineEditorTreeView );
        break;

    case PipelineCrop:
        m_cropImageContext->setupGUI( pipelineEditorTreeView );
        break;

    case PipelineExtractAndAnalyzeGraph:
        m_extractGraphContext->setupGUI( pipelineEditorTreeView );
        break;

    case PipelineSegmentNecroticRegion:
        m_segmentNecroticRegionContext->setupGUI( pipelineEditorTreeView );
        break;

    case PipelineSegmentNuclei20x:
        m_segmentNuclei20xContext->setupGUI( pipelineEditorTreeView );
        break;

    case PipelineSegmentNuclei60x:
        m_segmentNuclei60xContext->setupGUI( pipelineEditorTreeView );
        break;

    case PipelineSegmentSinusoidsAndBile20x:
        m_segmentSinBile20xContext->setupGUI( pipelineEditorTreeView );
        break;

    case PipelineSegmentSinusoidsAndBile60x:
        m_segmentSinBile60xContext->setupGUI( pipelineEditorTreeView );
        break;

    case PipelineSegmentVeins:
        m_segmentVeinsContext->setupGUI( pipelineEditorTreeView );
        seedPointSelectionPushButton->show();
        break;

#ifndef CS_TI_QUANT_ONLY

    case PipelineAddImages:
        m_AddImagesContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineAnalyzeNuclei:
        mAnalyzeNucleiContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineAnalyzeLobule:
        mAnalyzeLobuleContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineApproximateLobuleShape:
        m_estimateLobuleShapeContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineBackgroundElimination:
        m_backgroundEliminationContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineHConvexImageFilter:
        m_convexFilterContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineVisualize3D:
        m_visualize3DContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineAnalyzeHoechstDiffusion:
        mAnalyzeHoechstDiffContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineSegmentNucleiHoughTransformation:
        m_segmentNucleiWithHoughContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineSuperpixelPreparation:
        mpSuperpixelPreparationContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineObjectBasedSegmentation:
        mpObjectBasedSegmentationContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineTrainClassifier:
        mpTrainClassifierContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineSegmentCellMembrane:
        mpSegmentCellMembraneContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineObjectsToMatrix:
        mpObjectsToMatrixContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineCompareSegmentations:
        mpCompareSegmentationsContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineSegmentStellateCellCytoskeleton:
        m_segmentStellateCellsContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineAnalyzeStellateCells:
        m_analyzeStellateCellsContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineSkeletonization:
        m_skeletonizationContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineFileFormatConverter:
        mpFileFormatConverterContext->setupGUI( pipelineEditorTreeView );
        break;
    case PipelineTest:
        mpTestPipelineContext->setupGUI( pipelineEditorTreeView );       //there is no test context so far
        break;
#endif
    }

    pipelineEditorTreeView->resizeColumnToContents(0);
    pipelineEditorTreeView->resizeColumnToContents(1);
    pipelineEditorTreeView->setAlternatingRowColors(true);
    pipelineEditorTreeView->expandAll();

    connect( ((QAbstractItemModel*)(pipelineEditorTreeView->model())), SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)), this, SLOT(pipelineEditorTreeViewChanged(const QModelIndex&, const QModelIndex&)) );
}

void ImageProcessing::addJobButtonClicked()
{
    switch ( selectPipelineComboBox->currentIndex() )
    {
    case PipelineApproximateCellShape:
	    jobManager->AddJobs(JobManager::EstimateCellShapeJob, m_estimateCellShapeContext);
        break;

    case PipelineCLAHE:
	    jobManager->AddJobs(JobManager::CLAHEJob, m_CLAHEContext);
        break;

    case PipelineCrop:
	    jobManager->AddJobs(JobManager::CropDataSetJob, m_cropImageContext);
        break;

    case PipelineExtractAndAnalyzeGraph:
        jobManager->AddJobs(JobManager::ExtractAndAnalyzeGraphJob, m_extractGraphContext);
        break;

    case PipelineAnalyzeCells:
        jobManager->AddJobs(JobManager::AnalyzeCellsJob, mAnalyzeCellsContext);
        break;

    case PipelineSegmentNecroticRegion:
	    jobManager->AddJobs(JobManager::SegmentNecroticRegionJob, m_segmentNecroticRegionContext);
        break;

    case PipelineSegmentNuclei20x:
	    jobManager->AddJobs(JobManager::SegmentNuclei20xJob, m_segmentNuclei20xContext);
        break;

    case PipelineSegmentNuclei60x:
	    jobManager->AddJobs(JobManager::SegmentNuclei60xJob, m_segmentNuclei60xContext);
        break;

    case PipelineSegmentSinusoidsAndBile20x:
	    jobManager->AddJobs(JobManager::SegmentSinusoidsBile20xJob, m_segmentSinBile20xContext);
        break;

    case PipelineSegmentSinusoidsAndBile60x:
	    jobManager->AddJobs(JobManager::SegmentSinusoidsBile60xJob, m_segmentSinBile60xContext);
        break;
    case PipelineSegmentVeins:
    {
        QString warning;
        warning = QString::fromStdString("Warning: The Job Manager doesn't support this pipeline.");

        QMessageBox mb;
        mb.setText(warning);
        mb.exec();
        break;
    }

#ifndef CS_TI_QUANT_ONLY

    case PipelineAnalyzeNuclei:
        jobManager->AddJobs(JobManager::AnalyzeNucleiJob, mAnalyzeNucleiContext);
        break;

    case PipelineAnalyzeLobule:
        jobManager->AddJobs(JobManager::AnalyzeLobuleJob, mAnalyzeLobuleContext);
        break;

    case PipelineApproximateLobuleShape:
        jobManager->AddJobs(JobManager::EstimateLobuleShapeJob, m_estimateLobuleShapeContext);
        break;

    case PipelineBackgroundElimination:
	    jobManager->AddJobs(JobManager::BackgroundEliminationJob, m_backgroundEliminationContext);
        break;

    case PipelineHConvexImageFilter:
	    jobManager->AddJobs(JobManager::HConvexImageFilterJob, m_convexFilterContext);
        break;

    case PipelineSegmentNucleiHoughTransformation:
	    jobManager->AddJobs(JobManager::SegmentNucleiWithHoughJob, m_segmentNucleiWithHoughContext);
        break;

    case PipelineObjectBasedSegmentation:
	    jobManager->AddJobs(JobManager::ObjectBasedSegmentationJob, mpObjectBasedSegmentationContext);
        break;

    case PipelineSegmentStellateCellCytoskeleton:
	    jobManager->AddJobs(JobManager::SegmentAndClassifyStellateCellsJob, m_segmentStellateCellsContext);
        break;

    case PipelineSkeletonization:
	    jobManager->AddJobs(JobManager::SkeletonizationJob, m_skeletonizationContext);
        break;

    case PipelineSegmentCellMembrane:
        jobManager->AddJobs(JobManager::SegmentCellMembraneJob, mpSegmentCellMembraneContext);
        break;

    case PipelineFileFormatConverter:
        jobManager->AddJobs(JobManager::FileFormatConversionJob, mpFileFormatConverterContext);
        break;
#endif // CS_TI_QUANT_ONLY
    }
}


void ImageProcessing::startQueuedJobsButtonClicked()
{
    jobQueueProgressBar->show();
    mpJobManagerThread->start();
}


void ImageProcessing::loadJobsButtonClicked()
{
    std::string openDir = jobManager->GetJobFilePath();
    QString filename = QFileDialog::getOpenFileName(this,tr("Open Job File"),tr(openDir.c_str()),tr("*.txt"));
    if ( filename.length() )
        jobManager->InitWithJobQueueFile(filename.toStdString());
}


void ImageProcessing::pipelineEditorTreeViewChanged(const QModelIndex& topLeft, const QModelIndex& bottomRight)
{
    CSParameter* alteredParam = (CSParameter*)(topLeft.internalPointer());
    std::cout << "altered parameter name = " << alteredParam->name() << std::endl;
    CSParameterContext* paramContext = ( (QCSParameterModel*)(topLeft.model()) )->rootContext();
    std::cout << "parameter context name = " << paramContext->name() << std::endl;

    if(alteredParam->name().compare("Data dimensionality")==0) {
        std::string dataDimMode = ( (CSParameterChoice*)(paramContext->findParameter("Data dimensionality", 0)->dataPointer()) )->currentString();

        if(paramContext->findContext(PipelineName[PipelineCrop],0)!=NULL) {

            if(dataDimMode.compare("2D")==0) {
                paramContext->findParameter("z start")->setVisible(false);
                paramContext->findParameter("z end")->setVisible(false);
            }
            else if(dataDimMode.compare("3D")==0) {
                paramContext->findParameter("z start")->setVisible(true);
                paramContext->findParameter("z end")->setVisible(true);
            }
        }
        if(paramContext->findContext(PipelineName[PipelineSegmentNuclei60x],0)!=NULL) {

            if(dataDimMode.compare("2D")==0) {
                paramContext->findParameter("Voxel spacing z")->setVisible(false);
                paramContext->findParameter("Sample region size z")->setVisible(false);
                paramContext->findContext("1.3) Inverse Hole Filling on 1.2")->findParameter("Radius z")->setVisible(false);
                paramContext->findContext("1.5) Closing on 1.4")->findParameter("Kernel radius z")->setVisible(false);
                paramContext->findContext("1.6) Opening on 1.5")->findParameter("Kernel radius z")->setVisible(false);
            }
            else if(dataDimMode.compare("3D")==0) {
                paramContext->findParameter("Voxel spacing z")->setVisible(true);
                paramContext->findParameter("Sample region size z")->setVisible(true);
                paramContext->findContext("1.3) Inverse Hole Filling on 1.2")->findParameter("Radius z")->setVisible(true);
                paramContext->findContext("1.5) Closing on 1.4")->findParameter("Kernel radius z")->setVisible(true);
                paramContext->findContext("1.6) Opening on 1.5")->findParameter("Kernel radius z")->setVisible(true);
            }
        }
#ifndef CS_TI_QUANT_ONLY
        else if(paramContext->findContext(PipelineName[PipelineApproximateLobuleShape],0)!=NULL) {

            if(dataDimMode.compare("2D")==0)
                paramContext->findParameter("Voxel spacing z")->setVisible(false);
            else if(dataDimMode.compare("3D")==0)
                paramContext->findParameter("Voxel spacing z")->setVisible(true);
        }
        else if(paramContext->findContext(PipelineName[PipelineAnalyzeNuclei],0)!=NULL) {

            if(dataDimMode.compare("2D")==0)
                paramContext->findParameter("Voxel spacing z")->setVisible(false);
            else if(dataDimMode.compare("3D")==0)
                paramContext->findParameter("Voxel spacing z")->setVisible(true);
        }
        else if(paramContext->findContext(PipelineName[PipelineSuperpixelPreparation],0)!=NULL) {

            if(dataDimMode.compare("2D")==0) {
                paramContext->findParameter("Image voxel spacing z")->setVisible(false);
                paramContext->findParameter("Superpixel spacing z")->setVisible(false);
            }
            else if(dataDimMode.compare("3D")==0) {
                paramContext->findParameter("Image voxel spacing z")->setVisible(true);
                paramContext->findParameter("Superpixel spacing z")->setVisible(true);
            }
        }
        else if(paramContext->findContext(PipelineName[PipelineObjectBasedSegmentation],0)!=NULL) {

            if(dataDimMode.compare("2D")==0)
                paramContext->findParameter("Image voxel spacing z")->setVisible(false);
            else if(dataDimMode.compare("3D")==0)
                paramContext->findParameter("Image voxel spacing z")->setVisible(true);
        }
#endif
    }
    else if(alteredParam->name().compare("Input type")==0) {
        std::string inputType = ( (CSParameterChoice*)(paramContext->findParameter("Input type", 0)->dataPointer()) )->currentString();

        if(paramContext->findContext(PipelineName[PipelineExtractAndAnalyzeGraph],0)!=NULL) {

            if(inputType.compare("Skeleton Image")==0)
                paramContext->findContext("Graph Extraction Parameter", 0)->setVisible(true);
            else if(inputType.compare("Final Graph")==0)
                paramContext->findContext("Graph Extraction Parameter", 0)->setVisible(false);
        }
    }
#ifndef CS_TI_QUANT_ONLY
    else if(alteredParam->name().compare("Input file format")==0) {
        unsigned int currentIndex = ( (CSParameterChoice*)(paramContext->findParameter("Input file format", 0)->dataPointer()) )->currentIndex();

        if(paramContext->findContext(PipelineName[PipelineFileFormatConverter],0)!=NULL) {
            std::vector<FileFormatConverter::OutputData> outputDataFormats;
            std::vector<std::string> outputDataFormatNames;

            outputDataFormats = FileFormatConverter::AvailableConversions((FileFormatConverter::InputData)currentIndex);

            for(unsigned int i=0; i<outputDataFormats.size(); i++)
                outputDataFormatNames.push_back(FileFormatConverter::OutputDataName[outputDataFormats[i]]);
            mFileFormatConverterOutputChoice->setChoices(outputDataFormatNames, 0);
        }
    }
    else if(alteredParam->name().compare("Object graph filename")==0) {
        std::string graphFilename = *(std::string*)(paramContext->findParameter("Object graph filename", 0)->dataPointer());

        QFileInfo fileInfo;
        fileInfo.setFile(QString::fromStdString(graphFilename));

        if(paramContext->findContext(PipelineName[PipelineObjectBasedSegmentation],0)!=NULL) {
            if( fileInfo.exists() ) {
                vtkSmartPointer<vtkUndirectedGraph> graph = vtkSmartPointer<vtkUndirectedGraph>::New();

                vtkSmartPointer<vtkGraphReader> graphReader = vtkSmartPointer<vtkGraphReader>::New();
                graphReader->SetFileName(graphFilename.c_str());
                graphReader->Update();
                graphReader->GetOutput()->ToUndirectedGraph(graph);

                int i;
                for(i=0; i<LabelMapType::NumberUnaryFeatures; i++) {
                    bool hasFeature = LabelMapType::GraphHasFeature(graph, static_cast<UnaryFeatureType>(i));
                    paramContext->findParameter("With " + LabelMapType::UnaryFeatureName[i])->setVisible(hasFeature);
                    mObjectBasedSegmentationWithFeature[i] = hasFeature;
                }

                for(int j=0; j<LabelMapType::NumberBinaryFeatures; j++) {
                    bool hasFeature = LabelMapType::GraphHasFeature(graph, static_cast<BinaryFeatureType>(j));
                    paramContext->findParameter("With " + LabelMapType::BinaryFeatureName[j])->setVisible(hasFeature);
                    mObjectBasedSegmentationWithFeature[i+j] = hasFeature;
                }
            }
            else {
                int i;
                for(i=0; i<LabelMapType::NumberUnaryFeatures; i++) {
                    paramContext->findParameter("With " + LabelMapType::UnaryFeatureName[i])->setVisible(true);
                    mObjectBasedSegmentationWithFeature[i] = true;
                }

                for(int j=0; j<LabelMapType::NumberBinaryFeatures; j++) {
                    paramContext->findParameter("With " + LabelMapType::BinaryFeatureName[j])->setVisible(true);
                    mObjectBasedSegmentationWithFeature[i+j] = true;
                }
            }
        }
        else if(paramContext->findContext(PipelineName[PipelineTrainClassifier],0)!=NULL) {
            if( fileInfo.exists() ) {
                vtkSmartPointer<vtkUndirectedGraph> graph = vtkSmartPointer<vtkUndirectedGraph>::New();

                vtkSmartPointer<vtkGraphReader> graphReader = vtkSmartPointer<vtkGraphReader>::New();
                graphReader->SetFileName(graphFilename.c_str());
                graphReader->Update();
                graphReader->GetOutput()->ToUndirectedGraph(graph);

                int i;
                for(i=0; i<LabelMapType::NumberUnaryFeatures; i++) {
                    bool hasFeature = LabelMapType::GraphHasFeature(graph, static_cast<UnaryFeatureType>(i));
                    paramContext->findParameter("With " + LabelMapType::UnaryFeatureName[i])->setVisible(hasFeature);
                    mTrainClassifierWithFeature[i] = hasFeature;
                }

                for(int j=0; j<LabelMapType::NumberBinaryFeatures; j++) {
                    bool hasFeature = LabelMapType::GraphHasFeature(graph, static_cast<BinaryFeatureType>(j));
                    paramContext->findParameter("With " + LabelMapType::BinaryFeatureName[j])->setVisible(hasFeature);
                    mTrainClassifierWithFeature[i+j] = hasFeature;
                }
            }
            else {
                int i;
                for(i=0; i<LabelMapType::NumberUnaryFeatures; i++) {
                    paramContext->findParameter("With " + LabelMapType::UnaryFeatureName[i])->setVisible(true);
                    mTrainClassifierWithFeature[i] = true;
                }

                for(int j=0; j<LabelMapType::NumberBinaryFeatures; j++) {
                    paramContext->findParameter("With " + LabelMapType::BinaryFeatureName[j])->setVisible(true);
                    mTrainClassifierWithFeature[i+j] = true;
                }
            }
        }
    }
#endif
    else if(alteredParam->name().compare("Network type")==0) {
        std::string networkType = ( (CSParameterChoice*)(paramContext->findParameter("Network type", 0)->dataPointer()) )->currentString();

        if(paramContext->findContext(PipelineName[PipelineExtractAndAnalyzeGraph],0)!=NULL) {
            std::string folder = *(std::string*)(paramContext->findParameter("Dataset containing folder", 0)->dataPointer());
            std::string dataSetFileList = *(std::string*)(paramContext->findParameter("Dataset file list", 0)->dataPointer());
            std::string inputType = ( (CSParameterChoice*)(paramContext->findParameter("Input type", 0)->dataPointer()) )->currentString();

            if( networkType.compare("Bile Canaliculi Network")==0 ) {
                m_extractGraph_ResamplingFactor = 1;
                m_extractGraph_ResamplingMaxDist = 0.00;
                m_extractGraph_RemoveDeadEndsThreshold = 2.00;
                m_extractGraph_CollapseIsecNodes = 1.50;
                m_extractGraph_GeomPruning = 10.00;
                m_extractGraph_AnalysisFilePrefix = "Bile_";
                m_extractGraph_StaticRadius = 0.50;
                m_extractGraph_MaskDeadEnds = false;
            }
            else if( networkType.compare("Sinusoidal Network")==0 ) {
                m_extractGraph_ResamplingFactor = 1;
                m_extractGraph_ResamplingMaxDist = 0.00;
                m_extractGraph_RemoveDeadEndsThreshold = 7.00;
                m_extractGraph_CollapseIsecNodes = 3.00;
                m_extractGraph_GeomPruning = 10.00;
                m_extractGraph_AnalysisFilePrefix = "Sin_";
                m_extractGraph_StaticRadius = 2.50;
                m_extractGraph_MaskDeadEnds = true;
            }

            if( folder.compare("/path/to/dataset_containing_folder")!=0 ) {
                QString filename;
                QFileInfo fileInfo;

                if( networkType.compare("Bile Canaliculi Network")==0 ) {
                    m_extractGraph_GraphInputFile = folder + "Bile_graph0.txt";
                    m_extractGraph_SkeletonInputFile = folder + "bile_step6_skeleton.tif";          //60x case
                    m_extractGraph_NetworkSegmentationFile = folder + "bile_step6_bin.tif";         //60x case

                    fileInfo.setFile(QString::fromStdString(m_extractGraph_SkeletonInputFile));

                    if(!fileInfo.exists()) {
                        m_extractGraph_SkeletonInputFile = folder + "bile_step5_skeleton.tif";      //20x case
                        m_extractGraph_NetworkSegmentationFile = folder + "bile_step5_bin.tif";     //20x case
                    }
                }
                else {
                    m_extractGraph_GraphInputFile = folder + "Sin_graph0.txt";
                    m_extractGraph_SkeletonInputFile = folder + "sinus_step5_skeleton.tif";         //60x case
                    m_extractGraph_NetworkSegmentationFile = folder + "sinus_step5_bin.tif";        //60x case

                    fileInfo.setFile(QString::fromStdString(m_extractGraph_SkeletonInputFile));

                    if(!fileInfo.exists()) {
                        m_extractGraph_SkeletonInputFile = folder + "sinus_step4_skeleton.tif";      //20x case
                        m_extractGraph_NetworkSegmentationFile = folder + "sinus_step4_bin.tif";     //20x case
                    }
                }
                fileInfo.setFile(QString::fromStdString(m_extractGraph_SkeletonInputFile));
                if(!fileInfo.exists())  m_extractGraph_SkeletonInputFile = "/path/to/skeleton_image.tif";
                fileInfo.setFile(QString::fromStdString(m_extractGraph_GraphInputFile));
                if(!fileInfo.exists())  m_extractGraph_GraphInputFile = "/path/to/graph_0_file.txt";
                fileInfo.setFile(QString::fromStdString(m_extractGraph_NetworkSegmentationFile));
                if(!fileInfo.exists())  m_extractGraph_NetworkSegmentationFile = "/path/to/network_bin_file.tif";
            }

            if( dataSetFileList.compare("/path/to/file.ias")!=0 ) {
                if( networkType.compare("Bile Canaliculi Network")==0 ) {
                    m_extractGraph_SkeletonInputFile = ImageAnalysisSummaryFileIO::GetEntry(BileSkeleton, dataSetFileList);
                    m_extractGraph_GraphInputFile = ImageAnalysisSummaryFileIO::GetEntry(BileGraph, dataSetFileList);
                    m_extractGraph_NetworkSegmentationFile = ImageAnalysisSummaryFileIO::GetEntry(BileSegmentationBin, dataSetFileList);
                }
                else {
                    m_extractGraph_SkeletonInputFile = ImageAnalysisSummaryFileIO::GetEntry(SinusoidSkeleton, dataSetFileList);
                    m_extractGraph_GraphInputFile = ImageAnalysisSummaryFileIO::GetEntry(SinusoidGraph, dataSetFileList);
                    m_extractGraph_NetworkSegmentationFile = ImageAnalysisSummaryFileIO::GetEntry(SinusoidSegmentationBin, dataSetFileList);
                }
                if(m_extractGraph_SkeletonInputFile.empty()) m_extractGraph_SkeletonInputFile = "/path/to/skeleton_image.tif";
                if(m_extractGraph_GraphInputFile.empty()) m_extractGraph_GraphInputFile = "/path/to/graph_0_file.txt";
                if(m_extractGraph_NetworkSegmentationFile.empty()) m_extractGraph_NetworkSegmentationFile = "/path/to/network_bin_file.tif";
            }
        }
    }
    else if(alteredParam->name().compare("Dataset containing folder")==0) {
        std::string folder = *(std::string*)(paramContext->findParameter("Dataset containing folder", 0)->dataPointer());

        if(paramContext->findContext(PipelineName[PipelineExtractAndAnalyzeGraph],0)!=NULL) {
            std::string networkType = ( (CSParameterChoice*)(paramContext->findParameter("Network type", 0)->dataPointer()) )->currentString();

            QString filename;
            QFileInfo fileInfo;
            QDir dir;
            dir.setPath(QString::fromStdString(folder));

            if( dir.exists() ) {
                m_extractGraph_DataSetID = dir.dirName().toStdString();

                if( networkType.compare("Bile Canaliculi Network")==0 ) {
                    m_extractGraph_GraphInputFile = folder + "Bile_graph0.txt";
                    m_extractGraph_SkeletonInputFile = folder + "bile_step6_skeleton.tif";          //60x case
                    m_extractGraph_NetworkSegmentationFile = folder + "bile_step6_bin.tif";         //60x case

                    fileInfo.setFile(QString::fromStdString(m_extractGraph_SkeletonInputFile));

                    if(!fileInfo.exists()) {
                        m_extractGraph_SkeletonInputFile = folder + "bile_step5_skeleton.tif";      //20x case
                        m_extractGraph_NetworkSegmentationFile = folder + "bile_step5_bin.tif";     //20x case
                    }
                }
                else {
                    m_extractGraph_GraphInputFile = folder + "Sin_graph0.txt";
                    m_extractGraph_SkeletonInputFile = folder + "sinus_step5_skeleton.tif";         //60x case
                    m_extractGraph_NetworkSegmentationFile = folder + "sinus_step5_bin.tif";        //60x case

                    fileInfo.setFile(QString::fromStdString(m_extractGraph_SkeletonInputFile));

                    if(!fileInfo.exists()) {
                        m_extractGraph_SkeletonInputFile = folder + "sinus_step4_skeleton.tif";      //20x case
                        m_extractGraph_NetworkSegmentationFile = folder + "sinus_step4_bin.tif";     //20x case
                    }
                }
                fileInfo.setFile(QString::fromStdString(m_extractGraph_SkeletonInputFile));
                if(!fileInfo.exists())  m_extractGraph_SkeletonInputFile = "/path/to/skeleton_image.tif";
                fileInfo.setFile(QString::fromStdString(m_extractGraph_GraphInputFile));
                if(!fileInfo.exists())  m_extractGraph_GraphInputFile = "/path/to/graph_0_file.txt";
                fileInfo.setFile(QString::fromStdString(m_extractGraph_NetworkSegmentationFile));
                if(!fileInfo.exists())  m_extractGraph_NetworkSegmentationFile = "/path/to/network_bin_file.tif";

                m_extractGraph_CVSegmentationFile = folder + "vein_central_bin.tif";
                m_extractGraph_PVSegmentationFile = folder + "vein_portal_bin.tif";
                m_extractGraph_CellSegmentationFile = folder + "cellShape_step3_bin.tif";
                m_extractGraph_NecRegSegFile = folder + "necroticRegion_step7_bin.tif";

                fileInfo.setFile(QString::fromStdString(m_extractGraph_CVSegmentationFile));
                if(!fileInfo.exists())  m_extractGraph_CVSegmentationFile = "/path/to/vein_central_bin_file.tif";
                fileInfo.setFile(QString::fromStdString(m_extractGraph_PVSegmentationFile));
                if(!fileInfo.exists())  m_extractGraph_PVSegmentationFile = "/path/to/vein_portal_bin_file.tif";
                fileInfo.setFile(QString::fromStdString(m_extractGraph_CellSegmentationFile));
                if(!fileInfo.exists())  m_extractGraph_CellSegmentationFile = "/path/to/cell_bin_file.tif";
                fileInfo.setFile(QString::fromStdString(m_extractGraph_NecRegSegFile));
                if(!fileInfo.exists())  m_extractGraph_NecRegSegFile = "/path/to/necReg_bin_file.tif";
            }
        }
    }
    else if(alteredParam->name().compare("Dataset file list")==0) {
        std::string dataSetFileList = *(std::string*)(paramContext->findParameter("Dataset file list", 0)->dataPointer());

        if(paramContext->findContext(PipelineName[PipelineExtractAndAnalyzeGraph],0)!=NULL) {
            std::string networkType = ( (CSParameterChoice*)(paramContext->findParameter("Network type", 0)->dataPointer()) )->currentString();

            if( dataSetFileList.compare("/path/to/file.ias")!=0 ) {
                if( networkType.compare("Bile Canaliculi Network")==0 ) {
                    m_extractGraph_SkeletonInputFile = ImageAnalysisSummaryFileIO::GetEntry(BileSkeleton, dataSetFileList);
                    m_extractGraph_GraphInputFile = ImageAnalysisSummaryFileIO::GetEntry(BileGraph, dataSetFileList);
                    m_extractGraph_NetworkSegmentationFile = ImageAnalysisSummaryFileIO::GetEntry(BileSegmentationBin, dataSetFileList);
                }
                else {
                    m_extractGraph_SkeletonInputFile = ImageAnalysisSummaryFileIO::GetEntry(SinusoidSkeleton, dataSetFileList);
                    m_extractGraph_GraphInputFile = ImageAnalysisSummaryFileIO::GetEntry(SinusoidGraph, dataSetFileList);
                    m_extractGraph_NetworkSegmentationFile = ImageAnalysisSummaryFileIO::GetEntry(SinusoidSegmentationBin, dataSetFileList);
                }
                if(m_extractGraph_SkeletonInputFile.empty()) m_extractGraph_SkeletonInputFile = "/path/to/skeleton_image.tif";
                if(m_extractGraph_GraphInputFile.empty()) m_extractGraph_GraphInputFile = "/path/to/graph_0_file.txt";
                if(m_extractGraph_NetworkSegmentationFile.empty()) m_extractGraph_NetworkSegmentationFile = "/path/to/network_bin_file.tif";

                m_extractGraph_CVSegmentationFile = ImageAnalysisSummaryFileIO::GetEntry(CentralVeinSegmentationBin, dataSetFileList);
                if(m_extractGraph_CVSegmentationFile.empty()) m_extractGraph_CVSegmentationFile = "/path/to/vein_central_bin_file.tif";
                m_extractGraph_PVSegmentationFile = ImageAnalysisSummaryFileIO::GetEntry(PortalVeinSegmentationBin, dataSetFileList);
                if(m_extractGraph_PVSegmentationFile.empty()) m_extractGraph_PVSegmentationFile = "/path/to/vein_portal_bin_file.tif";
                m_extractGraph_CellSegmentationFile = ImageAnalysisSummaryFileIO::GetEntry(CellShapeBin, dataSetFileList);
                if(m_extractGraph_CellSegmentationFile.empty()) m_extractGraph_CellSegmentationFile = "/path/to/cell_bin_file.tif";
                m_extractGraph_NecRegSegFile = ImageAnalysisSummaryFileIO::GetEntry(NecroticRegionSegmentationBin, dataSetFileList);
                if(m_extractGraph_NecRegSegFile.empty()) m_extractGraph_NecRegSegFile = "/path/to/necReg_bin_file.tif";
            }
        }
        else if(paramContext->findContext(PipelineName[PipelineApproximateCellShape],0)!=NULL) {

            if( dataSetFileList.compare("/path/to/file.ias")!=0 ) {
                m_estimateCellShape_Nuclei_Filename = ImageAnalysisSummaryFileIO::GetEntry(HepNucleiSegmentationBin, dataSetFileList);
                if(m_estimateCellShape_Nuclei_Filename.empty()) m_estimateCellShape_Nuclei_Filename = "/path/to/hep_nuclei_bin_file.tif";
                m_estimateCellShape_Bile_Filename = ImageAnalysisSummaryFileIO::GetEntry(BileSegmentationBin, dataSetFileList);
                if(m_estimateCellShape_Bile_Filename.empty()) m_estimateCellShape_Bile_Filename = "/path/to/bile_bin_file.tif";
                m_estimateCellShape_Sinusoid_Filename = ImageAnalysisSummaryFileIO::GetEntry(SinusoidSegmentationBin, dataSetFileList);
                if(m_estimateCellShape_Sinusoid_Filename.empty()) m_estimateCellShape_Sinusoid_Filename = "/path/to/sin_bin_file.tif";
                m_estimateCellShape_CV_Filename = ImageAnalysisSummaryFileIO::GetEntry(CentralVeinSegmentationBin, dataSetFileList);
                if(m_estimateCellShape_CV_Filename.empty()) m_estimateCellShape_CV_Filename = "/path/to/CV_bin_file.tif --optional";
                m_estimateCellShape_PV_Filename = ImageAnalysisSummaryFileIO::GetEntry(PortalVeinSegmentationBin, dataSetFileList);
                if(m_estimateCellShape_PV_Filename.empty()) m_estimateCellShape_PV_Filename = "/path/to/PV_bin_file.tif --optional";
            }
        }
        else if(paramContext->findContext(PipelineName[PipelineAnalyzeCells],0)!=NULL) {

            if( dataSetFileList.compare("/path/to/file.ias")!=0 ) {
                mAnalyzeCellsCellShapeFilename = ImageAnalysisSummaryFileIO::GetEntry(CellShapeBin, dataSetFileList);
                if(mAnalyzeCellsCellShapeFilename.empty()) mAnalyzeCellsCellShapeFilename = "/path/to/cellShape_bin_file.tif";
                mAnalyzeCellsNucleiFilename = ImageAnalysisSummaryFileIO::GetEntry(HepNucleiSegmentationBin, dataSetFileList);
                if(mAnalyzeCellsNucleiFilename.empty()) mAnalyzeCellsNucleiFilename = "/path/to/nuclei_bin_file.tif";
                mAnalyzeCellsBileFilename = ImageAnalysisSummaryFileIO::GetEntry(BileSegmentationBin, dataSetFileList);
                if(mAnalyzeCellsBileFilename.empty()) mAnalyzeCellsBileFilename = "/path/to/bile_bin_file.tif";
                mAnalyzeCellsSinusoidFilename = ImageAnalysisSummaryFileIO::GetEntry(SinusoidSegmentationBin, dataSetFileList);
                if(mAnalyzeCellsSinusoidFilename.empty()) mAnalyzeCellsSinusoidFilename = "/path/to/sinusoid_bin_file.tif";
                mAnalyzeCellsCVFilename = ImageAnalysisSummaryFileIO::GetEntry(CentralVeinSegmentationBin, dataSetFileList);
                if(mAnalyzeCellsCVFilename.empty()) mAnalyzeCellsCVFilename = "/path/to/CV_bin_file.tif --optional";
                mAnalyzeCellsPVFilename = ImageAnalysisSummaryFileIO::GetEntry(PortalVeinSegmentationBin, dataSetFileList);
                if(mAnalyzeCellsPVFilename.empty()) mAnalyzeCellsPVFilename = "/path/to/PV_bin_file.tif --optional";
            }
        }
    }
}


void ImageProcessing::seedPointSelectionButtonClicked()
{
	QString dppivChannel = QString::fromStdString( *(std::string*)(m_segmentVeinsContext->findParameter("DPPIV channel", 0)->dataPointer()) );
	QFileInfo fileInfo;
	fileInfo.setFile(dppivChannel);

	if(!fileInfo.exists()) {
		QString messErr = QString("Please specify the dppiv channel file first!");
		
		QMessageBox mb;
		mb.setText(messErr);
		mb.exec();
	}
	else {
		m_CVSeedPoints.clear();
		m_PVSeedPoints.clear();
		m_numberCVSeedPoints = m_CVSeedPoints.size();
		m_numberPVSeedPoints = m_PVSeedPoints.size();

		displays["Seed Point Selection"] = new QCSVTKDisplay();
		displays["Seed Point Selection"]->setWindowTitle("Seed Point Selection on DPPIV Channel");
		displays["Seed Point Selection"]->setGeometry(QRect(0, 0, 800, 600));

		QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Preferred);
		displays["Seed Point Selection"]->setSizePolicy(sizePolicy1);

		QWidget* imageViewerSidebar = new QWidget(displays["Seed Point Selection"]);
		imageViewerSidebar->setObjectName(QString::fromUtf8("seedPointSelectionSidebar"));

		displays["Seed Point Selection"]->addControlWidget(imageViewerSidebar);

		QSizePolicy sizePolicy2(QSizePolicy::Fixed, QSizePolicy::Minimum);
		sizePolicy2.setHorizontalStretch(0);
		sizePolicy2.setVerticalStretch(0);
		sizePolicy2.setHeightForWidth(imageViewerSidebar->sizePolicy().hasHeightForWidth());
		imageViewerSidebar->setSizePolicy(sizePolicy2);

		QVBoxLayout* verticalLayout = new QVBoxLayout(imageViewerSidebar);
		verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));

		mp_finishCVSeedPointSelectionButton = new QPushButton(imageViewerSidebar);
		mp_finishCVSeedPointSelectionButton->setObjectName(QString::fromUtf8("finishCVSeedPointSelectionButton"));
		mp_finishCVSeedPointSelectionButton->setText(QApplication::translate("seedPointSelectionSidebar", "Finish CV Seed Point Selection", 0, QApplication::UnicodeUTF8));

		verticalLayout->addWidget(mp_finishCVSeedPointSelectionButton);
		connect( mp_finishCVSeedPointSelectionButton, SIGNAL(clicked(bool)), this, SLOT(finishCVSeedPointSelectionButtonClicked()) );

		mp_finishPVSeedPointSelectionButton = new QPushButton(imageViewerSidebar);
		mp_finishPVSeedPointSelectionButton->setObjectName(QString::fromUtf8("finishPVSeedPointSelectionButton"));
		mp_finishPVSeedPointSelectionButton->setText(QApplication::translate("seedPointSelectionSidebar", "Finish PV Seed Point Selection", 0, QApplication::UnicodeUTF8));
		mp_finishPVSeedPointSelectionButton->setEnabled(false);

		verticalLayout->addWidget(mp_finishPVSeedPointSelectionButton);
		connect( mp_finishPVSeedPointSelectionButton, SIGNAL(clicked(bool)), this, SLOT(finishPVSeedPointSelectionButtonClicked()) );

		QPushButton* finishSeedPointSelectionButton = new QPushButton(imageViewerSidebar);
		finishSeedPointSelectionButton->setObjectName(QString::fromUtf8("finishSeedPointSelectionButton"));
		finishSeedPointSelectionButton->setText(QApplication::translate("seedPointSelectionSidebar", "Finish Seed Point Selection", 0, QApplication::UnicodeUTF8));

		verticalLayout->addWidget(finishSeedPointSelectionButton);
		connect( finishSeedPointSelectionButton, SIGNAL(clicked(bool)), this, SLOT(finishSeedPointSelectionButtonClicked()) );
		//-----------------------------------------------------------

		//TODO: encapsualte all this in a imageData and imageDataIO wrapper class
		try {
			typedef unsigned char                                       PixelCharType;                              //workaround: itk-vtk-connector can not handle images templated over bool pixels
			typedef itk::Image<PixelCharType, 3>                        CScalar3DImageType;                         //therefore read the (expected!) bool image and recast it to something that the connector understands

			typedef itk::ImageFileReader<CScalar3DImageType>            CScalar3DReaderType;
			typedef itk::ImageToVTKImageFilter<CScalar3DImageType>      ITKVTKScalar3DConnectorType;

			if(BasePipeline<3>::GetNumberOfDimensions(dppivChannel.toStdString()) != 3)
			    throw std::string("Please specify a DPPIV channel file with three dimensions.");

			CScalar3DReaderType::Pointer reader = CScalar3DReaderType::New();
			reader->SetFileName( dppivChannel.toStdString() );
			reader->UseStreamingOff();
			reader->ReleaseDataFlagOn();

		#if (ITK_VERSION_MAJOR >= 4)
			reader->SetImageIO( itk::TIFFImageIO::New() );
		#endif

			ITKVTKScalar3DConnectorType::Pointer itkvtkConnectorInputImage = ITKVTKScalar3DConnectorType::New();
			itkvtkConnectorInputImage->SetInput(reader->GetOutput());
			itkvtkConnectorInputImage->ReleaseDataFlagOn();

			vtkSmartPointer<vtkImageFlip> flipYInputFilter = vtkSmartPointer<vtkImageFlip>::New();
			flipYInputFilter->SetFilteredAxis(1);                           // flip y axis
			flipYInputFilter->SetInput(itkvtkConnectorInputImage->GetOutput());
			flipYInputFilter->Update();
			//---------------------------------------------------------

			viewers["Seed Point Selection Viewer"] = ImageViewer2D::New();
			viewers["Seed Point Selection Viewer"]->SetInput(flipYInputFilter->GetOutput());
			viewers["Seed Point Selection Viewer"]->GetImageActor()->InterpolateOff();
			viewers["Seed Point Selection Viewer"]->SetColorLevel(255);
			viewers["Seed Point Selection Viewer"]->SetColorWindow(255);

			displays["Seed Point Selection"]->setRenderWindow(viewers["Seed Point Selection Viewer"]->GetRenderWindow());
			viewers["Seed Point Selection Viewer"]->SetupInteractor(displays["Seed Point Selection"]->getInteractor());
			displays["Seed Point Selection"]->show();
		}
		catch( std::string& err )
		{
		    std::stringstream s;
		    s << "Sorry. An error occured.\n" << err << std::endl;
		    QString messErr = QString::fromStdString(s.str());

		    QMessageBox mb;
		    mb.setText(messErr);
		    mb.exec();
		}
		catch( itk::ExceptionObject & err )
		{
			std::stringstream s;
			s << "Sorry. An error occured.\n" << err << std::endl;
			QString messErr = QString::fromStdString(s.str());
		
			QMessageBox mb;
			mb.setText(messErr);
			mb.exec();
		}
		catch( std::exception& err )
		{
			std::stringstream s;
			s << "Sorry. An error occured.\n" << err.what() << std::endl;
			QString messErr = QString::fromStdString(s.str());
		
			QMessageBox mb;
			mb.setText(messErr);
			mb.exec();
		}
	}
}

void ImageProcessing::finishCVSeedPointSelectionButtonClicked()
{
    m_CVSeedPoints.clear();

    viewers["Seed Point Selection Viewer"]->GetSeedPicker()->SetYAxisFlipCorrection(true);
    viewers["Seed Point Selection Viewer"]->GetSeedPicker()->GetSeedPoints(m_CVSeedPoints);

    std::cout << "CV seed points:" << std::endl;
    for(unsigned int i=0; i<m_CVSeedPoints.size(); i++) {
        std::cout << "index " << i << ": (" << m_CVSeedPoints[i][0] << ", " << m_CVSeedPoints[i][1] << ", " << m_CVSeedPoints[i][2] << ")" << std::endl;
    }
    viewers["Seed Point Selection Viewer"]->GetSeedPicker()->ClearSeedPoints();

    mp_finishCVSeedPointSelectionButton->setEnabled(false);
    mp_finishPVSeedPointSelectionButton->setEnabled(true);
    m_numberCVSeedPoints = m_CVSeedPoints.size();
    displays["Seed Point Selection"]->getVTKWidget()->setFocus();
}

void ImageProcessing::finishPVSeedPointSelectionButtonClicked()
{
    m_PVSeedPoints.clear();

    viewers["Seed Point Selection Viewer"]->GetSeedPicker()->SetYAxisFlipCorrection(true);
    viewers["Seed Point Selection Viewer"]->GetSeedPicker()->GetSeedPoints(m_PVSeedPoints);

    std::cout << "PV seed points:" << std::endl;
    for(unsigned int i=0; i<m_PVSeedPoints.size(); i++) {
        std::cout << "index " << i << ": (" << m_PVSeedPoints[i][0] << ", " << m_PVSeedPoints[i][1] << ", " << m_PVSeedPoints[i][2] << ")" << std::endl;
    }
    viewers["Seed Point Selection Viewer"]->GetSeedPicker()->ClearSeedPoints();

    mp_finishPVSeedPointSelectionButton->setEnabled(false);
    m_numberPVSeedPoints = m_PVSeedPoints.size();
}

void ImageProcessing::finishSeedPointSelectionButtonClicked()
{
	displays["Seed Point Selection"]->close();
	QCSVTKDisplay* displayToDelete = displays["Seed Point Selection"];

	displays.erase(displays.find("Seed Point Selection"));
	delete displayToDelete;
}


void ImageProcessing::moveObjectsWithSlider1(int newPos)
{
	std::cout << "slider 1 has new value: " << newPos << std::endl;

	std::map<std::string, vtkProp3D*>::iterator it;

	for(it=actors.begin(); it != actors.end(); it++) {
		double *pos = (*it).second->GetPosition();
		(*it).second->SetPosition((-1)*newPos, pos[1], pos[2]);
	}

#if !defined(__BUILD_MAC__)
	displays["Input image"]->getRenderWindow()->Render();
#endif
	displays["Input image"]->update();
}

void ImageProcessing::moveObjectsWithSlider2(int newPos)
{
	std::cout << "slider 2 has new value: " << newPos << std::endl;

	std::map<std::string, vtkProp3D*>::iterator it;

	for(it=actors.begin(); it != actors.end(); it++) {
		double *pos = (*it).second->GetPosition();
		(*it).second->SetPosition(pos[0], (-1)*newPos, pos[2]);
	}

#if !defined(__BUILD_MAC__)
	displays["Input image"]->getRenderWindow()->Render();
#endif
	displays["Input image"]->update();
}

void ImageProcessing::moveObjectsWithSlider3(int newPos)
{
	std::cout << "slider 3 has new value: " << newPos << std::endl;

	std::map<std::string, vtkProp3D*>::iterator it;

	for(it=actors.begin(); it != actors.end(); it++) {
		double *pos = (*it).second->GetPosition();
		(*it).second->SetPosition(pos[0], pos[1], (-1)*newPos);
	}

#if !defined(__BUILD_MAC__)
	displays["Input image"]->getRenderWindow()->Render();
#endif
	displays["Input image"]->update();
}


void ImageProcessing::toggleVisibilityBile(bool toggled)
{
    if(actors.count("Bile")>0)
        actors["Bile"]->SetVisibility( toggled? 0: 1 );

#if !defined(__BUILD_MAC__)
	displays["Input image"]->getRenderWindow()->Render();
#endif

	displays["Input image"]->update();
}

void ImageProcessing::toggleVisibilitySinusoids(bool toggled)
{
    if(actors.count("Sinusoids")>0)
        actors["Sinusoids"]->SetVisibility( toggled? 0: 1 );

#if !defined(__BUILD_MAC__)
	displays["Input image"]->getRenderWindow()->Render();
#endif

    displays["Input image"]->update();
}

void ImageProcessing::toggleVisibilityNuclei(bool toggled)
{
    if(actors.count("Nuclei")>0)
        actors["Nuclei"]->SetVisibility( toggled? 0: 1 );

#if !defined(__BUILD_MAC__)
	displays["Input image"]->getRenderWindow()->Render();
#endif

    displays["Input image"]->update();
}

void ImageProcessing::toggleVisibilityGraph(bool toggled)
{
    if(actors.count("GraphEdges")>0)
        actors["GraphEdges"]->SetVisibility( toggled? 0: 1 );
    if(actors.count("GraphNodes")>0)
        actors["GraphNodes"]->SetVisibility( toggled? 0: 1 );

#if !defined(__BUILD_MAC__)
	displays["Input image"]->getRenderWindow()->Render();
#endif

    displays["Input image"]->update();
}

void ImageProcessing::toggleVisibilityCellShape(bool toggled)
{
    if(actors.count("CellShape")>0)
        actors["CellShape"]->SetVisibility( toggled? 0: 1 );

#if !defined(__BUILD_MAC__)
    displays["Input image"]->getRenderWindow()->Render();
#endif

    displays["Input image"]->update();
}

void ImageProcessing::exportAsPOVButtonClicked()
{
	std::string filename = "../povScene.pov";

	vtkPOVExporter *povexp = vtkPOVExporter::New();
	povexp->SetRenderWindow(displays["Input image"]->getRenderWindow());
	povexp->SetFileName(filename.c_str());
	povexp->Write();
}

void ImageProcessing::deleteSelectedJobs()
{
    QAbstractItemModel * model = jobTableView->model();
    QModelIndexList indexList = jobTableView->selectionModel()->selectedRows();

    qSort(indexList);

    std::list<QModelIndex> selectedIndexes = indexList.toStdList();
    std::list<QModelIndex>::reverse_iterator indexRevIter = selectedIndexes.rbegin();
    for (; indexRevIter != selectedIndexes.rend(); ++indexRevIter )
        model->removeRow((*indexRevIter).row());

    toggleDeleteJobsButton();
}


void ImageProcessing::toggleDeleteJobsButton()
{
    QModelIndexList indexList = jobTableView->selectionModel()->selectedRows();

    deleteJobsButton->setEnabled(indexList.size());
}
