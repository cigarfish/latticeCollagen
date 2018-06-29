///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ImageProcessing.h                                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-11-09                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef IMAGEPROCESSING_H_
#define IMAGEPROCESSING_H_

#include <QCSVTKDisplay.h>

#include "itkIndex.h"

#include <vtkActor.h>

#include "../ImageViewer2D.h"
#include "../../images/filters/convertFilters/LabelImageToGraphFilter.h"
#include "../../images/filters/convertFilters/LabelImageToLabelMapGraphFilter.h"
#include "../../images/tools/LabelMapGraph.h"
#include "../../tools/JobManager.h"
#include "../../tools/parameters/CSParameterChoice.h"

#include "ui_ImageProcessing.h"

class CSParameterContext;

class ImageProcessingWorker;


class ImageProcessing : public QWidget, private Ui::ImageProcessing
{
	Q_OBJECT

    friend class ImageProcessingWorker;

public:

    enum Pipeline {
        PipelineAnalyzeCells,
        PipelineApproximateCellShape,
        PipelineCLAHE,
        PipelineCrop,
        PipelineExtractAndAnalyzeGraph,
        PipelineSegmentNecroticRegion,
        PipelineSegmentNuclei20x,
        PipelineSegmentNuclei60x,
        PipelineSegmentSinusoidsAndBile20x,
        PipelineSegmentSinusoidsAndBile60x,
        PipelineSegmentVeins,
#ifndef CS_TI_QUANT_ONLY
        PipelineAddImages,
        PipelineAnalyzeLobule,
        PipelineAnalyzeNuclei,
        PipelineApproximateLobuleShape,
        PipelineBackgroundElimination,
        PipelineHConvexImageFilter,
        PipelineAnalyzeHoechstDiffusion,
        PipelineSegmentNucleiHoughTransformation,
        PipelineVisualize3D,
        PipelineSuperpixelPreparation,
        PipelineObjectBasedSegmentation,
        PipelineTrainClassifier,
        PipelineSegmentCellMembrane,
        PipelineObjectsToMatrix,
        PipelineCompareSegmentations,
        PipelineSegmentStellateCellCytoskeleton,
        PipelineAnalyzeStellateCells,
        PipelineSkeletonization,
        PipelineFileFormatConverter,
        PipelineTest
#endif
    };

    enum DisplayMode {
        DisplayGraphWithoutAnalysis,
        DisplayGraphWithAnalysis
    };

    static const std::string PipelineName[];

	ImageProcessing(QWidget *parent=0);
	~ImageProcessing();

signals:
    //*********3d Visualization Widget*************************
	void setVisPageBileImagePath(QString);
	void setVisPageSinImagePath(QString);
	void setVisPageNucImagePath(QString);
	void setVisPageGraphPath(QString);
	void setVisPageCellShapeImagePath(QString);
	//*********************************************************

private slots:
    //*********Main Options Widget*****************************
    void selectionPipelineChanged(int);
    void startPipelineButtonClicked();
    void loadJobsButtonClicked();
	void addJobButtonClicked();
    void deleteSelectedJobs();
	void startQueuedJobsButtonClicked();
    void toggleDeleteJobsButton();
	//*********************************************************

	//*********Pipeline Editor slot****************************
	void pipelineEditorTreeViewChanged(const QModelIndex& topLeft, const QModelIndex& bottomRight);
	//*********************************************************

	//*********Seed Point Selection Widget*********************
	void seedPointSelectionButtonClicked();
	void finishCVSeedPointSelectionButtonClicked();
	void finishPVSeedPointSelectionButtonClicked();
	void finishSeedPointSelectionButtonClicked();
	//*********************************************************

	//*********3d Visualization Widget*************************
	void moveObjectsWithSlider1(int newPos);
	void moveObjectsWithSlider2(int newPos);
	void moveObjectsWithSlider3(int newPos);

    void toggleVisibilityBile(bool toggled);
    void toggleVisibilitySinusoids(bool toggled);
    void toggleVisibilityNuclei(bool toggled);
    void toggleVisibilityGraph(bool toggled);
    void toggleVisibilityCellShape(bool toggled);
	void exportAsPOVButtonClicked();
	//*********************************************************

    //*********Thread feedback*********************************
    void pipelineFinished();
    void processPipelineError(QString errorMsg);
    void displayResult( void * pipelineInput, int dm, bool displayGraphAlreadyDisplayed );
    //*********************************************************

protected:
    //Need to access LabelMapType:: members
    typedef itk::RGBPixel<unsigned char>                                            CRGBPixelType;
    typedef itk::Image<unsigned long int, 3>                                        LScalarImageType;
    typedef itk::Image<CRGBPixelType, 3>                                            CRGBImageType;
    typedef LabelImageToGraphFilter<3>                                              LabelImageToGraphFilterType;
    typedef itk::LabelImageToLabelMapGraphFilter<CRGBImageType, LScalarImageType>   LabelImageToLabelMapGraphFilterType;
    typedef typename LabelImageToLabelMapGraphFilterType::OutputImageType           LabelMapType;
    typedef typename LabelMapType::UnaryFeatures                                    UnaryFeatureType;
    typedef typename LabelMapType::BinaryFeatures                                   BinaryFeatureType;


    std::map<std::string, QCSVTKDisplay*> displays;     //TODO: deprecated with new view class system
	std::map<std::string, ImageViewer2D*> viewers;      //TODO: define base view class (with interaction & information widgets + different display widgets (ImageViewer2D, GraphLayoutView, Volume renderer, PolyData renderer, ..)
	std::map<std::string, vtkProp3D*> actors;

    std::vector< void * > smartPointerContainer;

	bool m_testOutput;
	bool m_isBusy;

	JobManager *jobManager;
    QThread    *mpJobManagerThread;

    ImageProcessingWorker * mpWorker;
    QThread               * mpWorkerThread;

	QSlider* slider1;
	QSlider* slider2;
	QSlider* slider3;

	QPushButton* mp_finishCVSeedPointSelectionButton;
	QPushButton* mp_finishPVSeedPointSelectionButton;


    CSParameterContext * mpParameterRoot;

	//Pipeline: Extract Graph+++++++++++++++++++++++++++++++++++
	CSParameterContext *m_extractGraphContext;

	CSParameterChoice *m_extractGraph_NetworkType;
	CSParameterChoice *m_extractGraph_InputType;

	std::string m_extractGraph_Folder;
	std::string m_extractGraph_DataSetFileList;
	std::string m_extractGraph_SkeletonInputFile;
	std::string m_extractGraph_GraphInputFile;

	int m_extractGraph_ResamplingFactor;
	double m_extractGraph_ResamplingMaxDist;
	double m_extractGraph_RemoveDeadEndsThreshold;
	bool m_extractGraph_RemoveAllDeadEnds;
	double m_extractGraph_CollapseIsecNodes;
	double m_extractGraph_GeomPruning;

	CSParameterChoice *m_extractGraph_AnalysisMode;
	std::string m_extractGraph_AnalysisFilePrefix;
	std::string m_extractGraph_DataSetID;
	std::string m_extractGraph_NetworkSegmentationFile;
	std::string m_extractGraph_CVSegmentationFile;
	std::string m_extractGraph_PVSegmentationFile;
	std::string m_extractGraph_CellSegmentationFile;
	bool m_extractGraph_MaskDeadEnds;
	std::string m_extractGraph_CVMaskFile;
	std::string m_extractGraph_PVMaskFile;
	std::string m_extractGraph_NecRegSegFile;
	double m_extractGraph_VoxSpacing[3];
	double m_extractGraph_StaticRadius;
	bool m_extractGraph_WriteToAnalysisFile;

	bool m_extractGraph_Display_VanillaGraph;
	bool m_extractGraph_Display_ResampledGraph;
	bool m_extractGraph_Display_PrunedGraph;
	bool m_extractGraph_Display_CollapsedGraph;
	bool m_extractGraph_Display_FinalGraph;

	std::string m_extractGraph_Filename;
	CSParameterChoice *m_extractGraph_Save_Mode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Sin & Bile 20x Segmentation+++++++++++++++++++++
	CSParameterContext *m_segmentSinBile20xContext;

	std::string m_segmentSinBile20x_DPPIV_Filename;
	std::string m_segmentSinBile20x_DMs_Filename;
	bool m_segmentSinBile20x_UseNecReg;
	std::string m_segmentSinBile20x_NecReg_Filename;
	std::string m_segmentSinBile20x_OtherInput_Filename;
	std::string m_segmentSinBile20x_CVSeg_Filename;
	std::string m_segmentSinBile20x_PVSeg_Filename;
	CSParameterChoice *m_segmentSin20x_DPPIV_Threshold_Mode;
	int m_segmentSin20x_DPPIV_AdapOtsuThreshold_RegionSize[3];
	int m_segmentSin20x_DPPIV_AdapOtsuThreshold_NumSamples;
	int m_segmentSin20x_DPPIV_ManualThreshold;
	CSParameterChoice *m_segmentSin20x_DMs_Threshold_Mode;
	int m_segmentSin20x_DMs_AdapOtsuThreshold_RegionSize[3];
	int m_segmentSin20x_DMs_AdapOtsuThreshold_NumSamples;
	int m_segmentSin20x_DMs_ManualThreshold;
	int m_segmentSin20x_HoleFilling_RegionSize[3];
	int m_segmentSin20x_HoleFilling_MajThreshold;
	int m_segmentSin20x_Closing_RadStrucElem;
	int m_segmentSin20x_Opening_RadStrucElem;
	int m_segmentSin20x_MaskVein_RegionSize[3];
	int m_segmentSin20x_MinObjSize;

	CSParameterChoice *m_segmentBile20x_DPPIV_Threshold_Mode;
	int m_segmentBile20x_DPPIV_AdapOtsuThreshold_RegionSize[3];
	int m_segmentBile20x_DPPIV_AdapOtsuThreshold_NumSamples;
	int m_segmentBile20x_DPPIV_ManualThreshold;
	int m_segmentBile20x_HoleFilling_RegionSize[3];
	int m_segmentBile20x_HoleFilling_MajThreshold;
	int m_segmentBile20x_InvHoleFilling_RegionSize[3];
	int m_segmentBile20x_InvHoleFilling_MajThreshold;
	int m_segmentBile20x_MaskVein_RegionSize[3];
	int m_segmentBile20x_MinObjSize;

	std::string m_segmentSinBile20x_Log_Filename;
	std::string m_segmentSin20x_Filename_Prefix;
	CSParameterChoice *m_segmentSin20x_Save_Mode;
	std::string m_segmentBile20x_Filename_Prefix;
	CSParameterChoice *m_segmentBile20x_Save_Mode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Sin & Bile 60x Segmentation+++++++++++++++++++++
	CSParameterContext *m_segmentSinBile60xContext;

	std::string m_segmentSinBile60x_DPPIV_Filename;
	std::string m_segmentSinBile60x_DMs_Filename;
	std::string m_segmentSinBile60x_CVSeg_Filename;
	std::string m_segmentSinBile60x_PVSeg_Filename;
	double m_segmentSinBile60x_VoxelSpacing[3];
	bool m_segmentSinBile60x_UseNecReg;
	std::string m_segmentSinBile60x_NecReg_Filename;
	std::string m_segmentSinBile60x_OtherInput_Filename;
	int m_segmentSin60x_Closing_RadStrucElem;
	CSParameterChoice *m_segmentSin60x_DPPIV_Threshold_Mode;
	int m_segmentSin60x_DPPIV_AdapOtsuThreshold_RegionSize[3];
	int m_segmentSin60x_DPPIV_AdapOtsuThreshold_NumSamples;
	int m_segmentSin60x_DPPIV_ManualThreshold;
	CSParameterChoice *m_segmentSin60x_DMs_Threshold_Mode;
	int m_segmentSin60x_DMs_AdapOtsuThreshold_RegionSize[3];
	int m_segmentSin60x_DMs_AdapOtsuThreshold_NumSamples;
	int m_segmentSin60x_DMs_ManualThreshold;
	int m_segmentSin60x_InvHoleFilling_RegionSize[3];
	int m_segmentSin60x_InvHoleFilling_MajThreshold;
	bool m_segmentSin60x_CavityFilling_WithRescaling;
	int m_segmentSin60x_CavityFilling_Radius;
	double m_segmentSin60x_CavityFilling_MinThreshold;
	double m_segmentSin60x_CavityFilling_MaxThreshold;
	int m_segmentSin60x_Closing_xRadStrucElem;
	int m_segmentSin60x_Closing_yRadStrucElem;
	int m_segmentSin60x_Closing_zRadStrucElem;
	int m_segmentSin60x_Opening_xRadStrucElem;
	int m_segmentSin60x_Opening_yRadStrucElem;
	int m_segmentSin60x_Opening_zRadStrucElem;
	int m_segmentSin60x_MaskVein_RegionSize[3];
	int m_segmentSin60x_MinObjSize;

	int m_segmentBile60x_Median_RadStrucElem;
	int m_segmentBile60x_GreyOpening_RadStrucElem;
	CSParameterChoice *m_segmentBile60x_DPPIV_Threshold_Mode;
	int m_segmentBile60x_DPPIV_AdapOtsuThreshold_RegionSize[3];
	int m_segmentBile60x_DPPIV_AdapOtsuThreshold_NumSamples;
	int m_segmentBile60x_DPPIV_ManualThreshold;
	int m_segmentBile60x_HoleFilling_RegionSize[3];
	int m_segmentBile60x_HoleFilling_MajThreshold;
	int m_segmentBile60x_BinOpening_RadStrucElem;
	int m_segmentBile60x_InvHoleFilling_RegionSize[3];
	int m_segmentBile60x_InvHoleFilling_MajThreshold;
	int m_segmentBile60x_MaskVein_RegionSize[3];
	int m_segmentBile60x_MinObjSize;

	std::string m_segmentSinBile60x_Log_Filename;
	std::string m_segmentSin60x_Filename_Prefix;
	CSParameterChoice *m_segmentSin60x_Save_Mode;
	std::string m_segmentBile60x_Filename_Prefix;
	CSParameterChoice *m_segmentBile60x_Save_Mode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Nuclei 20x Segmentation+++++++++++++++++++++++++
	CSParameterContext *m_segmentNuclei20xContext;

	std::string m_segmentNuclei20x_DAPI_Filename;
	std::string m_segmentNuclei20x_CVSeg_Filename;
	std::string m_segmentNuclei20x_PVSeg_Filename;
	bool m_segmentNuclei20x_UseNecReg;
	std::string m_segmentNuclei20x_NecReg_Filename;
	double m_segmentNuclei20x_VoxelSpacing[3];
	int m_segmentNuclei20x_Median_Radius;
	CSParameterChoice *m_segmentNuclei20x_DAPI_Threshold_Mode;
	int m_segmentNuclei20x_DAPI_AdapOtsuThreshold_RegionSize[3];
	int m_segmentNuclei20x_DAPI_AdapOtsuThreshold_NumSamples;
	int m_segmentNuclei20x_DAPI_ManualThreshold;
	int m_segmentNuclei20x_HoleFilling_RegionSize[3];
	int m_segmentNuclei20x_HoleFilling_MajThreshold;
	int m_segmentNuclei20x_Opening_RadStrucElem;
	double m_segmentNuclei20x_MinDiameterThreshold;
	double m_segmentNuclei20x_MaxDiameterThreshold;
	double m_segmentNuclei20x_WatershedFloodLevel;

	std::string m_segmentNuclei20x_Log_Filename;
	std::string m_segmentNuclei20x_Filename_Prefix;
	CSParameterChoice *m_segmentNuclei20x_Save_Mode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Nuclei 60x Segmentation+++++++++++++++++++++++++
	CSParameterContext *m_segmentNuclei60xContext;

	CSParameterChoice *m_segmentNuclei60x_DataDimMode;
	std::string m_segmentNuclei60x_DAPI_Filename;
	std::string m_segmentNuclei60x_CVSeg_Filename;
	std::string m_segmentNuclei60x_PVSeg_Filename;
	double m_segmentNuclei60x_VoxelSpacing[3];
	int m_segmentNuclei60x_Median_Radius;
	int m_segmentNuclei60x_Opening1_RadStrucElem;
	CSParameterChoice *m_segmentNuclei60x_DAPI_Threshold_Mode;
	int m_segmentNuclei60x_DAPI_AdapOtsuThreshold_RegionSize[3];
	int m_segmentNuclei60x_DAPI_AdapOtsuThreshold_NumSamples;
	int m_segmentNuclei60x_DAPI_ManualMinThreshold;
	int m_segmentNuclei60x_DAPI_ManualMaxThreshold;
	int m_segmentNuclei60x_InvHoleFilling_RegionSize[3];
	int m_segmentNuclei60x_InvHoleFilling_MajThreshold;
	bool m_segmentNuclei60x_CavityFilling_WithResampling;
	int m_segmentNuclei60x_CavityFilling_Radius;
	double m_segmentNuclei60x_CavityFilling_MinThreshold;
	double m_segmentNuclei60x_CavityFilling_MaxThreshold;
	int m_segmentNuclei60x_Closing_xRadStrucElem;
	int m_segmentNuclei60x_Closing_yRadStrucElem;
	int m_segmentNuclei60x_Closing_zRadStrucElem;
	int m_segmentNuclei60x_Opening_xRadStrucElem;
	int m_segmentNuclei60x_Opening_yRadStrucElem;
	int m_segmentNuclei60x_Opening_zRadStrucElem;
	double m_segmentNuclei60x_MinNonHepDiameterThreshold;
	double m_segmentNuclei60x_MaxNonHepDiameterThreshold;
	double m_segmentNuclei60x_MinHepDiameterThreshold;
	double m_segmentNuclei60x_MaxHepDiameterThreshold;
	double m_segmentNuclei60x_roundness;
	double m_segmentNuclei60x_WatershedFloodLevel;

	std::string m_segmentNuclei60x_Log_Filename;
	std::string m_segmentNuclei60x_Filename_Prefix;
	CSParameterChoice *m_segmentNuclei60x_Save_Mode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Nuclei Segmentation with Hough++++++++++++++++++
	CSParameterContext *m_segmentNucleiWithHoughContext;

	std::string m_segmentNucleiWithHough_DAPI_Filename;
	double m_segmentNucleiWithHough_VoxelSpacing[3];
	int m_segmentNucleiWithHough_Median_Radius;
	int m_segmentNucleiWithHough_Opening1_RadStrucElem;
	int m_segmentNucleiWithHough_Hough_NumSpheres;
	int m_segmentNucleiWithHough_Hough_Threshold;
	int m_segmentNucleiWithHough_Hough_GradientThreshold;
	double m_segmentNucleiWithHough_Hough_OutputThreshold;
	double m_segmentNucleiWithHough_Hough_MinRadius;
	double m_segmentNucleiWithHough_Hough_MaxRadius;
	double m_segmentNucleiWithHough_Hough_SigmaGradient;
	double m_segmentNucleiWithHough_Hough_Variance;
	double m_segmentNucleiWithHough_Hough_SphereRadiusRatio;
	double m_segmentNucleiWithHough_Hough_VotingRadiusRatio;
	double m_segmentNucleiWithHough_Hough_SamplingRatio;
	int m_segmentNucleiWithHough_SphereToImage_ResampleFactor;
	int m_segmentNucleiWithHough_SphereToImage_DilationRadius;

	double m_segmentNucleiWithHough_MinNonHepDiameterThreshold;
	double m_segmentNucleiWithHough_MinHepDiameterThreshold;
	double m_segmentNucleiWithHough_MaxHepDiameterThreshold;

	std::string m_segmentNucleiWithHough_Log_Filename;
	std::string m_segmentNucleiWithHough_Filename_Prefix;
	CSParameterChoice *m_segmentNucleiWithHough_Save_Mode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Estimate Cell Shape+++++++++++++++++++++++++++++
	CSParameterContext *m_estimateCellShapeContext;

	std::string m_estimateCellShape_DataSetFileList;
	std::string m_estimateCellShape_DPPIV_Filename;
	std::string m_estimateCellShape_Nuclei_Filename;
	std::string m_estimateCellShape_Bile_Filename;
	std::string m_estimateCellShape_Sinusoid_Filename;
	std::string m_estimateCellShape_CV_Filename;
	std::string m_estimateCellShape_PV_Filename;
	bool m_estimateCellShape_UseNecReg;
	std::string m_estimateCellShape_NecReg_Filename;
	double m_estimateCellShape_VoxelSpacing[3];
	double m_estimateCellShape_Bile_Weight;
	double m_estimateCellShape_MinCellDiameterThreshold;
	double m_estimateCellShape_WatershedFloodLevel;

	bool m_estimateCellShape_WithPositionFilter;
	int m_estimateCellShape_NumberCells;
	double m_estimateCellShape_Pos[3];

	std::string m_estimateCellShape_Log_Filename;
	std::string m_estimateCellShape_Filename_Prefix;
	CSParameterChoice *m_estimateCellShape_Save_Mode;
	bool m_estimateCellShape_WriteCellOutlineToFile;
	bool m_estimateCellShape_WriteSinusoidOutlineToFile;
	bool m_estimateCellShape_WriteSinusoidGraphToFile;
	std::string m_estimateCellShape_SinusoidSkeleton_Filename;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Estimate Lobule Shape+++++++++++++++++++++++++++++
	CSParameterContext *m_estimateLobuleShapeContext;

	CSParameterChoice *m_estimateLobuleShape_DataDimMode;
	std::string m_estimateLobuleShape_RawData_Filename;
	std::string m_estimateLobuleShape_CV_Filename;
	int m_estimateLobuleShape_CV_Intensity;
	std::string m_estimateLobuleShape_PV_Filename;
	int m_estimateLobuleShape_PV_Intensity;
	std::string m_estimateLobuleShape_Tissue_Filename;
	int m_estimateLobuleShape_Tissue_Intensity;
	int m_estimateLobuleShape_Void_Intensity;
	double m_estimateLobuleShape_VoxelSpacing[3];
	double m_estimateLobuleShape_PV_Weight;
	double m_estimateLobuleShape_MinLobuleDiameterThreshold;
	double m_estimateLobuleShape_MaxLobuleDiameterThreshold;
	double m_estimateLobuleShape_WatershedFloodLevel;
	CSParameterChoice *m_estimateLobuleShape_AnalysisMode;
	std::string m_estimateLobuleShape_DataSetID;

	std::string m_estimateLobuleShape_Log_Filename;
	std::string m_estimateLobuleShape_Filename_Prefix;
	CSParameterChoice *m_estimateLobuleShape_Save_Mode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Stellate Cell Segmentation++++++++++++++++++++++
	CSParameterContext *m_segmentStellateCellsContext;

	std::string m_segmentStellateCells_Desmin_Filename;
	std::string m_segmentStellateCells_DAPI_Filename;
	std::string m_segmentStellateCells_NonHepNuclei_Filename;
	std::string m_segmentStellateCells_HepNuclei_Filename;
	double m_segmentStellateCells_VoxelSpacing[3];
	int m_segmentStellateCells_Median_Radius;
	int m_segmentStellateCells_Opening_RadStrucElem;
	bool m_segmentStellateCells_WithConvexFilter;
	int m_segmentStellateCells_ConvexFilter_Height;
	CSParameterChoice *m_segmentStellateCells_Threshold_Mode;
	int m_segmentStellateCells_AdapOtsuThreshold_RegionSize[3];
	int m_segmentStellateCells_AdapOtsuThreshold_NumSamples;
	int m_segmentStellateCells_Desmin_ManualThreshold;

	int m_segmentStellateCells_InvHoleFilling_RegionSize[3];
	int m_segmentStellateCells_InvHoleFilling_MajThreshold;
	int m_segmentStellateCells_MinSizeThreshold;

	double m_segmentStellateCells_NucleiBodySuperposition;

	std::string m_segmentStellateCells_Log_Filename;
	std::string m_segmentStellateCells_Filename_Prefix;
	CSParameterChoice *m_segmentStellateCells_Save_Mode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Cell Analysis++++++++++++++++++++++++++++++++++
	CSParameterContext *mAnalyzeCellsContext;

	std::string mAnalyzeCellsDataSetID;
	std::string mAnalyzeCellsDataSetFileList;
	std::string mAnalyzeCellsCellShapeFilename;
	std::string mAnalyzeCellsNucleiFilename;
	std::string mAnalyzeCellsSinusoidFilename;
	std::string mAnalyzeCellsBileFilename;
	std::string mAnalyzeCellsCVFilename;
	std::string mAnalyzeCellsPVFilename;
	bool mAnalyzeCellsSaveAsGraph;

	double mAnalyzeCellsVoxSpacing[3];
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Nuclei Analysis+++++++++++++++++++++++++++++++++
	CSParameterContext *mAnalyzeNucleiContext;

	CSParameterChoice *mAnalyzeNucleiDataDimMode;
	std::string mAnalyzeNucleiDataSetID;
	std::string mAnalyzeNucleiHepNucleiFilename;
	std::string mAnalyzeNucleiNonHepNucleiFilename;
	int mAnalyzeNucleiNonProlifNucleiIntensity;
	int mAnalyzeNucleiProlifNucleiIntensity;
	std::string mAnalyzeNucleiCVFilename;
	int mAnalyzeNucleiCVIntensity;
	std::string mAnalyzeNucleiPVFilename;
	int mAnalyzeNucleiPVIntensity;
	std::string mAnalyzeNucleiTissueFilename;
	int mAnalyzeNucleiTissueIntensity;
	std::string mAnalyzeNucleiNecRegionFilename;
	int mAnalyzeNucleiNecRegionIntensity;
	int mAnalyzeNucleiVoidIntensity;
	std::string mAnalyzeNucleiLobuleFilename;
	int mAnalyzeNucleiLobulesIntensity;
	bool mAnalyzeNucleiSaveAsGraph;

	double mAnalyzeNucleiVoxSpacing[3];
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Lobule Analysis+++++++++++++++++++++++++++++++++
	CSParameterContext *mAnalyzeLobuleContext;

	std::string mAnalyzeLobuleCVFilename;
	std::string mAnalyzeLobulePVFilename;

	double mAnalyzeLobuleVoxSpacing[3];
	CSParameterChoice *mAnalyzeLobuleSaveMode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Hoechst Diffusion Analysis++++++++++++++++++++++
	CSParameterContext *mAnalyzeHoechstDiffContext;

	std::string mAnalyzeHoechstDiffDataSetID;
	std::string mAnalyzeHoechstDiffFilename;

	double mAnalyzeHoechstDiffTimeStep;
	double mAnalyzeHoechstDiffVoxSpacing[2];
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    //Pipeline: Stellate Cell Analysis++++++++++++++++++++++++++
    CSParameterContext *m_analyzeStellateCellsContext;

    std::string m_analyzeStellateCells_dataSetID;
    std::string m_analyzeStellateCells_HepSeg_Filename;
    std::string m_analyzeStellateCells_NonHepSeg_Filename;
    std::string m_analyzeStellateCells_StCSeg_Filename;
    std::string m_analyzeStellateCells_StCBodies_Filename;
    std::string m_analyzeStellateCells_Hepatocyte_Filename;
    std::string m_analyzeStellateCells_CV_Filename;
    std::string m_analyzeStellateCells_PV_Filename;

    double m_analyzeStellateCells_VoxSpacing[3];
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Necrotic Region Segmentation++++++++++++++++++++
	CSParameterContext *m_segmentNecroticRegionContext;

	std::string m_segmentNecroticRegion_DPPIV_Filename;
	std::string m_segmentNecroticRegion_DMs_Filename;
	CSParameterChoice *m_segmentNecroticRegion_DMs_Threshold_Mode;
	int m_segmentNecroticRegion_DMs_AdapOtsuThreshold_RegionSize[3];
	int m_segmentNecroticRegion_DMs_AdapOtsuThreshold_NumSamples;
	int m_segmentNecroticRegion_DMs_ManualThreshold;
	int m_segmentNecroticRegion_Erode1_RadStrucElem;
	int m_segmentNecroticRegion_Dilate1_RadStrucElem;
	int m_segmentNecroticRegion_Erode2_RadStrucElem;
	int m_segmentNecroticRegion_Dilate2_RadStrucElem;
	int m_segmentNecroticRegion_RemvObjSizeThreshold;

	std::string m_segmentNecroticRegion_Log_Filename;
	std::string m_segmentNecroticRegion_Filename_Prefix;
	CSParameterChoice *m_segmentNecroticRegion_Save_Mode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Visualization 3D++++++++++++++++++++++++++++++++
	CSParameterContext *m_visualize3DContext;

	std::string m_visualize3DBileFilename;
	std::string m_visualize3DSinusoidFilename;
	std::string m_visualize3DNucleiFilename;
	std::string m_visualize3DGraphFilename;
	std::string m_visualize3DCellShapeFilename;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Vein Segmentation+++++++++++++++++++++++++++++++
	CSParameterContext *m_segmentVeinsContext;

	std::string m_segmentVeins_DAPI_Filename;
	std::string m_segmentVeins_DPPIV_Filename;
	std::string m_segmentVeins_DMs_Filename;
	int m_segmentVeins_DAPI_Threshold;
	int m_segmentVeins_DPPIV_Threshold;
	int m_segmentVeins_DMs_Threshold;
	int m_segmentVeins_OpeningRadius;
	int m_segmentVeins_ClosingRadius;
	bool m_segmentVeins_With_Resampling;
	bool m_segmentVeins_With_Segmentation;
	bool m_segmentVeins_With_LobulePrediction;
	int m_numberCVSeedPoints;
	int m_numberPVSeedPoints;

    std::vector<itk::Index<3> > m_CVSeedPoints;
    std::vector<itk::Index<3> > m_PVSeedPoints;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    //Pipeline: Skeletonization+++++++++++++++++++++++++++++++++
    CSParameterContext *m_skeletonizationContext;

    std::string m_skeletonization_Filename;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: Crop image++++++++++++++++++++++++++++++++++++++
	CSParameterContext *m_cropImageContext;

	CSParameterChoice *m_cropImageContext_DataDimMode;
	std::string m_cropImage_Filename;
	int m_cropImage[3][2];
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	//Pipeline: CLAHE+++++++++++++++++++++++++++++++++++++++++++
	CSParameterContext *m_CLAHEContext;

	CSParameterChoice *m_CLAHE_DataDimMode;
	std::string m_CLAHE_Filename;
	unsigned int m_CLAHE_HistWindRad;
	unsigned int m_CLAHE_StepSize;
	double m_CLAHE_ClipLevel;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Background Elminiation++++++++++++++++++++++++++
	CSParameterContext *m_backgroundEliminationContext;

	CSParameterChoice *m_backgroundElimination_DataDimMode;
	std::string m_backgroundElimination_Filename;
	unsigned int m_backgroundElimination_MedianKernelRadius;
	unsigned int m_backgroundElimination_KernelRadius;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Convex Filter+++++++++++++++++++++++++++++++++++
	CSParameterContext *m_convexFilterContext;

	CSParameterChoice *m_convexFilter_DataDimMode;
	std::string m_convexFilter_Filename;
	unsigned int m_convexFilter_Height;
	bool m_convexFilter_FullyConnected;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Add images++++++++++++++++++++++++++++++++++++++
	CSParameterContext *m_AddImagesContext;

	CSParameterChoice *m_AddImages_DataDimMode;
	std::string m_AddImages_Filename1;
	std::string m_AddImages_Filename2;
	unsigned int m_AddImages_Intensity1;
	unsigned int m_AddImages_Intensity2;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Superpixel Preparation++++++++++++++++++++++++++
	CSParameterContext *mpSuperpixelPreparationContext;

	CSParameterChoice* mSuperpixelPreparationDataDimMode;
	CSParameterChoice* mSuperpixelPreparationStartingPoint;
	std::string mSuperpixelPreparationRawImageFilename;
	std::string mSuperpixelPreparationObjectImageFilename;
	std::string mSuperpixelPreparationObjectGraphFilename;
	CSParameterChoice* mSuperpixelPreparationSuperpixelMode;
	double mSuperpixelPreparationVoxSpacing[3];
	double mSuperpixelPreparationSuperpixelSpacing[3];
	unsigned int mSuperpixelPreparationSuperpixelNumber;
	bool* mSuperpixelPreparationWithFeature;

	std::string mSuperpixelPreparationLogFilename;
	std::string mSuperpixelPreparationFilenamePrefix;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Object-based Segmentation+++++++++++++++++++++++
	CSParameterContext *mpObjectBasedSegmentationContext;

	CSParameterChoice* mObjectBasedSegmentationDataDimMode;
	std::string mObjectBasedSegmentationRawImageFilename;
	std::string mObjectBasedSegmentationObjectImageFilename;
	std::string mObjectBasedSegmentationObjectGraphFilename;
	std::string mObjectBasedSegmentationTrainingDatabaseFilename;
	CSParameterChoice* mObjectBasedSegmentationSegmentationTargetMode;
	double mObjectBasedSegmentationVoxSpacing[3];
	CSParameterChoice* mObjectBasedSegmentationSegmentationMode;
	CSParameterChoice* mObjectBasedSegmentationClassifierChoice;
	bool mObjectBasedSegmentationKeepOnlyTargetClassObjects;
	bool* mObjectBasedSegmentationWithFeature;

	std::string mObjectBasedSegmentationLogFilename;
	std::string mObjectBasedSegmentationFilenamePrefix;
	CSParameterChoice* mpObjectBasedSegmentationSaveMode;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Pipeline: Train Classifier++++++++++++++++++++++++++++++++
    CSParameterContext *mpTrainClassifierContext;

    CSParameterChoice *mTrainDataDimMode;
    CSParameterChoice* mTrainActionChoice;
    CSParameterChoice* mTrainClassifierChoice;
    CSParameterChoice* mTrainClassLabelChoice;
    std::string mTrainClassifierOrigImageFilename;
    std::string mTrainClassifierObjectImageFilename;
    std::string mTrainClassifierObjectGraphFilename;
    std::string mTrainClassifierTrainingDatasetFilename;
    std::string mTrainClassifierTrainingDatabaseFilename;
    bool mTrainClassifierWithIncremental;
    bool* mTrainClassifierWithFeature;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Pipeline: Cell Membane Segmentation on BCat Filter++++++++
    CSParameterContext *mpSegmentCellMembraneContext;

    std::string mSegmentCellMembraneBCatFilename;
    std::string mSegmentCellMembraneSinFilename;
    std::string mSegmentCellMembraneNucFilename;
    std::string mSegmentCellMembraneCVFilename;
    std::string mSegmentCellMembranePVFilename;
    double mSegmentCellMembraneVoxelSpacing[3];
    CSParameterChoice *mSegmentCellMembraneBCatThresholdMode;
    int mSegmentCellMembraneBCatAdapOtsuThresholdRegionSize[3];
    int mSegmentCellMembraneBCatAdapOtsuThresholdNumSamples;
    int mSegmentCellMembraneBCatManualMinThreshold;
    int mSegmentCellMembraneBCatManualMaxThreshold;
    int mSegmentCellMembraneHoleFilling_RegionSize[3];
    int mSegmentCellMembraneHoleFilling_MajThreshold;
    int mSegmentCellMembraneInvHoleFilling_RegionSize[3];
    int mSegmentCellMembraneInvHoleFilling_MajThreshold;
    int mSegmentCellMembraneClosingRadStrucElem[3];
    int mSegmentCellMembraneOpeningRadStrucElem[3];
    double mSegmentCellMembraneMinCellDiameterThreshold;
    double mSegmentCellMembraneWatershedFloodLevel;
    std::string mSegmentCellMembraneLogFilename;
    std::string mSegmentCellMembraneFilenamePrefix;
    CSParameterChoice *mSegmentCellMembraneSaveMode;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Objects to Matrix Filter++++++++++++++++++++++++
	CSParameterContext *mpObjectsToMatrixContext;

	std::string mObjectsToMatrixImageFilename;
	std::string mObjectsToMatrixColoredImageFilename;
	double mObjectsToMatrixSpacing[3];
	unsigned int mObjectsToMatrixRows;
	unsigned int mObjectsToMatrixCols;
	unsigned int mObjectsToMatrixDim[3];
	bool mObjectsToMatrixFullyConnected;

	bool mObjectsToMatrixWithFilter;
	unsigned int mObjectsToMatrixFilterBySlide;
	bool mObjectsToMatrixWriteObjectIdentificationFile;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Compare Segmentations Filter+++++++++++++++++++++
	CSParameterContext *mpCompareSegmentationsContext;

	int mCompareSegmentationsDatasetID;
	std::string mCompareSegmentationsGoldImageFilename;
	std::string mCompareSegmentationsEvalImageFilename;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: File Format Converter Filter++++++++++++++++++++
	CSParameterContext *mpFileFormatConverterContext;

	std::string mFileFormatConverterFilename;
	CSParameterChoice* mFileFormatConverterInputChoice;
	CSParameterChoice* mFileFormatConverterOutputChoice;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Pipeline: Test Pipeline Dummy+++++++++++++++++++++++++++++
	CSParameterContext *mpTestPipelineContext;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
};

#endif /* IMAGEPROCESSING_H_ */
