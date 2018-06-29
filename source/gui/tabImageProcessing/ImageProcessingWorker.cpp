////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  ImageProcessingWorker.cpp                                     //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2014-05-30 18:28:51                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "ImageProcessingWorker.h"

#include "../../images/filters/analysisFilters/AnalyzeCellsFilter.h"
#include "../../images/filters/analysisFilters/AnalyzeLobuleFilter.h"
#include "../../images/filters/analysisFilters/AnalyzeNucleiFilter.h"
#include "../../images/filters/analysisFilters/AnalyzeStellateCellsFilter.h"
#include "../../images/filters/analysisFilters/CompareSegmentations.h"
#include "../../images/filters/graphFilters/AnalyzeBileNetworkFilter.h"
#include "../../images/pipelines/CropDataSet.h"
#include "../../images/pipelines/ExtractGraph.h"
#include "../../images/pipelines/EstimateCellShape.h"
#include "../../images/pipelines/EstimateLobuleShape.h"
#include "../../images/pipelines/ObjectBasedSegmentation.h"
#include "../../images/pipelines/ObjectsToMatrixFilter.h"
#include "../../images/pipelines/PreprocessDataSet.h"
#include "../../images/pipelines/ProcessHoechstDiffusionFilter.h"
#include "../../images/pipelines/SegmentBileNetworkOnDPPIV20x.h"
#include "../../images/pipelines/SegmentBileNetworkOnDPPIV60x.h"
#include "../../images/pipelines/SegmentCellMembraneOnBCat60x.h"
#include "../../images/pipelines/SegmentNecroticRegionOnDM.h"
#include "../../images/pipelines/SegmentNucleiOnDAPI20x.h"
#include "../../images/pipelines/SegmentNucleiOnDAPI60x.h"
#include "../../images/pipelines/SegmentNucleiOnDAPIWithHough.h"
#include "../../images/pipelines/SegmentAndClassifyStellateCells.h"
#include "../../images/pipelines/SegmentSinusoidalNetworkOnTwoChannels20x.h"
#include "../../images/pipelines/SegmentSinusoidalNetworkOnTwoChannels60x.h"
#include "../../images/pipelines/SegmentVeins.h"
#include "../../images/pipelines/Skeletonization.h"
#include "../../images/pipelines/SuperpixelPreparation.h"
#include "../../images/pipelines/TrainClassifier.h"
#include "../../images/tools/FileFormatConverter.h"

#include "GraphViewer.h"

#include "../../tools/parameters/CSParameterContextTemporary.h"


//TODO: deleting pipeline objects after execution: thorough testing!!
void
ImageProcessingWorker::execute()
{
    if (mpPipelineContextCopy)
    {
        delete mpPipelineContextCopy;
        mpPipelineContextCopy = NULL;
    }

    CSParameterContext * currentContext = mpParameterContextRoot->findContext( ImageProcessing::PipelineName[mPipelineID] );
    mpPipelineContextCopy = new CSParameterContextTemporary( currentContext );
    CSParameterContext * pipelineContext = (CSParameterContext *) mpPipelineContextCopy;

    unsigned int ImageDimension = 3;

    try{
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


        switch ( mPipelineID )
        {
        case ImageProcessing::PipelineAnalyzeCells:
        {
            AnalyzeCellsFilter* anaylzeCells = new AnalyzeCellsFilter();

            anaylzeCells->SetParameterContext(pipelineContext);
            anaylzeCells->Update();
        }
        break;

        case ImageProcessing::PipelineApproximateCellShape:
        {
            EstimateCellShape* estimateCellShape = new EstimateCellShape();

            estimateCellShape->SetParameterContext(pipelineContext);
            estimateCellShape->Update();

            delete estimateCellShape;
        }
        break;

        case ImageProcessing::PipelineSegmentVeins:
        {
            SegmentVeins* segmentVeins = new SegmentVeins();

            segmentVeins->SetParameterContext(pipelineContext);
            segmentVeins->SetCVSeedPoints(mpCaller->m_CVSeedPoints);
            segmentVeins->SetPVSeedPoints(mpCaller->m_PVSeedPoints);
            segmentVeins->Update();
        }
        break;

        case ImageProcessing::PipelineExtractAndAnalyzeGraph:
        {
            ExtractGraph *pipeline = new ExtractGraph();
            pipeline->SetParameterContext(pipelineContext);
            pipeline->Update();

            std::string analysisMode = ( (CSParameterChoice*)(pipelineContext->findParameter("Analysis mode", 0)->dataPointer()) )->currentString();
            bool withAnalysis = (analysisMode.compare("With Analysis")==0);
            AnalyzeBileNetworkFilter *analyzeBile = new AnalyzeBileNetworkFilter();
            if( withAnalysis ) {
                analyzeBile->SetParameterContext(pipelineContext);
                analyzeBile->SetGraphsToAnalyze(pipeline->GetFinalGraphs());
                analyzeBile->SaveGraphsWithArrays(true);
                analyzeBile->SetGraphFilename(pipeline->GetNameOfFirstGraph());
                analyzeBile->Update();
                emit displayResult( (void *) analyzeBile, ImageProcessing::DisplayGraphWithAnalysis, true );
                emit displayResult( (void *) pipeline, ImageProcessing::DisplayGraphWithoutAnalysis, true );
            }
            else
                emit displayResult( (void *) pipeline, ImageProcessing::DisplayGraphWithoutAnalysis, false );
        }
        break;

        case ImageProcessing::PipelineCrop:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                CropDataSet<2>* cropDataSet = new CropDataSet<2>();

                cropDataSet->SetParameterContext(pipelineContext);
                cropDataSet->Update();

                delete cropDataSet;
            }
            else if(ImageDimension == 3) {
                CropDataSet<3>* cropDataSet = new CropDataSet<3>();

                cropDataSet->SetParameterContext(pipelineContext);
                cropDataSet->Update();

                delete cropDataSet;
            }
        }
        break;

        case ImageProcessing::PipelineCLAHE:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                PreprocessDataSet<2>* preprocessDataSet = new PreprocessDataSet<2>();

                preprocessDataSet->SetAlgorithm(PreprocessDataSet<2>::CLAHE);
                preprocessDataSet->SetParameterContext(pipelineContext);
                preprocessDataSet->Update();

                delete preprocessDataSet;
            }
            else if(ImageDimension == 3) {                                   //template arguments have to be const :/
                PreprocessDataSet<3>* preprocessDataSet = new PreprocessDataSet<3>();

                preprocessDataSet->SetAlgorithm(PreprocessDataSet<3>::CLAHE);
                preprocessDataSet->SetParameterContext(pipelineContext);
                preprocessDataSet->Update();

                delete preprocessDataSet;
            }
        }
        break;

        case ImageProcessing::PipelineSegmentSinusoidsAndBile20x:
        {
            SegmentSinusoidalNetworkOnTwoChannels20x* extractSinusNetwork = new SegmentSinusoidalNetworkOnTwoChannels20x();
            SegmentBileNetworkOnDPPIV20x* extractBileNetwork = new SegmentBileNetworkOnDPPIV20x();

            extractSinusNetwork->SetParameterContext(pipelineContext);
            extractSinusNetwork->Update();

            extractBileNetwork->SetParameterContext(pipelineContext);
            extractBileNetwork->Update();
        }
        break;

        case ImageProcessing::PipelineSegmentSinusoidsAndBile60x:
        {
            SegmentSinusoidalNetworkOnTwoChannels60x* extractSinusNetwork = new SegmentSinusoidalNetworkOnTwoChannels60x();
            SegmentBileNetworkOnDPPIV60x* extractBileNetwork = new SegmentBileNetworkOnDPPIV60x();

            extractSinusNetwork->SetParameterContext(pipelineContext);
            extractSinusNetwork->Update();

            extractBileNetwork->SetParameterContext(pipelineContext);
            extractBileNetwork->Update();
        }
        break;

        case ImageProcessing::PipelineSegmentNecroticRegion:
        {
            SegmentNecroticRegionOnDM* extractNecroticRegion = new SegmentNecroticRegionOnDM();

            extractNecroticRegion->SetParameterContext(pipelineContext);
            extractNecroticRegion->Update();
        }
        break;

        case ImageProcessing::PipelineSegmentNuclei20x:
        {
            SegmentNucleiOnDAPI20x* segmentNucleiOnDAPI20x = new SegmentNucleiOnDAPI20x();

            segmentNucleiOnDAPI20x->SetParameterContext(pipelineContext);
            segmentNucleiOnDAPI20x->Update();
        }
        break;

        case ImageProcessing::PipelineSegmentNuclei60x:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                SegmentNucleiOnDAPI60x<2>* segmentNucleiOnDAPI60x = new SegmentNucleiOnDAPI60x<2>();

                segmentNucleiOnDAPI60x->SetParameterContext(pipelineContext);
                segmentNucleiOnDAPI60x->Update();

                delete segmentNucleiOnDAPI60x;
            }
            else if(ImageDimension == 3) {
                SegmentNucleiOnDAPI60x<3>* segmentNucleiOnDAPI60x = new SegmentNucleiOnDAPI60x<3>();

                segmentNucleiOnDAPI60x->SetParameterContext(pipelineContext);
                segmentNucleiOnDAPI60x->Update();

                delete segmentNucleiOnDAPI60x;
            }
        }
        break;

#ifndef CS_TI_QUANT_ONLY

        case ImageProcessing::PipelineAddImages:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                PreprocessDataSet<2>* preprocessDataSet = new PreprocessDataSet<2>();

                preprocessDataSet->SetAlgorithm(PreprocessDataSet<2>::AddImages);
                preprocessDataSet->SetParameterContext(pipelineContext);
                preprocessDataSet->Update();

                delete preprocessDataSet;
            }
            else if(ImageDimension == 3) {                                   //template arguments have to be const :/
                PreprocessDataSet<3>* preprocessDataSet = new PreprocessDataSet<3>();

                preprocessDataSet->SetAlgorithm(PreprocessDataSet<3>::AddImages);
                preprocessDataSet->SetParameterContext(pipelineContext);
                preprocessDataSet->Update();

                delete preprocessDataSet;
            }
        }
        break;

        case ImageProcessing::PipelineAnalyzeNuclei:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                AnalyzeNucleiFilter<2>* analyzeNuclei = new AnalyzeNucleiFilter<2>();

                analyzeNuclei->SetParameterContext(pipelineContext);
                analyzeNuclei->Update();

                delete analyzeNuclei;
            }
            else if(ImageDimension == 3) {                                   //template arguments have to be const :/
                AnalyzeNucleiFilter<3>* analyzeNuclei = new AnalyzeNucleiFilter<3>();

                analyzeNuclei->SetParameterContext(pipelineContext);
                analyzeNuclei->Update();

                delete analyzeNuclei;
            }
        }
        break;

        case ImageProcessing::PipelineAnalyzeLobule:
        {
            AnalyzeLobuleFilter* analyzeLobule = new AnalyzeLobuleFilter();

            analyzeLobule->SetParameterContext(pipelineContext);
            analyzeLobule->Update();

            delete analyzeLobule;
        }
        break;

        case ImageProcessing::PipelineApproximateLobuleShape:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                EstimateLobuleShape<2>* estimateLobuleShape = new EstimateLobuleShape<2>();

                estimateLobuleShape->SetParameterContext(pipelineContext);
                estimateLobuleShape->Update();

                delete estimateLobuleShape;
            }
            else if(ImageDimension == 3) {                                   //template arguments have to be const :/
                EstimateLobuleShape<3>* estimateLobuleShape = new EstimateLobuleShape<3>();

                estimateLobuleShape->SetParameterContext(pipelineContext);
                estimateLobuleShape->Update();

                delete estimateLobuleShape;
            }
        }
        break;

        case ImageProcessing::PipelineSegmentNucleiHoughTransformation:
        {
            SegmentNucleiOnDAPIWithHough* segmentNucleiOnDAPIWithHough = new SegmentNucleiOnDAPIWithHough();

            segmentNucleiOnDAPIWithHough->SetParameterContext(pipelineContext);
            segmentNucleiOnDAPIWithHough->Update();
        }
        break;

        case ImageProcessing::PipelineSegmentStellateCellCytoskeleton:
        {
            SegmentAndClassifyStellateCells* segmentStellateCells = new SegmentAndClassifyStellateCells();

            segmentStellateCells->SetParameterContext(pipelineContext);
            segmentStellateCells->Update();
        }
        break;

        case ImageProcessing::PipelineSuperpixelPreparation:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                SuperpixelPreparation<2>* superpixelPreparationPipeline = new SuperpixelPreparation<2>();

                superpixelPreparationPipeline->SetParameterContext(pipelineContext);
                superpixelPreparationPipeline->Update();

                delete superpixelPreparationPipeline;
            }
            else {
                SuperpixelPreparation<3>* superpixelPreparationPipeline = new SuperpixelPreparation<3>();

                superpixelPreparationPipeline->SetParameterContext(pipelineContext);
                superpixelPreparationPipeline->Update();

                delete superpixelPreparationPipeline;
            }

        }
        break;

        case ImageProcessing::PipelineObjectBasedSegmentation:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                ObjectBasedSegmentation<2>* objectBasedSegPipeline = new ObjectBasedSegmentation<2>();

                objectBasedSegPipeline->SetParameterContext(pipelineContext);
                objectBasedSegPipeline->Update();

                emit displayResult( (void *) objectBasedSegPipeline, ImageProcessing::DisplayGraphWithoutAnalysis, false );
            }
            else {
                ObjectBasedSegmentation<3>* objectBasedSegPipeline = new ObjectBasedSegmentation<3>();

                objectBasedSegPipeline->SetParameterContext(pipelineContext);
                objectBasedSegPipeline->Update();

                emit displayResult( (void *) objectBasedSegPipeline, ImageProcessing::DisplayGraphWithoutAnalysis, false );
            }

        }
        break;

        case ImageProcessing::PipelineTrainClassifier:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                TrainClassifier<2>* trainClassifierPipeline = new TrainClassifier<2>();

                trainClassifierPipeline->SetParameterContext(pipelineContext);
                trainClassifierPipeline->Update();

                delete trainClassifierPipeline;
            }
            else {
                TrainClassifier<3>* trainClassifierPipeline = new TrainClassifier<3>();

                trainClassifierPipeline->SetParameterContext(pipelineContext);
                trainClassifierPipeline->Update();

                delete trainClassifierPipeline;
            }
        }
        break;

        case ImageProcessing::PipelineSegmentCellMembrane:
        {
            SegmentCellMembraneOnBCat60x* segmentCellMembrane = new SegmentCellMembraneOnBCat60x();

            segmentCellMembrane->SetParameterContext(pipelineContext);
            segmentCellMembrane->Update();
        }
        break;

        case ImageProcessing::PipelineObjectsToMatrix:
        {
            ObjectsToMatrixFilter* objectsToMatrix = new ObjectsToMatrixFilter();

            objectsToMatrix->SetParameterContext(pipelineContext);
            objectsToMatrix->Update();
        }
        break;

        case ImageProcessing::PipelineCompareSegmentations:
        {
            CompareSegmentations* compareSegs = new CompareSegmentations();

            compareSegs->SetParameterContext(pipelineContext);
            compareSegs->Update();
        }
        break;

        case ImageProcessing::PipelineAnalyzeStellateCells:
        {
            AnalyzeStellateCellsFilter* anaylzeStellateCells = new AnalyzeStellateCellsFilter();

            anaylzeStellateCells->SetParameterContext(pipelineContext);
            anaylzeStellateCells->Update();
        }
        break;

        case ImageProcessing::PipelineAnalyzeHoechstDiffusion:
        {
            ProcessHoechstDiffusionFilter* processHoechstDiff = new ProcessHoechstDiffusionFilter();

            processHoechstDiff->SetParameterContext(pipelineContext);
            processHoechstDiff->Update();
        }
        break;

        case ImageProcessing::PipelineBackgroundElimination:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                PreprocessDataSet<2>* preprocessDataSet = new PreprocessDataSet<2>();

                preprocessDataSet->SetAlgorithm(PreprocessDataSet<2>::BackgroundElimination);
                preprocessDataSet->SetParameterContext(pipelineContext);
                preprocessDataSet->Update();

                delete preprocessDataSet;
            }
            else if(ImageDimension == 3) {                                   //template arguments have to be const :/
                PreprocessDataSet<3>* preprocessDataSet = new PreprocessDataSet<3>();

                preprocessDataSet->SetAlgorithm(PreprocessDataSet<3>::BackgroundElimination);
                preprocessDataSet->SetParameterContext(pipelineContext);
                preprocessDataSet->Update();

                delete preprocessDataSet;
            }
        }
        break;

        case ImageProcessing::PipelineHConvexImageFilter:
        {
            if(ImageDimension == 2) {                                   //template arguments have to be const :/
                PreprocessDataSet<2>* preprocessDataSet = new PreprocessDataSet<2>();

                preprocessDataSet->SetAlgorithm(PreprocessDataSet<2>::HConvexImageFilter);
                preprocessDataSet->SetParameterContext(pipelineContext);
                preprocessDataSet->Update();

                delete preprocessDataSet;
            }
            else if(ImageDimension == 3) {                                   //template arguments have to be const :/
                PreprocessDataSet<3>* preprocessDataSet = new PreprocessDataSet<3>();

                preprocessDataSet->SetAlgorithm(PreprocessDataSet<3>::HConvexImageFilter);
                preprocessDataSet->SetParameterContext(pipelineContext);
                preprocessDataSet->Update();

                delete preprocessDataSet;
            }
        }
        break;

        case ImageProcessing::PipelineSkeletonization:
        {
            SkeletonizationPipeline* skeletonizeDataSet = new SkeletonizationPipeline();

            skeletonizeDataSet->SetParameterContext(pipelineContext);
            skeletonizeDataSet->Update();
        }
        break;

        case ImageProcessing::PipelineFileFormatConverter:
        {
            FileFormatConverter* pFileFormatConverter = new FileFormatConverter();

            pFileFormatConverter->SetParameterContext(pipelineContext);
            pFileFormatConverter->Update();
        }
        break;

        case ImageProcessing::PipelineTest:
        {
//            FileFormatConverter* converter = new FileFormatConverter();
//            converter->ComposeGraphFromDresdenFormat();
////            converter->ReadMevisDataAndConvertToImageData();

//#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4 || ITK_VERSION_MAJOR > 4)
//            typedef unsigned char                                       PixelCharType;
//            typedef itk::RGBPixel<PixelCharType>                        RGBPixelType;
//
//            typedef itk::Image<PixelCharType, 2>                        CScalarVoImageType;
//            typedef itk::Image<RGBPixelType, 2>                         RGBVoImageType;
//
//            typedef itk::ImageFileReader<CScalarVoImageType>                                                CScalar3DReaderType;
//            typedef itk::ImageFileReader<RGBVoImageType>                                                    RGB3DReaderType;
//            typedef itk::StreamingImageFilter<CScalarVoImageType, CScalarVoImageType>                       StreamingFilter;
//            typedef itk::ImageFileWriter<CScalarVoImageType>                                                CScalar3DWriterType;
//            typedef itk::ImageFileWriter<RGBVoImageType>                                                    RGBVoWriterType;
//
//            CScalar3DReaderType::Pointer reader = CScalar3DReaderType::New();
//            reader->SetFileName("/scratch/friebel/TestData/LoadWholeSlideScans/test.ndpi");
//
//            itk::SCIFIOImageIO::Pointer io = itk::SCIFIOImageIO::New();
//            io->DebugOn();
//
//            reader->SetImageIO( io );
//
//            StreamingFilter::Pointer streamer = StreamingFilter::New();
//            streamer->SetInput( reader->GetOutput() );
//            streamer->SetNumberOfStreamDivisions( 300 );
//
//            itk::SCIFIOImageIO::Pointer ioOut = itk::SCIFIOImageIO::New();
//            ioOut->DebugOn();
//
//            CScalar3DWriterType::Pointer writer = CScalar3DWriterType::New();
//            writer->SetFileName("/scratch/friebel/TestData/LoadWholeSlideScans/test.tiff");
//            writer->SetInput(streamer->GetOutput());
//            writer->SetNumberOfStreamDivisions( 300 );
////            writer->SetImageIO( ioOut );
//#if (ITK_VERSION_MAJOR >= 4)
//            writer->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//
////            reader->UpdateOutputInformation();
////            std::cout << "Get series counts " << io->GetSeriesCount() << std::endl;
////            io->SetSeries(3);
////            reader->Modified();
//
//            try
//            {
//                writer->Update();
//            }
//            catch (itk::ExceptionObject &e)
//            {
//                std::cerr << e << std::endl;
//            }
//#endif
        }
        break;

        case ImageProcessing::PipelineVisualize3D:
        {
            //TODO: encapsulate this in volume viewer widget class, similar to graph viewer

//            typedef unsigned char                                       PixelCharType;
//            typedef itk::Image<PixelCharType, 3>                        CScalar3DImageType;
//
//            typedef itk::ImageFileReader<CScalar3DImageType>            CScalar3DReaderType;
//            typedef itk::ImageToVTKImageFilter<CScalar3DImageType>      ITKVTKScalar3DConnectorType;
//
//            std::string bileImagePath, sinusoidImagePath, nucleiImagePath, graphFilePath, cellShapeImagePath;
//            bileImagePath = *(std::string*)(pipelineContext->findParameter("Input Bile Image:", 0)->dataPointer());
//            sinusoidImagePath = *(std::string*)(pipelineContext->findParameter("Input Sinusoid Image:", 0)->dataPointer());
//            nucleiImagePath = *(std::string*)(pipelineContext->findParameter("Input Nuclei Image:", 0)->dataPointer());
//            graphFilePath = *(std::string*)(pipelineContext->findParameter("Input Graph Input:", 0)->dataPointer());
//            cellShapeImagePath = *(std::string*)(pipelineContext->findParameter("Input Cell Shape Image:", 0)->dataPointer());
//
//            vtkSmartPointer<vtkPlane> plane1, plane2, plane3, plane4, plane5, plane6;
//
//            plane1 = vtkSmartPointer<vtkPlane>::New();
//            plane1->SetOrigin(0.0, 150.0, 0.0);
//            plane1->SetNormal(150.0, 0.0, 0.0);
//
//            plane2 = vtkSmartPointer<vtkPlane>::New();
//            plane2->SetOrigin(150.0, 0.0, 0.0);
//            plane2->SetNormal(0.0, 150.0, 0.0);
//
//            plane3 = vtkSmartPointer<vtkPlane>::New();
//            plane3->SetOrigin(150.0, 300.0, 0.0);
//            plane3->SetNormal(0.0, -150.0, 0.0);
//
//            plane4 = vtkSmartPointer<vtkPlane>::New();
//            plane4->SetOrigin(300.0, 150.0, 0.0);
//            plane4->SetNormal(-150.0, 0.0, 0.0);
//
//            plane5 = vtkSmartPointer<vtkPlane>::New();
//            plane5->SetOrigin(150.0, 150.0, 0.0);
//            plane5->SetNormal(0.0, 0.0, 15.0);
//
//            plane6 = vtkSmartPointer<vtkPlane>::New();
//            plane6->SetOrigin(150.0, 150.0, 30.0);
//            plane6->SetNormal(0.0, 0.0, -15.0);
//
//            //CREATE-DISPLAY-AND-GUI-TOOLS-----------------------------------------------------------------------------
//            displays["Input image"] = new QCSVTKDisplay();
//            displays["Input image"]->setGeometry(QRect(0, 0, 800, 600));
//            displays["Input image"]->setWindowTitle("Test Page Display");
//
//            QWidget* volVizSidebar = new QWidget(displays["Input image"]);
//            volVizSidebar->setObjectName(QString::fromUtf8("volVizSidebar"));
//            volVizSidebar->setFixedWidth(200);
//            volVizSidebar->setGeometry(QRect(600, 0, 200, 200));
//            volVizSidebar->setSizePolicy( QSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding) );
//
//            displays["Input image"]->addControlWidget(volVizSidebar);
//
//            QGridLayout* controlLayout = new QGridLayout(volVizSidebar);
//            controlLayout->setObjectName(QString::fromUtf8("controlLayout"));
//            volVizSidebar->setLayout(controlLayout);
//
//            QPalette *p = new QPalette;
//            p->setColor(QPalette::WindowText, QColor(0, 0, 0));
//
//            slider1 = new QSlider(Qt::Vertical, volVizSidebar);
//            slider1->setRange(0, 1024);
//            slider1->setValue(0);
//            controlLayout->addWidget(slider1, 0, 0, 2, 1, Qt::AlignCenter);
//
//            connect( slider1, SIGNAL(valueChanged(int)), this, SLOT(moveObjectsWithSlider1(int)) );
//
//            QLabel * slider1Label = new QLabel("Move in X direction", volVizSidebar);
//            controlLayout->addWidget(slider1Label, 0, 1, 1, 1, Qt::AlignLeft );
//
//            QLCDNumber* coordAxis1 = new QLCDNumber(volVizSidebar);
//            coordAxis1->setPalette(*p);
//            coordAxis1->setSegmentStyle(QLCDNumber::Filled);
//            controlLayout->addWidget(coordAxis1, 1, 1, 1, 1, Qt::AlignCenter);
//
//            connect( slider1, SIGNAL(valueChanged(int)), coordAxis1, SLOT(display(int)) );
//
//
//            slider2 = new QSlider(Qt::Vertical, volVizSidebar);
//            slider2->setRange(0, 1024);
//            slider2->setValue(0);
//            controlLayout->addWidget(slider2, 2, 0, 2, 1, Qt::AlignCenter);
//
//            connect( slider2, SIGNAL(valueChanged(int)), this, SLOT(moveObjectsWithSlider2(int)) );
//
//            QLabel * slider2Label = new QLabel("Move in Y direction", volVizSidebar);
//            controlLayout->addWidget(slider2Label, 2, 1, 1, 1, Qt::AlignLeft );
//
//            QLCDNumber* coordAxis2 = new QLCDNumber(volVizSidebar);
//            coordAxis2->setPalette(*p);
//            coordAxis2->setSegmentStyle(QLCDNumber::Filled);
//            controlLayout->addWidget(coordAxis2, 3, 1, 1, 1, Qt::AlignCenter);
//
//            connect( slider2, SIGNAL(valueChanged(int)), coordAxis2, SLOT(display(int)) );
//
//
//            slider3 = new QSlider(Qt::Vertical, volVizSidebar);
//            slider3->setRange(0, 150);
//            slider3->setValue(0);
//            controlLayout->addWidget(slider3, 4, 0, 2, 1, Qt::AlignCenter);
//
//            connect( slider3, SIGNAL(valueChanged(int)), this, SLOT(moveObjectsWithSlider3(int)) );
//
//            QLabel * slider3Label = new QLabel("Move in Z direction", volVizSidebar);
//            controlLayout->addWidget(slider3Label, 4, 1, 1, 1, Qt::AlignLeft );
//
//            QLCDNumber* coordAxis3 = new QLCDNumber(volVizSidebar);
//            coordAxis3->setPalette(*p);
//            coordAxis3->setSegmentStyle(QLCDNumber::Filled);
//            controlLayout->addWidget(coordAxis3, 5, 1, 1, 1, Qt::AlignCenter);
//
//            connect( slider3, SIGNAL(valueChanged(int)), coordAxis3, SLOT(display(int)) );
//
//
//            QPushButton * visibilityBileButton = new QPushButton(QString("Hide bile ducts"), volVizSidebar);
//            visibilityBileButton->setCheckable(true);
//            controlLayout->addWidget( visibilityBileButton, 6, 0, 1, 2 );
//
//            connect( visibilityBileButton, SIGNAL(toggled(bool)), this, SLOT(toggleVisibilityBile(bool)) );
//
//            QPushButton * visibilitySinusoidsButton = new QPushButton(QString("Hide sinusoids"), volVizSidebar);
//            visibilitySinusoidsButton->setCheckable(true);
//            controlLayout->addWidget( visibilitySinusoidsButton, 7, 0, 1, 2 );
//
//            connect( visibilitySinusoidsButton, SIGNAL(toggled(bool)), this, SLOT(toggleVisibilitySinusoids(bool)) );
//
//            QPushButton * visibilityNucleiButton = new QPushButton(QString("Hide nuclei"), volVizSidebar);
//            visibilityNucleiButton->setCheckable(true);
//            controlLayout->addWidget( visibilityNucleiButton, 8, 0, 1, 2 );
//
//            connect( visibilityNucleiButton, SIGNAL(toggled(bool)), this, SLOT(toggleVisibilityNuclei(bool)) );
//
//            QPushButton * visibilityGraphButton = new QPushButton(QString("Hide graph"), volVizSidebar);
//            visibilityGraphButton->setCheckable(true);
//            controlLayout->addWidget( visibilityGraphButton, 9, 0, 1, 2 );
//
//            connect( visibilityGraphButton, SIGNAL(toggled(bool)), this, SLOT(toggleVisibilityGraph(bool)) );
//
//            QPushButton * visibilityCellShapeButton = new QPushButton(QString("Hide cell shape"), volVizSidebar);
//            visibilityCellShapeButton->setCheckable(true);
//            controlLayout->addWidget( visibilityCellShapeButton, 10, 0, 1, 2 );
//
//            connect( visibilityCellShapeButton, SIGNAL(toggled(bool)), this, SLOT(toggleVisibilityCellShape(bool)) );
//
//            QPushButton * exportAsPOVButton = new QPushButton(QString("Export as .pov"), volVizSidebar);
//            exportAsPOVButton->setCheckable(true);
//            controlLayout->addWidget( exportAsPOVButton, 11, 0, 1, 2 );
//
//            connect( exportAsPOVButton, SIGNAL(clicked(bool)), this, SLOT(exportAsPOVButtonClicked()) );
//
//
//            vtkRenderWindow * renWin = displays["Input image"]->getRenderWindow();
//
//            vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
//            ren1->SetBackground(1,1,1);
//
//            renWin->AddRenderer(ren1);
//            //---------------------------------------------------------------------------------------------------------
//
//            CScalar3DReaderType::Pointer reader = CScalar3DReaderType::New();
//            reader->SetFileName(bileImagePath);
//            reader->ReleaseDataFlagOn();
//
//#if (ITK_VERSION_MAJOR >= 4)
//            reader->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//
//            ITKVTKScalar3DConnectorType::Pointer itkvtkConnectorInputImage = ITKVTKScalar3DConnectorType::New();
//            itkvtkConnectorInputImage->SetInput(reader->GetOutput());
//            itkvtkConnectorInputImage->ReleaseDataFlagOn();
//
//            smartPointerContainer.push_back((void *) itkvtkConnectorInputImage);
//
//            vtkSmartPointer<vtkContourFilter> contour = vtkSmartPointer<vtkContourFilter>::New();
//            contour->SetValue(0, 150);
//            contour->SetInput(itkvtkConnectorInputImage->GetOutput());
//            contour->ReleaseDataFlagOn();
//            contour->Update();
//
//            //          vtkSmartPointer<vtkTriangleFilter> tris = vtkTriangleFilter::New();
//            //          tris->SetInputConnection(contour->GetOutputPort());
//            //          tris->GetOutput()->ReleaseDataFlagOn();
//            //
//            //          vtkSmartPointer<vtkStripper> strip = vtkStripper::New();
//            //          strip->SetInputConnection(tris->GetOutputPort());
//            ////            strip->GetOutput()->ReleaseDataFlagOn();
//
//            //          std::cout << "Before: There are " << contour->GetOutput()->GetNumberOfPolys() << " polygons." << std::endl;
//            //
//            //          vtkSmartPointer<vtkQuadricClustering> decimate = vtkSmartPointer<vtkQuadricClustering>::New();
//            //          decimate->SetInputConnection(contour->GetOutputPort());
//            ////            decimate->AutoAdjustNumberOfDivisionsOff();
//            ////            decimate->SetNumberOfDivisions(300, 300, 200);
//            ////            decimate->Update();
//            //
//            //          std::cout << "After: There are " << decimate->GetOutput()->GetNumberOfPolys() << " polygons." << std::endl;
//
//            vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//            mapper->SetInputConnection(contour->GetOutputPort());
//            mapper->ScalarVisibilityOff();
//            //          mapper->ImmediateModeRenderingOn();
//            mapper->AddClippingPlane(plane1);
//            mapper->AddClippingPlane(plane2);
//            mapper->AddClippingPlane(plane3);
//            mapper->AddClippingPlane(plane4);
//            mapper->AddClippingPlane(plane5);
//            mapper->AddClippingPlane(plane6);
//
//            // Create an actor for the polygonal data.
//            vtkSmartPointer<vtkLODActor> actor = vtkSmartPointer<vtkLODActor>::New();
//            actor->SetMapper(mapper);
//            actor->SetNumberOfCloudPoints(500);
//            actor->GetProperty()->LightingOn();
//            actor->GetProperty()->SetColor(0, 1, 0);
//            actor->GetProperty()->SetOpacity(0.6);
//
//            actors["Bile"] = (vtkProp3D::SafeDownCast(actor));
//
//            ren1->AddActor(actor);
//
//            //          // Create transfer mapping scalar value to opacity.
//            //          vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
//            //          opacityTransferFunction->AddPoint(80.0,  0.0);
//            //          opacityTransferFunction->AddPoint(120.0, 0.2);
//            //          opacityTransferFunction->AddPoint(255.0, 0.2);
//            //
//            //          // Create transfer mapping scalar value to color.
//            //          vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
//            //          colorTransferFunction->AddRGBPoint(80.0,  0.0, 0.0, 0.0);
//            //          colorTransferFunction->AddRGBPoint(120.0, 0.0, 0.0, 1.0);
//            //          colorTransferFunction->AddRGBPoint(160.0, 1.0, 0.0, 0.0);
//            //          colorTransferFunction->AddRGBPoint(200.0, 0.0, 1.0, 0.0);
//            //          colorTransferFunction->AddRGBPoint(255.0, 0.0, 1.0, 1.0);
//            //
//            //          // The property describes how the data will look.
//            //          vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
//            //          volumeProperty->SetColor(colorTransferFunction);
//            //          volumeProperty->SetScalarOpacity(opacityTransferFunction);
//            //          volumeProperty->ShadeOff();
//            //          volumeProperty->SetInterpolationTypeToLinear();
//            //
//            //          // The mapper that renders the volume data.
//            //          vtkSmartPointer<vtkSmartVolumeMapper> volumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
//            //          volumeMapper->SetInput(itkvtkConnectorInputImage->GetOutput());
//            //          volumeMapper->AddClippingPlane(plane1);
//            //          volumeMapper->AddClippingPlane(plane2);
//            //          volumeMapper->AddClippingPlane(plane3);
//            //          volumeMapper->AddClippingPlane(plane4);
//            //          volumeMapper->AddClippingPlane(plane5);
//            //          volumeMapper->AddClippingPlane(plane6);
//            //
//            //          vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
//            //          volume->SetMapper(volumeMapper);
//            //          volume->SetProperty(volumeProperty);
//            //
//            //          actors["Bile"] = (vtkProp3D::SafeDownCast(volume));
//            //
//            //          ren1->AddViewProp(volume);
//
//
//            CScalar3DReaderType::Pointer reader2 = CScalar3DReaderType::New();
//            reader2->SetFileName(sinusoidImagePath);
//            reader2->ReleaseDataFlagOn();
//
//#if (ITK_VERSION_MAJOR >= 4)
//            reader2->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//
//            ITKVTKScalar3DConnectorType::Pointer itkvtkConnectorInputImage2 = ITKVTKScalar3DConnectorType::New();
//            itkvtkConnectorInputImage2->SetInput(reader2->GetOutput());
//            itkvtkConnectorInputImage2->ReleaseDataFlagOn();
//
//            smartPointerContainer.push_back((void *) itkvtkConnectorInputImage2);
//
//            vtkSmartPointer<vtkContourFilter> contour2 = vtkSmartPointer<vtkContourFilter>::New();
//            contour2->SetValue(0, 150);
//            contour2->SetInput(itkvtkConnectorInputImage2->GetOutput());
//
//            vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//            mapper2->SetInputConnection(contour2->GetOutputPort());
//            mapper2->ScalarVisibilityOff();
//            //              mapper2->ImmediateModeRenderingOn();
//            mapper2->AddClippingPlane(plane1);
//            mapper2->AddClippingPlane(plane2);
//            mapper2->AddClippingPlane(plane3);
//            mapper2->AddClippingPlane(plane4);
//            mapper2->AddClippingPlane(plane5);
//            mapper2->AddClippingPlane(plane6);
//
//            vtkSmartPointer<vtkLODActor> actor2 = vtkSmartPointer<vtkLODActor>::New();
//            actor2->SetMapper(mapper2);
//            actor2->SetNumberOfCloudPoints(500);
//            actor2->GetProperty()->LightingOn();
//            actor2->GetProperty()->SetColor(1, 0, 0);
//            actor2->GetProperty()->SetOpacity(0.6);
//
//            actors["Sinusoids"] = (vtkProp3D::SafeDownCast(actor2));
//
//            ren1->AddActor(actor2);
//
//
//            CScalar3DReaderType::Pointer reader3 = CScalar3DReaderType::New();
//            reader3->SetFileName(nucleiImagePath);
//            reader3->ReleaseDataFlagOn();
//
//#if (ITK_VERSION_MAJOR >= 4)
//            reader3->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//
//            ITKVTKScalar3DConnectorType::Pointer itkvtkConnectorInputImage3 = ITKVTKScalar3DConnectorType::New();
//            itkvtkConnectorInputImage3->SetInput(reader3->GetOutput());
//            itkvtkConnectorInputImage3->ReleaseDataFlagOn();
//
//            smartPointerContainer.push_back((void *) itkvtkConnectorInputImage3);
//
//            vtkSmartPointer<vtkContourFilter> contour3 = vtkSmartPointer<vtkContourFilter>::New();
//            contour3->SetValue(0, 150);
//            contour3->SetInput(itkvtkConnectorInputImage3->GetOutput());
//
//            vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
//            mapper3->SetInputConnection(contour3->GetOutputPort());
//            mapper3->ScalarVisibilityOff();
//            //              mapper3->ImmediateModeRenderingOn();
//            mapper3->AddClippingPlane(plane1);
//            mapper3->AddClippingPlane(plane2);
//            mapper3->AddClippingPlane(plane3);
//            mapper3->AddClippingPlane(plane4);
//            mapper3->AddClippingPlane(plane5);
//            mapper3->AddClippingPlane(plane6);
//
//            vtkSmartPointer<vtkLODActor> actor3 = vtkSmartPointer<vtkLODActor>::New();
//            actor3->SetMapper(mapper3);
//            actor3->SetNumberOfCloudPoints(500);
//            actor3->GetProperty()->LightingOn();
//            actor3->GetProperty()->SetColor(0, 0, 1);
//            actor3->GetProperty()->SetOpacity(1.0);
//
//            actors["Nuclei"] = (vtkProp3D::SafeDownCast(actor3));
//
//            ren1->AddActor(actor3);
//
//
//            CScalar3DReaderType::Pointer reader4 = CScalar3DReaderType::New();
//            reader4->SetFileName(cellShapeImagePath);
//            reader4->ReleaseDataFlagOn();
//
//#if (ITK_VERSION_MAJOR >= 4)
//            reader4->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//
//            ITKVTKScalar3DConnectorType::Pointer itkvtkConnectorInputImage4 = ITKVTKScalar3DConnectorType::New();
//            itkvtkConnectorInputImage4->SetInput(reader4->GetOutput());
//            itkvtkConnectorInputImage4->ReleaseDataFlagOn();
//
//            smartPointerContainer.push_back((void *) itkvtkConnectorInputImage3);
//
//            vtkSmartPointer<vtkContourFilter> contour4 = vtkSmartPointer<vtkContourFilter>::New();
//            contour4->SetValue(0, 150);
//            contour4->SetInput(itkvtkConnectorInputImage4->GetOutput());
//
//            vtkSmartPointer<vtkPolyDataMapper> mapper4 = vtkSmartPointer<vtkPolyDataMapper>::New();
//            mapper4->SetInputConnection(contour4->GetOutputPort());
//            mapper4->ScalarVisibilityOff();
//            //              mapper3->ImmediateModeRenderingOn();
//            mapper4->AddClippingPlane(plane1);
//            mapper4->AddClippingPlane(plane2);
//            mapper4->AddClippingPlane(plane3);
//            mapper4->AddClippingPlane(plane4);
//            mapper4->AddClippingPlane(plane5);
//            mapper4->AddClippingPlane(plane6);
//
//            vtkSmartPointer<vtkLODActor> actor4 = vtkSmartPointer<vtkLODActor>::New();
//            actor4->SetMapper(mapper4);
//            actor4->SetNumberOfCloudPoints(500);
//            actor4->GetProperty()->LightingOn();
//            actor4->GetProperty()->SetColor(0.6, 0.6, 0.6);
//            actor4->GetProperty()->SetOpacity(1.0);
//
//            actors["CellShape"] = (vtkProp3D::SafeDownCast(actor4));
//
//            ren1->AddActor(actor4);
//
//
//            std::string graphFilename = graphFilePath;
//            if(!graphFilename.empty()) {
//                vtkSmartPointer<vtkGraphReader> readerGraph = vtkSmartPointer<vtkGraphReader>::New();
//                readerGraph->ReleaseDataFlagOn();
//
//                vtkSmartPointer<vtkGraphToPolyData> graphToPoly = vtkSmartPointer<vtkGraphToPolyData>::New();
//                graphToPoly->SetInputConnection(readerGraph->GetOutputPort());
//                graphToPoly->ReleaseDataFlagOn();
//
//                vtkSmartPointer<vtkGlyphSource2D> glyphSource = vtkGlyphSource2D::New();
//                glyphSource->SetGlyphTypeToVertex();
//
//                vtkSmartPointer<vtkGlyph3D> vertGlyph = vtkGlyph3D::New();
//                vertGlyph->SetInputConnection(0, graphToPoly->GetOutputPort());
//                vertGlyph->SetInputConnection(1, glyphSource->GetOutputPort());
//                vertGlyph->ReleaseDataFlagOn();
//
//                vtkSmartPointer<vtkMultiBlockDataSet> graphEdges = vtkSmartPointer<vtkMultiBlockDataSet>::New();
//                vtkSmartPointer<vtkMultiBlockDataSet> graphVertices = vtkSmartPointer<vtkMultiBlockDataSet>::New();
//
//
//                size_t lastSlash, lastDot;
//
//                lastSlash = graphFilename.find_last_of("/\\");
//                std::string graphPath = graphFilename.substr(0, lastSlash+1);
//
//                lastDot = graphFilename.find_last_of(".");
//                std::string graphFileExtension = graphFilename.substr(lastDot);
//
//                int nameLength = graphFilename.length() - graphPath.length() - graphFileExtension.length();
//                graphFilename = graphFilename.substr(lastSlash+1, nameLength);
//                graphFilename = graphFilename.substr(0, graphFilename.length()-1);
//
//                int i=0;
//                while(true) {
//                    std::stringstream name;
//                    name << graphPath << graphFilename << i << ".txt";
//                    std::cout << "load graph file " << name.str() << std::endl;
//
//                    readerGraph->SetFileName(name.str().c_str());
//                    int s = readerGraph->OpenVTKFile();
//
//                    if(s != 0) {
//                        graphToPoly->Update();
//                        vertGlyph->Update();
//
//                        vtkPolyData* pdE = vtkPolyData::New();
//                        vtkPolyData* pdG = vtkPolyData::New();
//
//                        pdE->ShallowCopy(graphToPoly->GetOutput());
//                        pdG->ShallowCopy(vertGlyph->GetOutput());
//
//                        graphVertices->SetBlock(i, pdG);
//                        graphEdges->SetBlock(i, pdE);
//                        pdG->Delete();
//                        pdE->Delete();
//
//                        readerGraph->CloseVTKFile();
//                    }
//                    else
//                        break;
//
//                    i++;
//                }
//
//                vtkSmartPointer<vtkCompositeDataGeometryFilter> edgesFilter = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();
//                edgesFilter->SetInput(graphEdges);
//                edgesFilter->Update();
//
//                vtkSmartPointer<vtkCompositeDataGeometryFilter> vertexFilter = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();
//                vertexFilter->SetInput(graphVertices);
//                vertexFilter->Update();
//
//                vtkSmartPointer<vtkPolyDataMapper> edgeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//                edgeMapper->SetInputConnection(0, edgesFilter->GetOutputPort(0));
//                edgeMapper->ImmediateModeRenderingOn();
//                edgeMapper->AddClippingPlane(plane1);
//                edgeMapper->AddClippingPlane(plane2);
//                edgeMapper->AddClippingPlane(plane3);
//                edgeMapper->AddClippingPlane(plane4);
//                edgeMapper->AddClippingPlane(plane5);
//                edgeMapper->AddClippingPlane(plane6);
//
//                vtkActor* edgeActor = vtkActor::New();
//                edgeActor->SetMapper(edgeMapper);
//                edgeActor->GetProperty()->SetColor(0.0, 0.0, 0.0);
//                edgeActor->GetProperty()->SetLineWidth(5.0);
//
//                vtkPolyDataMapper* vertMapper = vtkPolyDataMapper::New();
//                vertMapper->SetInputConnection(0, vertexFilter->GetOutputPort(0));
//                vertMapper->ImmediateModeRenderingOn();
//                vertMapper->AddClippingPlane(plane1);
//                vertMapper->AddClippingPlane(plane2);
//                vertMapper->AddClippingPlane(plane3);
//                vertMapper->AddClippingPlane(plane4);
//                vertMapper->AddClippingPlane(plane5);
//                vertMapper->AddClippingPlane(plane6);
//
//                vtkActor* vertActor = vtkActor::New();
//                vertActor->SetMapper(vertMapper);
//                vertActor->GetProperty()->SetPointSize(5);
//
//                actors["BileGraphEdges"] = edgeActor;
//                actors["BileGraphNodes"] = vertActor;
//
//                ren1->AddActor(edgeActor);
//                ren1->AddActor(vertActor);
//            }
//
//            ren1->ResetCamera();
//
//            //          vtkBoxWidget *boxWidget = vtkBoxWidget::New();
//            //          boxWidget->SetInteractor(ren1->GetRenderWindow()->GetInteractor());
//            //          boxWidget->SetPlaceFactor(1.25);
//            //          boxWidget->SetProp3D(actor);
//            //          boxWidget->PlaceWidget();
//            //          boxWidget->On();
//
//            //          // 3D texture mode. For coverage.
//            //          volumeMapper->SetRequestedRenderModeToRayCastAndTexture();
//
//            // Software mode, for coverage. It also makes sure we will get the same
//            // regression image on all platforms.
//            //          volumeMapper->SetRequestedRenderModeToRayCast();
//
//            //          parentWidget->resize(800,600);
//
//            //---------------------------------------------------------------------------------------------------------
//            ren1->GetActiveCamera()->SetPosition(150, 150, -800);
//            ren1->GetActiveCamera()->SetFocalPoint(150, 150, 0);
//            ren1->GetActiveCamera()->SetRoll(180);
//
//
//            vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();                                          //TODO: should go into QCSVTKDisplay
//            axes->SetScale(0.5);
//            axes->SetXAxisLabelText("X");
//            axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0,0,0);
//            axes->SetYAxisLabelText("Y");
//            axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0,0,0);
//            axes->SetZAxisLabelText("Z");
//            axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0,0,0);
//
//            vtkOrientationMarkerWidget* widget = vtkOrientationMarkerWidget::New();
//            widget->SetOutlineColor(0.9300, 0.5700, 0.1300);
//            widget->SetOrientationMarker(axes);
//            widget->SetInteractor(displays["Input image"]->getInteractor());
//            widget->SetViewport(0.0, 0.0, 0.2, 0.2);
//            widget->SetEnabled(1);
//            widget->InteractiveOn();
//            //---------------------------------------------------------------------------------------------------------
//
//
//            displays["Input image"]->show();
        }
        break;
#endif // CS_TI_QUANT_ONLY
        }
    }
    catch( std::string& err )
    {
        std::stringstream s;
        s << "Sorry. An error occured.\n" << err << std::endl;
        QString messErr = QString::fromStdString(s.str());

        emit aborted( messErr );
    }
    catch( itk::ExceptionObject & err )
	{
		std::stringstream s;
		s << "Sorry. An error occured.\n" << err << std::endl;
		QString messErr = QString::fromStdString(s.str());

        emit aborted( messErr );
	}
    catch( std::exception& err )
	{
		std::stringstream s;
		s << "Sorry. An error occured.\n" << err.what() << std::endl;
		QString messErr = QString::fromStdString(s.str());

        emit aborted( messErr );
	}

    emit finished();
}
