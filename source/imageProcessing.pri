
DEFINES += CS_BUILD_IMAGEPROCESSING

include(vtk.pri)
include(itk.pri)


DEPENDPATH += images \
              images/tools \
              images/pipelines \
              images/filters/convertFilters \
              images/filters/graphFilters \
              images/filters/imageFilters

INCLUDEPATH += images \
               images/filters/imageFilters \
               images/pipelines \
               images/filters/convertFilters \
               images/filters/graphFilters \
               images/tools

HEADERS += images/CSImage.h \
           images/pipelines/BasePipeline.h \
           images/pipelines/BasePipeline.tpp \
           images/pipelines/CropDataSet.h \
           images/pipelines/CropDataSet.tpp \
           images/pipelines/EstimateCellShape.h \
           images/pipelines/EstimateLobuleShape.h \
           images/pipelines/EstimateLobuleShape.tpp \
           images/pipelines/ExtractGraph.h \
           images/pipelines/LabelMapGraphBasePipeline.h \
           images/pipelines/LabelMapGraphBasePipeline.tpp \
           images/pipelines/ObjectBasedSegmentation.h \
           images/pipelines/ObjectBasedSegmentation.tpp \
           images/pipelines/ObjectsToMatrixFilter.h \
           images/pipelines/PreprocessDataSet.h \
           images/pipelines/PreprocessDataSet.tpp \
           images/pipelines/ResampleNetwork.h \
           images/pipelines/ProcessHoechstDiffusionFilter.h \
           images/pipelines/SegmentAndClassifyStellateCells.h \
           images/pipelines/SegmentBileNetworkOnDPPIV20x.h \
           images/pipelines/SegmentBileNetworkOnDPPIV60x.h \
           images/pipelines/SegmentCellMembraneOnBCat60x.h \
           images/pipelines/SegmentNecroticRegionOnDM.h \
           images/pipelines/SegmentNucleiOnDAPI20x.h \
           images/pipelines/SegmentNucleiOnDAPI60x.h \
           images/pipelines/SegmentNucleiOnDAPI60x.tpp \
           images/pipelines/SegmentNucleiOnDAPIWithHough.h \
           images/pipelines/SegmentSinusoidalNetworkOnTwoChannels20x.h \
           images/pipelines/SegmentSinusoidalNetworkOnTwoChannels60x.h \
           images/pipelines/SegmentVeins.h \
           images/pipelines/Skeletonization.h \
           images/pipelines/SuperpixelPreparation.h \
           images/pipelines/SuperpixelPreparation.tpp \
           images/pipelines/TrainClassifier.h \
           images/pipelines/TrainClassifier.tpp \
           images/tools/FeatureLabelObject.h \
           images/tools/FeatureLabelObject.tpp \
           images/tools/FileFormatConverter.h \
           images/tools/GraphAnnotationHelper.h \
           images/tools/ImageAnalysisSummaryFileIO.h \
           images/tools/LabelMapGraph.h \
           images/tools/LabelMapGraph.tpp \
           images/tools/LineConstIteratorDerivative.h \
           images/tools/LineConstIteratorDerivative.tpp \
           images/filters/analysisFilters/AnalyzeCellsFilter.h \
           images/filters/analysisFilters/AnalyzeLobuleFilter.h \
           images/filters/analysisFilters/AnalyzeNucleiFilter.h \
           images/filters/analysisFilters/AnalyzeNucleiFilter.tpp \
           images/filters/analysisFilters/AnalyzeStellateCellsFilter.h \
           images/filters/analysisFilters/CompareSegmentations.h \
           images/filters/analysisFilters/Lobule2DHandlingFilter.h \
           images/filters/classifiers/GraphBaseClassifier.h \
           images/filters/classifiers/GraphBaseClassifier.tpp \
           images/filters/classifiers/GraphEdgeBaseClassifier.h \
           images/filters/classifiers/GraphEdgeBaseClassifier.tpp \
           images/filters/classifiers/GraphEdgeKNNClassifier.h \
           images/filters/classifiers/GraphEdgeKNNClassifier.tpp \
           images/filters/classifiers/GraphEdgePythonSVMClassifier.h \
           images/filters/classifiers/GraphEdgePythonSVMClassifier.tpp \
           images/filters/classifiers/GraphEdgeTestClassifier.h \
           images/filters/classifiers/GraphEdgeTestClassifier.tpp \
           images/filters/classifiers/GraphVertexBaseClassifier.h \
           images/filters/classifiers/GraphVertexBaseClassifier.tpp \
           images/filters/classifiers/GraphVertexPythonSVMClassifier.h \
           images/filters/classifiers/GraphVertexPythonSVMClassifier.tpp \
           images/filters/convertFilters/GraphToImageFilter.h \
           images/filters/convertFilters/LabelImageToLabelMapGraphFilter.h \
           images/filters/convertFilters/LabelImageToLabelMapGraphFilter.tpp \
           images/filters/convertFilters/LabelImageToGraphFilter.h \
           images/filters/convertFilters/LabelImageToGraphFilter.tpp \
           images/filters/convertFilters/LabelMapGraphToLabelImageFilter.h \
           images/filters/convertFilters/LabelMapGraphToLabelImageFilter.tpp \
           images/filters/convertFilters/SkeletonImageToGraphFilter.h \
           images/filters/graphFilters/AnalyzeBileNetworkFilter.h \
           images/filters/graphFilters/AnalyzeGraphFilter.h \
           images/filters/graphFilters/CollapseIntersectionNodesFilter.h \
           images/filters/graphFilters/GeometricalThresholdUndirectedGraphFilter.h \
           images/filters/graphFilters/GraphHandlingHelper.h \
           images/filters/graphFilters/ResampleUndirectedGraphFilter.h \
           images/filters/graphFilters/TopologicalPruneGraphFilter.h \
           images/filters/imageFilters/CavityFillingImageFilter.h \
           images/filters/imageFilters/CavityFillingImageFilter.tpp \
           images/filters/imageFilters/CSImageFilter.h \
           images/filters/imageFilters/FastCLAHEImageFilter.h \
           images/filters/imageFilters/FastCLAHEImageFilter.tpp \
           images/filters/imageFilters/IterativeCavityFillingImageFilter.h \
           images/filters/imageFilters/IterativeCavityFillingImageFilter.tpp \
           images/filters/imageFilters/itkAdaptiveOtsuThresholdingImageFilter.h \
           images/filters/imageFilters/itkAdaptiveOtsuThresholdingImageFilter.tpp \
           images/filters/imageFilters/itkBinaryThinningImageFilter3D.h \
           images/filters/imageFilters/itkBinaryThinningImageFilter3D.tpp \
           images/filters/imageFilters/itkHoughTransformRadialVotingImageFilter.h \
           images/filters/imageFilters/itkHoughTransformRadialVotingImageFilter.tpp \
           images/filters/imageFilters/itkMaskedSpatialObjectToImageFilter.h \
           images/filters/imageFilters/itkMaskedSpatialObjectToImageFilter.tpp \
           images/filters/imageFilters/LabelShapeKeepNObjectsNextToPosImageFilter.h \
           images/filters/imageFilters/LabelShapeKeepNObjectsNextToPosImageFilter.tpp \
           images/filters/imageFilters/MultiChannelBinaryThresholdImageFilter.h \
           images/filters/imageFilters/MultiChannelBinaryThresholdImageFilter.tpp \
           images/filters/imageFilters/MultiChannelBinaryThresholdImageFunction.h \
           images/filters/imageFilters/MultiChannelBinaryThresholdImageFunction.tpp \
           images/filters/imageFilters/MultiChannelConnectedThresholdImageFilter.h \
           images/filters/imageFilters/MultiChannelConnectedThresholdImageFilter.tpp \
           images/filters/imageFilters/SLICImageFilter.h \
           images/filters/imageFilters/SLICImageFilter.tpp \
           tools/JobManager.h 
           

contains( TARGET, CSCore ) {
SOURCES += tools/JobManager.cpp \
           images/CSImage.cpp \
           images/pipelines/EstimateCellShape.cpp \
           images/pipelines/ExtractGraph.cpp \
           images/pipelines/ObjectsToMatrixFilter.cpp \
           images/pipelines/ResampleNetwork.cpp \
           images/pipelines/ProcessHoechstDiffusionFilter.cpp \
           images/pipelines/SegmentAndClassifyStellateCells.cpp \
           images/pipelines/SegmentBileNetworkOnDPPIV20x.cpp \
           images/pipelines/SegmentBileNetworkOnDPPIV60x.cpp \
           images/pipelines/SegmentCellMembraneOnBCat60x.cpp \
           images/pipelines/SegmentNecroticRegionOnDM.cpp \
           images/pipelines/SegmentNucleiOnDAPI20x.cpp \
           images/pipelines/SegmentNucleiOnDAPIWithHough.cpp \
           images/pipelines/SegmentSinusoidalNetworkOnTwoChannels20x.cpp \
           images/pipelines/SegmentSinusoidalNetworkOnTwoChannels60x.cpp \
           images/pipelines/SegmentVeins.cpp \
           images/pipelines/Skeletonization.cpp \
           images/tools/FileFormatConverter.cpp \
           images/tools/GraphAnnotationHelper.cpp \
           images/tools/ImageAnalysisSummaryFileIO.cpp \
           images/filters/analysisFilters/AnalyzeCellsFilter.cpp \
           images/filters/analysisFilters/AnalyzeLobuleFilter.cpp \
           images/filters/analysisFilters/AnalyzeStellateCellsFilter.cpp \
           images/filters/analysisFilters/CompareSegmentations.cpp \
           images/filters/analysisFilters/Lobule2DHandlingFilter.cpp \
           images/filters/convertFilters/GraphToImageFilter.cpp \
           images/filters/convertFilters/SkeletonImageToGraphFilter.cpp \
           images/filters/graphFilters/AnalyzeBileNetworkFilter.cpp \
           images/filters/graphFilters/AnalyzeGraphFilter.cpp \
           images/filters/graphFilters/CollapseIntersectionNodesFilter.cpp \
           images/filters/graphFilters/GeometricalThresholdUndirectedGraphFilter.cpp \
           images/filters/graphFilters/GraphHandlingHelper.cpp \
           images/filters/graphFilters/ResampleUndirectedGraphFilter.cpp \
           images/filters/graphFilters/TopologicalPruneGraphFilter.cpp
}

contains( TARGET, CSGUI ) {
   message("CSGUI with image processing")
   DEPENDPATH += gui/tabImageProcessing \
                 gui/tabToolsWithItkVtk

   INCLUDEPATH += gui/tabImageProcessing \
                  gui/tabToolsWithItkVtk

   HEADERS += gui/tabImageProcessing/GraphViewer.h \
              gui/tabImageProcessing/ImageProcessing.h \
              gui/tabImageProcessing/JobManagerModel.h \
              gui/tabImageProcessing/QCSVTKDisplay.h \
              gui/tabImageProcessing/QCSVTKDisplayRender.h \
              gui/tabImageProcessing/ImageProcessingWorker.h \
              gui/ImageViewer2D.h \
              gui/SeedPicker.h \
              gui/tabToolsWithItkVtk/toolsWithItkVtk.h \
              gui/tabToolsWithItkVtk/vtkPOVExporterWireFrame.h

   FORMS += gui/tabImageProcessing/ImageProcessing.ui \
            gui/tabToolsWithItkVtk/toolsWithItkVtk.ui

   SOURCES += gui/tabImageProcessing/GraphViewer.cpp \
              gui/tabImageProcessing/ImageProcessing.cpp \
              gui/tabImageProcessing/JobManagerModel.cpp \
              gui/tabImageProcessing/QCSVTKDisplay.cpp \
              gui/tabImageProcessing/QCSVTKDisplayRender.cpp \
              gui/tabImageProcessing/ImageProcessingWorker.cpp \
              gui/ImageViewer2D.cpp \
              gui/SeedPicker.cpp \
              tools/JobManager.cpp \
              gui/tabToolsWithItkVtk/toolsWithItkVtk.cpp \
              gui/tabToolsWithItkVtk/vtkPOVExporterWireFrame.cpp

}
