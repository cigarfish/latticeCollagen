///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ExtractGraph.h                                                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-12-09                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef EXTRACTGRAPH_H_
#define EXTRACTGRAPH_H_

#include "BasePipeline.h"

#include <string>
#include <vector>

#include <itkNeighborhoodIterator.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkUndirectedGraph.h>

#include "../tools/GraphAnnotationHelper.h"

class ExtractGraph : public BasePipeline<3>
{
    static const unsigned int   dim2 = 2;
    static const unsigned int   dim3 = 3;

    typedef bool                                    PixelBoolType;

    typedef itk::Image<PixelBoolType, dim2>         BScalar2DImageType;
    typedef itk::Image<PixelBoolType, dim3>         BScalar3DImageType;

    typedef itk::ImageFileReader<BScalar2DImageType>    BScalar2DReaderType;
    typedef itk::ImageFileWriter<BScalar2DImageType>    BScalar2DWriterType;

    typedef itk::ImageFileReader<BScalar3DImageType>    BScalar3DReaderType;
    typedef itk::ImageFileWriter<BScalar3DImageType>    BScalar3DWriterType;

public:
    ExtractGraph();

	std::string GetNameOfFirstGraph() { return m_nameOfFirstGraph.toStdString(); };
    vtkSmartPointer<vtkImageData> GetInputImageData();
    std::vector< vtkSmartPointer<vtkUndirectedGraph> > GetVanillaGraphs()   { return m_graphs; };
    std::vector< vtkSmartPointer<vtkUndirectedGraph> > GetResampledGraphs() { return m_resampledGraphs; };
    std::vector< vtkSmartPointer<vtkUndirectedGraph> > GetPrunedGraphs()    { return m_prunedGraphs; };
    std::vector< vtkSmartPointer<vtkUndirectedGraph> > GetCollapsedGraphs() { return m_collapsedGraphs; };
    std::vector< vtkSmartPointer<vtkUndirectedGraph> > GetFinalGraphs()     { return m_geomThresholdedGraphs; };
    bool Is3DData() { return m_3DMode; };

    void Update();
private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteFinalGraphsInNicksXMLFormat();
    void WriteDataSetSummary();

    void SaveGraphs(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, int filenameIndex);
    void PrintBasicGraphInfo(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, std::string graphName);


    std::vector< vtkSmartPointer<vtkUndirectedGraph> > m_graphs;
    std::vector< vtkSmartPointer<vtkUndirectedGraph> > m_resampledGraphs;
    std::vector< vtkSmartPointer<vtkUndirectedGraph> > m_prunedGraphs;
    std::vector< vtkSmartPointer<vtkUndirectedGraph> > m_collapsedGraphs;
    std::vector< vtkSmartPointer<vtkUndirectedGraph> > m_geomThresholdedGraphs;

	QString m_fullFilename;
    QFileInfo m_infoFullFilename;
    std::string m_path;
    std::string m_filename;
    std::string m_fileExtension;

    QString m_networkFilename;
    QString m_cvFilename;
    QString m_pvFilename;

    CScalarVoImageType::Pointer mp_networkBin;
    FScalarVoImageType::Pointer mp_networkDistMap;
    bool mHasCV;
    bool mHasPV;
    FScalarVoImageType::Pointer mp_cvDistMap;
    FScalarVoImageType::Pointer mp_pvDistMap;
    double m_maxCVeinConnectednessDistance;
    double m_maxPVeinConnectednessDistance;

    unsigned int m_entryPoint;

    unsigned int m_resamplingFactor;
    double m_maxResamplingDist;
    double m_deadEndPruneThreshold;
    bool m_deadEndPruneAll;
    double m_collapseThreshold;
    double m_pruningAngle;

    double m_defaultRadius;

    unsigned int m_num_dim;
    CScalarVoImageType::SpacingType mSpacing;

    bool m_3DMode;
    bool m_testOutput;

    bool m_saveEverything;                   //1 - save everything, 0 - save only essentials
    std::string m_savePrefix;
    std::string m_filenameSave[3];
    QString m_nameOfFirstGraph;
};

#endif /* EXTRACT2DGRAPH_H_ */
