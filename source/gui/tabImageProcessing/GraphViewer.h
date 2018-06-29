///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SLICImageFilter.h                                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-07-09                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHVIEWER_H_
#define GRAPHVIEWER_H_

#include "QCSVTKDisplay.h"

#include <QFileInfo>
#include <QModelIndex>
#include <QString>
#include <QtGui/QTreeView>

#include <vtkCubeAxesActor.h>
#include <vtkGraphLayoutView.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkPoints.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkUndirectedGraph.h>
#include <vtkViewTheme.h>

#include "../../tools/parameters/CSParameterContext.h"
#include "../../tools/parameters/CSParameter.h"


class AnalyzeBileNetworkFilter;

class GraphViewer : public QCSVTKDisplay
{
    Q_OBJECT

public:
    GraphViewer();
    virtual ~GraphViewer();

    void Enable2DDisplayMode() { mIs3DMode = false; };     //default 3d mode
    void Enable3DDisplayMode() { mIs3DMode = true; };     //default 3d mode

    std::string GetName() { return mName; };

    void SetGraph(vtkSmartPointer<vtkUndirectedGraph> graph) { mpDisplayGraph = graph; mDisplayGraphReady = true; };
    void SetGraphs(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs) { mGraphs = graphs; mDisplayGraphReady = false; };
	void SetGraphFilePath(std::string graphFile) { mGraphPath = QString::fromStdString(graphFile); mLoadGraphFromDisk = true; mDisplayGraphReady = false; };
    void SetName(std::string name) { mName = name; };
    void SetScaling(double x, double y, double z) { mpScaling[0] = x; mpScaling[1] = y; mpScaling[2] = z; };

    void Show(bool withGraphPreparation = false);

private slots:
    void AdjustDisplayOptions(const QModelIndex& topLeft, const QModelIndex& bottomRight);
    void PrintGraphInformation(void);

private:
    void LoadGraphsFromDisk();
    void SetUpDisplayGraph();
    void BuildGraphAnnotations();
    void BuildGraphAnnotationSelectionMenus();

    void ToggleAxes();
    void ToggleScaling(bool startupCall);
    void SetupGUI();

    std::string mName;				//displayed widget name

    QString mGraphPath;
    bool mLoadGraphFromDisk;

    std::vector< vtkSmartPointer<vtkUndirectedGraph> > mGraphs;
    vtkSmartPointer<vtkUndirectedGraph> mpDisplayGraph;

    bool mDisplayGraphReady;

    AnalyzeBileNetworkFilter* mpAnalyzeBile;

    bool mIs3DMode;

    vtkSmartPointer<vtkGraphLayoutView> mpGraphLayoutView;
    vtkSmartPointer<vtkViewTheme> mTheme;
    QWidget* mpGraphViewSidebar;
    QTreeView* mpGraphVisOptTreeView;

    vtkSmartPointer<vtkCubeAxesActor> mpAxes;
    const static std::string mAxesUnits[];

    vtkSmartPointer<vtkScalarBarActor> mpVertexColorLegend;
    vtkSmartPointer<vtkScalarBarActor> mpEdgeColorLegend;

    bool mIsScaled;
    double* mpScaling;

    //Parameter Context for Graph Visualization+++++++++++++++++
    CSParameterContext* mpParameters;

    CSParameterContext* mpDisplayOptionsContext;
    bool                mWithAxes;
    CSParameterContext* mpVertexOptionsContext;
    CSParameterChoice*  mpVertexColoringChoice;
    CSParameterChoice*  mpVertexScalingChoice;
    CSParameterChoice*  mpVertexLabelChoice;
    CSParameterContext* mpEdgeOptionsContext;
    CSParameterChoice*  mpEdgeColoringChoice;
    CSParameterChoice*  mpEdgeLabelChoice;
    int                 mEdgeWidth;

    std::vector<std::string> mDisplayModeVertexLabel;
    std::vector<std::string> mDisplayModeVertexColoring;
    std::vector<std::string> mDisplayModeVertexScaling;
    std::vector<std::string> mDisplayModeEdgeLabel;
    std::vector<std::string> mDisplayModeEdgeColoring;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
};

#endif /* GRAPHVIEWER_H_ */
