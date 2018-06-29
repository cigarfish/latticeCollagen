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

#include <sstream>

#include <QString>
#include <QtGui/QGridLayout>
#include <QtGui/QWidget>

#include <vtkAnnotationLink.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkDataSetAttributes.h>
#include <vtkDiscretizableColorTransferFunction.h>
#include <vtkDoubleArray.h>
#include <vtkEdgeListIterator.h>
#include <vtkEventQtSlotConnect.h>
#include <vtkGraphReader.h>
#include <vtkGraphToGlyphs.h>
#include <vtkInteractorStyleRubberBand2D.h>
#include <vtkInteractorStyleRubberBand3D.h>
#include <vtkLookupTable.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkVertexListIterator.h>

#include <vtkIdTypeArray.h>
#include <vtkRenderedGraphRepresentation.h>

#include "../QCSParameterModel.h"
#include "../QCSParameterDelegate.h"
#include "../../images/tools/GraphAnnotationHelper.h"
#include "../../images/filters/graphFilters/AnalyzeBileNetworkFilter.h"

#include "GraphViewer.h"


const std::string GraphViewer::mAxesUnits[] = {"pixel", "micron"};



GraphViewer::GraphViewer()
{
    mName = "default";
    mLoadGraphFromDisk = false;
    mIs3DMode = true;

    mpDisplayGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
    mpDisplayGraph = NULL;

    mDisplayGraphReady = false;

    //Setup graph visualization options context***********************************************************************************
    mWithAxes = true;

    mIsScaled = true;
    mpScaling = new double[3];
    mpScaling[0] = 1.; mpScaling[1] = 1.; mpScaling[2] = 1.;

    mDisplayModeVertexColoring.push_back("no selection");
    mpVertexColoringChoice = new CSParameterChoice(mDisplayModeVertexColoring, 0);

    mDisplayModeVertexScaling.push_back("no selection");
    mpVertexScalingChoice = new CSParameterChoice(mDisplayModeVertexScaling, 0);

    mDisplayModeVertexLabel.push_back("no selection");
    mpVertexLabelChoice = new CSParameterChoice(mDisplayModeVertexLabel, 0);

    mDisplayModeEdgeColoring.push_back("no selection");
    mpEdgeColoringChoice = new CSParameterChoice(mDisplayModeEdgeColoring, 0);

    mDisplayModeEdgeLabel.push_back("no selection");
    mpEdgeLabelChoice = new CSParameterChoice(mDisplayModeEdgeLabel, 0);

    mEdgeWidth = 1;


    mpParameters = new CSParameterContext("Graph Visualization Options");
    mpParameters->setDebug(true);

    mpDisplayOptionsContext = mpParameters->addContext("General Display Options");
#ifndef CS_TI_QUANT_ONLY
    mpDisplayOptionsContext->addParameter("Show axes", CSParameter::Bool, &mWithAxes, "");
    mpDisplayOptionsContext->addParameter("With scaling", CSParameter::Bool, &mIsScaled, "");
#endif
    mpDisplayOptionsContext->addParameter("x-spacing", CSParameter::Double, &(mpScaling[0]), "micron");
    mpDisplayOptionsContext->addParameter("y-spacing", CSParameter::Double, &(mpScaling[1]), "micron");
    mpDisplayOptionsContext->addParameter("z-spacing", CSParameter::Double, &(mpScaling[2]), "micron");
#ifndef CS_TI_QUANT_ONLY
    mpDisplayOptionsContext->setVisible(false);
#endif
    mpVertexOptionsContext = mpParameters->addContext("Vertex Display Options");
    mpVertexOptionsContext->addParameter("Vertex color", CSParameter::Choice, mpVertexColoringChoice, "");
    mpVertexOptionsContext->addParameter("Vertex scale", CSParameter::Choice, mpVertexScalingChoice, "");
    mpVertexOptionsContext->addParameter("Vertex label", CSParameter::Choice, mpVertexLabelChoice, "");
    mpEdgeOptionsContext = mpParameters->addContext("Edge Display Options");
    mpEdgeOptionsContext->addParameter("Edge color", CSParameter::Choice, mpEdgeColoringChoice, "");
    mpEdgeOptionsContext->addParameter("Edge label", CSParameter::Choice, mpEdgeLabelChoice, "");
    mpEdgeOptionsContext->addParameter("Edge width", CSParameter::Int, &mEdgeWidth, "");
    //****************************************************************************************************************************
}


GraphViewer::~GraphViewer()
{
    delete mpGraphViewSidebar;
    delete mpGraphVisOptTreeView;

    delete mpParameters;

    delete mpVertexColoringChoice;
    delete mpVertexScalingChoice;
    delete mpVertexLabelChoice;
    delete mpEdgeColoringChoice;
    delete mpEdgeLabelChoice;

    delete mpAnalyzeBile;
}


void GraphViewer::LoadGraphsFromDisk()
{
    std::cout << "GraphViewer: Load graphs from disk..." << std::endl;

    vtkSmartPointer<vtkGraphReader> reader = vtkSmartPointer<vtkGraphReader>::New();

	QFileInfo fileInfo;
    std::string path, filename, fileExt;

	fileInfo.setFile(mGraphPath);
	std::cout << "fileInfo.suffix() = " << fileInfo.suffix().toStdString() << std::endl;

	if(!fileInfo.exists() || fileInfo.suffix()!=QString("txt")) {
	    std::cout << "Error: throw() !!" << std::endl;
		throw std::string("Please specify a graph file (*.txt)");
	}

	path = (fileInfo.path() + QString("/")).toStdString();
	filename = fileInfo.baseName().toStdString();
	fileExt = (QString(".") + fileInfo.suffix()).toStdString();

    std::string prefix;
    int startIndex = 0;
    bool isSeq  = FilenameParser::IsSequence(filename, prefix, startIndex);

	std::cout << "Load graphs from disk under " << mGraphPath.toStdString() << std::endl;
    int i=startIndex;
    do {
        std::stringstream name;
        if(isSeq)  name << path << prefix << i << ".txt";
        else            name << path << prefix << ".txt";
        std::cout << "load graph file " << name.str() << std::endl;

		QString outPath = QString::fromStdString(name.str());
		fileInfo.setFile(outPath);

        if(fileInfo.exists()) {
			reader->SetFileName(name.str().c_str());
			bool isVTKFile = reader->OpenVTKFile();
			bool hasHeader = reader->ReadHeader();

			if(!isVTKFile || !hasHeader)
			    throw std::string("The specified file is not a VTK graph file (*.txt)");

            reader->Update();

            vtkSmartPointer<vtkUndirectedGraph> undirectedGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
            reader->GetOutput()->ToUndirectedGraph(undirectedGraph);

            mGraphs.push_back(undirectedGraph);

            reader->CloseVTKFile();
        }
        else
            break;

        i++;
    } while(isSeq);
    std::cout << "Loaded " << mGraphs.size() << " graphs from disk." << std::endl;

    std::cout << "GraphViewer: In " << mName << "-Display number of single graphs is " << mGraphs.size() << std::endl;
}


void GraphViewer::SetUpDisplayGraph()
{
    std::cout << "GraphViewer: Setup display graph..." << std::endl;

    vtkSmartPointer<vtkMutableUndirectedGraph> displayGraphHelper = vtkSmartPointer<vtkMutableUndirectedGraph>::New();;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    std::map<vtkIdType, vtkIdType> *vertexIDMap = new std::map<vtkIdType, vtkIdType>[mGraphs.size()];

    for(int k=0; k<mGraphs[0]->GetVertexData()->GetNumberOfArrays(); k++)
        displayGraphHelper->GetVertexData()->DeepCopy(mGraphs[0]->GetVertexData());
    for(int k=0; k<mGraphs[0]->GetEdgeData()->GetNumberOfArrays(); k++)
        displayGraphHelper->GetEdgeData()->DeepCopy(mGraphs[0]->GetEdgeData());
    displayGraphHelper->Squeeze();

    for(unsigned int i=0; i<mGraphs.size(); i++) {
        for(vtkIdType j=0; j<mGraphs[i]->GetNumberOfVertices(); j++) {
            vtkIdType v = displayGraphHelper->AddVertex();
            points->InsertPoint(v, mGraphs[i]->GetPoint(j));

            displayGraphHelper->GetVertexData()->GetPedigreeIds()->InsertVariantValue(v, mGraphs[i]->GetVertexData()->GetPedigreeIds()->GetVariantValue(j));
            for(int k=0; k<mGraphs[i]->GetVertexData()->GetNumberOfArrays(); k++)
                displayGraphHelper->GetVertexData()->GetArray(k)->InsertVariantValue(v, mGraphs[i]->GetVertexData()->GetArray(k)->GetVariantValue(j));

            vertexIDMap[i].insert( std::pair<vtkIdType, vtkIdType>(j,v) );
        }

        vtkSmartPointer<vtkEdgeListIterator> eIt = vtkSmartPointer<vtkEdgeListIterator>::New();
        mGraphs[i]->GetEdges(eIt);

        while(eIt->HasNext()) {
            vtkEdgeType edge = eIt->Next();
            vtkEdgeType e = displayGraphHelper->AddEdge(vertexIDMap[i][edge.Source], vertexIDMap[i][edge.Target]);

            for(int k=0; k<mGraphs[i]->GetEdgeData()->GetNumberOfArrays(); k++)
                displayGraphHelper->GetEdgeData()->GetArray(k)->InsertVariantValue(e.Id, mGraphs[i]->GetEdgeData()->GetArray(k)->GetVariantValue(edge.Id));
        }
    }
    displayGraphHelper->SetPoints(points);

    mpDisplayGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
    mpDisplayGraph->CheckedDeepCopy(displayGraphHelper);

    mDisplayGraphReady = true;

    std::cout << "GraphViewer: Graph has " << mpDisplayGraph->GetNumberOfVertices() << " vertices " << std::endl;
}


void GraphViewer::BuildGraphAnnotations()
{
    std::cout << "GraphViewer: Build graph annotations..." << std::endl;

    std::map<vtkIdType, int> edgeIdToBranchId;
    std::map<vtkIdType, int> edgeIdToSecOrderBranchId;
    std::map<vtkIdType, int> vertexIdToDEndState;
    std::map<vtkIdType, int> vertexIdToIsecState;
    std::map<vtkIdType, int> vertexIdToSecOState;

    GraphAnnotationHelper anno;
    anno.EnableVertexIDAnnotation();
    anno.EnableVertexDegreeAnnotation();
    anno.EnableEdgeIDAnnotation();
    anno.EnableEdgeLengthAnnotation(mpScaling[0], mpScaling[1], mpScaling[2]);
    anno.AddPredefinedAnnotations(mpDisplayGraph);

    std::vector< vtkSmartPointer<vtkUndirectedGraph> >  graphWrapperForAnalysis;
    graphWrapperForAnalysis.push_back(mpDisplayGraph);

//    mpAnalyzeBile = new AnalyzeBileNetworkFilter();
//    mpAnalyzeBile->SetGraphsToAnalyze(graphWrapperForAnalysis);
//    mpAnalyzeBile->AssembleFirstOrderBranches(0, false, edgeIdToBranchId, vertexIdToDEndState, vertexIdToIsecState);
//    mpAnalyzeBile->AssembleSecondOrderBranches(0, vertexIdToSecOState, edgeIdToSecOrderBranchId);
//
//    anno.AddCustomEdgeAnnotation(mpDisplayGraph, "BranchIDs", edgeIdToBranchId, -1);                               //and finally add the branch labels to the edges
//    anno.AddCustomEdgeAnnotation(mpDisplayGraph, "SecondOrderBranchIDs", edgeIdToSecOrderBranchId, -1);
//    anno.AddCustomVertexAnnotation(mpDisplayGraph, "DeadEndBranchFlag", vertexIdToDEndState, -1);
//    anno.AddCustomVertexAnnotation(mpDisplayGraph, "IntersectionBranchFlag", vertexIdToIsecState, -1);
//    anno.AddCustomVertexAnnotation(mpDisplayGraph, "SecondOrderBranchFlag", vertexIdToSecOState, -1);
}


void GraphViewer::BuildGraphAnnotationSelectionMenus()
{
    vtkDataSetAttributes* vd = mpDisplayGraph->GetVertexData();
    for(int i = 0; i < vd->GetNumberOfArrays(); ++i) {
        if(vd->GetArrayName(i) && vd->GetArrayName(i) == '\0') {
            mDisplayModeVertexColoring.push_back( "unnamed array" );
            mDisplayModeVertexScaling.push_back( "unnamed array" );
            mDisplayModeVertexLabel.push_back( "unnamed array" );
        }
        else {
            mDisplayModeVertexColoring.push_back( vd->GetArrayName(i) );
            mDisplayModeVertexScaling.push_back( vd->GetArrayName(i) );
            mDisplayModeVertexLabel.push_back( vd->GetArrayName(i) );
        }
    }
    if(vd->GetNumberOfArrays()==0) {
        mDisplayModeVertexLabel.push_back( "no data available" );
        mDisplayModeVertexColoring.push_back( "no data available" );
        mDisplayModeVertexScaling.push_back( "no data available" );
    }
    mpVertexLabelChoice->setChoices(mDisplayModeVertexLabel, 0);
    mpVertexColoringChoice->setChoices(mDisplayModeVertexColoring, 0);
    mpVertexScalingChoice->setChoices(mDisplayModeVertexScaling, 0);


    vtkDataSetAttributes* ed = mpDisplayGraph->GetEdgeData();
    for(int i = 0; i < ed->GetNumberOfArrays(); ++i) {
        if(ed->GetArrayName(i) && ed->GetArrayName(i) == '\0') {
            mDisplayModeEdgeLabel.push_back( "unnamed array" );
            mDisplayModeEdgeColoring.push_back( "unnamed array" );
        }
        else {
            mDisplayModeEdgeLabel.push_back( ed->GetArrayName(i) );
            mDisplayModeEdgeColoring.push_back( ed->GetArrayName(i) );
        }
    }
    if(ed->GetNumberOfArrays()==0) {
        mDisplayModeEdgeLabel.push_back( "no data available" );
        mDisplayModeEdgeColoring.push_back( "no data available" );
    }

    mpEdgeLabelChoice->setChoices(mDisplayModeEdgeLabel, 0);
    mpEdgeColoringChoice->setChoices(mDisplayModeEdgeColoring, 0);
}


void GraphViewer::AdjustDisplayOptions(const QModelIndex& topLeft, const QModelIndex& bottomRight)
{
    CSParameter* alteredParam = (CSParameter*)(topLeft.internalPointer());

    if(alteredParam->name().compare("Show axes")==0)
        ToggleAxes();
    if(alteredParam->name().compare("With scaling")==0 || alteredParam->name().compare("x-spacing")==0 || alteredParam->name().compare("y-spacing")==0 || alteredParam->name().compare("z-spacing")==0)
        ToggleScaling(false);
    else if(alteredParam->name().compare("Vertex color")==0) {
        std::string colorChoice = ( (CSParameterChoice*)(alteredParam->dataPointer()) )->currentString();

        if(colorChoice.compare("no selection")==0 || colorChoice.compare("no data available")==0) {
            mpGraphLayoutView->ColorVerticesOff();

            mpVertexColorLegend->SetVisibility(0);
        }
        else {
            vtkSmartPointer<vtkDataArray> dataArray = vtkDataArray::SafeDownCast(mpDisplayGraph->GetVertexData()->GetArray(colorChoice.c_str()));

            if(dataArray != NULL) {
                vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
                lut->SetRange(dataArray->GetRange());
                lut->Build();

                mTheme->SetScalePointLookupTable(true);
                mTheme->SetPointLookupTable(lut);

                mpGraphLayoutView->ApplyViewTheme(mTheme);
                mpGraphLayoutView->SetVertexColorArrayName(colorChoice.c_str());
                mpGraphLayoutView->ColorVerticesOn();

                mpVertexColorLegend->SetLookupTable(lut);
                mpVertexColorLegend->SetVisibility(true);
            }
            else {
                mpGraphLayoutView->ColorVerticesOff();

                mpVertexColorLegend->SetVisibility(false);
            }
        }
    }
    else if(alteredParam->name().compare("Vertex scale")==0) {
        std::string scaleChoice = ( (CSParameterChoice*)(alteredParam->dataPointer()) )->currentString();

        std::cout << "scaleChoice = " << scaleChoice << std::endl;

        if(scaleChoice.compare("no selection")==0 || scaleChoice.compare("no data available")==0) {
            mpGraphLayoutView->ScaledGlyphsOff();
        }
        else {
            vtkSmartPointer<vtkDataArray> dataArray = vtkDataArray::SafeDownCast(mpDisplayGraph->GetVertexData()->GetArray(scaleChoice.c_str()));

            if(dataArray != NULL) {
                std::cout << "I am right here in the machine room" << std::endl;

                vtkRenderedGraphRepresentation::SafeDownCast(mpGraphLayoutView->GetRepresentation())->SetGlyphType(vtkGraphToGlyphs::CIRCLE);

                mpGraphLayoutView->ScaledGlyphsOn();
                mpGraphLayoutView->SetScalingArrayName(scaleChoice.c_str());
            }
            else {
                mpGraphLayoutView->ScaledGlyphsOff();
            }
        }
    }
    else if(alteredParam->name().compare("Vertex label")==0) {
        std::string labelChoice = ( (CSParameterChoice*)(alteredParam->dataPointer()) )->currentString();

        if(labelChoice.compare("no selection")==0 || labelChoice.compare("no data available")==0) {
            mpGraphLayoutView->SetVertexLabelVisibility(false);
        }
        else {
            mpGraphLayoutView->SetVertexLabelVisibility(true);
            mpGraphLayoutView->SetVertexLabelArrayName(labelChoice.c_str());
        }
    }
    else if(alteredParam->name().compare("Edge color")==0) {
        std::string colorChoice = ( (CSParameterChoice*)(alteredParam->dataPointer()) )->currentString();

        if(colorChoice.compare("no selection")==0 || colorChoice.compare("no data available")==0) {
            mpGraphLayoutView->ColorEdgesOff();
            mpEdgeColorLegend->SetVisibility(0);
        }
        else {
            vtkSmartPointer<vtkDataArray> dataArray = vtkDataArray::SafeDownCast(mpDisplayGraph->GetEdgeData()->GetArray(colorChoice.c_str()));

            if(dataArray != NULL) {
                vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
                lut->SetRange(dataArray->GetRange());
                lut->Build();

                mTheme->SetScaleCellLookupTable(true);
                mTheme->SetCellLookupTable(lut);

                mpGraphLayoutView->ApplyViewTheme(mTheme);
                mpGraphLayoutView->SetEdgeColorArrayName(colorChoice.c_str());
                mpGraphLayoutView->ColorEdgesOn();

                mpEdgeColorLegend->SetLookupTable(lut);
                mpEdgeColorLegend->SetVisibility(true);
            }
            else {
                mpGraphLayoutView->ColorEdgesOff();

                mpEdgeColorLegend->SetVisibility(false);
            }
        }
    }
    else if(alteredParam->name().compare("Edge label")==0) {
        std::string labelChoice = ( (CSParameterChoice*)(alteredParam->dataPointer()) )->currentString();

        if(labelChoice.compare("no selection")==0 || labelChoice.compare("no data available")==0)
            mpGraphLayoutView->SetEdgeLabelVisibility(false);
        else {
            mpGraphLayoutView->SetEdgeLabelVisibility(true);
            mpGraphLayoutView->SetEdgeLabelArrayName(labelChoice.c_str());
        }
    }
    else if(alteredParam->name().compare("Edge width")==0) {
        mEdgeWidth = *(int*)(alteredParam->dataPointer());

        mTheme->SetLineWidth(mEdgeWidth);

        mpGraphLayoutView->ApplyViewTheme(mTheme);
    }
    this->update();
}


void GraphViewer::PrintGraphInformation(void)
{
    // Forward events
    vtkSelection* selection = mpGraphLayoutView->GetRepresentation()->GetAnnotationLink()->GetCurrentSelection();
    vtkSelectionNode* vertices = NULL;
    vtkSelectionNode* edges = NULL;

    //Prepare selection data
    if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX) {
        vertices = selection->GetNode(0);
    }
    else if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::EDGE) {
        edges = selection->GetNode(0);
    }

    if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX) {
        vertices = selection->GetNode(1);
    }
    else if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::EDGE) {
        edges = selection->GetNode(1);
    }

    vtkIdTypeArray* vertexList = vtkIdTypeArray::SafeDownCast(vertices->GetSelectionList());
    vtkIdTypeArray* edgeList = vtkIdTypeArray::SafeDownCast(edges->GetSelectionList());

    //Print general information
    std::stringstream outputString;

    outputString << vertexList->GetNumberOfTuples() << " vertices selected." << std::endl;
    outputString << edgeList->GetNumberOfTuples() << " edges selected." << std::endl;
    outputString << std::endl;

    //Print out vertex data
    if(vertexList->GetNumberOfTuples()>0) {
        outputString << "Vertex Ids: ";

        for(vtkIdType i=0; i<vertexList->GetNumberOfTuples(); i++) {
            vtkIdType vertexId = vertexList->GetValue(i);
            double *pos = mpDisplayGraph->GetPoint(vertexId);

            outputString << vertexId << "(" << pos[0] << "," << pos[1] << "," << pos[2] << "); ";
        }
    }
    outputString << std::endl;
    outputString << std::endl;

    //Print out edge data
    double accEdgeLength = 0;

    vtkDoubleArray* edgeLengths;
    if(mpDisplayGraph->GetEdgeData()->HasArray("EdgeLengths"))
        edgeLengths = vtkDoubleArray::SafeDownCast( mpDisplayGraph->GetEdgeData()->GetArray("EdgeLengths") );

    if(edgeList->GetNumberOfTuples()>0) {
        outputString << "Edge Ids: ";

        for(vtkIdType i=0; i<edgeList->GetNumberOfTuples(); i++) {
            vtkIdType edgeId = edgeList->GetValue(i);

            if(mpDisplayGraph->GetEdgeData()->HasArray("EdgeLengths")) {
                double edgeLength = edgeLengths->GetValue(edgeId);
                accEdgeLength += edgeLength;

                outputString << edgeId << ", length = " << edgeLength << "; ";
            }
            else
                outputString << edgeId << "; ";
        }
    }
    outputString << std::endl;

    if(mpDisplayGraph->GetEdgeData()->HasArray("EdgeLengths"))
        outputString << "total edge length = " << accEdgeLength;
    else
        outputString << "total edge length unknown, no edge length array attached to graph";


    mpConsole->setText(QString(outputString.str().c_str()));
}


void GraphViewer::ToggleAxes()
{
    mpGraphLayoutView->GetRenderer()->RemoveViewProp(mpAxes.GetPointer());

    if(mWithAxes) {
        mpAxes = vtkSmartPointer<vtkCubeAxesActor>::New();
        mpAxes->SetBounds(mpDisplayGraph->GetBounds());
        mpAxes->SetCamera(mpGraphLayoutView->GetRenderer()->GetActiveCamera());
        mpAxes->SetFlyModeToClosestTriad();
        if(mIsScaled) {
            mpAxes->SetXUnits(mAxesUnits[1].c_str());
            mpAxes->SetYUnits(mAxesUnits[1].c_str());
            mpAxes->SetZUnits(mAxesUnits[1].c_str());
        }
        else {
            mpAxes->SetXUnits(mAxesUnits[0].c_str());
            mpAxes->SetYUnits(mAxesUnits[0].c_str());
            mpAxes->SetZUnits(mAxesUnits[0].c_str());
        }
        mpAxes->XAxisMinorTickVisibilityOff();
        mpAxes->YAxisMinorTickVisibilityOff();
        mpAxes->ZAxisMinorTickVisibilityOff();

        mpGraphLayoutView->GetRenderer()->AddViewProp(mpAxes.GetPointer());
    }
    else
        mpGraphLayoutView->GetRenderer()->RemoveViewProp(mpAxes.GetPointer());
}


void GraphViewer::ToggleScaling(bool startupCall)
{
    vtkSmartPointer<vtkPoints> oldPoints = mpDisplayGraph->GetPoints();
    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();

    if(mIsScaled)
        for(unsigned int i=0; i<oldPoints->GetNumberOfPoints(); i++)
            newPoints->InsertNextPoint(oldPoints->GetPoint(i)[0]*mpScaling[0], oldPoints->GetPoint(i)[1]*mpScaling[1], oldPoints->GetPoint(i)[2]*mpScaling[2]);
    else
        for(unsigned int i=0; i<oldPoints->GetNumberOfPoints(); i++)
            newPoints->InsertNextPoint(oldPoints->GetPoint(i)[0]/mpScaling[0], oldPoints->GetPoint(i)[1]/mpScaling[1], oldPoints->GetPoint(i)[2]/mpScaling[2]);

    mpDisplayGraph->SetPoints(newPoints);

    mpGraphLayoutView->ResetCamera();

    if(!startupCall) {
        ToggleAxes();
        ToggleAxes();
        mpGraphLayoutView->ResetCamera();
    }
}


void GraphViewer::SetupGUI()
{
    BuildGraphAnnotationSelectionMenus();

    //Create QCSVTKWidget
    this->setWindowTitle(mName.c_str());
    this->resize(800,600);

    //Create vtk graph viewer
    mpGraphLayoutView = vtkSmartPointer<vtkGraphLayoutView>::New();
    mpGraphLayoutView->AddRepresentationFromInput(mpDisplayGraph);
    mpGraphLayoutView->SetLayoutStrategy("Pass Through");

    if(mIs3DMode) {
        vtkSmartPointer<vtkInteractorStyleRubberBand3D> printSelectedInteractor = vtkSmartPointer<vtkInteractorStyleRubberBand3D>::New();
        mpGraphLayoutView->SetInteractionModeTo3D();
        mpGraphLayoutView->SetInteractorStyle(printSelectedInteractor);

        vtkEventQtSlotConnect* vtk_qt_connector = vtkEventQtSlotConnect::New();
        vtk_qt_connector->Connect(mpGraphLayoutView, vtkCommand::SelectionChangedEvent, (GraphViewer*)this, SLOT(PrintGraphInformation(void)), 0, 1.0);
    }
    else {
        vtkSmartPointer<vtkInteractorStyleRubberBand2D> printSelectedInteractor = vtkSmartPointer<vtkInteractorStyleRubberBand2D>::New();
        mpGraphLayoutView->SetInteractionModeTo2D();
        mpGraphLayoutView->SetInteractorStyle(printSelectedInteractor);

        vtkEventQtSlotConnect* vtk_qt_connector = vtkEventQtSlotConnect::New();
        vtk_qt_connector->Connect(mpGraphLayoutView, vtkCommand::SelectionChangedEvent, (GraphViewer*)this, SLOT(PrintGraphInformation(void)), 0, 1.0);
    }

    mpGraphLayoutView->ResetCamera();
    mpGraphLayoutView->SetRenderWindow(this->getRenderWindow());                                    //add vtk graph viewer to QCSVTKWidget

    mTheme = vtkSmartPointer<vtkViewTheme>::New();

    //Create sidebar for QCSVTKWidget
    mpGraphViewSidebar = new QWidget(this);
    mpGraphViewSidebar->setObjectName(QString::fromUtf8("graphViewSidebar"));
    mpGraphViewSidebar->setMinimumWidth(200);
    mpGraphViewSidebar->setMaximumWidth(500);
    mpGraphViewSidebar->setGeometry(QRect(600, 0, 200, 200));
    mpGraphViewSidebar->setSizePolicy( QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding) );

    this->addControlWidget(mpGraphViewSidebar);                                                     //add control sidebar to QCSVTKWidget

    QSizePolicy sizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(mpGraphViewSidebar->sizePolicy().hasHeightForWidth());

    //Create layout for sidebar
    QGridLayout* controlLayout = new QGridLayout(mpGraphViewSidebar);
    controlLayout->setObjectName(QString::fromUtf8("controlLayout"));
    mpGraphViewSidebar->setLayout(controlLayout);

    //Create tree view element for sidebar
    mpGraphVisOptTreeView = new QTreeView(mpGraphViewSidebar);
    mpGraphVisOptTreeView->setObjectName(QString::fromUtf8("mpGraphVisOptTreeView"));
    sizePolicy.setHeightForWidth(mpGraphVisOptTreeView->sizePolicy().hasHeightForWidth());
    mpGraphVisOptTreeView->setSizePolicy(sizePolicy);

    controlLayout->addWidget(mpGraphVisOptTreeView);

    mpParameters->setupGUI( mpGraphVisOptTreeView );

    mpGraphVisOptTreeView->resizeColumnToContents(0);
    mpGraphVisOptTreeView->resizeColumnToContents(1);
    mpGraphVisOptTreeView->setAlternatingRowColors(true);
    mpGraphVisOptTreeView->expandAll();

    connect( ((QAbstractItemModel*)(mpGraphVisOptTreeView->model())), SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)), (GraphViewer*)this, SLOT(AdjustDisplayOptions(const QModelIndex&, const QModelIndex&)) );

    ToggleScaling(true);
    ToggleAxes();

    mpVertexColorLegend = vtkSmartPointer<vtkScalarBarActor>::New();
    mpVertexColorLegend->SetVisibility(0);
    mpVertexColorLegend->SetTitle("vertex color legend");
    mpVertexColorLegend->SetMaximumHeightInPixels(300);
    mpVertexColorLegend->SetMaximumWidthInPixels(100);
    mpVertexColorLegend->SetPosition(0.82, 0.1);
    mpGraphLayoutView->GetRenderer()->AddActor(mpVertexColorLegend);

    mpEdgeColorLegend = vtkSmartPointer<vtkScalarBarActor>::New();
    mpEdgeColorLegend->SetVisibility(0);
    mpEdgeColorLegend->SetTitle("edge color legend");
    mpEdgeColorLegend->SetMaximumHeightInPixels(300);
    mpEdgeColorLegend->SetMaximumWidthInPixels(100);
    mpEdgeColorLegend->SetPosition(0.82, 0.63);
    mpGraphLayoutView->GetRenderer()->AddActor(mpEdgeColorLegend);
}


void GraphViewer::Show(bool overwriteGraphAnnotations)
{
    //TODO: for coloring some kind of color table selection in GUI necessary
    //TODO: present some kind of selection to other classes regarding the graph annotations (to prevent overwriting by GraphViewer)

    //TODO: GUI element that holds a plain graph viewer (detour using the bile/sinusoid graph analysis tab is not user-friendly)

    if(mLoadGraphFromDisk) LoadGraphsFromDisk();
    if(!mDisplayGraphReady) SetUpDisplayGraph();
    if(overwriteGraphAnnotations) BuildGraphAnnotations();
    SetupGUI();

    this->show();
}
