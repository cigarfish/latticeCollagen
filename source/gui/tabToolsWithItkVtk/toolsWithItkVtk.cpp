///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  toolsWithItkVtk.cpp                                                  //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes@neitsch.de>                               //
//    Created:  2013-08-13                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "toolsWithItkVtk.h"

#include <map>

#include <QVector3D>

#include "../QCSGLDisplay.h"
#include "../CSGLArena.h"

#include "../../model/BasicDatatypes/Color.h"

#if defined( CS_BUILD_IMAGEPROCESSING )
  #include "itkImage.h"
  #include "itkImageFileWriter.h"
  #if (ITK_VERSION_MAJOR >= 4)
    #include <itkTIFFImageIO.h>
  #endif
  
  #include <vtkContourFilter.h>
  #include <vtkDataSetAttributes.h>
  #include <vtkDecimatePro.h>
  #include <vtkDoubleArray.h>
  #include <vtkImageFlip.h>  
  #include <vtkGraphReader.h>
  #include <vtkLODActor.h>
  #include <vtkMutableUndirectedGraph.h>
  #include <vtkPoints.h>
  #include <vtkProperty.h>
  #include <vtkSmartPointer.h>
  #include <vtkQuadricClustering.h>

  #include <vtkSmartPointer.h>
  #include <vtkPoints.h>
  #include <vtkMutableUndirectedGraph.h>
  #include <vtkGraphToPolyData.h>
  #include <vtkPolyDataMapper.h>
  #include <vtkActor.h>
  #include <vtkRenderWindow.h>
  #include <vtkRenderer.h>
  #include <vtkRenderWindowInteractor.h>

  #include <../tabImageProcessing/QCSVTKDisplayRender.h>
  #include <../../../images/filters/convertFilters/GraphToImageFilter.h>
  #include "../../../images/tools/GraphAnnotationHelper.h"
  
  #include <itkImageFileReader.h>
  #include <itkImageToVTKImageFilter.h>

  typedef unsigned char                                       PixelCharType;
  
  typedef itk::Image<PixelCharType, 3>                        CScalar3DImageType;
  typedef itk::ImageFileReader<CScalar3DImageType>            CScalar3DReaderType;
  typedef itk::ImageToVTKImageFilter<CScalar3DImageType>      ITKVTKScalar3DConnectorType;

#endif


// Constructor
toolsWithItkVtk::toolsWithItkVtk(QWidget *parent) : QWidget(parent)
{

  // Setting up the Qt Designer code
  setupUi(this);

  connect( loadButton, SIGNAL(clicked(bool)), this, SLOT(load()) );

}
toolsWithItkVtk::~toolsWithItkVtk()
{}

void toolsWithItkVtk::load(){
  typedef unsigned char                                   CScalarPixelType;

  typedef itk::Image<CScalarPixelType, 3>                 CScalarVoImageType;

  typedef itk::ImageFileWriter<CScalarVoImageType>        ScalarVoWriterType;

  GraphAnnotationHelper anno;
  
  vtkSmartPointer<vtkGraphReader> reader = vtkSmartPointer<vtkGraphReader>::New();
  reader->SetFileName("../../input/sin_graph1.txt");
  reader->Update();
  
    std::cout << "load graph file " << "../../input/sin_graph1.txt" << std::endl;

  vtkSmartPointer<vtkUndirectedGraph> graph = vtkSmartPointer<vtkUndirectedGraph>::New();
  reader->GetOutput()->ToUndirectedGraph(graph);
  //attend radius!

  std::map<vtkIdType,float> radiusMap;
  
  //for(unsigned int i=0; i<graph->GetNumberOfEdges(); i++)
  //  radiusMap.insert(std::pair<vtkIdType,double>(i, 2.5));
  
  anno.AddCustomEdgeAnnotation(graph, "radius", radiusMap, 2.5);
  
  std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphVector;
  //graphVector.push_back(graph);



  //create a new graph (example)
  vtkSmartPointer<vtkMutableUndirectedGraph> g = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
  // Add 4 vertices to the graph
  vtkIdType v0 = g->AddVertex();
  vtkIdType v1 = g->AddVertex();
  vtkIdType v2 = g->AddVertex();
  vtkIdType v3 = g->AddVertex();
  vtkIdType v4 = g->AddVertex();
  vtkIdType v5 = g->AddVertex();
  vtkIdType v6 = g->AddVertex();
  vtkIdType v7 = g->AddVertex();
  vtkIdType v8 = g->AddVertex();
  vtkIdType v9 = g->AddVertex();
  vtkIdType v10 = g->AddVertex();
  vtkIdType v11 = g->AddVertex();
  vtkIdType v12 = g->AddVertex();
  vtkIdType v13 = g->AddVertex();
  vtkIdType v14 = g->AddVertex();
  vtkIdType v15 = g->AddVertex();
  vtkIdType v16 = g->AddVertex();
  vtkIdType v17 = g->AddVertex();
  vtkIdType v18 = g->AddVertex();
  vtkIdType v19 = g->AddVertex();
  vtkIdType v20 = g->AddVertex();
  vtkIdType v21 = g->AddVertex();




  // Add 3 edges to the graph
  g->AddEdge ( v0, v3 );
  g->AddEdge ( v3, v5 );
  g->AddEdge ( v5, v7 );
  g->AddEdge ( v7, v9 );
  g->AddEdge ( v9, v11 );
  g->AddEdge ( v11, v13 );
  g->AddEdge ( v13, v15 );
  g->AddEdge ( v15, v17 );
  g->AddEdge ( v17, v19 );
  g->AddEdge ( v19, v1 );
  g->AddEdge ( v0, v2 );
  g->AddEdge ( v3, v4 );
  g->AddEdge ( v5, v6 );
  g->AddEdge ( v7, v8 );
  g->AddEdge ( v9, v10 );
  g->AddEdge ( v11, v12 );
  g->AddEdge ( v13, v14 );
  g->AddEdge ( v15, v16 );
  g->AddEdge ( v17, v18 );
  g->AddEdge ( v19, v20 );
  g->AddEdge ( v1, v21 );


  // Create 4 points - one for each vertex
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
 double radius = 10;
 double height = 10;

  points->InsertNextPoint(0.0, 0.0, -5);
  points->InsertNextPoint(0.0, 0.0, 5);
  points->InsertNextPoint( -radius*0.5 ,radius*sqrt(3.)*0.5, -5);
  points->InsertNextPoint(0.0, 0.0, -4);
  points->InsertNextPoint(-radius, 0, -4);
  points->InsertNextPoint(0.0, 0.0, -3);
  points->InsertNextPoint(radius ,0, -3);
  points->InsertNextPoint(0.0, 0.0, -2);
  points->InsertNextPoint(radius*0.5 ,radius*sqrt(3.)*0.5, -2);
  points->InsertNextPoint(0.0, 0.0, -1);
  points->InsertNextPoint(-radius*0.5 ,-radius*sqrt(3.)*0.5, -1);
  points->InsertNextPoint(0.0, 0.0, 0);
  points->InsertNextPoint(radius*0.5 ,-radius*sqrt(3.)*0.5, 0);
  points->InsertNextPoint(0.0, 0.0, 1);
  points->InsertNextPoint(-radius ,0 ,1);
  points->InsertNextPoint(0.0, 0.0, 2);
  points->InsertNextPoint(radius*0.5, radius*sqrt(3.)*0.5, 2);
  points->InsertNextPoint(0.0, 0.0, 3);
  points->InsertNextPoint(radius*0.5 ,-radius*sqrt(3.)*0.5, 3);
  points->InsertNextPoint(0.0, 0.0, 4);
  points->InsertNextPoint(-radius*0.5 ,radius*sqrt(3.)*0.5 ,4);
  points->InsertNextPoint(-radius*0.5 ,-radius*sqrt(3.)*0.5 ,5 );

 // Add the coordinates of the points to the graph
  g->SetPoints(points);

  vtkSmartPointer<vtkUndirectedGraph> g2 = vtkSmartPointer<vtkUndirectedGraph>::New();
  //g2->ToUndirectedGraph(graph);
  g2->DeepCopy(g);
  //add radius
  anno.AddCustomEdgeAnnotation(g2, "radius", radiusMap, 7./20.);//7

  graphVector.push_back(g2);//graph


  // visualization
  /*
  // Convert the graph to a polydata
  vtkSmartPointer<vtkGraphToPolyData> graphToPolyData = vtkSmartPointer<vtkGraphToPolyData>::New();
  graphToPolyData->SetInput(graph);
  graphToPolyData->Update();
  // Create a mapper and actor
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(graphToPolyData->GetOutputPort());
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  // Add the actor to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(.3, .6, .3); // Background color green
  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
  */




  CScalarVoImageType::Pointer image = CScalarVoImageType::New();

  GraphToImageFilter::GraphToImage(graphVector, image, "radius", 5., 5.);
//GraphToImage(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, itk::SmartPointer< itk::Image<unsigned char, 3> > &image, double sampleFactor = 1., double frameFactor = 50.);

  ScalarVoWriterType::Pointer writer = ScalarVoWriterType::New();
  writer->SetFileName("../../input/image.tif");
  writer->SetInput(image);
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
  writer->Update();



  ITKVTKScalar3DConnectorType::Pointer itkvtkConnectorInputImage = ITKVTKScalar3DConnectorType::New();
  itkvtkConnectorInputImage->SetInput(image);
  itkvtkConnectorInputImage->ReleaseDataFlagOn();

  vtkSmartPointer<vtkImageFlip> flipYInputFilter = vtkSmartPointer<vtkImageFlip>::New();
  flipYInputFilter->SetFilteredAxis(1);                           // flip y axis
  flipYInputFilter->SetInput(itkvtkConnectorInputImage->GetOutput());
  flipYInputFilter->Update();

  vtkSmartPointer<vtkContourFilter> contour = vtkSmartPointer<vtkContourFilter>::New();
   contour->SetValue(0, 150);
   contour->SetInput(itkvtkConnectorInputImage->GetOutput());
   contour->ReleaseDataFlagOn();
   contour->Update();

   
  vtkQuadricClustering *decimate = vtkQuadricClustering::New();
   decimate->SetNumberOfXDivisions(390);
   decimate->SetNumberOfYDivisions(390);
   decimate->SetNumberOfZDivisions(390);
   decimate->SetInput(contour->GetOutput());
   
   /*
   vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
  input->ShallowCopy(contour->GetOutput());

   vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
#if VTK_MAJOR_VERSION <= 5
  decimate->SetInputConnection(input->GetProducerPort());
#else
  decimate->SetInputData(input);
#endif
  //decimate->SetTargetReduction(.99); //99% reduction (if there was 100 triangles, now there will be 1)
  decimate->SetTargetReduction(.99); //10% reduction (if there was 100 triangles, now there will be 90)
  decimate->Update();
  */

//  vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
//  extractEdges->SetInputConnection(contour->GetOutputPort());
//  extractEdges->Update();
  
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//  mapper->SetInputConnection(extractEdges->GetOutputPort());
  mapper->SetInputConnection(decimate->GetOutputPort());
  mapper->ScalarVisibilityOff();
  //              mapper2->ImmediateModeRenderingOn();

  vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper2->SetInputConnection(decimate->GetOutputPort());
  mapper2->ScalarVisibilityOff();

  vtkSmartPointer<vtkLODActor> actor = vtkSmartPointer<vtkLODActor>::New();
  actor->SetMapper(mapper);
  actor->SetNumberOfCloudPoints(500);
  actor->GetProperty()->LightingOn();
  actor->GetProperty()->SetColor(1, 0, 0);
  actor->GetProperty()->SetOpacity(1.0);
  actor->GetProperty()->SetSpecular(0.9);


  vtkSmartPointer<vtkLODActor> actor2 = vtkSmartPointer<vtkLODActor>::New();
  actor2->SetMapper(mapper2);
  actor2->SetNumberOfCloudPoints(500);
  actor2->GetProperty()->LightingOn();
  actor2->GetProperty()->SetColor(1, 1, 0);
  actor2->GetProperty()->SetOpacity(1.0);//transparency
  actor2->GetProperty()->SetRepresentationToWireframe();


  QCSVTKDisplayRender* displays = new QCSVTKDisplayRender();
  displays->setGeometry(QRect(0, 0, 800, 600));
  displays->setWindowTitle("Test Page Display");

  displays->renWin = displays->getRenderWindow();

  vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
  ren1->SetBackground(1,1,1);
  ren1->AddActor(actor);
  ren1->AddActor(actor2);

  displays->renWin->AddRenderer(ren1);
  
  ren1->ResetCamera();

  displays->show();
}

