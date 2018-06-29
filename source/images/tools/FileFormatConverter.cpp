///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  MevisFileFormatConverter.h                                            //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-08-01                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "FileFormatConverter.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <math.h>
#include <QDateTime>
#include <QFileInfo>

#if (ITK_VERSION_MAJOR >= 4)
#include <itkNrrdImageIO.h>
#include <itkTIFFImageIO.h>
#endif
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkEdgeListIterator.h"
#include "vtkFloatArray.h"
#include "vtkGraphReader.h"
#include "vtkGraphWriter.h"
#include "vtkImageExport.h"
#include "vtkImageStencil.h"
#include "vtkMergeGraphs.h"
#include "vtkMutableGraphHelper.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkSTLReader.h"
#include "vtkUnsignedLongArray.h"

#include "GraphAnnotationHelper.h"
#include "../filters/graphFilters/GraphHandlingHelper.h"
#include "../../tools/input/FilenameParser.h"
#include "../../gui/tabImageProcessing/GraphViewer.h"



const std::string FileFormatConverter::InputDataName[] =
{
        "Graph in Mevis Format",
        "Graph in Dresden Format",
        "Triangualation in Dresden Format",
        "8-bit RGB 2d image",
        "8-bit RGB 3d image",
        ""
};

const std::string FileFormatConverter::OutputDataName[] =
{
        "Graph in TiQuant Format",
        "Binary image",
        "8-bit greyscale 2d image",
        "8-bit greyscale 3d image",
        ""
};


const std::vector<FileFormatConverter::OutputData> FileFormatConverter::AvailableConversions(FileFormatConverter::InputData input)
{
    std::vector<FileFormatConverter::OutputData> availableOutputData;

    switch(input)
    {
    case Graph_Mevis:
    {
        availableOutputData.push_back(Graph_TiQuant);
    }
    break;
    case Graph_Dresden:
    {
        availableOutputData.push_back(Graph_TiQuant);
    }
    break;
    case Triangualation_Dresden:
    {
        availableOutputData.push_back(Binary_Image);
    }
    break;
    case Bit8_RGB_2D_Image:
    {
        availableOutputData.push_back(Bit8_Greyscale_2D_Image);
    }
    break;
    case Bit8_RGB_3D_Image:
    {
        availableOutputData.push_back(Bit8_Greyscale_3D_Image);
    }
    break;
    default:
        std::cout << "FileFormatConverter::AvailableConversions():  no known conversion for \""
                  << InputDataName[input] <<"\".\n";
        break;
    }
    return availableOutputData;
}


FileFormatConverter::FileFormatConverter()
{
}


FileFormatConverter::~FileFormatConverter()
{
}


void FileFormatConverter::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    this->m_paramContext->dump(parameter);

    file.open((mPath + "log" + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-file-format-converter-------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-file-format-converter---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


void FileFormatConverter::ParseParameterContext()
{
    if(this->m_paramContext->findContext("File Format Converter",0)==NULL) {
        std::cout << "Error: FileFormatConverter: Invalid parameter context" << std::endl;
        return;
    }

    mInputFilename = *(std::string*)(this->m_paramContext->findParameter("Input filename", 0)->dataPointer());

    mInfoFullFilename.setFile(QString::fromStdString( mInputFilename) );
    bool fileExists = mInfoFullFilename.exists();

    if(!fileExists)
        throw std::string("Please specify input file.");

    mPath = (mInfoFullFilename.path() + QString("/")).toStdString();
    mFilename = mInfoFullFilename.completeBaseName().toStdString();
    mFilenameExtension = (QString(".") + mInfoFullFilename.suffix()).toStdString();

    std::string inputFormat = ( (CSParameterChoice*)(m_paramContext->findParameter("Input file format", 0)->dataPointer()) )->currentString();
    for(int i=0; i<NumberInputFileFormats; i++) {
        if( inputFormat.compare(InputDataName[i])==0 ) {
            mInputFileFormat = (InputData)i;
            break;
        }
    }

    std::string outputFormat = ( (CSParameterChoice*)(m_paramContext->findParameter("Output file format", 0)->dataPointer()) )->currentString();
    for(int i=0; i<NumberOutputFileFormats; i++) {
        if( outputFormat.compare(OutputDataName[i])==0 ) {
            mOutputFileFormat = (OutputData)i;
            break;
        }
    }
}


void FileFormatConverter::ComputeBounds(vtkSmartPointer<vtkGraph> graph, double*& bounds, double*& offset)
{
    double rim[3];

    graph->ComputeBounds();
    graph->GetBounds(bounds);

    std::cout << "Original bounds: bounds[0] = " << bounds[0] << " bounds[1] = " << bounds[1] << " bounds[2] = " << bounds[2] << " bounds[3] = " << bounds[3] << " bounds[4] = " << bounds[4] <<
            " bounds[5] = " << bounds[5] << std::endl;

    rim[0] = (bounds[1]-bounds[0])/50.;
    rim[1] = (bounds[3]-bounds[2])/50.;
    rim[2] = (bounds[5]-bounds[4])/50.;

    if(bounds[0]-rim[0]<0)
        offset[0] = 0;
    else
        offset[0] = bounds[0]-rim[0];

    if(bounds[2]-rim[1]<0)
        offset[1] = 0;
    else
        offset[1] = bounds[2]-rim[1];

    if(bounds[4]-rim[2]<0)
        offset[2] = 0;
    else
        offset[2] = bounds[4]-rim[2];

    bounds[0] -= offset[0];
    bounds[1] -= offset[0];
    bounds[2] -= offset[1];
    bounds[3] -= offset[1];
    bounds[4] -= offset[2];
    bounds[5] -= offset[2];

    std::cout << "Corrected bounds: bounds[0] = " << bounds[0] << " bounds[1] = " << bounds[1] << " bounds[2] = " << bounds[2] << " bounds[3] = " << bounds[3] << " bounds[4] = " << bounds[4] <<
                " bounds[5] = " << bounds[5] << std::endl;
}


void FileFormatConverter::ComputeBounds(vtkSmartPointer<vtkGraph> graph1, vtkSmartPointer<vtkGraph> graph2, double*& boundsBoth, double*& offset)
{
    double boundsG1[6], boundsG2[6];
    double rim[3];

    graph1->ComputeBounds();
    graph1->GetBounds(boundsG1);

    graph2->ComputeBounds();
    graph2->GetBounds(boundsG2);

    boundsBoth[0] = std::min(boundsG1[0], boundsG2[0]);
    boundsBoth[2] = std::min(boundsG1[2], boundsG2[2]);
    boundsBoth[4] = std::min(boundsG1[4], boundsG2[4]);
    boundsBoth[1] = std::max(boundsG1[1], boundsG2[1]);
    boundsBoth[3] = std::max(boundsG1[3], boundsG2[3]);
    boundsBoth[5] = std::max(boundsG1[5], boundsG2[5]);

    std::cout << "Original bounds: bounds[0] = " << boundsBoth[0] << " bounds[1] = " << boundsBoth[1] << " bounds[2] = " << boundsBoth[2] << " bounds[3] = " << boundsBoth[3] << " bounds[4] = " << boundsBoth[4] <<
            " bounds[5] = " << boundsBoth[5] << std::endl;

    rim[0] = (boundsBoth[1]-boundsBoth[0])/50.;
    rim[1] = (boundsBoth[3]-boundsBoth[2])/50.;
    rim[2] = (boundsBoth[5]-boundsBoth[4])/50.;

    if(boundsBoth[0]-rim[0]<0)
        offset[0] = 0;
    else
        offset[0] = boundsBoth[0]-rim[0];

    if(boundsBoth[2]-rim[1]<0)
        offset[1] = 0;
    else
        offset[1] = boundsBoth[2]-rim[1];

    if(boundsBoth[4]-rim[2]<0)
        offset[2] = 0;
    else
        offset[2] = boundsBoth[4]-rim[2];

    boundsBoth[0] -= offset[0];
    boundsBoth[1] -= offset[0];
    boundsBoth[2] -= offset[1];
    boundsBoth[3] -= offset[1];
    boundsBoth[4] -= offset[2];
    boundsBoth[5] -= offset[2];

    std::cout << "Corrected bounds: bounds[0] = " << boundsBoth[0] << " bounds[1] = " << boundsBoth[1] << " bounds[2] = " << boundsBoth[2] << " bounds[3] = " << boundsBoth[3] << " bounds[4] = " << boundsBoth[4] <<
                " bounds[5] = " << boundsBoth[5] << std::endl;
}


vtkSmartPointer<vtkMutableDirectedGraph> FileFormatConverter::ReadAndConvertMevisGraphToVTKGraph(std::string path)
{
    vtkSmartPointer<vtkMutableDirectedGraph> graph = vtkSmartPointer<vtkMutableDirectedGraph>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    std::ifstream file;
    std::stringstream data;
    std::map<int, int> pIDToTID;

    file.open(path.c_str());
    data << file.rdbuf();
    file.close();

    std::string word;
    data >> word;

    int numNodes, numEdges;

    vtkSmartPointer<vtkDoubleArray> edgeRadiusData = vtkSmartPointer<vtkDoubleArray>::New();
    edgeRadiusData->SetNumberOfComponents(1);
    edgeRadiusData->SetName("radius");

    while(!data.eof()) {
        if(word.compare("#") == 0) {
            data >> word;

            if(word.compare("numNodes") == 0) {
                data >> numNodes;

                data >> word;    //= "#"
                data >> word;    //= "id"
                data >> word;    //= "xCoord"
                data >> word;    //= "yCoord"
                data >> word;    //= "zCoord"

                for(int i=0; i<numNodes; i++)
                {
                    int id;
                    long double x,y,z;

                    data >> id;
                    data >> x;   data >> y;   data >> z;

                    pIDToTID[id] = graph->AddVertex();
                    points->InsertNextPoint(x, y, z);

                    std::cout << "Node " << id << " added" << std::endl;
                }

            }
            else if(word.compare("numEdges") == 0) {
                data >> numEdges;
                std::cout << "numEdges " << numEdges << std::endl;

                data >> word;    //= "#"
                data >> word;    //= "id"
                data >> word;    //= "initNode"
                data >> word;    //= "termNode"
                data >> word;    //= "parentEdge"
                data >> word;    //= "nDau"
                data >> word;    //= "dauEdge"
                data >> word;    //= "..."
                data >> word;    //= "radiusRelativeToParent"
                data >> word;    //= "flowThroughEdge"

                for(int i=0; i<numEdges; i++)
                {
                    int id, initNode, termNode, parEdge, dauEdge1, dauEdge2;
                    long double relRad, absRad, flow, length;

                    data >> id;
                    data >> initNode;    data >> termNode;
                    data >> parEdge;
                    data >> word;    //ndau
                    data >> dauEdge1;    data >> dauEdge2;
                    data >> relRad;
                    data >> word;    //ATTENTION: normally here flow; but in some files this is +nan
                    data >> word;    //= "#"
                    data >> word;    //= "rad"
                    data >> word;    //= "="
                    data >> absRad;
                    data >> word;    //= ","
                    data >> word;    //= "len"
                    data >> word;    //= "="
                    data >> length;

                    graph->AddEdge(pIDToTID[initNode], pIDToTID[termNode]);
                    edgeRadiusData->InsertNextValue(absRad);

                    std::cout << "Edge " << id << " added" << std::endl;
                }
            }
            else {
                cout << "File parsing error: at " << word << ", 'numNodes' or 'numEdges' expected" << endl;
                exit(0);
            }
        }
        else {
            cout << "File parsing error: at " << word << ", '#' expected" << endl;
            exit(0);
        }
        data >> word;
    }
    graph->SetPoints(points);
    graph->GetEdgeData()->AddArray(edgeRadiusData);

    return graph;
}


std::vector< vtkSmartPointer<vtkMutableUndirectedGraph> > FileFormatConverter::ReadAndConvertDresdenGraphToVTKGraph(std::string path)
{
    std::vector< vtkSmartPointer<vtkMutableUndirectedGraph> > graphs;

    vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIds = vtkSmartPointer<vtkUnsignedLongArray>::New();
    vtkSmartPointer<vtkDoubleArray> radiusData = vtkSmartPointer<vtkDoubleArray>::New();
    radiusData->SetNumberOfComponents(1);
    radiusData->SetName("radius");

    std::ifstream file;
    std::stringstream data;
    std::multimap<int,int> edges_pedIdToPedId;

    file.open(path.c_str());
    data << file.rdbuf();
    file.close();

    std::string header;
    std::getline(data, header);   //skip header

    double lastGraphId_d = 0;

    int debugOutput = 0;
    while(!data.eof()) {
        std::string line;
        std::getline(data, line);
        if(line.empty())
            break;
        std::istringstream lineStream(line);

        std::string graphId_s, nodeId_s, x_s, y_s, z_s, radius_s, neighbor_s;
        std::vector<std::string> neighbors_s;

        std::getline(lineStream, graphId_s, ',');
        std::getline(lineStream, nodeId_s, ',');
        std::getline(lineStream, x_s, ',');
        std::getline(lineStream, y_s, ',');
        std::getline(lineStream, z_s, ',');
        std::getline(lineStream, radius_s, ',');
        while(std::getline(lineStream, neighbor_s, ','))
            neighbors_s.push_back(neighbor_s);

//        if(debugOutput < 10 || debugOutput>5555) {
//            std::cout << "line = " << line << " was cut to:" << std::endl;
//        }

        double graphId_d, nodeId_d, x_d, y_d, z_d, radius_d;
        std::vector<int> neighbors_d;

        graphId_d = ::atof(graphId_s.c_str());
        nodeId_d = ::atof(nodeId_s.c_str());
        x_d = ::atof(x_s.c_str());
        y_d = ::atof(y_s.c_str());
        z_d = ::atof(z_s.c_str());
        radius_d = ::atof(radius_s.c_str());
        for(unsigned int i=0; i<neighbors_s.size(); i++) {
            int neighbor_d = ::atoi(neighbors_s[i].c_str());
            neighbors_d.push_back(neighbor_d);
        }

//        if(debugOutput < 10 || debugOutput>5555) {
//            std::cout << "graphId_d = " << graphId_d << " nodeId_d = " << nodeId_d << " x_d = " << x_d << " y_d = " << y_d << " z_d  = " << z_d << " radius_d = " << radius_d << " numNeighbors = " <<  neighbors_d.size() << std::endl;
//            for(unsigned int i=0; i<neighbors_d.size(); i++)
//                std::cout << " neighbors_d[" << i << "] = " << neighbors_d[i];
//            std::cout << std::endl;
//        }

        if(lastGraphId_d != graphId_d) {
            lastGraphId_d = graphId_d;

            graph->SetPoints(points);
            graph->GetVertexData()->SetPedigreeIds(pedigreeIds);
            graph->GetVertexData()->AddArray(radiusData);

            for(std::multimap<int, int>::iterator it=edges_pedIdToPedId.begin(); it!=edges_pedIdToPedId.end(); ++it)
                if(graph->GetEdgeId(graph->FindVertex((*it).first), graph->FindVertex((*it).second)) == -1) //dresden format mentions edges at both nodes, we don't want two edges between the same two nodes
                    graph->AddEdge(graph->FindVertex((*it).first), graph->FindVertex((*it).second));

            graphs.push_back(graph);

            graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
            points = vtkSmartPointer<vtkPoints>::New();
            pedigreeIds = vtkSmartPointer<vtkUnsignedLongArray>::New();
            radiusData = vtkSmartPointer<vtkDoubleArray>::New();
            radiusData->SetNumberOfComponents(1);
            radiusData->SetName("radius");

            edges_pedIdToPedId.clear();
        }

        graph->AddVertex();
        pedigreeIds->InsertNextValue(nodeId_d);
        points->InsertNextPoint(x_d, y_d, z_d);
        radiusData->InsertNextValue(radius_d);

        for(unsigned int i=0; i<neighbors_d.size(); i++)
            edges_pedIdToPedId.insert(std::pair<int,int>(nodeId_d,neighbors_d[i]));

        debugOutput++;
    }
    graph->SetPoints(points);
    graph->GetVertexData()->SetPedigreeIds(pedigreeIds);
    graph->GetVertexData()->AddArray(radiusData);

    for(std::multimap<int, int>::iterator it=edges_pedIdToPedId.begin(); it!=edges_pedIdToPedId.end(); ++it)
        if(graph->GetEdgeId(graph->FindVertex((*it).first), graph->FindVertex((*it).second)) == -1) //dresden format mentions edges at both nodes, we don't want two edges between the same two nodes
            graph->AddEdge(graph->FindVertex((*it).first), graph->FindVertex((*it).second));

    graphs.push_back(graph);

    return graphs;
}


vtkSmartPointer<vtkPolyData> FileFormatConverter::ReadSTLFile(std::string path)
{
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(path.c_str());
    reader->Update();

    return reader->GetOutput();
}


vtkSmartPointer<vtkImageData> FileFormatConverter::ConvertTriangulationToImageData(vtkSmartPointer<vtkPolyData> polyData)
{
    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
    double bounds[6];
    polyData->GetBounds(bounds);

    double spacing[3]; // desired volume spacing
    spacing[0] = 0.5;
    spacing[1] = 0.5;
    spacing[2] = 0.5;
    whiteImage->SetSpacing(spacing);

    // compute dimensions
    int dim[3];
    for(int i=0; i<3; i++)
        dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));

    whiteImage->SetDimensions(dim);
    whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

    double origin[3];
    origin[0] = bounds[0] + spacing[0] / 2;
    origin[1] = bounds[2] + spacing[1] / 2;
    origin[2] = bounds[4] + spacing[2] / 2;
    whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
    whiteImage->SetScalarTypeToUnsignedChar();
    whiteImage->AllocateScalars();
#else
    whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif
    // fill the image with foreground voxels:
    unsigned char inval = 255;
    unsigned char outval = 0;

    vtkIdType count = whiteImage->GetNumberOfPoints();
    for(vtkIdType i=0; i<count; ++i)
        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);

    // polygonal data --> image stencil:
    vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
    pol2stenc->SetInput(polyData);
#else
    pol2stenc->SetInputData(polyData);
#endif
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
    pol2stenc->Update();

    // cut the corresponding white image and set the background:
    vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
    imgstenc->SetInput(whiteImage);
    imgstenc->SetStencil(pol2stenc->GetOutput());
#else
    imgstenc->SetInputData(whiteImage);
    imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(outval);
    imgstenc->Update();

    return imgstenc->GetOutput();
}


void FileFormatConverter::ConvertVtkToItkImageData(vtkSmartPointer<vtkImageData> vtkImageData, C3DImageType::Pointer &itkImageData)
{
    vtkSmartPointer<vtkImageExport> inputImageExporter = vtkImageExport::New();
    inputImageExporter->SetInput(vtkImageData);

    typename VTKImageImportType::Pointer inputImageImporter = VTKImageImportType::New();
//    ConnectPipelines(inputImageExporter.GetPointer(), inputImageImporter);

    inputImageImporter->SetUpdateInformationCallback(inputImageExporter->GetUpdateInformationCallback());

    inputImageImporter->SetPipelineModifiedCallback(inputImageExporter->GetPipelineModifiedCallback());
    inputImageImporter->SetWholeExtentCallback(inputImageExporter->GetWholeExtentCallback());
    inputImageImporter->SetSpacingCallback(inputImageExporter->GetSpacingCallback());
    inputImageImporter->SetOriginCallback(inputImageExporter->GetOriginCallback());
    inputImageImporter->SetScalarTypeCallback(inputImageExporter->GetScalarTypeCallback());

    inputImageImporter->SetNumberOfComponentsCallback(inputImageExporter->GetNumberOfComponentsCallback());

    inputImageImporter->SetPropagateUpdateExtentCallback(inputImageExporter->GetPropagateUpdateExtentCallback());
    inputImageImporter->SetUpdateDataCallback(inputImageExporter->GetUpdateDataCallback());
    inputImageImporter->SetDataExtentCallback(inputImageExporter->GetDataExtentCallback());
    inputImageImporter->SetBufferPointerCallback(inputImageExporter->GetBufferPointerCallback());
    inputImageImporter->SetCallbackUserData(inputImageExporter->GetCallbackUserData());

    itkImageData = const_cast<C3DImageType*>(inputImageImporter->GetOutput());
    itkImageData->Update();

    if(itkImageData.IsNull())
        throw std::string("ConvertVtkToItkImageData failed: ITK image is invalid.");
}


void FileFormatConverter::ConvertRGBToGreyscaleImageData(RGBC2DImageType::Pointer &rgbImageData, C2DImageType::Pointer &scalarImageData, int mode)
{
    itk::ImageRegionConstIterator<RGBC2DImageType> itRGB(rgbImageData, rgbImageData->GetLargestPossibleRegion());
    itk::ImageRegionIterator<C2DImageType> itGrey(scalarImageData, scalarImageData->GetLargestPossibleRegion());

    while(!itRGB.IsAtEnd())
    {
        if(mode==0) {
            int grey = (double)(itRGB.Get().GetRed() + itRGB.Get().GetGreen() + itRGB.Get().GetBlue()) / 3.;
            itGrey.Set(grey);
        }
        else if(mode==1)
            itGrey.Set(itRGB.Get().GetRed());
        else if(mode==2)
            itGrey.Set(itRGB.Get().GetGreen());
        else if(mode==3)
            itGrey.Set(itRGB.Get().GetBlue());

        ++itRGB;
        ++itGrey;
    }

    if(scalarImageData.IsNull())
        throw std::string("ConvertRGBToGreyscaleImageData failed: scalar image is invalid.");
}


void FileFormatConverter::ConvertRGBToGreyscaleImageData(RGBC3DImageType::Pointer &rgbImageData, C3DImageType::Pointer &scalarImageData, int mode)
{
    itk::ImageRegionConstIterator<RGBC3DImageType> itRGB(rgbImageData, rgbImageData->GetLargestPossibleRegion());
    itk::ImageRegionIterator<C3DImageType> itGrey(scalarImageData, scalarImageData->GetLargestPossibleRegion());

    while(!itRGB.IsAtEnd())
    {
        if(mode==0) {
            int grey = (double)(itRGB.Get().GetRed() + itRGB.Get().GetGreen() + itRGB.Get().GetBlue()) / 3.;
            itGrey.Set(grey);
        }
        else if(mode==1)
            itGrey.Set(itRGB.Get().GetRed());
        else if(mode==2)
            itGrey.Set(itRGB.Get().GetGreen());
        else if(mode==3)
            itGrey.Set(itRGB.Get().GetBlue());

        ++itRGB;
        ++itGrey;
    }

    if(scalarImageData.IsNull())
        throw std::string("ConvertRGBToGreyscaleImageData failed: scalar image is invalid.");
}


void FileFormatConverter::ReadGraph(std::string fullFilename, std::vector< vtkSmartPointer<vtkMutableUndirectedGraph> > &graphs)
{
    QFileInfo infoFilename;
    infoFilename.setFile(QString::fromStdString(fullFilename));
    std::string path = (infoFilename.path() + QString("/")).toStdString();
    std::string filename = infoFilename.baseName().toStdString();
    std::string fileExtension = (QString(".") + infoFilename.suffix()).toStdString();

    vtkSmartPointer<vtkGraphReader> reader = vtkSmartPointer<vtkGraphReader>::New();

    std::string prefix;
    int startIndex = 0;
    bool isSeq  = FilenameParser::IsSequence(filename, prefix, startIndex);

    std::cout << "Load graphs from disk under " << path << std::endl;
    int i=startIndex;
    do {
        std::stringstream name;
        if(isSeq)   name << path << prefix << i << fileExtension;
        else        name << path << prefix << fileExtension;
        std::cout << "load graph file " << name.str() << std::endl;

        QString graphFile = QString::fromStdString(name.str());
        QFileInfo infoFile;
        infoFile.setFile(graphFile);

        if(infoFile.exists()) {
            reader->SetFileName(name.str().c_str());
            bool isVTKFile = reader->OpenVTKFile();
            bool hasHeader = reader->ReadHeader();

            if(!isVTKFile || !hasHeader)
                throw std::string("The specified file is not a VTK graph file (*.txt)");

            reader->Update();

            vtkSmartPointer<vtkMutableUndirectedGraph> undirectedGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
            reader->GetOutput()->ToUndirectedGraph(undirectedGraph);

            graphs.push_back(undirectedGraph);

            reader->CloseVTKFile();
        }
        else
            break;

        i++;
    } while(isSeq);
}


void FileFormatConverter::SaveGraph(std::string path, vtkSmartPointer<vtkMutableDirectedGraph> graph)
{
    vtkSmartPointer<vtkGraphWriter> writer = vtkSmartPointer<vtkGraphWriter>::New();
    writer->SetFileName(path.c_str());
    writer->SetInput(graph);
    writer->Update();
}


void FileFormatConverter::SaveGraph(std::string path, vtkSmartPointer<vtkMutableUndirectedGraph> graph, double spacingX, double spacingY, double spacingZ)
{
    vtkSmartPointer<vtkPoints> points = graph->GetPoints();
    for(unsigned int i=0; i<graph->GetNumberOfVertices(); i++) {
        points->SetPoint(i, points->GetPoint(i)[0]*spacingX, points->GetPoint(i)[1]*spacingY, points->GetPoint(i)[2]*spacingZ);
    }
    graph->SetPoints(points);

    vtkSmartPointer<vtkGraphWriter> writer = vtkSmartPointer<vtkGraphWriter>::New();
    writer->SetFileName(path.c_str());
    writer->SetInput(graph);
    writer->Update();
}


void FileFormatConverter::ReadImage(std::string path, RGBC2DImageType::Pointer &image)
{
//    if(BasePipeline::GetNumberOfDimensions(path) != 2)
//        throw std::string("Data " + path + " has different dimensionality than expected by pipeline.");

    typename RGBC2DReaderType::Pointer reader = RGBC2DReaderType::New();
    reader->SetFileName(path);
    reader->ReleaseDataBeforeUpdateFlagOn();
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    reader->Update();

    image = reader->GetOutput();
    image->DisconnectPipeline();
}


void FileFormatConverter::SaveImage(std::string path, C2DImageType::Pointer image)
{
    typename C2DWriterType::Pointer writer = C2DWriterType::New();
    writer->SetFileName(path);
    writer->SetInput(image);
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();
}


void FileFormatConverter::SaveImage(std::string path, C3DImageType::Pointer image)
{
    typename C3DWriterType::Pointer writer = C3DWriterType::New();
    writer->SetFileName(path);
    writer->SetInput(image);
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();
}


//Deprecated: not supported by new input -> output scheme
void FileFormatConverter::ReadMevisDataAndConvertToImageData()
{
    std::string inputMevisFilenameCV = "/scratch/friebel/TestData/LobulePrediction/DPH-077-TreeA_02.tree";
    std::string inputMevisFilenamePV = "/scratch/friebel/TestData/LobulePrediction/DPH-077-TreeB_02.tree";
    std::string outputMevisFilenameCVEndSections = "/scratch/friebel/TestData/LobulePrediction/cvEndSections.tif";
    std::string outputMevisFilenameCVRestSections = "/scratch/friebel/TestData/LobulePrediction/cvRestSections.tif";
    std::string outputMevisFilenamePVEndSections = "/scratch/friebel/TestData/LobulePrediction/pvEndSections.tif";
    std::string outputMevisFilenamePVRestSections = "/scratch/friebel/TestData/LobulePrediction/pvRestSections.tif";

    vtkSmartPointer<vtkMutableDirectedGraph> graphCV = ReadAndConvertMevisGraphToVTKGraph(inputMevisFilenameCV);
    vtkSmartPointer<vtkMutableDirectedGraph> graphPV = ReadAndConvertMevisGraphToVTKGraph(inputMevisFilenamePV);

    std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs;
    vtkSmartPointer<vtkUndirectedGraph> undirectedGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
    graphCV->ToUndirectedGraph(undirectedGraph);
    graphs.push_back(undirectedGraph);
    graphPV->ToUndirectedGraph(undirectedGraph);
    graphs.push_back(undirectedGraph);

    double scaling [] = {1,1,1};

    GraphViewer *graphViewer = new GraphViewer;
    graphViewer->SetName("CV & PV Tree");
    graphViewer->SetScaling(scaling[0], scaling[1], scaling[2]);
    graphViewer->SetGraphs(graphs);
    graphViewer->Show(true);

    vtkSmartPointer<vtkGraphWriter> cvWriter = vtkSmartPointer<vtkGraphWriter>::New();
    cvWriter->SetFileName("/scratch/friebel/TestData/LobulePrediction/cvTree.txt");
    cvWriter->SetInput(graphCV);
    cvWriter->Update();

    vtkSmartPointer<vtkGraphWriter> pvWriter = vtkSmartPointer<vtkGraphWriter>::New();
    pvWriter->SetFileName("/scratch/friebel/TestData/LobulePrediction/pvTree.txt");
    pvWriter->SetInput(graphPV);
    pvWriter->Update();

    double* bounds = new double[6];
    double* offset = new double[3];

    ComputeBounds(graphCV, graphPV, bounds, offset);

    ConvertTreeToImage(graphCV, outputMevisFilenameCVEndSections, outputMevisFilenameCVRestSections, bounds, offset);
    ConvertTreeToImage(graphPV, outputMevisFilenamePVEndSections, outputMevisFilenamePVRestSections, bounds, offset);

    delete [] bounds;
    delete [] offset;
}


//Deprecated: not supported by new input -> output scheme
void FileFormatConverter::ConvertTreeToImage(vtkSmartPointer<vtkMutableDirectedGraph> graph, std::string pathEndSections, std::string pathRestSections, double *bounds, double *offset)
{
    double sampleFac = 2.;

    C3DImageType::SizeType size;
    size[0] = bounds[1]*sampleFac;
    size[1] = bounds[3]*sampleFac;
    size[2] = bounds[5]*sampleFac;

    C3DImageType::SpacingType spacing;
    spacing[0] =  1. / sampleFac;
    spacing[1] =  1. / sampleFac;
    spacing[2] =  1. / sampleFac;

    std::cout << "size = " << size << " spacing = " << spacing << std::endl;

    const CPixelType airHounsfieldUnits  = 0;
    const CPixelType boneHounsfieldUnits = 255;

    GroupType::Pointer endSectionsGroup = GroupType::New();
    GroupType::Pointer restSectionsGroup = GroupType::New();

    vtkDoubleArray* edgeRadius = vtkDoubleArray::SafeDownCast( graph->GetEdgeData()->GetArray("radius") );


    vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
    graph->GetEdges(it);

    int numEdges = 0;

    while(it->HasNext()) {
        vtkEdgeType e = it->Next();

        TubeType::PointListType list;

        TubePointType source, target;
        source.SetPosition(graph->GetPoint(e.Source)[0]-offset[0], graph->GetPoint(e.Source)[1]-offset[1], graph->GetPoint(e.Source)[2]-offset[2]);
        source.SetRadius(edgeRadius->GetValue(e.Id));
        list.push_back(source);
//        std::cout << "tube " << numEdges << " starting point = (" << (graph->GetPoint(e.Source)[0]-mOffset[0])*sampleFac << ", " << (graph->GetPoint(e.Source)[1]-mOffset[1])*sampleFac
//                << ", " << (graph->GetPoint(e.Source)[2]-mOffset[2])*sampleFac << ")" << std::endl;

        target.SetPosition(graph->GetPoint(e.Target)[0]-offset[0], graph->GetPoint(e.Target)[1]-offset[1], graph->GetPoint(e.Target)[2]-offset[2]);
        target.SetRadius(edgeRadius->GetValue(e.Id));
        list.push_back(target);
//        std::cout << "tube " << numEdges << " end point = (" << (graph->GetPoint(e.Target)[0]-mOffset[0])*sampleFac << ", " << (graph->GetPoint(e.Target)[1]-mOffset[1])*sampleFac
//                << ", " << (graph->GetPoint(e.Target)[2]-mOffset[2])*sampleFac << ")" << std::endl;
//        std::cout << "tube " << numEdges << " radius = " << edgeRadius->GetValue(e.Id) << "(*" << sampleFac << ")" << std::endl;

        TubeType::Pointer tube = TubeType::New();
        tube->SetPoints(list);

        if(graph->GetDegree(e.Target) == 1) {
//            std::cout << "Add tube number " << numEdges << " to endSectionsGroup" << std::endl;
            endSectionsGroup->AddSpatialObject(tube);
        }
        else {
//            std::cout << "Add tube number " << numEdges << " to restSectionsGroup" << std::endl;
            restSectionsGroup->AddSpatialObject(tube);
        }

        numEdges++;
    }

    SpatialObjectToImageFilterType::Pointer endSectionConverter = SpatialObjectToImageFilterType::New();
    endSectionConverter->SetInput(endSectionsGroup);
    endSectionConverter->SetInsideValue(boneHounsfieldUnits);
    endSectionConverter->SetOutsideValue(airHounsfieldUnits);
    endSectionConverter->SetSize(size);
    endSectionConverter->SetSpacing(spacing);
    endSectionConverter->SetMaskResampleFactor(4);
    endSectionConverter->SetMaskDilationSize(2);
    endSectionConverter->Update();

    C3DWriterType::Pointer endSectionWriter = C3DWriterType::New();
    endSectionWriter->SetFileName(pathEndSections);
    endSectionWriter->SetInput(endSectionConverter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    endSectionWriter->SetImageIO( itk::TIFFImageIO::New() );
#endif
    endSectionWriter->Update();

    SpatialObjectToImageFilterType::Pointer restSectionConverter = SpatialObjectToImageFilterType::New();
    restSectionConverter->SetInput(restSectionsGroup);
    restSectionConverter->SetInsideValue(boneHounsfieldUnits);
    restSectionConverter->SetOutsideValue(airHounsfieldUnits);
    restSectionConverter->SetSize(size);
    restSectionConverter->SetSpacing(spacing);
    restSectionConverter->SetMaskResampleFactor(4);
    restSectionConverter->SetMaskDilationSize(2);
    restSectionConverter->Update();

    C3DWriterType::Pointer restSectionWriter = C3DWriterType::New();
    restSectionWriter->SetFileName(pathRestSections);
    restSectionWriter->SetInput(restSectionConverter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    restSectionWriter->SetImageIO( itk::TIFFImageIO::New() );
#endif
    restSectionWriter->Update();
}


//Deprecated: not supported by new input -> output scheme
void FileFormatConverter::ComposeSegmentationFilesForInterfaceToStefan()
{
    std::string bileSegFile = "/home/ls1/friebel/Workspace/Data/Confocal/Bile_Canaliculi_Analysis/Deconvoluted/122/bile_step6_bin.tif";
    std::string sinSegFile = "/home/ls1/friebel/Workspace/Data/Confocal/Bile_Canaliculi_Analysis/Deconvoluted/122/sinus_step5_bin.tif";
    std::string cveinSegFile = "/home/ls1/friebel/Workspace/Data/Confocal/Bile_Canaliculi_Analysis/Deconvoluted/122/vein_central_bin.tif";
    std::string pveinSegFile = "/home/ls1/friebel/Workspace/Data/Confocal/Bile_Canaliculi_Analysis/Deconvoluted/122/vein_portal_bin.tif";
    std::string cellSegFile = "/home/ls1/friebel/Workspace/Data/Confocal/Bile_Canaliculi_Analysis/Deconvoluted/122/cellShape_step3_bin.tif";


    C3DReaderType::Pointer bileReader = C3DReaderType::New();
    bileReader->SetFileName(bileSegFile);
#if (ITK_VERSION_MAJOR >= 4)
    bileReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    bileReader->Update();

    C3DReaderType::Pointer sinReader = C3DReaderType::New();
    sinReader->SetFileName(sinSegFile);
#if (ITK_VERSION_MAJOR >= 4)
    sinReader->SetImageIO( itk::TIFFImageIO::New() );
#endif

    C3DReaderType::Pointer cvReader = C3DReaderType::New();
    cvReader->SetFileName(cveinSegFile);
#if (ITK_VERSION_MAJOR >= 4)
    cvReader->SetImageIO( itk::TIFFImageIO::New() );
#endif

    C3DReaderType::Pointer pvReader = C3DReaderType::New();
    pvReader->SetFileName(pveinSegFile);
#if (ITK_VERSION_MAJOR >= 4)
    pvReader->SetImageIO( itk::TIFFImageIO::New() );
#endif

    C3DReaderType::Pointer cellReader = C3DReaderType::New();
    cellReader->SetFileName(cellSegFile);
#if (ITK_VERSION_MAJOR >= 4)
    cellReader->SetImageIO( itk::TIFFImageIO::New() );
#endif

    AddImageFilterType::Pointer add1Filter = AddImageFilterType::New ();
    add1Filter->SetInput1(sinReader->GetOutput());
    add1Filter->SetInput2(cvReader->GetOutput());

    AddImageFilterType::Pointer add2Filter = AddImageFilterType::New ();
    add2Filter->SetInput1(add1Filter->GetOutput());
    add2Filter->SetInput2(pvReader->GetOutput());

    MaskImageFilterType::Pointer mask1Filter = MaskImageFilterType::New ();
    mask1Filter->SetInput(cellReader->GetOutput());
    mask1Filter->SetMaskImage(add2Filter->GetOutput());
    mask1Filter->SetOutsideValue(100);
    mask1Filter->Update();

    C3DImageType::Pointer image = mask1Filter->GetOutput();
    image->DisconnectPipeline();

    C3DWriterType::Pointer writer = C3DWriterType::New();
    writer->SetFileName("/home/ls1/friebel/Workspace/Data/Confocal/Bile_Canaliculi_Analysis/Deconvoluted/122/aComposition1.tif");
    writer->SetInput(image);
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();

    ImageToShapeLabelMapFilterType::Pointer cellImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    cellImageToShaLabMapFilter->SetInput(cellReader->GetOutput());
    cellImageToShaLabMapFilter->SetFullyConnected(false);

    LabelMapToLabelImageFilterType::Pointer cellLabelMapToImage = LabelMapToLabelImageFilterType::New();
    cellLabelMapToImage->SetInput(cellImageToShaLabMapFilter->GetOutput());
    cellLabelMapToImage->Update();

    LImageType::Pointer cellLabelImage = cellLabelMapToImage->GetOutput();
    cellLabelImage->DisconnectPipeline();

    CastLIImageFilterType::Pointer castFilter = CastLIImageFilterType::New();
    castFilter->SetInput(cellLabelImage);

    IWriterType::Pointer iwriter = IWriterType::New();
    iwriter->SetFileName("/home/ls1/friebel/Workspace/Data/Confocal/Bile_Canaliculi_Analysis/Deconvoluted/122/aComposition0.nrrd");
    iwriter->SetInput(castFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    iwriter->SetImageIO( itk::NrrdImageIO::New() );
#endif
    iwriter->Update();

    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);

    BoundaryConditionType boundaryCondition;
    boundaryCondition.SetConstant(0);

    NeighborhoodIteratorType it(radius, image, image->GetLargestPossibleRegion());
    it.SetBoundaryCondition(boundaryCondition);

    const unsigned int neighborhoodSize = 6;

    NeighborhoodIteratorType::OffsetType neighbor[neighborhoodSize];
    neighbor[0][0] = 1;     neighbor[1][0] = -1;
    neighbor[0][1] = 0;     neighbor[1][1] = 0;
    neighbor[0][2] = 0;     neighbor[1][2] = 0;

    neighbor[2][0] = 0;     neighbor[3][0] = 0;
    neighbor[2][1] = 1;     neighbor[3][1] = -1;
    neighbor[2][2] = 0;     neighbor[3][2] = 0;

    neighbor[4][0] = 0;     neighbor[5][0] = 0;
    neighbor[4][1] = 0;     neighbor[5][1] = 0;
    neighbor[4][2] = 1;     neighbor[5][2] = -1;

    bool changes = true;

    while(changes) {
        changes = false;

        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
            if(it.GetCenterPixel() == 0) {
                std::set<unsigned long> cellHits;        //encountered cell label
                bool inContactWithSinusoid = false;

                for(unsigned int i=0; i<neighborhoodSize; i++) {
                    if( image->GetLargestPossibleRegion().IsInside(it.GetIndex()+neighbor[i]) ) {
                        unsigned char pI = image->GetPixel(it.GetIndex()+neighbor[i]);
                        unsigned long pC = cellLabelImage->GetPixel(it.GetIndex()+neighbor[i]);

                        if(pI==100 || pI==175)
                            inContactWithSinusoid = true;
                        if(pC!=0 && cellHits.count(pC)==0)
                            cellHits.insert(pC);
                    }
                }
                if(inContactWithSinusoid && cellHits.size()<2) {
                    it.SetCenterPixel(175);
                    changes = true;
                }
            }
        }
    }

    writer->SetFileName("/home/ls1/friebel/Workspace/Data/Confocal/Bile_Canaliculi_Analysis/Deconvoluted/122/aComposition2.tif");
    writer->SetInput(image);
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();
}


//Deprecated: not supported by new input -> output scheme
void FileFormatConverter::ComposeSegmentationFilesForInterfaceToStefan2()
{
    std::string sinGraphFile = "/scratch/friebel/TestData/connectNetworkAndVeins/Sin_graph0.txt";
    std::string cveinSegFile = "/scratch/friebel/TestData/connectNetworkAndVeins/vein_central_bin.tif";
    std::string pveinSegFile = "/scratch/friebel/TestData/connectNetworkAndVeins/vein_portal_bin.tif";

    std::vector< vtkSmartPointer<vtkMutableUndirectedGraph> > sinGraphs;
    ReadGraph(sinGraphFile, sinGraphs);

    GraphAnnotationHelper *anno = new GraphAnnotationHelper();
    anno->EnableVertexIDAnnotation();
    anno->EnableEdgeIDAnnotation();
    anno->AddPredefinedAnnotations(sinGraphs[0]);

    vtkSmartPointer<vtkPoints> sinPoints = vtkSmartPointer<vtkPoints>::New();
    sinPoints = sinGraphs[0]->GetPoints();

    vtkSmartPointer<vtkMutableGraphHelper> graphHelper = vtkMutableGraphHelper::New();
    graphHelper->SetGraph(sinGraphs[0]);

    vtkSmartPointer<vtkMergeGraphs> mergeGraphsFilter = vtkMergeGraphs::New();

    for(unsigned int i=0; i<sinGraphs.size(); i++) {
        unsigned long numVert = graphHelper->GetGraph()->GetNumberOfVertices();

        std::map<vtkIdType, int> vertIdToState;
        for(unsigned int j=0; j<sinGraphs[i]->GetNumberOfVertices(); j++)
            vertIdToState.insert(std::pair<vtkIdType, int>(j,0));
        anno->AddCustomVertexAnnotation(sinGraphs[i], "Network type", vertIdToState, -1);

        if(i!=0) {
            anno->SetVertexIDOffset(numVert);
            anno->AddPredefinedAnnotations(sinGraphs[i]);

            vtkSmartPointer<vtkPoints> sinTempPoints = vtkSmartPointer<vtkPoints>::New();
            sinTempPoints = sinGraphs[i]->GetPoints();
            for(unsigned int j=0; j<sinTempPoints->GetNumberOfPoints(); j++)
                sinPoints->InsertNextPoint(sinTempPoints->GetPoint(j));

            mergeGraphsFilter->ExtendGraph(graphHelper, sinGraphs[i]);
        }
    }

    C3DReaderType::Pointer cvReader = C3DReaderType::New();
    cvReader->SetFileName(cveinSegFile);
#if (ITK_VERSION_MAJOR >= 4)
    cvReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    cvReader->Update();

    ImageToShapeLabelMapFilterType::Pointer cvImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    cvImageToShaLabMapFilter->SetInput(cvReader->GetOutput());
    cvImageToShaLabMapFilter->SetFullyConnected(false);
    cvImageToShaLabMapFilter->Update();

    LabelMapToLabelImageFilterType::Pointer cvLabelMapToImage = LabelMapToLabelImageFilterType::New();
    cvLabelMapToImage->SetInput(cvImageToShaLabMapFilter->GetOutput());
    cvLabelMapToImage->Update();

    C3DReaderType::Pointer pvReader = C3DReaderType::New();
    pvReader->SetFileName(pveinSegFile);
#if (ITK_VERSION_MAJOR >= 4)
    pvReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    pvReader->Update();

    ImageToShapeLabelMapFilterType::Pointer pvImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    pvImageToShaLabMapFilter->SetInput(pvReader->GetOutput());
    pvImageToShaLabMapFilter->SetFullyConnected(false);
    pvImageToShaLabMapFilter->Update();

    LabelMapToLabelImageFilterType::Pointer pvLabelMapToImage = LabelMapToLabelImageFilterType::New();
    pvLabelMapToImage->SetInput(pvImageToShaLabMapFilter->GetOutput());
    pvLabelMapToImage->Update();


    vtkSmartPointer<vtkMutableUndirectedGraph> veinGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    std::map<vtkIdType, int> vertIdToState;


    std::vector<int> slicesToVisit;
    std::map<int, std::map<unsigned long, unsigned long*> > slicesOfVeins;       //slice nr // vein id // central index

    double zSize = cvReader->GetOutput()->GetLargestPossibleRegion().GetSize(2);
    double stepWidthIdeal = 10.;
    int stepWidth = zSize / (double)ceil(zSize / stepWidthIdeal);
    std::cout << "stepWidth = " << stepWidth << std::endl;

    slicesToVisit.push_back(0);
    int j=1;
    while(slicesToVisit[slicesToVisit.size()-1] <= zSize-1) {
        slicesToVisit.push_back((stepWidth*j)-1);
        j++;
    }
    slicesToVisit.pop_back();
    if(slicesToVisit[slicesToVisit.size()-1] != zSize-1)
        slicesToVisit.push_back(zSize-1);

    for(unsigned int i=0; i<slicesToVisit.size(); i++)
        std::cout << "slicesToVisit[" << i << "] = " << slicesToVisit[i] << std::endl;


    for(unsigned int i=0; i<slicesToVisit.size(); i++) {
        std::map<unsigned long, unsigned long*> veinCentrals;
        for(unsigned int j=0; j<cvImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); j++) {
            unsigned long *idx = new unsigned long[3];
            idx[0] = 0; idx[1] = 0; idx[2] = 0;
            veinCentrals.insert(std::pair<unsigned long, unsigned long*>(cvImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(j)->GetLabel(), idx));
        }
        slicesOfVeins.insert(std::pair<int, std::map<unsigned long, unsigned long*> >(slicesToVisit[i], veinCentrals));
    }

    std::map<unsigned long, vtkIdType> veinsLastVert;
    for(unsigned int i=0; i<slicesToVisit.size(); i++) {
        std::vector<unsigned long> veinIds;
        std::map<unsigned long, unsigned long> veinPixelNumber;
        for(unsigned int j=0; j<cvImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); j++) {
            veinIds.push_back(cvImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(j)->GetLabel());
            veinPixelNumber.insert(std::pair<unsigned long, unsigned long>(cvImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(j)->GetLabel(), 0));
        }
        itk::Index<3> sliceOrigin;
        sliceOrigin[0] = 0;
        sliceOrigin[1] = 0;
        sliceOrigin[2] = slicesToVisit[i];

        itk::Size<3> sliceSize;
        sliceSize[0] = cvLabelMapToImage->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
        sliceSize[1] = cvLabelMapToImage->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
        sliceSize[2] = 1;

        typedef itk::ImageRegion<3> ImageRegionType;
        ImageRegionType imageSlice;
        imageSlice.SetIndex(sliceOrigin);
        imageSlice.SetSize(sliceSize);

        itk::ImageRegionConstIterator<LImageType> iter(cvLabelMapToImage->GetOutput(), imageSlice);
        for(iter = iter.Begin(); !iter.IsAtEnd(); ++iter) {
            if(iter.Value() != 0) {
                unsigned long labelId = iter.Value();
                typename LImageType::IndexType idx = iter.GetIndex();
                slicesOfVeins[slicesToVisit[i]][labelId][0] += idx[0];
                slicesOfVeins[slicesToVisit[i]][labelId][1] += idx[1];
                slicesOfVeins[slicesToVisit[i]][labelId][2] += idx[2];
                veinPixelNumber[labelId] += 1;
            }
        }
        for(unsigned int j=0; j<veinIds.size(); j++) {
            slicesOfVeins[slicesToVisit[i]][veinIds[j]][0] = (double)slicesOfVeins[slicesToVisit[i]][veinIds[j]][0] / (double)veinPixelNumber[veinIds[j]];
            slicesOfVeins[slicesToVisit[i]][veinIds[j]][1] = (double)slicesOfVeins[slicesToVisit[i]][veinIds[j]][1] / (double)veinPixelNumber[veinIds[j]];
            slicesOfVeins[slicesToVisit[i]][veinIds[j]][2] = (double)slicesOfVeins[slicesToVisit[i]][veinIds[j]][2] / (double)veinPixelNumber[veinIds[j]];

            vtkIdType newVertId = veinGraph->AddVertex();
            points->InsertNextPoint(slicesOfVeins[slicesToVisit[i]][veinIds[j]][0], slicesOfVeins[slicesToVisit[i]][veinIds[j]][1], slicesOfVeins[slicesToVisit[i]][veinIds[j]][2]);
            vertIdToState.insert(std::pair<vtkIdType, int>(newVertId,1));

            if(i==0)
                veinsLastVert.insert(std::pair<unsigned long, vtkIdType>(j, newVertId));
            else {
                veinGraph->AddEdge(veinsLastVert[j], newVertId);
                veinsLastVert[j] = newVertId;
            }

            std::cout << "center of central vein " << j << " in slice " << slicesToVisit[i] << " is [" << slicesOfVeins[slicesToVisit[i]][veinIds[j]][0] << ", " << slicesOfVeins[slicesToVisit[i]][veinIds[j]][1] << ", " << slicesOfVeins[slicesToVisit[i]][veinIds[j]][2] << "]" << std::endl;
            veinPixelNumber[veinIds[j]] = 0;
        }
    }


    veinsLastVert.clear();
    for(unsigned int i=0; i<slicesToVisit.size(); i++) {
        for(unsigned int j=0; j<pvImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); j++) {
            unsigned long veinId = pvImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(j)->GetLabel();
            slicesOfVeins[slicesToVisit[i]][veinId][0] = 0;
            slicesOfVeins[slicesToVisit[i]][veinId][1] = 0;
            slicesOfVeins[slicesToVisit[i]][veinId][2] = 0;
        }
    }

    for(unsigned int i=0; i<slicesToVisit.size(); i++) {
        std::vector<unsigned long> veinIds;
        std::map<unsigned long, unsigned long> veinPixelNumber;
        for(unsigned int j=0; j<pvImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); j++) {
            veinIds.push_back(pvImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(j)->GetLabel());
            veinPixelNumber.insert(std::pair<unsigned long, unsigned long>(pvImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(j)->GetLabel(), 0));
        }
        itk::Index<3> sliceOrigin;
        sliceOrigin[0] = 0;
        sliceOrigin[1] = 0;
        sliceOrigin[2] = slicesToVisit[i];

        itk::Size<3> sliceSize;
        sliceSize[0] = pvLabelMapToImage->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
        sliceSize[1] = pvLabelMapToImage->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
        sliceSize[2] = 1;

        typedef itk::ImageRegion<3> ImageRegionType;
        ImageRegionType imageSlice;
        imageSlice.SetIndex(sliceOrigin);
        imageSlice.SetSize(sliceSize);

        itk::ImageRegionConstIterator<LImageType> iter(pvLabelMapToImage->GetOutput(), imageSlice);
        for(iter = iter.Begin(); !iter.IsAtEnd(); ++iter) {
            if(iter.Value() != 0) {
                unsigned long labelId = iter.Value();
                typename LImageType::IndexType idx = iter.GetIndex();
                slicesOfVeins[slicesToVisit[i]][labelId][0] += idx[0];
                slicesOfVeins[slicesToVisit[i]][labelId][1] += idx[1];
                slicesOfVeins[slicesToVisit[i]][labelId][2] += idx[2];
                veinPixelNumber[labelId] += 1;
            }
        }
        for(unsigned int j=0; j<veinIds.size(); j++) {
            slicesOfVeins[slicesToVisit[i]][veinIds[j]][0] = (double)slicesOfVeins[slicesToVisit[i]][veinIds[j]][0] / (double)veinPixelNumber[veinIds[j]];
            slicesOfVeins[slicesToVisit[i]][veinIds[j]][1] = (double)slicesOfVeins[slicesToVisit[i]][veinIds[j]][1] / (double)veinPixelNumber[veinIds[j]];
            slicesOfVeins[slicesToVisit[i]][veinIds[j]][2] = (double)slicesOfVeins[slicesToVisit[i]][veinIds[j]][2] / (double)veinPixelNumber[veinIds[j]];

            vtkIdType newVertId = veinGraph->AddVertex();
            points->InsertNextPoint(slicesOfVeins[slicesToVisit[i]][veinIds[j]][0], slicesOfVeins[slicesToVisit[i]][veinIds[j]][1], slicesOfVeins[slicesToVisit[i]][veinIds[j]][2]);
            vertIdToState.insert(std::pair<vtkIdType, int>(newVertId,2));

            if(i==0)
                veinsLastVert.insert(std::pair<unsigned long, vtkIdType>(j, newVertId));
            else {
                veinGraph->AddEdge(veinsLastVert[j], newVertId);
                veinsLastVert[j] = newVertId;
            }

            std::cout << "center of portal vein " << j << " in slice " << slicesToVisit[i] << " is [" << slicesOfVeins[slicesToVisit[i]][veinIds[j]][0] << ", " << slicesOfVeins[slicesToVisit[i]][veinIds[j]][1] << ", " << slicesOfVeins[slicesToVisit[i]][veinIds[j]][2] << "]" << std::endl;
            veinPixelNumber[veinIds[j]] = 0;
        }
    }

    for(unsigned int j=0; j<points->GetNumberOfPoints(); j++)
        sinPoints->InsertNextPoint(points->GetPoint(j));

    veinGraph->SetPoints(points);

    int veinPedIdFirst = graphHelper->GetGraph()->GetNumberOfVertices();
    int numVeinVertices = veinGraph->GetNumberOfVertices();

    anno->AddCustomVertexAnnotation(veinGraph, "Network type", vertIdToState, -1);
    anno->SetVertexIDOffset(graphHelper->GetGraph()->GetNumberOfVertices());
    anno->AddPredefinedAnnotations(veinGraph);

    mergeGraphsFilter->ExtendGraph(graphHelper, veinGraph);
    graphHelper->GetGraph()->SetPoints(sinPoints);

    vtkSmartPointer<vtkMutableUndirectedGraph> allInGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    graphHelper->GetGraph()->ToUndirectedGraph(allInGraph);

    vtkIntArray* networkTypeArray = vtkIntArray::SafeDownCast( allInGraph->GetVertexData()->GetArray("Network type") );
    vtkFloatArray* distToCVArray = vtkFloatArray::SafeDownCast( allInGraph->GetVertexData()->GetArray("distance to CV") );
    vtkFloatArray* distToPVArray = vtkFloatArray::SafeDownCast( allInGraph->GetVertexData()->GetArray("distance to PV") );

    std::vector<int> cvPedIds, pvPedIds;
    for(unsigned int i=veinPedIdFirst; i<veinPedIdFirst+numVeinVertices; i++) {
        if(networkTypeArray->GetValue(allInGraph->FindVertex(i)) == 1)
            cvPedIds.push_back(i);
        if(networkTypeArray->GetValue(allInGraph->FindVertex(i)) == 2)
            pvPedIds.push_back(i);
    }

    for(int j=0; j<allInGraph->GetNumberOfVertices(); j++) {
        if(networkTypeArray->GetValue(j) == 0 && allInGraph->GetDegree(j) == 1) {
            if(distToCVArray->GetValue(j) < 20.) {
                double distMin = 10000;
                int distMinVert = -1;

                double* deadEndPt = allInGraph->GetPoint(j);
                double a2 = deadEndPt[0]; double b2 = deadEndPt[1]; double c2 = deadEndPt[2];       //no clue why this is necessary
                for(unsigned int i=0; i<cvPedIds.size(); i++) {
                    double* veinPt = allInGraph->GetPoint(allInGraph->FindVertex(cvPedIds[i]));
                    double a1 = veinPt[0]; double b1 = veinPt[1]; double c1 = veinPt[2];

                    double dist = sqrt(pow(a1-a2, 2) + pow(b1-b2, 2) + pow(c1-c2, 2));
                    if(dist < distMin) {
                        distMin = dist;
                        distMinVert = allInGraph->FindVertex(cvPedIds[i]);
                    }
                }
                allInGraph->AddEdge(j, distMinVert);
            }
            else if(distToPVArray->GetValue(j) < 20.) {
                double distMin = 10000;
                int distMinVert = -1;

                double* deadEndPt = allInGraph->GetPoint(j);
                double a2 = deadEndPt[0]; double b2 = deadEndPt[1]; double c2 = deadEndPt[2];       //no clue why this is necessary
                for(unsigned int i=0; i<pvPedIds.size(); i++) {
                    double* veinPt = allInGraph->GetPoint(allInGraph->FindVertex(pvPedIds[i]));
                    double a1 = veinPt[0]; double b1 = veinPt[1]; double c1 = veinPt[2];

                    double dist = sqrt(pow(a1-a2, 2) + pow(b1-b2, 2) + pow(c1-c2, 2));
                    if(dist < distMin) {
                        distMin = dist;
                        distMinVert = allInGraph->FindVertex(pvPedIds[i]);
                    }
                }
                allInGraph->AddEdge(j, distMinVert);
            }
        }
    }

    SaveGraph("/scratch/friebel/TestData/connectNetworkAndVeins/veinGraph.txt", veinGraph);
    SaveGraph("/scratch/friebel/TestData/connectNetworkAndVeins/sinAndVeinGraph.txt", allInGraph, 0.621, 0.621, 0.54);
}


//Deprecated: not supported by new input -> output scheme
void FileFormatConverter::ComposeGraphFromDresdenFormat()
{
    std::string sinGraphFile = "/home/ls1/friebel/Workspace/Data/Confocal/Kalaidzi/kalaidzi_results/Sinusoids_converted0.txt";
    std::string bileGraphFile = "/home/ls1/friebel/Workspace/Data/Confocal/Kalaidzi/kalaidzi_results/Bile_canaliculi_converted0.txt";
    std::string cveinSegFile = "/home/ls1/friebel/Workspace/Data/Confocal/Kalaidzi/raw_data/vein_central_bin_mir.tif";
    std::string pveinSegFile = "/home/ls1/friebel/Workspace/Data/Confocal/Kalaidzi/raw_data/vein_portalField_bin_mir.tif";

    std::vector< vtkSmartPointer<vtkMutableUndirectedGraph> > sinGraphs;
    ReadGraph(sinGraphFile, sinGraphs);

    std::vector< vtkSmartPointer<vtkMutableUndirectedGraph> > bileGraphs;
    ReadGraph(bileGraphFile, bileGraphs);


    C3DImageType::SpacingType imageSpacing;
    imageSpacing.Fill(0.3);

    C3DReaderType::Pointer cvReader = C3DReaderType::New();
    cvReader->SetFileName(cveinSegFile);
#if (ITK_VERSION_MAJOR >= 4)
    cvReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    cvReader->Update();

    C3DImageType::Pointer cvMask = cvReader->GetOutput();
    cvMask->DisconnectPipeline();
    cvMask->SetSpacing(imageSpacing);

    SignedMaurerDistanceMapImageFilterType::Pointer cvDistMap = SignedMaurerDistanceMapImageFilterType::New();
    cvDistMap->SetInput(cvMask);
    cvDistMap->UseImageSpacingOn();
    cvDistMap->SquaredDistanceOff();
    cvDistMap->Update();

    C3DReaderType::Pointer pvReader = C3DReaderType::New();
    pvReader->SetFileName(pveinSegFile);
#if (ITK_VERSION_MAJOR >= 4)
    pvReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    pvReader->Update();

    C3DImageType::Pointer pvMask = pvReader->GetOutput();
    pvMask->DisconnectPipeline();
    pvMask->SetSpacing(imageSpacing);

    SignedMaurerDistanceMapImageFilterType::Pointer pvDistMap = SignedMaurerDistanceMapImageFilterType::New();
    pvDistMap->SetInput(pvMask);
    pvDistMap->UseImageSpacingOn();
    pvDistMap->SquaredDistanceOff();
    pvDistMap->Update();


    GraphAnnotationHelper *anno = new GraphAnnotationHelper();

    for(unsigned int i=0; i<sinGraphs.size(); i++) {
        vtkSmartPointer<vtkPoints> points = sinGraphs[i]->GetPoints();
        std::map<vtkIdType, float> vertIdToCVDist, vertIdToPVDist;
        std::map<vtkIdType,int> vertexIdToVeinState;
        std::stringstream name;

        for(unsigned int j=0; j<sinGraphs[i]->GetNumberOfVertices(); j++) {
            double* point = points->GetPoint(j);
            itk::Index<3> idx;
            for(unsigned int k=0; k<3; k++)
                idx[k] = point[k];

            vertIdToCVDist.insert(std::pair<vtkIdType, float>(j, cvDistMap->GetOutput()->GetPixel(idx)));
            vertIdToPVDist.insert(std::pair<vtkIdType, float>(j, pvDistMap->GetOutput()->GetPixel(idx)));
            vertexIdToVeinState.insert(std::pair<vtkIdType,int>(j, GraphHandlingHelper::GetVeinConnectednessStatus(sinGraphs[i], j, cvDistMap->GetOutput(), pvDistMap->GetOutput(),
                    30., 30., true, true)));
        }
        anno->EnableEdgeLengthAnnotation(.3, .3, .3);
        anno->AddPredefinedAnnotations(sinGraphs[i]);
        anno->AddCustomVertexAnnotation(sinGraphs[i], "distance to CV", vertIdToCVDist, -1);
        anno->AddCustomVertexAnnotation(sinGraphs[i], "distance to PV", vertIdToPVDist, -1);
        anno->AddCustomVertexAnnotation(sinGraphs[i], "vein state", vertexIdToVeinState, -100);

        name << "/home/ls1/friebel/Workspace/Data/Confocal/Kalaidzi/kalaidzi_results/" << "Sinusoids_converted" << i << ".txt";
        SaveGraph(name.str(), sinGraphs[i]);
    }

    for(unsigned int i=0; i<bileGraphs.size(); i++) {
        vtkSmartPointer<vtkPoints> points = bileGraphs[i]->GetPoints();
        std::map<vtkIdType, float> vertIdToCVDist, vertIdToPVDist;
        std::map<vtkIdType,int> vertexIdToVeinState;
        std::stringstream name;

        for(unsigned int j=0; j<bileGraphs[i]->GetNumberOfVertices(); j++) {
            double* point = points->GetPoint(j);
            itk::Index<3> idx;
            for(unsigned int k=0; k<3; k++)
                idx[k] = point[k];

            vertIdToCVDist.insert(std::pair<vtkIdType, float>(j, cvDistMap->GetOutput()->GetPixel(idx)));
            vertIdToPVDist.insert(std::pair<vtkIdType, float>(j, pvDistMap->GetOutput()->GetPixel(idx)));
            vertexIdToVeinState.insert(std::pair<vtkIdType,int>(j, GraphHandlingHelper::GetVeinConnectednessStatus(bileGraphs[i], j, cvDistMap->GetOutput(), pvDistMap->GetOutput(),
                    -1000., 15., true, true)));
        }
        anno->EnableEdgeLengthAnnotation(.3, .3, .3);
        anno->AddPredefinedAnnotations(bileGraphs[i]);
        anno->AddCustomVertexAnnotation(bileGraphs[i], "distance to CV", vertIdToCVDist, -1);
        anno->AddCustomVertexAnnotation(bileGraphs[i], "distance to PV", vertIdToPVDist, -1);
        anno->AddCustomVertexAnnotation(bileGraphs[i], "vein state", vertexIdToVeinState, -100);

        name << "/home/ls1/friebel/Workspace/Data/Confocal/Kalaidzi/kalaidzi_results/" << "Bile_canaliculi_converted" << i << ".txt";
        SaveGraph(name.str(), bileGraphs[i]);
    }
}


void FileFormatConverter::Update()
{
    ParseParameterContext();

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();


    std::cout << "FileFormatConverter: " << std::endl;
    std::cout << " dir: " << mPath << std::endl;
    std::cout << " filename: " << mFilename << std::endl;
    std::cout << " ext: " << mFilenameExtension << std::endl;


    switch(mInputFileFormat)
    {
    case Graph_Mevis:
    {
        switch(mOutputFileFormat)
        {
        case Graph_TiQuant:
        {
            vtkSmartPointer<vtkMutableDirectedGraph> graph = ReadAndConvertMevisGraphToVTKGraph(mPath+mFilename+mFilenameExtension);
            SaveGraph(mPath+mFilename+"_converted.txt", graph);
        }
        break;
        default:
        {
            std::cout << "FileFormatConverter::Update(): Conversion from " << InputDataName[mInputFileFormat] << " to " << OutputDataName[mOutputFileFormat] << " not known!" << std::endl;
        }
        break;
        }
    }
    break;
    case Graph_Dresden:
    {
        switch(mOutputFileFormat)
        {
        case Graph_TiQuant:
        {
            std::vector< vtkSmartPointer<vtkMutableUndirectedGraph> > graphs = ReadAndConvertDresdenGraphToVTKGraph(mPath+mFilename+mFilenameExtension);
            for(unsigned int i=0; i<graphs.size(); i++) {
                std::stringstream name;
                name << mPath << mFilename << "_converted" << i << ".txt";
                SaveGraph(name.str(), graphs[i]);
            }
        }
        break;
        default:
        {
            std::cout << "FileFormatConverter::Update(): Conversion from " << InputDataName[mInputFileFormat] << " to " << OutputDataName[mOutputFileFormat] << " not known!" << std::endl;
        }
        break;
        }
    }
    break;
    case Triangualation_Dresden:
    {
        switch(mOutputFileFormat)
        {
        case Binary_Image:
        {
            vtkSmartPointer<vtkPolyData> polyData = ReadSTLFile(mPath+mFilename+mFilenameExtension);
            vtkSmartPointer<vtkImageData> vtkImageData = ConvertTriangulationToImageData(polyData);
            C3DImageType::Pointer itkImageData;
            ConvertVtkToItkImageData(vtkImageData, itkImageData);
            SaveImage(mPath+mFilename+"_converted.tif", itkImageData);
        }
        break;
        default:
        {
            std::cout << "FileFormatConverter::Update(): Conversion from " << InputDataName[mInputFileFormat] << " to " << OutputDataName[mOutputFileFormat] << " not known!" << std::endl;
        }
        break;
        }
    }
    break;
    case Bit8_RGB_2D_Image:
    {
        switch(mOutputFileFormat)
        {
        case Bit8_Greyscale_2D_Image:
        {
            RGBC2DImageType::Pointer rgbImage = RGBC2DImageType::New();
            ReadImage(mInputFilename, rgbImage);
            C2DImageType::Pointer scalarImage = C2DImageType::New();
            scalarImage->SetRegions(rgbImage->GetLargestPossibleRegion());
            scalarImage->Allocate();

            ConvertRGBToGreyscaleImageData(rgbImage, scalarImage, 0);
            SaveImage(mPath+mFilename+"_converted.tif", scalarImage);
        }
        break;
        default:
        {
            std::cout << "FileFormatConverter::Update(): Conversion from " << InputDataName[mInputFileFormat] << " to " << OutputDataName[mOutputFileFormat] << " not known!" << std::endl;
        }
        break;
        }
    }
    break;
    case Bit8_RGB_3D_Image:
    {
        switch(mOutputFileFormat)
        {
        case Bit8_Greyscale_3D_Image:
        {
            RGBC3DImageType::Pointer rgbImage = RGBC3DImageType::New();
            RGBC3DImageType::SpacingType spacing;
            spacing.Fill(1.);
            BasePipeline::ReadImage(mInputFilename, rgbImage, spacing);
            C3DImageType::Pointer scalarImage = C3DImageType::New();
            scalarImage->SetRegions(rgbImage->GetLargestPossibleRegion());
            scalarImage->SetSpacing(spacing);
            scalarImage->Allocate();

            ConvertRGBToGreyscaleImageData(rgbImage, scalarImage, 0);
            SaveImage(mPath+mInputFilename+"_converted.tif", scalarImage);
        }
        break;
        default:
        {
            std::cout << "FileFormatConverter::Update(): Conversion from " << InputDataName[mInputFileFormat] << " to " << OutputDataName[mOutputFileFormat] << " not known!" << std::endl;
        }
        break;
        }
    }
    break;
    default:
    {
        std::cout << "FileFormatConverter::Update(): Input format " << InputDataName[mInputFileFormat] << std::endl;
    }
    break;
    }

    WriteLogFile(timeStamp);
}


