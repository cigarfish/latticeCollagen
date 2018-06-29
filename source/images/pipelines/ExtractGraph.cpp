///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ExtractGraph.cpp                                                     //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-12-09                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "ExtractGraph.h"

#include <QDateTime>

#include "itkImageToVTKImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#if (ITK_VERSION_MAJOR >= 4)
#include "itkTIFFImageIOFactory.h"
#endif

#include <vtkGraphReader.h>
#include <vtkGraphWriter.h>
#include <vtkImageFlip.h>
#include <vtkVertexListIterator.h>

#include "../filters/convertFilters/SkeletonImageToGraphFilter.h"
#include "../filters/graphFilters/CollapseIntersectionNodesFilter.h"
#include "../filters/graphFilters/GeometricalThresholdUndirectedGraphFilter.h"
#include "../filters/graphFilters/ResampleUndirectedGraphFilter.h"
#include "../filters/graphFilters/TopologicalPruneGraphFilter.h"
#include "../tools/ImageAnalysisSummaryFileIO.h"
#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"


#define dataForNick 0
#define buildSixNeighbourhoodDataForNick 1
#define pathToDatForNick "/home/ls1/friebel/Workspace/Data/Confocal/TestData/ResampleNetwork/3D/041/spacing_5micron/"
#define radDataForNick 2.5


ExtractGraph::ExtractGraph()
{
	m_testOutput = false;

	m_filenameSave[0] = "van_graph";
	m_filenameSave[1] = "pruned_graph";
	m_filenameSave[2] = "graph";

#if (ITK_VERSION_MAJOR >= 4)
	itk::TIFFImageIOFactory::RegisterOneFactory();
#endif
}


vtkSmartPointer<vtkImageData> ExtractGraph::GetInputImageData()
{
	if(m_entryPoint==0) {

		itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO((m_path + m_filename + m_fileExtension).c_str(), itk::ImageIOFactory::ReadMode);
		if(imageIO->CanReadFile((m_path + m_filename + m_fileExtension).c_str())) {
			imageIO->ReadImageInformation();

			std::cout << "numDimensions: " << imageIO->GetNumberOfDimensions() << std::endl;
			if(imageIO->GetNumberOfDimensions()==2) 		m_3DMode = false;
			else if(imageIO->GetNumberOfDimensions()==3) 	m_3DMode = true;
			else {											std::cout << "Error: Dimension number of input data invalid!" << std::endl; return NULL; }

			if(!m_3DMode) {
				typedef unsigned char																PixelCharType;								//workaround: itk-vtk-connector can not handle images templated over bool pixels
				typedef itk::Image<PixelCharType, dim2>												CScalar2DImageType;							//therefore read the (expected!) bool image and recast it to something that the connector understands

				typedef itk::ImageToVTKImageFilter<CScalar2DImageType>								ITKVTKScalar2DConnectorType;
				typedef itk::RescaleIntensityImageFilter<BScalar2DImageType, CScalar2DImageType> 	RescalIntensity2DFilterType;

				BScalar2DReaderType::Pointer skeletonReader = BScalar2DReaderType::New();
				skeletonReader->SetFileName(m_path + m_filename + m_fileExtension);

				RescalIntensity2DFilterType::Pointer rescaleToCScalarImage = RescalIntensity2DFilterType::New();
				rescaleToCScalarImage->SetInput(skeletonReader->GetOutput());

				ITKVTKScalar2DConnectorType::Pointer itkvtkConnectorInputImage = ITKVTKScalar2DConnectorType::New();
				itkvtkConnectorInputImage->SetInput(rescaleToCScalarImage->GetOutput());

				vtkSmartPointer<vtkImageFlip> flipYInputFilter = vtkSmartPointer<vtkImageFlip>::New();
				flipYInputFilter->SetFilteredAxis(1); 							// flip y axis
				flipYInputFilter->SetInput(itkvtkConnectorInputImage->GetOutput());
				flipYInputFilter->Update();

				return flipYInputFilter->GetOutput();
			}
			else {
				typedef unsigned char																PixelCharType;
				typedef itk::Image<PixelCharType, dim3>												CScalar3DImageType;

				typedef itk::ImageToVTKImageFilter<CScalar3DImageType>								ITKVTKScalar3DConnectorType;
				typedef itk::RescaleIntensityImageFilter<BScalar3DImageType, CScalar3DImageType> 	RescalIntensity3DFilterType;

				BScalar3DReaderType::Pointer skeletonReader = BScalar3DReaderType::New();
				skeletonReader->SetFileName(m_path + m_filename + m_fileExtension);

				RescalIntensity3DFilterType::Pointer rescaleToCScalarImage = RescalIntensity3DFilterType::New();
				rescaleToCScalarImage->SetInput(skeletonReader->GetOutput());

				ITKVTKScalar3DConnectorType::Pointer itkvtkConnectorInputImage = ITKVTKScalar3DConnectorType::New();
				itkvtkConnectorInputImage->SetInput(rescaleToCScalarImage->GetOutput());

				vtkSmartPointer<vtkImageFlip> flipYInputFilter = vtkSmartPointer<vtkImageFlip>::New();
				flipYInputFilter->SetFilteredAxis(1); 							// flip y axis
				flipYInputFilter->SetInput(itkvtkConnectorInputImage->GetOutput());
				flipYInputFilter->Update();

				return flipYInputFilter->GetOutput();
			}
		}
		else {
			std::cout << "Error: No valid file for PruneGraph::Update! (specified filename=" << m_filename << ")" << std::endl;
		}
	}

	return NULL;
}


void ExtractGraph::PrintBasicGraphInfo(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, std::string graphName)
{
    int count_total = 0, count_end=0, count_regular=0, count_cross=0;
    for(unsigned int i=0; i<graphs.size(); i++) {
        count_total += graphs[i]->GetNumberOfVertices();

        for(vtkIdType j=0; j<graphs[i]->GetNumberOfVertices(); j++) {
            if(graphs[i]->GetDegree(j)==1)      count_end++;
            else if(graphs[i]->GetDegree(j)==2) count_regular++;
            else                                count_cross++;
        }
    }
    std::cout << "number vertices in " << graphName << ": " << count_total << std::endl;
    std::cout << "with " << count_end << " end nodes, "<< count_regular << " regular nodes and " << count_cross << " cross nodes." << std::endl;
}


void ExtractGraph::SaveGraphs(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, int filenameIndex)
{
    for(unsigned int i=0; i<graphs.size(); i++) {
        std::stringstream name;
        name << m_path << m_savePrefix << m_filenameSave[filenameIndex] << i << ".txt";
		if(filenameIndex == 2 && i==0)  m_nameOfFirstGraph = QString::fromStdString(name.str());
//        std::cout << "save " << m_savePrefix << m_filenameSave[filenameIndex] << " to " << name.str() << std::endl;

        vtkSmartPointer<vtkGraphWriter> writer = vtkSmartPointer<vtkGraphWriter>::New();
        writer->SetFileName(name.str().c_str());
        writer->SetInput(graphs[i]);
        writer->Update();
    }
}


void ExtractGraph::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    file.open((m_path + "log_" + m_savePrefix + "graph" + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-extract-graph-------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << "file = " << m_path + m_filename + m_fileExtension << "\n";

    if(m_entryPoint == 0) {
        file << "voxel size (" << mSpacing[0] << ", " << mSpacing[1] <<  ", " << mSpacing[2] << ")\n";
        file << "step 1: resampling factor = " << m_resamplingFactor << " with maximal resampling distance = " << m_maxResamplingDist << "\n";
        file << "step 2: prune dead ends shorter than = " << m_deadEndPruneThreshold << "\n";
        file << "step 3: collapse intersection nodes distant less than = " << m_collapseThreshold << "\n";
        file << "step 4: remove regular nodes that introduce deviation from straight line smaller than " << m_pruningAngle << " degrees \n";
    }

    file << "End-extract-graph---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


//TODO: deprecated,should be moved to FileFormatConverter and should use AnalyzeBileNetwork::SetRadiusMeasurement & SetVeinConnectednessMeasurement
void ExtractGraph::WriteFinalGraphsInNicksXMLFormat()
{
    typedef itk::Image<unsigned char, 3>                CScalarVoImageType;
    typedef itk::ImageFileReader<CScalarVoImageType>    ScalarVoReaderType;

    CScalarVoImageType::Pointer deadEndMasksCV;
    CScalarVoImageType::Pointer deadEndMasksPV;

    std::string full_filename, path, filename, ext;

    full_filename = *(std::string*)(m_paramContext->findParameter("Specify central vein mask file", 0)->dataPointer());
    bool f0 = FilenameParser::ParseFilename(full_filename, path, filename, ext);

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO((path + filename + ext).c_str(), itk::ImageIOFactory::ReadMode);
    imageIO->ReadImageInformation();

    ScalarVoReaderType::Pointer reader = ScalarVoReaderType::New();
    reader->SetFileName(path + filename + ext);
    reader->ReleaseDataBeforeUpdateFlagOn();
    reader->Update();

    CScalarVoImageType::Pointer maskImage = reader->GetOutput();
    maskImage->DisconnectPipeline();

    deadEndMasksCV = maskImage.GetPointer();

    full_filename = *(std::string*)(m_paramContext->findParameter("Specify portal vein mask file", 0)->dataPointer());
    bool f1 = FilenameParser::ParseFilename(full_filename, path, filename, ext);

    itk::ImageIOBase::Pointer imageIO1 = itk::ImageIOFactory::CreateImageIO((path + filename + ext).c_str(), itk::ImageIOFactory::ReadMode);
    imageIO1->ReadImageInformation();

    ScalarVoReaderType::Pointer reader1 = ScalarVoReaderType::New();
    reader1->SetFileName(path + filename + ext);
    reader1->ReleaseDataBeforeUpdateFlagOn();
    reader1->Update();

    CScalarVoImageType::Pointer maskImage1 = reader1->GetOutput();
    maskImage1->DisconnectPipeline();

    deadEndMasksPV = maskImage1.GetPointer();


    for(unsigned int i=0; i<m_prunedGraphs.size(); i++) {
        std::fstream file;

        std::stringstream filename;
        filename << m_path << m_savePrefix << m_filenameSave[2] << i << ".xml";

        file.open(filename.str().c_str(), std::ios::out);

        file << "<?xml version=\"1.0\"?>\n";
        file << "<Model>\n";
		file << "<VesselGraph dimensions=\"" << m_num_dim << "\" x=\"" << maskImage->GetLargestPossibleRegion().GetSize()[0] << "\" y=\"" << maskImage->GetLargestPossibleRegion().GetSize()[1] << "\" z=\"" << maskImage->GetLargestPossibleRegion().GetSize()[2] << "\">\n";

        for(vtkIdType v=0; v<m_prunedGraphs[i]->GetNumberOfVertices(); ++v) {
            double pos[3];
            m_prunedGraphs[i]->GetPoint(v, pos);

            if(m_prunedGraphs[i]->GetDegree(v)==1) {
                CScalarVoImageType::IndexType idx;
                idx[0] = m_prunedGraphs[i]->GetPoint(v)[0];
                idx[1] = m_prunedGraphs[i]->GetPoint(v)[1];
                idx[2] = m_prunedGraphs[i]->GetPoint(v)[2];

                if(deadEndMasksCV->GetPixel(idx)==255) {    //in vicninity of cv
                    std::cout << "so this guy is near a cv" << std::endl;
                    file << "<Node id=\"" << v << "\" x=\"" << pos[0]*mSpacing[0]+mSpacing[0]/2 << "\" y=\"" << pos[1]*mSpacing[1]+mSpacing[1]/2 << "\" z=\"" << pos[2]*mSpacing[2]+mSpacing[2]/2 << "\" type=\"" << 1 << "\" pressure=\"" << 0.0 << "\"/>\n";
                }
                else if(deadEndMasksPV->GetPixel(idx)==255) {       //in vicninity of pv
                    std::cout << "so this guy is near a pv" << std::endl;
                    file << "<Node id=\"" << v << "\" x=\"" << pos[0]*mSpacing[0]+mSpacing[0]/2 << "\" y=\"" << pos[1]*mSpacing[1]+mSpacing[1]/2 << "\" z=\"" << pos[2]*mSpacing[2]+mSpacing[2]/2 << "\" type=\"" << 1 << "\" pressure=\"" << 10.0 << "\"/>\n";
                }
                else   //not near a vein
                    file << "<Node id=\"" << v << "\" x=\"" << pos[0]*mSpacing[0]+mSpacing[0]/2 << "\" y=\"" << pos[1]*mSpacing[1]+mSpacing[1]/2 << "\" z=\"" << pos[2]*mSpacing[2]+mSpacing[2]/2 << "\" type=\"" << 1 << "\" pressure=\"\"/>\n";
            }
            else
                file << "<Node id=\"" << v << "\" x=\"" << pos[0]*mSpacing[0]+mSpacing[0]/2 << "\" y=\"" << pos[1]*mSpacing[1]+mSpacing[1]/2 << "\" z=\"" << pos[2]*mSpacing[2]+mSpacing[2]/2 << "\" type=\"" << 2 << "\" pressure=\"\"/>\n";
        }

        for(vtkIdType e=0; e<m_prunedGraphs[i]->GetNumberOfEdges(); ++e)
            file << "<Segment id=\"" << e << "\" node1=\"" << m_prunedGraphs[i]->GetSourceVertex(e) << "\" node2=\"" << m_prunedGraphs[i]->GetTargetVertex(e) << "\" radiusStatic=\"" << radDataForNick << "\"/>\n";


        file << "</VesselGraph>\n";
        file << "</Model>\n";
        file.close();
    }
}


void ExtractGraph::WriteDataSetSummary()
{
    if(m_entryPoint==0) {
        std::string networkType = ( (CSParameterChoice*)(m_paramContext->findParameter("Network type", 0)->dataPointer()) )->currentString();
        if( networkType.compare("Bile Canaliculi Network")==0 )
			ImageAnalysisSummaryFileIO::AddEntry(BileGraph, m_path, m_nameOfFirstGraph.toStdString());
        else if( networkType.compare("Sinusoidal Network")==0 )
            ImageAnalysisSummaryFileIO::AddEntry(SinusoidGraph, m_path, m_nameOfFirstGraph.toStdString());
    }
}


void ExtractGraph::ParseParameterContext()
{
    if(m_paramContext->findContext("Extract and Analyze Graph",0)==NULL) {
        std::cout << "Error: ExtractGraph: Invalid parameter context: " << std::endl;
        return;
    }

    std::string readMode = ( (CSParameterChoice*)(m_paramContext->findContext("Reader", 0)->findParameter("Input type", 0)->dataPointer()) )->currentString();
    if( readMode.compare("Skeleton Image")==0 ) {
        m_entryPoint = 0;
        
		m_fullFilename = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Skeleton image", 0)->dataPointer()) );
		m_infoFullFilename.setFile(m_fullFilename);

		if(!m_infoFullFilename.exists())
			throw std::string("Please specify a skeleton file");

        m_networkFilename = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Specify network segmentation file", 0)->dataPointer()) );
		m_infoFullFilename.setFile(m_networkFilename);

		if(!m_infoFullFilename.exists())
			throw std::string("Please specify the network segmentation mask");
    }
    else if( readMode.compare("Final Graph")==0 ) {
        m_entryPoint = 1;
        m_fullFilename = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Graph file", 0)->dataPointer()) );
        m_infoFullFilename.setFile(m_fullFilename);

		if(!m_infoFullFilename.exists() || m_infoFullFilename.suffix()!=QString("txt"))
			throw std::string("Please specify a graph file (*.txt)");

		m_nameOfFirstGraph = m_fullFilename;
    }

    m_infoFullFilename.setFile(m_fullFilename);
	m_path = (m_infoFullFilename.path() + QString("/")).toStdString();
	m_filename = m_infoFullFilename.baseName().toStdString();
	m_fileExtension = (QString(".") + m_infoFullFilename.suffix()).toStdString();

    m_resamplingFactor = *(int*)(m_paramContext->findParameter("Resampling factor", 0)->dataPointer());
    m_maxResamplingDist = *(double*)(m_paramContext->findParameter("Maximal resampling distance", 0)->dataPointer());
    m_deadEndPruneThreshold = *(double*)(m_paramContext->findParameter("Length threshold", 0)->dataPointer());
    m_collapseThreshold = *(double*)(m_paramContext->findParameter("Distance threshold", 0)->dataPointer());
    m_pruningAngle = *(double*)(m_paramContext->findParameter("Drop angle", 0)->dataPointer());


    std::string networkType = ( (CSParameterChoice*)(m_paramContext->findParameter("Network type", 0)->dataPointer()) )->currentString();
#ifndef CS_TI_QUANT_ONLY
    m_deadEndPruneAll = *(bool*)(m_paramContext->findParameter("Remove all dead ends (regardless of length)", 0)->dataPointer());

    m_cvFilename = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Specify central vein segmentation file", 0)->dataPointer()) );
    m_infoFullFilename.setFile(m_cvFilename);
    mHasCV = m_infoFullFilename.exists();

    m_pvFilename = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Specify portal vein segmentation file", 0)->dataPointer()) );
    m_infoFullFilename.setFile(m_pvFilename);
    mHasPV = m_infoFullFilename.exists();

    if( networkType.compare("Bile Canaliculi Network")==0 ) {
        m_maxCVeinConnectednessDistance = -1000;
        m_maxPVeinConnectednessDistance = 2.5;  //micron
    }
    else if( networkType.compare("Sinusoidal Network")==0 ) {
        m_maxCVeinConnectednessDistance = 10.0;  //micron
        m_maxPVeinConnectednessDistance = 10.0;  //micron
    }
#else
    m_deadEndPruneAll = false;
#endif

    m_defaultRadius = *(double*)(m_paramContext->findParameter("Assumed average radius of branch", 0)->dataPointer());

    mSpacing[0] = *(double*)(m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    mSpacing[1] = *(double*)(m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    mSpacing[2] = *(double*)(m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    std::string saveMode = ( (CSParameterChoice*)(m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    m_savePrefix = *(std::string*)(m_paramContext->findParameter("Output file prefix", 0)->dataPointer());
    if( saveMode.compare("Save everything")==0 )
        m_saveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        m_saveEverything = 0;
}


void ExtractGraph::Update()
{
	ParseParameterContext();

	std::string timeStamp;

	QDateTime time = QDateTime::currentDateTime();
	timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

	std::cout << "Extract graph: " << std::endl;
	std::cout << " dir: " << m_path << std::endl;
	std::cout << " file: " << m_filename << std::endl;
	std::cout << " ext: " << m_fileExtension << std::endl;


	if(m_entryPoint==0) {
		m_num_dim = GetNumberOfDimensions(m_path + m_filename + m_fileExtension);

		if(m_num_dim != 0) {
			std::cout << "numDimensions: " << m_num_dim << std::endl;

			if(m_num_dim==2) {
				m_3DMode = false;
				throw std::string("Please specify a network skeleton file with three dimensions.");
			}
			else if(m_num_dim==3) 	
				m_3DMode = true;
			else					
				throw std::string("Please specify a network skeleton file with three dimensions.");

			if(!m_3DMode) {
				BScalar2DReaderType::Pointer skeletonReader = BScalar2DReaderType::New();
				skeletonReader->SetFileName(m_path + m_filename + m_fileExtension);
				skeletonReader->Update();

				m_graphs = SkeletonImageToGraphFilter::skeletonToGraph2D(skeletonReader->GetOutput());
			}
			else {
				BScalar3DReaderType::Pointer skeletonReader = BScalar3DReaderType::New();
				skeletonReader->SetFileName(m_path + m_filename + m_fileExtension);
				skeletonReader->Update();

#if (!dataForNick)
				m_graphs = SkeletonImageToGraphFilter::skeletonToGraph3D(skeletonReader->GetOutput(), true);
#else
#if (buildSixNeighbourhoodDataForNick)
				m_graphs = SkeletonImageToGraphFilter::skeletonToGraph3D(skeletonReader->GetOutput(), false);
#else
				m_graphs = SkeletonImageToGraphFilter::skeletonToGraph3D(skeletonReader->GetOutput(), true);
#endif
#endif
			}
			if(m_testOutput) std::cout << "number graphs = "<< m_graphs.size() << std::endl;

			GraphAnnotationHelper anno1;
			anno1.EnableVertexIDAnnotation();
			anno1.EnableVertexDegreeAnnotation();
			anno1.EnableEdgeIDAnnotation();
			anno1.EnableEdgeLengthAnnotation(mSpacing[0], mSpacing[1], mSpacing[2]);
			for(unsigned int i=0; i<m_graphs.size(); i++)   anno1.AddPredefinedAnnotations(m_graphs[i]);

			if(m_saveEverything)    SaveGraphs(m_graphs, 0);
			PrintBasicGraphInfo(m_graphs, "vanilla graphs");
			//------------------------------------------------------------------------------..........................--------------------------------------------------------------------------------

			for(unsigned int i=0; i<m_graphs.size(); i++) {
				vtkSmartPointer<ResampleUndirectedGraphFilter> resampleGraph = vtkSmartPointer<ResampleUndirectedGraphFilter>::New();
				resampleGraph->SetInput(m_graphs[i]);
				resampleGraph->SetResamplingFactor(m_resamplingFactor);
				resampleGraph->SetMaxResamplingDist(m_maxResamplingDist);
				resampleGraph->SetTestOutput(m_testOutput);
				resampleGraph->Update();

				vtkSmartPointer<vtkUndirectedGraph> resampledUndirectedGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
				resampleGraph->GetOutput()->ToUndirectedGraph(resampledUndirectedGraph);
				if(resampledUndirectedGraph->GetNumberOfVertices()>0)
				    m_resampledGraphs.push_back(resampledUndirectedGraph);
			}

			GraphAnnotationHelper anno2;
			anno2.EnableEdgeIDAnnotation();
			anno2.EnableEdgeLengthAnnotation(mSpacing[0], mSpacing[1], mSpacing[2]);
			for(unsigned int i=0; i<m_resampledGraphs.size(); i++)  anno2.AddPredefinedAnnotations(m_resampledGraphs[i]);

			PrintBasicGraphInfo(m_resampledGraphs, "resampled graphs");
			//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

			if(GetNumberOfDimensions(m_cvFilename.toStdString()) != 3)
			    throw std::string("Please specify a central vein segmentation file with three dimensions.");
			if(GetNumberOfDimensions(m_pvFilename.toStdString()) != 3)
			    throw std::string("Please specify a portal vein segmentation file with three dimensions.");

			CScalarVoImageType::Pointer cvBin, pvBin;
			if(mHasCV) {
			    ReadImage(m_cvFilename.toStdString(), cvBin, mSpacing);
			    BuildNonZeroDistanceMap(cvBin, mp_cvDistMap, m_maxCVeinConnectednessDistance, mSpacing, true);
			}
			else
			    BasePipeline::CreateImage(mp_cvDistMap, mp_networkBin->GetLargestPossibleRegion(), m_maxCVeinConnectednessDistance, mSpacing);
			if(mHasPV) {
			    ReadImage(m_pvFilename.toStdString(), pvBin, mSpacing);
			    BuildNonZeroDistanceMap(pvBin, mp_pvDistMap, m_maxPVeinConnectednessDistance, mSpacing, true);
			}
			else
			    BasePipeline::CreateImage(mp_pvDistMap, mp_networkBin->GetLargestPossibleRegion(), m_maxPVeinConnectednessDistance, mSpacing);

			for(unsigned int i=0; i<m_resampledGraphs.size(); i++) {
				vtkSmartPointer<TopologicalPruneGraphFilter> pruneGraphFilter = vtkSmartPointer<TopologicalPruneGraphFilter>::New();
				pruneGraphFilter->SetInput(m_resampledGraphs[i]);
				pruneGraphFilter->SetDeadEndLengthThreshold(m_deadEndPruneThreshold);
				pruneGraphFilter->SetPruneAll(m_deadEndPruneAll);
				pruneGraphFilter->SetVeinDistanceMaps(mp_cvDistMap, mp_pvDistMap);
				pruneGraphFilter->SetImageSpacing(mSpacing);
				pruneGraphFilter->SetVeinConnectednessDistances(m_maxCVeinConnectednessDistance, m_maxPVeinConnectednessDistance);
				pruneGraphFilter->SetAssumedStandardRadius(m_defaultRadius);
				pruneGraphFilter->SetTestOutput(m_testOutput);
				pruneGraphFilter->Update();

				vtkSmartPointer<vtkUndirectedGraph> prunedGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
				pruneGraphFilter->GetOutput()->ToUndirectedGraph(prunedGraph);
				if(prunedGraph->GetNumberOfVertices()>0)
				    m_prunedGraphs.push_back(prunedGraph);
			}
			GraphAnnotationHelper anno2a;
			anno2a.EnableEdgeIDAnnotation();
			anno2a.EnableEdgeLengthAnnotation(mSpacing[0], mSpacing[1], mSpacing[2]);
			for(unsigned int i=0; i<m_prunedGraphs.size(); i++)  anno2a.AddPredefinedAnnotations(m_prunedGraphs[i]);

			PrintBasicGraphInfo(m_prunedGraphs, "pruned graphs");

#if (dataForNick)
			WriteFinalGraphsInNicksXMLFormat();
			SaveGraphs(m_prunedGraphs, 1);
#endif

			if(GetNumberOfDimensions(m_networkFilename.toStdString()) != 3)
				throw std::string("Please specify a network segmentation file with three dimensions.");

			ReadImage(m_networkFilename.toStdString(), mp_networkBin, mSpacing);
			BuildDistanceMap(mp_networkBin, mp_networkDistMap, mSpacing, false);

			for(unsigned int i=0; i<m_prunedGraphs.size(); i++) {
				vtkSmartPointer<CollapseIntersectionNodesFilter> collapseGraphFilter = vtkSmartPointer<CollapseIntersectionNodesFilter>::New();
				collapseGraphFilter->SetInput(m_prunedGraphs[i]);
				collapseGraphFilter->SetIsecDistanceThreshold(m_collapseThreshold);
				collapseGraphFilter->SetNetworkDistanceMap(mp_networkDistMap);
				collapseGraphFilter->SetImageSpacing(mSpacing);
				collapseGraphFilter->SetTestOutput(m_testOutput);
				collapseGraphFilter->Update();

				vtkSmartPointer<vtkUndirectedGraph> collapsedGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
				collapseGraphFilter->GetOutput()->ToUndirectedGraph(collapsedGraph);
				if(collapsedGraph->GetNumberOfVertices()>0)
				    m_collapsedGraphs.push_back(collapsedGraph);
			}

			GraphAnnotationHelper anno4;
			anno4.EnableVertexIDAnnotation();
			anno4.EnableVertexDegreeAnnotation();
			anno4.EnableEdgeIDAnnotation();
			anno4.EnableEdgeLengthAnnotation(mSpacing[0], mSpacing[1], mSpacing[2]);
			for(unsigned int i=0; i<m_collapsedGraphs.size(); i++)  anno4.AddPredefinedAnnotations(m_collapsedGraphs[i]);

			PrintBasicGraphInfo(m_collapsedGraphs, "collapsed graphs");
			//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

			for(unsigned int i=0; i<m_collapsedGraphs.size(); i++) {
				vtkSmartPointer<GeometricalThresholdUndirectedGraphFilter> geomThresholdGraph = vtkSmartPointer<GeometricalThresholdUndirectedGraphFilter>::New();
				geomThresholdGraph->SetInput(m_collapsedGraphs[i]);
				geomThresholdGraph->SetAngleThreshold(m_pruningAngle);
				geomThresholdGraph->SetTestOutput(m_testOutput);
				geomThresholdGraph->SetWithAngleAnnotation(true);
				geomThresholdGraph->Update();

				vtkSmartPointer<vtkUndirectedGraph> geomThresholdedUndirectedGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
				geomThresholdGraph->GetOutput()->ToUndirectedGraph(geomThresholdedUndirectedGraph);
				if(geomThresholdedUndirectedGraph->GetNumberOfVertices()>0)
				    m_geomThresholdedGraphs.push_back(geomThresholdedUndirectedGraph);
			}

			GraphAnnotationHelper anno3;
			anno3.EnableVertexIDAnnotation();
			anno3.EnableVertexDegreeAnnotation();
			anno3.EnableEdgeIDAnnotation();
			anno3.EnableEdgeLengthAnnotation(mSpacing[0], mSpacing[1], mSpacing[2]);
			for(unsigned int i=0; i<m_geomThresholdedGraphs.size(); i++)    anno3.AddPredefinedAnnotations(m_geomThresholdedGraphs[i]);

			SaveGraphs(m_geomThresholdedGraphs, 2);
			PrintBasicGraphInfo(m_geomThresholdedGraphs, "geometric thresholded graphs");
		}
		else
			throw std::string("Please specify a 2D or 3D skeleton file.");

	    WriteLogFile(timeStamp);
	    WriteDataSetSummary();
	}
	if(m_entryPoint==1) {
		vtkSmartPointer<vtkGraphReader> reader = vtkSmartPointer<vtkGraphReader>::New();
		m_3DMode = true;                                                                //TODO: find out if 3d mode appropriated

		std::string prefix;
		int startIndex = 0;
		bool isSeq  = FilenameParser::IsSequence(m_filename, prefix, startIndex);

		std::cout << "Load graphs from disk under " << m_path << std::endl;
		int i=startIndex;
		do {
		    std::stringstream name;
		    if(isSeq)   name << m_path << prefix << i << m_fileExtension;
		    else        name << m_path << prefix << m_fileExtension;
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

		        vtkSmartPointer<vtkUndirectedGraph> undirectedGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
		        reader->GetOutput()->ToUndirectedGraph(undirectedGraph);

		        m_geomThresholdedGraphs.push_back(undirectedGraph);

		        reader->CloseVTKFile();
		    }
		    else
		        break;

		    i++;
		} while(isSeq);
	}
}

