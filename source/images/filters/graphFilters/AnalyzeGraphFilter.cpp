///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeGraphFilter.cpp                                               //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-01-17                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "AnalyzeGraphFilter.h"

#include <math.h>

#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>

#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIOFactory.h>
#include <itkTIFFImageIO.h>
#endif // (ITK_VERSION_MAJOR >= 4)

#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkEdgeListIterator.h>
#include <vtkFloatArray.h>
#include <vtkGraphWriter.h>
#include <vtkIntArray.h>
#include <vtkVertexListIterator.h>

#include "../tools/GraphAnnotationHelper.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"


AnalyzeGraphFilter::AnalyzeGraphFilter()
{
    m_dataSetID = "";
    m_networkDataSetName = "";

    m_withRadiusMeasurement = false;

    m_saveGraphsWithArrays = false;
    m_graphDataSetPath = "";

	m_vox_spacing[0] = 1;
	m_vox_spacing[1] = 1;
	m_vox_spacing[2] = 1;
	m_voxel_volume = 1;
	m_dataset_volume = 1;
	m_effective_dataset_volume = 1;

    m_rad = 1;

    m_withDeadEndMasking = false;

	m_num_nodes = 0;
	m_num_intersec_nodes = 0;
	m_num_reg_nodes = 0;
	m_num_total_dead_end_nodes = 0;
	m_num_within_dead_end_nodes = 0;
	m_num_zero_nodes = 0;
	m_num_edges = 0;

	m_mean_edges_per_intersec_node = 0;
	m_vari_edges_per_intersec_node = 0;
	m_total_edge_length = 0;
	m_mean_radius = 0;
	m_num_neg_radius = 0;

    m_network_vol_calc = 0;
    m_network_vol_calc_per_vol = 0;
    m_network_vol_count = 0;
    m_network_vol_count_per_vol = 0;

    m_mean_isecBranch_minAngle = 0;
    m_mean_isecBranch_minAvgAngle = 0;

#if (ITK_VERSION_MAJOR >= 4)
    itk::TIFFImageIOFactory::RegisterOneFactory();
#endif // (ITK_VERSION_MAJOR >= 4)
}


AnalyzeGraphFilter::~AnalyzeGraphFilter()
{
	// TODO Auto-generated destructor stub
}


void AnalyzeGraphFilter::ParseParameterContext()
{
    if(m_paramContext->findContext("Extract and Analyze Graph",0)==NULL) {
        std::cout << "Error: AnalyzeGraph: Invalid parameter context: " << std::endl;
        return;
    }

    m_dataSetID = *(std::string*)(m_paramContext->findParameter("Specify data set name", 0)->dataPointer());

    m_networkDataSetFullFilename = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Specify network segmentation file", 0)->dataPointer()) );
	m_infoFile.setFile(m_networkDataSetFullFilename);

	if(!m_infoFile.exists())
		throw std::string("Please specify a network segmentation file.");
	else
		m_withRadiusMeasurement = true;

	m_networkDataSetPath = (m_infoFile.path() + QString("/")).toStdString();
    m_networkDataSetName = m_infoFile.baseName().toStdString();
    m_networkDataSetFileExtension = (QString(".") + m_infoFile.suffix()).toStdString();

    m_saveAnalysisFile = *(bool*)(m_paramContext->findParameter("Write results to analysis file", 0)->dataPointer());

    m_vox_spacing[0] = *(double*)(m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    m_vox_spacing[1] = *(double*)(m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    m_vox_spacing[2] = *(double*)(m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    m_rad = *(double*)(m_paramContext->findParameter("Assumed average radius of branch", 0)->dataPointer());

//    m_withDeadEndMasking = *(bool*)(m_paramContext->findParameter("Mask dead ends", 0)->dataPointer());
//
//    if(m_withDeadEndMasking) {
//        AddDeadEndMaskImage( *(std::string*)(m_paramContext->findParameter("Specify central vein mask file", 0)->dataPointer()) );
//        AddDeadEndMaskImage( *(std::string*)(m_paramContext->findParameter("Specify portal vein mask file", 0)->dataPointer()) );
//    }
}


void AnalyzeGraphFilter::AddDeadEndMaskImage(std::string full_filename)
{
    std::string path, filename, ext;

    bool f0 = FilenameParser::ParseFilename(full_filename, path, filename, ext);

    if(f0) {
        ifstream file((path + filename + ext).c_str());
        if(file.good())
        {
            itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO((path + filename + ext).c_str(), itk::ImageIOFactory::ReadMode);
            imageIO->ReadImageInformation();

            CScalarVoImageType::Pointer maskImage;
            BasePipeline<3>::ReadImage(path + filename + ext, maskImage, m_vox_spacing);

            m_deadEndMasks.push_back(maskImage.GetPointer());
        }
    }
}


double AnalyzeGraphFilter::MeasureDistance(vtkGraph *input, vtkIdType vertex, FScalarVoImageType::Pointer image)
{
    double pos[3];
    input->GetPoint(vertex, pos);

    CScalarVoImageType::IndexType idx;
    idx[0] = pos[0];
    idx[1] = pos[1];
    idx[2] = pos[2];

    return image->GetPixel(idx);
}


void AnalyzeGraphFilter::SaveGraphs(std::string path, std::string filename , std::string ext)
{
    for(unsigned int i=0; i<m_graphs.size(); i++) {
        std::stringstream name;
        name << path << filename << i << ext;
        if(i==0)    std::cout << "save analysis graph to " << name.str() << std::endl;

        vtkSmartPointer<vtkGraphWriter> writer = vtkSmartPointer<vtkGraphWriter>::New();
        writer->SetFileName(name.str().c_str());
        writer->SetInput(m_graphs[i]);
        writer->Update();
    }
}


void AnalyzeGraphFilter::CollectBasicImageInformation()
{
	itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO((m_networkDataSetPath + m_networkDataSetName + m_networkDataSetFileExtension).c_str(), itk::ImageIOFactory::ReadMode);
	if(imageIO->CanReadFile((m_networkDataSetPath + m_networkDataSetName + m_networkDataSetFileExtension).c_str())) {
		imageIO->ReadImageInformation();

		m_num_dim = imageIO->GetNumberOfDimensions();

		m_dim[0] = imageIO->GetDimensions(0);
		m_dim[1] = imageIO->GetDimensions(1);
		m_dim[2] = imageIO->GetDimensions(2);

		m_voxel_volume = m_vox_spacing[0] * m_vox_spacing[1] * m_vox_spacing[2] / 1000000000;       //in mm³

		m_dataset_volume *= ((long double)m_dim[0]*m_vox_spacing[0]);
		m_dataset_volume *= ((long double)m_dim[1]*m_vox_spacing[1]);
		m_dataset_volume *= ((long double)m_dim[2]*m_vox_spacing[2]);
		m_dataset_volume /= 1000000000;									                    //in mm³ (=1000³)

		m_effective_dataset_volume = m_dataset_volume;
	}
}


void AnalyzeGraphFilter::CollectBasicNodeEdgeInformation()
{
	long int num_rim_dead_end_nodes = 0;
	GraphAnnotationHelper anno;

	for(unsigned int i=0; i<m_graphs.size(); i++) {
		m_num_nodes += m_graphs[i]->GetNumberOfVertices();
		m_num_edges += m_graphs[i]->GetNumberOfEdges();

		vtkIntArray* degreeArray = vtkIntArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray(anno.GetArrayNameVertexDegree().c_str()) );
		vtkFloatArray* radiusArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("radius") );
		vtkFloatArray* minAngleArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("minimal branching angle (absolute)") );
		vtkFloatArray* minAvgAngleArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("minimal branching angle (averaged)") );
		vtkIntArray* deadEndStateArray = vtkIntArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("dead end branches") );

		for(int j=0; j<m_graphs[i]->GetNumberOfVertices(); j++) {
			if(degreeArray->GetValue(j) == 0) 			m_num_zero_nodes++;
			else if(degreeArray->GetValue(j) == 2) 		m_num_reg_nodes++;
			else if(degreeArray->GetValue(j) > 2) {
				m_num_intersec_nodes++;
				m_mean_edges_per_intersec_node += degreeArray->GetValue(j);
				m_mean_isecBranch_minAngle += minAngleArray->GetValue(j);
				m_mean_isecBranch_minAvgAngle += minAvgAngleArray->GetValue(j);
			}
			else if(degreeArray->GetValue(j) == 1) {
				m_num_total_dead_end_nodes++;

				if(deadEndStateArray->GetValue(j) != 1)
				    num_rim_dead_end_nodes++;
			}
			m_mean_radius += radiusArray->GetValue(j);
			if(radiusArray->GetValue(j)<0) m_num_neg_radius++;
		}
	}
	m_num_within_dead_end_nodes				= m_num_total_dead_end_nodes - num_rim_dead_end_nodes;
	m_num_nodes_per_vol 					= (long double)m_num_nodes / (long double)m_effective_dataset_volume;
	m_num_intersec_nodes_per_vol 			= (long double)m_num_intersec_nodes / (long double)m_effective_dataset_volume;
	m_num_within_dead_end_nodes_per_vol 	= (long double)m_num_within_dead_end_nodes / (long double)m_effective_dataset_volume;
	m_mean_edges_per_intersec_node 			= m_mean_edges_per_intersec_node / (long double)m_num_intersec_nodes;
	m_mean_isecBranch_minAngle              = m_mean_isecBranch_minAngle / (long double)m_num_intersec_nodes;
	m_mean_isecBranch_minAvgAngle           = m_mean_isecBranch_minAvgAngle / (long double)m_num_intersec_nodes;
	m_mean_radius                           = m_mean_radius / (long double)m_num_nodes;


	for(unsigned int i=0; i<m_graphs.size(); i++) {
		vtkIntArray* degreeArray = vtkIntArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray(anno.GetArrayNameVertexDegree().c_str()) );
		vtkDoubleArray* edgeLengthArray = vtkDoubleArray::SafeDownCast( m_graphs[i]->GetEdgeData()->GetArray(anno.GetArrayNameEdgeLength().c_str()) );

		for(int j=0; j<m_graphs[i]->GetNumberOfVertices(); j++)
			if(degreeArray->GetValue(j) > 2)
				m_vari_edges_per_intersec_node += pow((degreeArray->GetValue(j)-m_mean_edges_per_intersec_node), 2);

		vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
		m_graphs[i]->GetEdges(it);
		while(it->HasNext()) {
			vtkEdgeType e = it->Next();
			m_total_edge_length += edgeLengthArray->GetValue(e.Id);
		}
	}
	m_total_edge_length /= 1000;        //in mm

	m_vari_edges_per_intersec_node 	= m_vari_edges_per_intersec_node / (long double) m_num_intersec_nodes;
	m_total_edge_length_per_vol 	= m_total_edge_length / (long double)m_effective_dataset_volume;

    m_network_vol_calc = (m_total_edge_length * (m_rad * m_rad * itk::Math::pi) / 1000000); //edge_length is already in mm
    m_network_vol_calc_per_vol = m_network_vol_calc / (long double)m_effective_dataset_volume;
}


void AnalyzeGraphFilter::Update()
{
    ParseParameterContext();

	CollectBasicImageInformation();

	BasePipeline<3>::ReadImage(m_networkDataSetPath + m_networkDataSetName + m_networkDataSetFileExtension, m_networkBin, m_vox_spacing);
	if(m_withRadiusMeasurement) BasePipeline<3>::BuildDistanceMap(m_networkBin, m_networkDistMap, m_vox_spacing, false);

	float maxVectorLengthForAngleComputation = 0.;
	for(unsigned int i=0; i<3; i++)
	    if(m_vox_spacing[i]>maxVectorLengthForAngleComputation)
	        maxVectorLengthForAngleComputation = m_vox_spacing[i];
	maxVectorLengthForAngleComputation *= 10.;

	m_network_vol_count = BasePipeline<3>::ComputeVolumeOfRegions(m_networkBin, m_voxel_volume);
	m_network_vol_count_per_vol = m_network_vol_count / m_effective_dataset_volume;

	mpAnno = new GraphAnnotationHelper();
	for(unsigned int i=0; i<m_graphs.size(); i++) {
	    std::map<vtkIdType,float> vertexIdToRadius;
	    std::map<vtkIdType,float> vertexIdToMinAngle;
	    std::map<vtkIdType,float> vertexIdToAvgMinAngle;

	    vtkSmartPointer<vtkVertexListIterator> it = vtkSmartPointer<vtkVertexListIterator>::New();
	    m_graphs[i]->GetVertices(it);

	    while(it->HasNext()) {
	        vtkIdType v = it->Next();

	        std::vector<float> angles;
	        GraphHandlingHelper::ComputeMinAngles(m_graphs[i], v, maxVectorLengthForAngleComputation, angles);
	        float minAngle = 360;
	        float avgMinAngle = 0;
	        for(unsigned int i=0; i<angles.size(); i++) {
	            if(minAngle > angles[i])    minAngle = angles[i];
	            avgMinAngle += angles[i];
	        }
	        avgMinAngle /= (float)angles.size();
	        vertexIdToMinAngle.insert(std::pair<vtkIdType,float>(v, minAngle));
	        vertexIdToAvgMinAngle.insert(std::pair<vtkIdType,float>(v, avgMinAngle));

	        if(m_withRadiusMeasurement) vertexIdToRadius.insert(std::pair<vtkIdType,float>(v, MeasureDistance(m_graphs[i], v, m_networkDistMap)));
	    }
	    mpAnno->AddCustomVertexAnnotation(m_graphs[i], "minimal branching angle (absolute)", vertexIdToMinAngle, -1.);
	    mpAnno->AddCustomVertexAnnotation(m_graphs[i], "minimal branching angle (averaged)", vertexIdToAvgMinAngle, -1.);
	    if(m_withRadiusMeasurement) mpAnno->AddCustomVertexAnnotation(m_graphs[i], "radius", vertexIdToRadius, -1);
	}
	CollectBasicNodeEdgeInformation();

	std::cout << "m_dataSetID " << m_dataSetID << "; m_dataSetName " << m_networkDataSetName << "; numDimensions " << m_num_dim << "; dataSetDimensions ";
	for(int i=0; i<m_num_dim; i++)
		std::cout << m_dim[i] << ", ";
	std::cout << " data set volume " << m_dataset_volume << "; effective data set volume " << m_effective_dataset_volume << "; num_nodes " << m_num_nodes << "; num_edges " << m_num_edges << std::endl;

	std::cout << "num_zero_nodes " << m_num_zero_nodes << "; num_total_dead_end_nodes " << m_num_total_dead_end_nodes << "; num_within_dead_end_nodes " << m_num_within_dead_end_nodes << "; num_reg_nodes " << m_num_reg_nodes << "; num_intersec_nodes " << m_num_intersec_nodes << std::endl;
	std::cout << "num_nodes_per_vol " << m_num_nodes_per_vol << "; num_dead_end_nodes_per_volume " << m_num_within_dead_end_nodes_per_vol << "; num_intersec_nodes_per_vol " << m_num_intersec_nodes_per_vol << std::endl;
	std::cout << "mean_edges_per_intersec_node " << m_mean_edges_per_intersec_node << "; vari_edges_per_intersec_node " << m_vari_edges_per_intersec_node << std::endl;
	std::cout << "total_edge_length (mm)" << m_total_edge_length << "; total_edge_length_per_volume " << m_total_edge_length_per_vol << std::endl;
	std::cout << "mean radius (micron)" << m_mean_radius << "; m_num_neg_radius " << m_num_neg_radius << std::endl;

	m_networkBin->ReleaseData();
	m_networkBin = NULL;
	m_networkDistMap->ReleaseData();
	m_networkDistMap = NULL;

	for(unsigned int i=0; i<m_deadEndMasks.size(); i++) {
	    m_deadEndMasks[i]->ReleaseData();
	    m_deadEndMasks[i] = NULL;
	}
}



