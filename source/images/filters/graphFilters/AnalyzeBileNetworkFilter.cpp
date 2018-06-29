///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeBileNetworkFilter.cpp                                         //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-01-17                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "AnalyzeBileNetworkFilter.h"

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>

#include <vtkArrayWriter.h>
#include <vtkArrayData.h>
#include <vtkDataSetAttributes.h>
#include <vtkDelimitedTextReader.h>
#include <vtkDelimitedTextWriter.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkMath.h>
#include <vtkTable.h>
#include <vtkVariant.h>
#include <vtkVariantArray.h>
#include <vtkVertexListIterator.h>

#include "GraphHandlingHelper.h"
#include "../../tools/ImageAnalysisSummaryFileIO.h"
#include "../../../tools/parameters/CSParameter.h"
#include "../../../tools/parameters/CSParameterContext.h"
#include "../../../tools/Tools.h"


AnalyzeBileNetworkFilter::AnalyzeBileNetworkFilter()
{
    m_necRegDataSetName = "";
    m_cvDataSetName = "";
    m_pvDataSetName = "";

    m_withVeinConnectednessMeasurement = false;

    m_maxCVeinConnectednessDistance = 10.;  //micron
    m_maxPVeinConnectednessDistance = 10.;  //micron

    m_maxAdjacentCells = 0;
    m_maxMinAngles = 0;

    m_output_prefix = "";

    m_cv_volume = 0;
    m_pv_volume = 0;

	m_num_branches = 0;
	m_num_branches_with_dead_end = 0;
	m_num_branches_between_intersections = 0;
	m_num_branches_rest = 0;
	m_num_secOrderBranches = 0;

	m_mean_branch_length = 0;
	m_mean_branch_with_dead_end_length = 0;
	m_mean_branch_between_intersections_length = 0;
	m_mean_secOrderBranch_length = 0;
	m_mean_minSecOrderBranch_length = 0;
	m_mean_maxSecOrderBranch_length = 0;

	m_mean_branch_between_intersections_radius = 0;
	m_mean_branch_with_dead_end_radius = 0;
	m_mean_secOrderBranch_radius = 0;
}


AnalyzeBileNetworkFilter::~AnalyzeBileNetworkFilter()
{
	// TODO Auto-generated destructor stub
}


void AnalyzeBileNetworkFilter::ParseParameterContext()
{
    AnalyzeGraphFilter::ParseParameterContext();

    m_output_prefix = *(std::string*)(m_paramContext->findParameter("Output file prefix", 0)->dataPointer());
    
	m_necRegDataSetFullFilename = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Specify necrotic region segmentation file", 0)->dataPointer()) );
	m_infoNecRegDataSetFullFilename.setFile(m_necRegDataSetFullFilename);
	m_withNecReg = m_infoNecRegDataSetFullFilename.exists();

	if(m_infoNecRegDataSetFullFilename.exists()) {
		m_necRegDataSetPath = (m_infoNecRegDataSetFullFilename.path() + QString("/")).toStdString();
		m_necRegDataSetName = m_infoNecRegDataSetFullFilename.baseName().toStdString();
		m_necRegDataSetFileExtension = (QString(".") + m_infoNecRegDataSetFullFilename.suffix()).toStdString();
	}

	m_CellDataSetFullFilename = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Specify cell segmentation file", 0)->dataPointer()) );
	m_infoCellDataSetFullFilename.setFile(m_CellDataSetFullFilename);
	m_withCell = m_infoCellDataSetFullFilename.exists();

	if(m_infoCellDataSetFullFilename.exists()) {
	    m_cellDataSetPath = (m_infoCellDataSetFullFilename.path() + QString("/")).toStdString();
	    m_cellDataSetName = m_infoCellDataSetFullFilename.baseName().toStdString();
	    m_cellDataSetFileExtension = (QString(".") + m_infoCellDataSetFullFilename.suffix()).toStdString();
	}

    m_CVDataSetFullFilename = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Specify central vein segmentation file", 0)->dataPointer()) );
	m_infoCVDataSetFullFilename.setFile(m_CVDataSetFullFilename);
	m_withCV = m_infoCVDataSetFullFilename.exists();

	if(m_infoCVDataSetFullFilename.exists()) {
		m_cvDataSetPath = (m_infoCVDataSetFullFilename.path() + QString("/")).toStdString();
		m_cvDataSetName = m_infoCVDataSetFullFilename.baseName().toStdString();
		m_cvDataSetFileExtension = (QString(".") + m_infoCVDataSetFullFilename.suffix()).toStdString();
	}

    m_PVDataSetFullFilename = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Specify portal vein segmentation file", 0)->dataPointer()) );
	m_infoPVDataSetFullFilename.setFile(m_PVDataSetFullFilename);
	m_withPV = m_infoPVDataSetFullFilename.exists();

	if(m_infoPVDataSetFullFilename.exists()) {
		m_pvDataSetPath = (m_infoPVDataSetFullFilename.path() + QString("/")).toStdString();
		m_pvDataSetName = m_infoPVDataSetFullFilename.baseName().toStdString();
		m_pvDataSetFileExtension = (QString(".") + m_infoPVDataSetFullFilename.suffix()).toStdString();
	}

	if(m_infoCVDataSetFullFilename.exists() && m_infoCVDataSetFullFilename.exists()) 
		m_withVeinConnectednessMeasurement = true;

    std::string networkType = ( (CSParameterChoice*)(m_paramContext->findParameter("Network type", 0)->dataPointer()) )->currentString();
    if( networkType.compare("Bile Canaliculi Network")==0 ) {
        m_maxCVeinConnectednessDistance = -1000;
        m_maxPVeinConnectednessDistance = 10.0;  //micron
    }
    else if( networkType.compare("Sinusoidal Network")==0 ) {
        m_maxCVeinConnectednessDistance = 10.0;  //micron
        m_maxPVeinConnectednessDistance = 10.0;  //micron
    }
}


void AnalyzeBileNetworkFilter::MeasureBranchVeinDistances(vtkGraph *input, FirstOrderBranch &branch)
{
    long double branchCentroid[3];
    branchCentroid[0] = 0;
    branchCentroid[1] = 0;
    branchCentroid[2] = 0;

    for(unsigned int i=0; i<branch.nodes.size(); i++) {
        double pos[3];
        input->GetPoint(branch.nodes[i], pos);

        branchCentroid[0] += pos[0];
        branchCentroid[1] += pos[1];
        branchCentroid[2] += pos[2];
    }

    CScalarVoImageType::IndexType idx;
    idx[0] = branchCentroid[0] / (long double)branch.nodes.size();
    idx[1] = branchCentroid[1] / (long double)branch.nodes.size();
    idx[2] = branchCentroid[2] / (long double)branch.nodes.size();

    branch.distanceToCV = m_cvDistMap->GetPixel(idx);
    branch.distanceToPV = m_pvDistMap->GetPixel(idx);
}


std::set<int> AnalyzeBileNetworkFilter::FindAdjacentCells(vtkGraph *input, vtkIdType vertex, float radius)
{
    std::set<int> adjacentCells;

    double pos[3];
    input->GetPoint(vertex, pos);

    CScalarVoImageType::IndexType idx;
    idx[0] = pos[0];
    idx[1] = pos[1];
    idx[2] = pos[2];

    IteratorType2::RadiusType elementRadius;
    elementRadius[0] = ceil(radius / m_vox_spacing[0]);
    elementRadius[1] = ceil(radius / m_vox_spacing[1]);
    elementRadius[2] = ceil(radius / m_vox_spacing[2]);

//    Solution with sphere as search area, too slow for practical use
//    StructuringElementType structuringElement;
//    structuringElement.SetRadius(elementRadius);
//    structuringElement.CreateStructuringElement();
//
//    IteratorType siterator(structuringElement.GetRadius(), m_cellLabelImage, m_cellLabelImage->GetLargestPossibleRegion());
//    siterator.CreateActiveListFromNeighborhood(structuringElement);
//    siterator.NeedToUseBoundaryConditionOn();
//    siterator.SetLocation(idx);
//
//    for(IteratorType::ConstIterator i = siterator.Begin(); !i.IsAtEnd(); ++i) {
//        if(i.Get()!=0 && adjacentCells.count(i.Get())==0)
//            adjacentCells.insert(i.Get());
//    }

    //Instead use standard (matrix shaped) neighborhood iterator
    BoundaryConditionType cb;
    cb.SetConstant(0);

    IteratorType2 siterator(elementRadius, m_cellLabelImage, m_cellLabelImage->GetLargestPossibleRegion());
    siterator.SetLocation(idx);
    siterator.SetBoundaryCondition(cb);

    int numElements = (2*elementRadius[0]+1)*(2*elementRadius[1]+1)*(2*elementRadius[2]+1);

    for(unsigned int i=0; i<numElements; i++) {
        LScalarPixelType pxl = siterator.GetPixel(i);
        if(pxl!=0 && adjacentCells.count(pxl)==0)
            adjacentCells.insert(pxl);
    }

    return adjacentCells;
}


void AnalyzeBileNetworkFilter::MeasureBranchVeinDistances(SecondOrderBranch &branch)
{
    std::set<vtkIdType> vertexEncountered;

    long double branchCentroid[3];
    branchCentroid[0] = 0;
    branchCentroid[1] = 0;
    branchCentroid[2] = 0;

    for(unsigned int i=0; i<branch.firstOrderBranches.size(); i++) {
        for(unsigned int j=0; j<m_isecBranches[branch.graphNr][branch.firstOrderBranches[i]].nodes.size(); j++) {
            if(vertexEncountered.count(m_isecBranches[branch.graphNr][branch.firstOrderBranches[i]].nodes[j])==0) {
                double pos[3];
                m_graphs[branch.graphNr]->GetPoint(m_isecBranches[branch.graphNr][branch.firstOrderBranches[i]].nodes[j], pos);

                branchCentroid[0] += pos[0];
                branchCentroid[1] += pos[1];
                branchCentroid[2] += pos[2];

                vertexEncountered.insert(m_isecBranches[branch.graphNr][branch.firstOrderBranches[i]].nodes[j]);
            }
        }
    }

    CScalarVoImageType::IndexType idx;
    idx[0] = branchCentroid[0] / (long double)vertexEncountered.size();
    idx[1] = branchCentroid[1] / (long double)vertexEncountered.size();
    idx[2] = branchCentroid[2] / (long double)vertexEncountered.size();

    branch.distanceToCV = m_cvDistMap->GetPixel(idx);
    branch.distanceToPV = m_pvDistMap->GetPixel(idx);
}


void AnalyzeBileNetworkFilter::BuildBranchMapsAndPopulateGraphBranchArrays()
{
    m_isecBranches.clear();
    m_innerDeadEndBranches.clear();
    m_maskedDeadEndBranches.clear();
    m_secondOrderBranches.clear();

    for(unsigned int i=0; i<m_graphs.size(); i++) {
        std::vector<FirstOrderBranch> deadEndBranches, intersectionBranches, innerDeadEndBranches, maskedDeadEndBranches;
        std::vector<SecondOrderBranch> secondOrderBranches;

        std::map<vtkIdType, int> edgeIdToBranchId = GraphHandlingHelper::AssembleFirstOrderBranches(m_graphs[i], deadEndBranches, intersectionBranches);
        AnalyzeBranchesAndPopulateGraphArrays(m_graphs[i], intersectionBranches, deadEndBranches, edgeIdToBranchId);
        IdentifyAnalysableDeadEnds(m_graphs[i], deadEndBranches, innerDeadEndBranches, maskedDeadEndBranches);
        AddBranchStateGraphArrays(m_graphs[i], intersectionBranches, innerDeadEndBranches);

        m_isecBranches.insert(std::pair<int, std::vector<FirstOrderBranch> >(i, intersectionBranches));
        m_innerDeadEndBranches.insert(std::pair<int, std::vector<FirstOrderBranch> >(i, innerDeadEndBranches));
        m_maskedDeadEndBranches.insert(std::pair<int, std::vector<FirstOrderBranch> >(i, maskedDeadEndBranches));

        std::map<vtkIdType, int> edgeIdToSecOBranchId, vertexIdToSecOState;
        GraphHandlingHelper::AssembleSecondOrderBranches(m_graphs[i], secondOrderBranches, intersectionBranches, innerDeadEndBranches, maskedDeadEndBranches, vertexIdToSecOState, edgeIdToSecOBranchId);
        for(unsigned int j=0; j<secondOrderBranches.size(); j++)
            secondOrderBranches[j].graphNr = i;
        AnalyzeBranchesAndPopulateGraphArrays(m_graphs[i], secondOrderBranches, vertexIdToSecOState, edgeIdToSecOBranchId);

        m_secondOrderBranches.insert(std::pair<int, std::vector<SecondOrderBranch> >(i, secondOrderBranches));

        m_num_branches_between_intersections += intersectionBranches.size();
        m_num_branches_with_dead_end += innerDeadEndBranches.size();
        m_num_branches_rest += maskedDeadEndBranches.size();
        m_num_secOrderBranches += secondOrderBranches.size();
    }
//    mpAnno->AddCustomVertexAnnotation(m_graphs, "MinAngles", m_minAnglesPerGraph, m_maxMinAngles, -1);
    mpAnno->AddCustomVertexAnnotation(m_graphs, "Cell", m_adjacentCellsPerGraph, m_maxAdjacentCells, -1);

    m_num_branches = m_num_branches_between_intersections + m_num_branches_with_dead_end + m_num_branches_rest;
}


void AnalyzeBileNetworkFilter::IdentifyAnalysableDeadEnds(vtkGraph *input, std::vector<FirstOrderBranch> &deadEndBranches, std::vector<FirstOrderBranch> &analysableDeadEndBranches,
        std::vector<FirstOrderBranch> &nonAnalysableDeadEndBranches)
{
    vtkFloatArray* radiusArray;
    if(m_withRadiusMeasurement)
        radiusArray = vtkFloatArray::SafeDownCast( input->GetVertexData()->GetArray("radius") );

    for(unsigned int i=0; i<deadEndBranches.size(); i++) {
        bool isVertexAtDatasetBorder = false;
        bool isVertexConnectedToVein = false;

        if(input->GetDegree(deadEndBranches[i].firstVert) == 1) {
            if(m_withRadiusMeasurement)
                isVertexAtDatasetBorder |= GraphHandlingHelper::IsVertexAtDatasetBorder(input, deadEndBranches[i].firstVert, (2.0*radiusArray->GetValue(deadEndBranches[i].firstVert)), m_vox_spacing, m_dim);
            else
                isVertexAtDatasetBorder |= GraphHandlingHelper::IsVertexAtDatasetBorder(input, deadEndBranches[i].firstVert, (2.0*m_rad), m_vox_spacing, m_dim);

            int veinStatus = GraphHandlingHelper::GetVeinConnectednessStatus(input, deadEndBranches[i].firstVert, m_cvDistMap, m_pvDistMap,
                    m_maxCVeinConnectednessDistance, m_maxPVeinConnectednessDistance, m_withCV, m_withPV);
            if(veinStatus!=0)  isVertexConnectedToVein = true;
        }
        if(input->GetDegree(deadEndBranches[i].lastVert) == 1) {
            if(m_withRadiusMeasurement)
                isVertexAtDatasetBorder |= GraphHandlingHelper::IsVertexAtDatasetBorder(input, deadEndBranches[i].lastVert, (2.0*radiusArray->GetValue(deadEndBranches[i].lastVert)), m_vox_spacing, m_dim);
            else
                isVertexAtDatasetBorder |= GraphHandlingHelper::IsVertexAtDatasetBorder(input, deadEndBranches[i].lastVert, (2.0*m_rad), m_vox_spacing, m_dim);

            int veinStatus = GraphHandlingHelper::GetVeinConnectednessStatus(input, deadEndBranches[i].lastVert, m_cvDistMap, m_pvDistMap,
                    m_maxCVeinConnectednessDistance, m_maxPVeinConnectednessDistance, m_withCV, m_withPV);
            if(veinStatus!=0)  isVertexConnectedToVein = true;
        }

        if(isVertexAtDatasetBorder || isVertexConnectedToVein) {
            deadEndBranches[i].id = nonAnalysableDeadEndBranches.size();
            nonAnalysableDeadEndBranches.push_back(deadEndBranches[i]);
        }
        else {
            deadEndBranches[i].id = analysableDeadEndBranches.size();
            analysableDeadEndBranches.push_back(deadEndBranches[i]);
        }
    }
}


void AnalyzeBileNetworkFilter::AnalyzeBranchesAndPopulateGraphArrays(vtkGraph *input, std::vector<FirstOrderBranch> &isecBranches,
        std::vector<FirstOrderBranch> &deadEndBranches, std::map<vtkIdType, int> edgeIdToBranchId)
{
    //Collect and add branch independent data***************************************************

    std::map<vtkIdType,float> vertexIdToRadius;
    std::map<vtkIdType,float> vertexIdToMinAngle;
    std::map<vtkIdType,float> vertexIdToAvgMinAngle;
    std::map<vtkIdType,float> vertexIdToDistToCV;
    std::map<vtkIdType,float> vertexIdToDistToPV;
    std::map<vtkIdType,int> vertexIdToVeinState;
    std::multimap<vtkIdType,int> vertexIdToAdjacentCells;
//    std::multimap<vtkIdType,float> vertexIdToMinAngles;

    vtkSmartPointer<vtkVertexListIterator> it = vtkSmartPointer<vtkVertexListIterator>::New();
    input->GetVertices(it);

    while(it->HasNext()) {
        vtkIdType v = it->Next();

        vertexIdToDistToCV.insert(std::pair<vtkIdType,float>(v, MeasureDistance(input, v, m_cvDistMap)));
        vertexIdToDistToPV.insert(std::pair<vtkIdType,float>(v, MeasureDistance(input, v, m_pvDistMap)));

        std::vector<float> angles;
        GraphHandlingHelper::ComputeMinAngles(input, v, 4., angles);
        float minAngle = 360;
        float avgMinAngle = 0;
        for(unsigned int i=0; i<angles.size(); i++) {
            if(minAngle > angles[i])    minAngle = angles[i];
            avgMinAngle += angles[i];
//            vertexIdToMinAngles.insert(std::pair<vtkIdType,float>(v, angles[i]));
        }
        avgMinAngle /= (float)angles.size();
//        if(angles.size() > m_maxMinAngles)  m_maxMinAngles = angles.size();
        vertexIdToMinAngle.insert(std::pair<vtkIdType,float>(v, minAngle));
        vertexIdToAvgMinAngle.insert(std::pair<vtkIdType,float>(v, avgMinAngle));

        float radius;
        if(m_withRadiusMeasurement) {
            radius = MeasureDistance(input, v, m_networkDistMap);
            if(radius<0.0001) radius = 0.0;
            vertexIdToRadius.insert(std::pair<vtkIdType,float>(v, radius));
        }

        if(m_withVeinConnectednessMeasurement)
            vertexIdToVeinState.insert(std::pair<vtkIdType,int>(v, GraphHandlingHelper::GetVeinConnectednessStatus(input, v, m_cvDistMap, m_pvDistMap,
                    m_maxCVeinConnectednessDistance, m_maxPVeinConnectednessDistance, m_withCV, m_withPV)));

        if(m_withCell) {
            std::set<int> cells;
            if(m_withRadiusMeasurement) cells = FindAdjacentCells(input, v, radius + 1.2);
            else                        cells = FindAdjacentCells(input, v, m_rad + 1.2);

            for(std::set<int>::iterator ic = cells.begin(); ic != cells.end(); ++ic)
                vertexIdToAdjacentCells.insert(std::pair<vtkIdType,int>(v, *ic));

            if(cells.size() > m_maxAdjacentCells)
                m_maxAdjacentCells = cells.size();
        }
    }
    mpAnno->AddCustomVertexAnnotation(input, "distance to CV", vertexIdToDistToCV, -100);
    mpAnno->AddCustomVertexAnnotation(input, "distance to PV", vertexIdToDistToPV, -100);
    mpAnno->AddCustomVertexAnnotation(input, "minimal branching angle (absolute)", vertexIdToMinAngle, -1.);
    mpAnno->AddCustomVertexAnnotation(input, "minimal branching angle (averaged)", vertexIdToAvgMinAngle, -1.);
//    m_minAnglesPerGraph.push_back(vertexIdToMinAngles);
    if(m_withRadiusMeasurement)             mpAnno->AddCustomVertexAnnotation(input, "radius", vertexIdToRadius, -100);
    if(m_withVeinConnectednessMeasurement)  mpAnno->AddCustomVertexAnnotation(input, "vein state", vertexIdToVeinState, -100);
    if(m_withCell)                          m_adjacentCellsPerGraph.push_back(vertexIdToAdjacentCells);
    //Collect and add branch independent data****finished***************************************


    //Collect and add branch dependent data*****************************************************

    std::map<vtkIdType, float> vertexIdToBranchDistToCV;
    std::map<vtkIdType, float> vertexIdToBranchDistToPV;

    vtkFloatArray* radiusArray;
    if(m_withRadiusMeasurement) radiusArray = vtkFloatArray::SafeDownCast( input->GetVertexData()->GetArray("radius") );

    for(unsigned int i=0; i<isecBranches.size(); i++) {
        MeasureBranchVeinDistances(input, isecBranches[i]);                                                      // sets isecBranches[i].distanceToCV &&  isecBranches[i].distanceToPV

        double branch_radius = 0;

        for(unsigned int j=0; j<isecBranches[i].nodes.size(); j++) {
            vertexIdToBranchDistToCV[ isecBranches[i].nodes[j] ] = isecBranches[i].distanceToCV;
            vertexIdToBranchDistToPV[ isecBranches[i].nodes[j] ] = isecBranches[i].distanceToPV;

            if(m_withRadiusMeasurement)
                branch_radius += radiusArray->GetValue(isecBranches[i].nodes[j]);
        }
        if(m_withRadiusMeasurement)
            isecBranches[i].radiusMean = branch_radius / (float)isecBranches[i].nodes.size();
    }
    for(unsigned int i=0; i<deadEndBranches.size(); i++) {
        MeasureBranchVeinDistances(input, deadEndBranches[i]);                                         // sets isecBranches[i].distanceToCV &&  isecBranches[i].distanceToPV

        double branch_radius = 0;

        for(unsigned int j=0; j<deadEndBranches[i].nodes.size(); j++) {
            vertexIdToBranchDistToCV[ deadEndBranches[i].nodes[j] ] = deadEndBranches[i].distanceToCV;
            vertexIdToBranchDistToPV[ deadEndBranches[i].nodes[j] ] = deadEndBranches[i].distanceToPV;

            if(m_withRadiusMeasurement)
                branch_radius += radiusArray->GetValue(deadEndBranches[i].nodes[j]);
        }
        if(m_withRadiusMeasurement)
            deadEndBranches[i].radiusMean = branch_radius / (float)deadEndBranches[i].nodes.size();
    }
    mpAnno->AddCustomEdgeAnnotation(input, "branch Id", edgeIdToBranchId, -100);
    //For test visualization
    mpAnno->AddCustomVertexAnnotation(input, "branch CV Dist", vertexIdToBranchDistToCV, -100);
    mpAnno->AddCustomVertexAnnotation(input, "branch PV Dist", vertexIdToBranchDistToPV, -100);

    input->Squeeze();
    //Collect and add branch dependent data****finished*****************************************
}


void AnalyzeBileNetworkFilter::AnalyzeBranchesAndPopulateGraphArrays(vtkGraph *input, std::vector<SecondOrderBranch> &secOrderBranches,
        std::map<vtkIdType, int> &vertexIdToSecOState, std::map<vtkIdType, int> &edgeIdToSecOBranchId)
{
    std::map<vtkIdType, float> vertIDToMinSecBranch;
    std::map<vtkIdType, float> vertIDToMaxSecBranch;

    for(unsigned int i=0; i<secOrderBranches.size(); i++) {
        //Measure minimal and maximal second order branch length starting at second order vertex
        if(vertIDToMinSecBranch.count(secOrderBranches[i].firstVert) == 0)
            vertIDToMinSecBranch.insert(std::pair<vtkIdType, float>(secOrderBranches[i].firstVert, secOrderBranches[i].length));
        else {
            if(vertIDToMinSecBranch[secOrderBranches[i].firstVert]>secOrderBranches[i].length)
                vertIDToMinSecBranch[secOrderBranches[i].firstVert] = secOrderBranches[i].length;
        }
        if(vertIDToMinSecBranch.count(secOrderBranches[i].lastVert) == 0)
            vertIDToMinSecBranch.insert(std::pair<vtkIdType, float>(secOrderBranches[i].lastVert, secOrderBranches[i].length));
        else {
            if(vertIDToMinSecBranch[secOrderBranches[i].lastVert]>secOrderBranches[i].length)
                vertIDToMinSecBranch[secOrderBranches[i].lastVert] = secOrderBranches[i].length;
        }
        if(vertIDToMaxSecBranch.count(secOrderBranches[i].firstVert) == 0)
            vertIDToMaxSecBranch.insert(std::pair<vtkIdType, float>(secOrderBranches[i].firstVert, secOrderBranches[i].length));
        else {
            if(vertIDToMaxSecBranch[secOrderBranches[i].firstVert]<secOrderBranches[i].length)
                vertIDToMaxSecBranch[secOrderBranches[i].firstVert] = secOrderBranches[i].length;
        }
        if(vertIDToMaxSecBranch.count(secOrderBranches[i].lastVert) == 0)
            vertIDToMaxSecBranch.insert(std::pair<vtkIdType, float>(secOrderBranches[i].lastVert, secOrderBranches[i].length));
        else {
            if(vertIDToMaxSecBranch[secOrderBranches[i].lastVert]<secOrderBranches[i].length)
                vertIDToMaxSecBranch[secOrderBranches[i].lastVert] = secOrderBranches[i].length;
        }

        //Measure second order branch radii
        vtkFloatArray* radiusArray;
        if(m_withRadiusMeasurement) {
            radiusArray = vtkFloatArray::SafeDownCast( input->GetVertexData()->GetArray("radius") );
            secOrderBranches[i].radiusMean = 0.;
        }
        std::set<vtkIdType> vertexEncountered;

        for(unsigned int j=0; j<secOrderBranches[i].firstOrderBranches.size(); j++) {
            for(unsigned int k=0; k<m_isecBranches[secOrderBranches[i].graphNr][ secOrderBranches[i].firstOrderBranches[j] ].nodes.size(); k++) {

                if(m_withRadiusMeasurement && vertexEncountered.count( m_isecBranches[secOrderBranches[i].graphNr][ secOrderBranches[i].firstOrderBranches[j] ].nodes[k] ) == 0) {
                    vertexEncountered.insert(m_isecBranches[secOrderBranches[i].graphNr][secOrderBranches[i].firstOrderBranches[j]].nodes[k]);
                    double rad = radiusArray->GetValue(m_isecBranches[secOrderBranches[i].graphNr][ secOrderBranches[i].firstOrderBranches[j] ].nodes[k]);
                    secOrderBranches[i].radiusMean += rad;
                }
            }
        }
        if(m_withRadiusMeasurement)
            secOrderBranches[i].radiusMean /= (long double)vertexEncountered.size();
        MeasureBranchVeinDistances(secOrderBranches[i]);
    }
    mpAnno->AddCustomEdgeAnnotation(input, "second order branch Id", edgeIdToSecOBranchId, -3);
    mpAnno->AddCustomVertexAnnotation(input, "second order branches", vertexIdToSecOState, -100);
    mpAnno->AddCustomVertexAnnotation(input, "minimal second order branches", vertIDToMinSecBranch, -1);
    mpAnno->AddCustomVertexAnnotation(input, "maximal second order branches", vertIDToMaxSecBranch, -1);

    input->Squeeze();
}


void AnalyzeBileNetworkFilter::AddBranchStateGraphArrays(vtkGraph *input, std::vector<FirstOrderBranch> &isecBranches, std::vector<FirstOrderBranch> &innerDeadEndBranches)
{
    std::map<vtkIdType, int> vertexIdToDEndState;
    std::map<vtkIdType, int> vertexIdToIsecState;

    for(unsigned int i=0; i<isecBranches.size(); i++)
        for(unsigned int j=0; j<isecBranches[i].nodes.size(); j++)
            vertexIdToIsecState[ isecBranches[i].nodes[j] ] = 1;

    for(unsigned int i=0; i<innerDeadEndBranches.size(); i++)
        for(unsigned int j=0; j<innerDeadEndBranches[i].nodes.size(); j++)
            vertexIdToDEndState[ innerDeadEndBranches[i].nodes[j] ] = 1;

    mpAnno->AddCustomVertexAnnotation(input, "dead end branches", vertexIdToDEndState, -100);
    mpAnno->AddCustomVertexAnnotation(input, "intersection branches", vertexIdToIsecState, -100);

    input->Squeeze();
}


//std::vector< std::vector<vtkIdType> > AnalyzeBileNetworkFilter::AssembleCycles(int subgraph, vtkIdType node, int max_cycle_size)
//{
//    std::vector< std::vector<vtkIdType> > cycleExplorer;
//    std::vector< std::vector<vtkIdType> > cycles;
//
//    std::vector<vtkIdType> seed(1, node);
//    cycleExplorer.push_back( seed );
//
//    while(cycleExplorer.size()>0) {
////        for(unsigned int z=0; z<cycleExplorer.size(); z++) {
////            std::cout << "cycleExplorer " << z << " composed of nodes: ";
////            for(unsigned int y=0; y<cycleExplorer[z].size(); y++)
////                std::cout << cycleExplorer[z][y] << ", ";
////            std::cout << std::endl;
////        }
//
//        std::vector<vtkIdType> series = cycleExplorer[0];
//        cycleExplorer.erase(cycleExplorer.begin());
//
//        vtkIdType v = series[series.size()-1];
//
//        for(int j=0; j<m_graphs[subgraph]->GetDegree(v); j++) {
//            vtkIdType d = m_graphs[subgraph]->GetOutEdge(v, j).Target;
//
//            if(series.size()==1 || (series.size()>1 && d!=series[series.size()-2])) {
//                bool isElem = false;
//
//                if(d==series[0]) {
//                    cycles.push_back(series);
//                    isElem = true;
//                }
//
//                for(unsigned int j=1; j<series.size(); j++) {
//                    if(d==series[j]) {
//                        isElem = true;
//                        break;
//                    }
//                }
//
//                //test fuer maximal entfernung -- nur beschleunigende wirkung
//
//                if(!isElem && series.size()+1 <= max_cycle_size) {
//                    std::vector<vtkIdType> temp = series;
//                    temp.push_back(d);
//                    cycleExplorer.push_back(temp);
//                }
//            }
//        }
//    }
//
////    for(unsigned int i=0; i<cycles.size(); i++) {
////        std::cout << "cycle " << i << " composed of nodes: ";
////        for(unsigned int j=0; j<cycles[i].size(); j++)
////            std::cout << cycles[i][j] << ", ";
////        std::cout << std::endl;
////    }
//
//    std::vector< std::vector<vtkIdType> > orderedCycles;
//    for(unsigned int i=0; i<cycles.size(); i++) {
//        orderedCycles.push_back(cycles[i]);
//        std::sort( orderedCycles[i].begin(), orderedCycles[i].end() );
//
//        for(unsigned int j=0; j<orderedCycles.size()-1; j++) {
//            if(orderedCycles[orderedCycles.size()-1]==orderedCycles[j]) {
//                cycles.erase(cycles.begin()+i);
//                orderedCycles.erase(orderedCycles.begin()+i);
//                i--;
//            }
//        }
//    }
//
////    for(unsigned int i=0; i<cycles.size(); i++) {
////        std::cout << "cycle " << i << " composed of nodes: ";
////        for(unsigned int j=0; j<cycles[i].size(); j++)
////            std::cout << cycles[i][j] << ", ";
////        std::cout << std::endl;
////    }
//
//    return cycles;
//}


void AnalyzeBileNetworkFilter::CollectFirstOrderBranchInformation()
{
	m_num_edges_per_branch = m_num_edges / (long double)m_num_branches;
	m_num_branches_per_vol = m_num_branches / (long double)m_effective_dataset_volume;
	m_num_branches_with_dead_end_per_vol = m_num_branches_with_dead_end / (long double)m_effective_dataset_volume;
	m_num_branches_between_intersections_per_vol = m_num_branches_between_intersections / (long double)m_effective_dataset_volume;

	for(unsigned int i=0; i<m_isecBranches.size(); i++) {
	    for(unsigned int j=0; j<m_isecBranches[i].size(); j++) {
	        m_mean_branch_between_intersections_length += m_isecBranches[i][j].length;
	        if(m_withRadiusMeasurement) m_mean_branch_between_intersections_radius += m_isecBranches[i][j].radiusMean;
	    }
	}

	for(unsigned int i=0; i<m_innerDeadEndBranches.size(); i++) {
	    for(unsigned int j=0; j<m_innerDeadEndBranches[i].size(); j++) {
	        m_mean_branch_with_dead_end_length += m_innerDeadEndBranches[i][j].length;
	        if(m_withRadiusMeasurement) m_mean_branch_with_dead_end_radius += m_innerDeadEndBranches[i][j].radiusMean;
	    }
	}

	long double mean_rest_branch_length = 0, mean_rest_branch_radius = 0;
	for(unsigned int i=0; i<m_maskedDeadEndBranches.size(); i++) {
	    for(unsigned int j=0; j<m_maskedDeadEndBranches[i].size(); j++) {
	        mean_rest_branch_length += m_maskedDeadEndBranches[i][j].length;
	        if(m_withRadiusMeasurement) mean_rest_branch_radius += m_maskedDeadEndBranches[i][j].radiusMean;
	    }
	}

	m_mean_branch_length = m_mean_branch_between_intersections_length + m_mean_branch_with_dead_end_length;
	long number_real_branches = m_num_branches_with_dead_end + m_num_branches_between_intersections;

	m_mean_branch_length /= (long double)number_real_branches;
	m_mean_branch_with_dead_end_length /= (long double)m_num_branches_with_dead_end;
	m_mean_branch_between_intersections_length /= (long double)m_num_branches_between_intersections;
	if(m_withRadiusMeasurement) m_mean_branch_with_dead_end_radius /= (long double)m_num_branches_with_dead_end;
	if(m_withRadiusMeasurement) m_mean_branch_between_intersections_radius /= (long double)m_num_branches_between_intersections;
}


void AnalyzeBileNetworkFilter::CollectSecondOrderBranchInformation()
{
    m_num_secOrderBranches_per_vol = m_num_secOrderBranches / (long double)m_effective_dataset_volume;

    m_num_isecBranches_per_secOrderBranch = 0;
    m_mean_secOrderBranch_length = 0;
    m_mean_secOrderBranch_radius = 0;

    for(unsigned int i=0; i<m_secondOrderBranches.size(); i++) {
        for(unsigned int j=0; j<m_secondOrderBranches[i].size(); j++) {
            m_num_isecBranches_per_secOrderBranch += m_secondOrderBranches[i][j].firstOrderBranches.size();
            m_mean_secOrderBranch_length += m_secondOrderBranches[i][j].length;
            if(m_withRadiusMeasurement) m_mean_secOrderBranch_radius += (long double)m_secondOrderBranches[i][j].radiusMean;
        }
    }
    m_num_isecBranches_per_secOrderBranch /= (long double)m_num_secOrderBranches;
    m_mean_secOrderBranch_length /= (long double)m_num_secOrderBranches;
    if(m_withRadiusMeasurement) m_mean_secOrderBranch_radius = m_mean_secOrderBranch_radius / (long double)m_num_secOrderBranches;

    int numMinSecOBranchLengthMeasurements = 0;
    int numMaxSecOBranchLengthMeasurements = 0;
    for(unsigned int i=0; i<m_graphs.size(); i++) {
        vtkFloatArray* minSecOrderLengthArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("minimal second order branches") );
        vtkFloatArray* maxSecOrderLengthArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("maximal second order branches") );

        for(int j=0; j<m_graphs[i]->GetNumberOfVertices(); j++) {
            if(minSecOrderLengthArray->GetValue(j) != -1) {
                m_mean_minSecOrderBranch_length += minSecOrderLengthArray->GetValue(j);
                numMinSecOBranchLengthMeasurements++;
            }
            if(maxSecOrderLengthArray->GetValue(j) != -1) {
                m_mean_maxSecOrderBranch_length += maxSecOrderLengthArray->GetValue(j);
                numMaxSecOBranchLengthMeasurements++;
            }
        }
    }
    m_mean_minSecOrderBranch_length /= (double)numMinSecOBranchLengthMeasurements;
    m_mean_maxSecOrderBranch_length /= (double)numMaxSecOBranchLengthMeasurements;
}


void AnalyzeBileNetworkFilter::CollectCycleInformation()
{
    //Stub to be filled out
}


void AnalyzeBileNetworkFilter::WriteToFile()
{
	WriteGeneralInfoToFile();
	WriteNodeInfoToFile();
	WriteBranchInfoToFile();
	WriteSecOrderBranchInfoToFile();
//	WriteCycleInfoToFile();
}

void AnalyzeBileNetworkFilter::WriteGeneralInfoToFile()
{
	QFileInfo fileInfo;
	fileInfo.setFile(QString::fromStdString( (m_networkDataSetPath + "../" + m_output_prefix + "Output_file_1.txt").c_str() ));

	vtkSmartPointer<vtkTable> readTable;
	if(fileInfo.exists()) {
		vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
		reader->SetFileName((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_1.txt").c_str());
		reader->SetHaveHeaders(true);
		reader->DetectNumericColumnsOn();
		reader->SetFieldDelimiterCharacters(" ");
		reader->Update();
		
		readTable = reader->GetOutput();
	}
	else
		readTable = vtkTable::New();

	vtkSmartPointer<vtkVariantArray> head =  vtkVariantArray::SafeDownCast(readTable->GetColumnByName("header"));
	if(head==NULL) {
		vtkSmartPointer<vtkVariantArray> header = vtkSmartPointer<vtkVariantArray>::New();
		header->InsertNextValue( vtkVariant("#dim") );
		header->InsertNextValue( vtkVariant("dimX") );
		header->InsertNextValue( vtkVariant("dimY") );
		header->InsertNextValue( vtkVariant("dimZ") );
		header->InsertNextValue( vtkVariant("dataSetVolume") );
		header->InsertNextValue( vtkVariant("veinVolume") );
		header->InsertNextValue( vtkVariant("effectiveVolume") );
		header->InsertNextValue( vtkVariant("#nodes") );
		header->InsertNextValue( vtkVariant("#nodes/effVol") );
		header->InsertNextValue( vtkVariant("#edges") );
		header->InsertNextValue( vtkVariant("length_of_network(mm)") );
		header->InsertNextValue( vtkVariant("length_of_network(mm)/effVol") );
		header->InsertNextValue( vtkVariant("avg_radius(micron)") );
		header->InsertNextValue( vtkVariant("volCalc_of_network(mm³)") );
		header->InsertNextValue( vtkVariant("volCalc_of_network(mm³)/effVol") );
		header->InsertNextValue( vtkVariant("volCount_of_network(mm³)") );
		header->InsertNextValue( vtkVariant("volCount_of_network(mm³)/effVol") );
		header->InsertNextValue( vtkVariant("#fOBranches") );
		header->InsertNextValue( vtkVariant("#fOBranches/effVol") );
		header->InsertNextValue( vtkVariant("avg_#edges/fOBranch") );
		header->InsertNextValue( vtkVariant("#dend_fOBranches") );
		header->InsertNextValue( vtkVariant("#dend_fOBranches/effVol") );
		header->InsertNextValue( vtkVariant("#isec_fOBranches") );
		header->InsertNextValue( vtkVariant("#isec_fOBranches/effVol") );
		header->InsertNextValue( vtkVariant("#intersec_nodes") );
		header->InsertNextValue( vtkVariant("#intersec_nodes/effVol") );
		header->InsertNextValue( vtkVariant("avg_#branches/intersec_node") );
		header->InsertNextValue( vtkVariant("var_#branches/intersec_node") );
		header->InsertNextValue( vtkVariant("#sOBranches") );
		header->InsertNextValue( vtkVariant("#sOBranches/effVol") );
		header->InsertNextValue( vtkVariant("avg_#fOBranches/sOBranch") );
		header->InsertNextValue( vtkVariant("avg_fOBranch_length") );
		header->InsertNextValue( vtkVariant("avg_dEnd_fOBranch_length") );
		header->InsertNextValue( vtkVariant("avg_isec_fOBranch_length") );
		header->InsertNextValue( vtkVariant("avg_sOBranch_length") );
		header->InsertNextValue( vtkVariant("avg_minSOBranch_length") );
		header->InsertNextValue( vtkVariant("avg_maxSOBranch_length") );
		header->InsertNextValue( vtkVariant("avg_dEnd_fOBranch_radius") );
		header->InsertNextValue( vtkVariant("avg_isec_fOBranch_radius") );
		header->InsertNextValue( vtkVariant("avg_sOBranch_radius") );
		header->InsertNextValue( vtkVariant("avg_isecBranchMinAngle") );
		header->InsertNextValue( vtkVariant("avg_isecBranchAvgMinAngle") );
		header->SetName("header");

		readTable->AddColumn(header);
	}

	std::cout << "in table writer: m_mean_secOrderBranch_radius = " << m_mean_secOrderBranch_radius << std::endl;

	vtkSmartPointer<vtkVariantArray> dataSet = vtkSmartPointer<vtkVariantArray>::New();
	dataSet->InsertNextValue( vtkVariant(m_num_dim) );
	dataSet->InsertNextValue( vtkVariant(m_dim[0]) );
	dataSet->InsertNextValue( vtkVariant(m_dim[1]) );
	dataSet->InsertNextValue( vtkVariant(m_dim[2]) );
	dataSet->InsertNextValue( vtkVariant((double)m_dataset_volume) );
	dataSet->InsertNextValue( vtkVariant((double)(m_cv_volume + m_pv_volume)) );
	dataSet->InsertNextValue( vtkVariant((double)m_effective_dataset_volume) );
	dataSet->InsertNextValue( vtkVariant(m_num_nodes) );
	dataSet->InsertNextValue( vtkVariant((double)m_num_nodes_per_vol) );
	dataSet->InsertNextValue( vtkVariant(m_num_edges) );
	dataSet->InsertNextValue( vtkVariant((double)m_total_edge_length) );
	dataSet->InsertNextValue( vtkVariant((double)m_total_edge_length_per_vol) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_radius) );
	dataSet->InsertNextValue( vtkVariant((double)m_network_vol_calc) );
	dataSet->InsertNextValue( vtkVariant((double)m_network_vol_calc_per_vol) );
	dataSet->InsertNextValue( vtkVariant((double)m_network_vol_count) );
	dataSet->InsertNextValue( vtkVariant((double)m_network_vol_count_per_vol) );
	dataSet->InsertNextValue( vtkVariant(m_num_branches) );
	dataSet->InsertNextValue( vtkVariant((double)m_num_branches_per_vol) );
	dataSet->InsertNextValue( vtkVariant((double)m_num_edges_per_branch) );
	dataSet->InsertNextValue( vtkVariant(m_num_branches_with_dead_end) );
	dataSet->InsertNextValue( vtkVariant((double)m_num_branches_with_dead_end_per_vol) );
	dataSet->InsertNextValue( vtkVariant(m_num_branches_between_intersections) );
	dataSet->InsertNextValue( vtkVariant((double)m_num_branches_between_intersections_per_vol) );
	dataSet->InsertNextValue( vtkVariant(m_num_intersec_nodes) );
	dataSet->InsertNextValue( vtkVariant((double)m_num_intersec_nodes_per_vol) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_edges_per_intersec_node) );
	dataSet->InsertNextValue( vtkVariant((double)m_vari_edges_per_intersec_node) );
	dataSet->InsertNextValue( vtkVariant(m_num_secOrderBranches) );
	dataSet->InsertNextValue( vtkVariant((double)m_num_secOrderBranches_per_vol) );
	dataSet->InsertNextValue( vtkVariant((double)m_num_isecBranches_per_secOrderBranch) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_branch_length) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_branch_with_dead_end_length) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_branch_between_intersections_length) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_secOrderBranch_length) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_minSecOrderBranch_length) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_maxSecOrderBranch_length) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_branch_with_dead_end_radius) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_branch_between_intersections_radius) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_secOrderBranch_radius) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_isecBranch_minAngle) );
	dataSet->InsertNextValue( vtkVariant((double)m_mean_isecBranch_minAvgAngle) );
	dataSet->SetName(m_dataSetID.c_str());

	readTable->AddColumn(dataSet);

	readTable->Dump(25);

	vtkSmartPointer<vtkDelimitedTextWriter> writer = vtkSmartPointer<vtkDelimitedTextWriter>::New();
	writer->SetFileName((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_1.txt").c_str());
	writer->SetFieldDelimiter(" ");
	writer->SetInput(readTable);
	writer->Update();
}


void AnalyzeBileNetworkFilter::WriteNodeInfoToFile()
{
	std::fstream file1, file2, tempfile;
	bool hasHeader = false;

	file1.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_2.txt").c_str(), fstream::in);

	std::string line;
	getline(file1, line);

	if(line.find("dataSet")!=std::string::npos) {
		hasHeader = true;

		tempfile.open((m_networkDataSetPath + "../" + m_output_prefix + "tempfile_2.txt").c_str(), fstream::out);
		tempfile << line;

		while(!file1.eof()) {
			getline(file1, line);
			if(line.find(m_dataSetID)!=0) {
				tempfile << std::endl;
				tempfile << line;
			}
		}
		file1.close();
		tempfile.close();

		std::remove((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_2.txt").c_str());
		std::rename((m_networkDataSetPath + "../" + m_output_prefix + "tempfile_2.txt").c_str(), (m_networkDataSetPath + "../" + m_output_prefix + "Output_file_2.txt").c_str());
	}
	else
		file1.close();

	if(!hasHeader)
		file2.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_2.txt").c_str(), fstream::out);
	else
		file2.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_2.txt").c_str(), fstream::out | fstream::app);

	file2.flags(fstream::left | fstream::scientific);
	if(!hasHeader) {
		file2.width(40);
		file2 << "dataSet";
		file2.width(10);
		file2 << "isIsec";
		file2.width(10);
		file2 << "isDead";
		file2.width(10);
		file2 << "#edges";
		file2.width(20);
		file2 << "radius";
		file2.width(20);
		file2 << "minAngle";
		file2.width(20);
		file2 << "avgMinAngle";
		file2.width(20);
		file2 << "distToCV";
		file2.width(20);
		file2 << "distToPV";
		file2.width(20);
		file2 << "minSecOrderLength";
		file2.width(20);
		file2 << "maxSecOrderLength" << std::endl;
	}

	for(unsigned int i=0; i<m_graphs.size(); i++) {
	    vtkFloatArray* radiusArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("radius") );
	    vtkFloatArray* minAngleArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("minimal branching angle (absolute)") );
	    vtkFloatArray* avgMinAngleArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("minimal branching angle (averaged)") );
	    vtkFloatArray* distToCVArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("distance to CV") );
	    vtkFloatArray* distToPVArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("distance to PV") );
	    vtkIntArray* deadEndStateArray = vtkIntArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("dead end branches") );
	    vtkFloatArray* minSecOrderLengthArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("minimal second order branches") );
	    vtkFloatArray* maxSecOrderLengthArray = vtkFloatArray::SafeDownCast( m_graphs[i]->GetVertexData()->GetArray("maximal second order branches") );


		for(int j=0; j<m_graphs[i]->GetNumberOfVertices(); j++) {
			int degree = m_graphs[i]->GetDegree(j);

			if(degree!=0) {

			    if(degree==1) {
			        int deadEndState = deadEndStateArray->GetValue(j);

			        if(deadEndState == 1) {
			            file2.width(40);
			            file2 << m_dataSetID;
			            file2.width(10);
			            file2 << 0;
			            file2.width(10);
			            file2 << 1;
			            file2.width(10);
			            file2 << degree;
			            file2.width(20);
			            file2 << radiusArray->GetValue(j);
			            file2.width(20);
			            file2 << -1;
			            file2.width(20);
			            file2 << -1;
			            file2.width(20);
			            file2 << distToCVArray->GetValue(j);
			            file2.width(20);
			            file2 << distToPVArray->GetValue(j);
			            file2.width(20);
			            file2 << minSecOrderLengthArray->GetValue(j);
			            file2.width(20);
			            file2 << maxSecOrderLengthArray->GetValue(j) << std::endl;
			        }
				}
				else if(degree==2) {
				    file2.width(40);
				    file2 << m_dataSetID;
				    file2.width(10);
				    file2 << 0;
				    file2.width(10);
				    file2 << 0;
				    file2.width(10);
				    file2 << degree;
				    file2.width(20);
				    file2 << radiusArray->GetValue(j);
				    file2.width(20);
				    file2 << -1;
				    file2.width(20);
				    file2 << -1;
				    file2.width(20);
				    file2 << distToCVArray->GetValue(j);
				    file2.width(20);
				    file2 << distToPVArray->GetValue(j);
				    file2.width(20);
				    file2 << minSecOrderLengthArray->GetValue(j);
				    file2.width(20);
				    file2 << maxSecOrderLengthArray->GetValue(j) << std::endl;
				}
				else if(degree>2) {
					file2.width(40);
					file2 << m_dataSetID;
					file2.width(10);
					file2 << 1;
					file2.width(10);
					file2 << 0;
					file2.width(10);
					file2 << degree;
                    file2.width(20);
                    file2 << radiusArray->GetValue(j);
                    file2.width(20);
                    file2 << minAngleArray->GetValue(j);
                    file2.width(20);
                    file2 << avgMinAngleArray->GetValue(j);
                    file2.width(20);
                    file2 << distToCVArray->GetValue(j);
                    file2.width(20);
                    file2 << distToPVArray->GetValue(j);
                    file2.width(20);
                    file2 << minSecOrderLengthArray->GetValue(j);
                    file2.width(20);
                    file2 << maxSecOrderLengthArray->GetValue(j) << std::endl;
				}
			}
		}
	}
	
	file2.close();
}


void AnalyzeBileNetworkFilter::WriteBranchInfoToFile()
{
	std::fstream file1, file2, tempfile;
	bool hasHeader = false;

	file1.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_3.txt").c_str(), fstream::in);

	std::string line;
	getline(file1, line);

	if(line.find("dataSet")!=std::string::npos) {
		hasHeader = true;

		tempfile.open((m_networkDataSetPath + "../" + m_output_prefix + "tempfile_3.txt").c_str(), fstream::out);
		tempfile << line;

		while(!file1.eof()) {
			getline(file1, line);
			if(line.find(m_dataSetID)!=0) {
				tempfile << std::endl;
				tempfile << line;
			}
		}
		file1.close();
		tempfile.close();

		std::remove((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_3.txt").c_str());
		std::rename((m_networkDataSetPath + "../" + m_output_prefix + "tempfile_3.txt").c_str(), (m_networkDataSetPath + "../" + m_output_prefix + "Output_file_3.txt").c_str());
	}
	else
		file1.close();

	if(!hasHeader)
		file2.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_3.txt").c_str(), fstream::out);
	else
		file2.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_3.txt").c_str(), fstream::out | fstream::app);

	file2.flags(fstream::left | fstream::scientific);
	if(!hasHeader) {
		file2.width(40);
		file2 << "dataSet";
		file2.width(10);
		file2 << "isIsec";
		file2.width(10);
		file2 << "isDead";
		file2.width(20);
		file2 << "length";
		file2.width(20);
		file2 << "radius";
		file2.width(20);
		file2 << "distToCV";
		file2.width(20);
		file2 << "distToPV" << std::endl;
	}

	for(unsigned int i=0; i<m_isecBranches.size(); i++) {
	    for(unsigned int j=0; j<m_isecBranches[i].size(); j++) {
	        file2.width(40);
	        file2 << m_dataSetID;
	        file2.width(10);
	        file2 << 1;
	        file2.width(10);
	        file2 << 0;
	        file2.width(20);
	        file2 << m_isecBranches[i][j].length;
	        file2.width(20);
	        file2 << m_isecBranches[i][j].radiusMean;
	        file2.width(20);
	        file2 << m_isecBranches[i][j].distanceToCV;
	        file2.width(20);
	        file2 << m_isecBranches[i][j].distanceToPV << std::endl;
	    }
	}

	for(unsigned int i=0; i<m_innerDeadEndBranches.size(); i++) {
	    for(unsigned int j=0; j<m_innerDeadEndBranches[i].size(); j++) {
	        file2.width(40);
	        file2 << m_dataSetID;
	        file2.width(10);
	        file2 << 0;
	        file2.width(10);
	        file2 << 1;
	        file2.width(20);
	        file2 << m_innerDeadEndBranches[i][j].length;
	        file2.width(20);
	        file2 << m_innerDeadEndBranches[i][j].radiusMean;
	        file2.width(20);
	        file2 << m_innerDeadEndBranches[i][j].distanceToCV;
	        file2.width(20);
	        file2 << m_innerDeadEndBranches[i][j].distanceToPV << std::endl;
	    }
	}

//	for(unsigned int i=0; i<m_maskedDeadEndBranches.size(); i++) {
//	    for(unsigned int j=0; j<m_maskedDeadEndBranches[i].size(); j++) {
//	        file2.width(40);
//	        file2 << m_dataSetID;
//	        file2.width(10);
//	        file2 << 0;
//	        file2.width(10);
//	        file2 << 0;
//	        file2.width(10);
//	        file2 << m_maskedDeadEndBranches[i][j].length << std::endl;
//	    }
//	}

	file2.close();
}


void AnalyzeBileNetworkFilter::WriteSecOrderBranchInfoToFile()
{
    std::fstream file1, file2, tempfile;
    bool hasHeader = false;

    file1.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_5.txt").c_str(), fstream::in);

    std::string line;
    getline(file1, line);

    if(line.find("dataSet")!=std::string::npos) {
        hasHeader = true;

        tempfile.open((m_networkDataSetPath + "../" + m_output_prefix + "tempfile_5.txt").c_str(), fstream::out);
        tempfile << line;

        while(!file1.eof()) {
            getline(file1, line);
            if(line.find(m_dataSetID)!=0) {
                tempfile << std::endl;
                tempfile << line;
            }
        }
        file1.close();
        tempfile.close();

        std::remove((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_5.txt").c_str());
        std::rename((m_networkDataSetPath + "../" + m_output_prefix + "tempfile_5.txt").c_str(), (m_networkDataSetPath + "../" + m_output_prefix + "Output_file_5.txt").c_str());
    }
    else
        file1.close();

    if(!hasHeader)
        file2.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_5.txt").c_str(), fstream::out);
    else
        file2.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_5.txt").c_str(), fstream::out | fstream::app);

    file2.flags(fstream::left | fstream::scientific);
    if(!hasHeader) {
        file2.width(40);
        file2 << "dataSet";
        file2.width(30);
        file2 << "#subSegments";
        file2.width(20);
        file2 << "length";
        file2.width(20);
        file2 << "radius";
        file2.width(20);
        file2 << "distToCV";
        file2.width(20);
        file2 << "distToPV";
        file2.width(10);
        file2 << "partitioning" << std::endl;
    }

    for(unsigned int i=0; i<m_secondOrderBranches.size(); i++) {
        for(unsigned int j=0; j<m_secondOrderBranches[i].size(); j++) {
            file2.width(40);
            file2 << m_dataSetID;
            file2.width(30);
            file2 << m_secondOrderBranches[i][j].firstOrderBranches.size();
            file2.width(20);
            file2 << m_secondOrderBranches[i][j].length;
            file2.width(20);
            file2 << m_secondOrderBranches[i][j].radiusMean;
            file2.width(20);
            file2 << m_secondOrderBranches[i][j].distanceToCV;
            file2.width(20);
            file2 << m_secondOrderBranches[i][j].distanceToPV;
            file2.width(10);
            for(unsigned int k=0; k<m_secondOrderBranches[i][j].firstOrderBranches.size(); k++)
            {
                double percentOfSecOrdBranch = m_isecBranches[i][m_secondOrderBranches[i][j].firstOrderBranches[k]].length / m_secondOrderBranches[i][j].length;
                file2.precision(5);
                file2 << percentOfSecOrdBranch << "; ";
            }
            file2 << std::endl;
        }
    }
    //test

    file2.close();
}


void AnalyzeBileNetworkFilter::WriteCycleInfoToFile()
{
    std::vector<unsigned int> cyclesGraph;
    std::vector< std::vector<vtkIdType> > cycles;

    unsigned int number_cycles = 50;
    unsigned int max_draws = 500;
    unsigned int max_cycle_size = 11;

    Tools *tools = new Tools();

//    for(unsigned int i=0; i<m_graphs.size(); i++)
//        std::cout << "m_graphs[" << i << "] has " << m_graphs[i]->GetNumberOfVertices() << " vertices" << std::endl;

    while(max_draws>0 && cycles.size()<number_cycles) {
        vtkIdType node = tools->random->GetRandomUniformInt(0, m_num_nodes-1);
//        std::cout << "node " << node << " drawn: ";

        unsigned int num_graph = 0;
        while(num_graph < m_graphs.size()) {
            if(node >= m_graphs[num_graph]->GetNumberOfVertices())
                node -=  m_graphs[num_graph]->GetNumberOfVertices();
            else
                break;

            num_graph++;
        }
//        std::cout << " it is node " << node << " in graph " << num_graph << std::endl;

        std::vector< std::vector<vtkIdType> > cs;// = AssembleCycles(num_graph, node, max_cycle_size);

//        std::cout << "For this node " << cs.size() << " cycles were found" << std::endl;

        for(unsigned int i=0; i<cs.size(); i++) {
            cycles.push_back(cs[i]);
            cyclesGraph.push_back(num_graph);
        }
        max_draws--;
    }

//    std::cout << "The following " << cycles.size() << " cycles were found: " << std::endl;
//    for(unsigned int i=0; i<cycles.size(); i++) {
//        std::cout << "cycle " << i << " composed of nodes: ";
//        for(unsigned int j=0; j<cycles[i].size(); j++)
//            std::cout << cycles[i][j] << ", ";
//        std::cout << std::endl;
//    }

    //TODO this into CollectCycleInfo
    //evaluate cycle properties
    std::vector<double> lengths;
    std::vector<int> num_isec_nodes;
    std::vector<double> minDistsToCenter;
    std::vector<double> maxDistsToCenter;

    for(unsigned int i=0; i<cycles.size(); i++) {
        vtkDoubleArray* edge_lengths = vtkDoubleArray::SafeDownCast( m_graphs[cyclesGraph[i]]->GetEdgeData()->GetArray(mpAnno->GetArrayNameEdgeLength().c_str()) );

        double length = 0.0;
        int num_isec = 0;
        std::vector< std::vector<double> > scaledPoints;
        double center[3] = {0, 0, 0};

//        std::cout << "cycle " << i << " composed of nodes: ";
//        for(unsigned int j=0; j<cycles[i].size(); j++)
//            std::cout << cycles[i][j] << ", ";

        for(unsigned int j=0; j<cycles[i].size(); j++) {
            if(j>0)
                length += edge_lengths->GetValue( m_graphs[cyclesGraph[i]]->GetEdgeId(cycles[i][j], cycles[i][j-1]));

            if(m_graphs[cyclesGraph[i]]->GetDegree(cycles[i][j]) > 2)
                num_isec++;

            double a[3];
            m_graphs[cyclesGraph[i]]->GetPoint(cycles[i][j], a);

            a[0] *= m_vox_spacing[0];
            a[1] *= m_vox_spacing[1];
            a[2] *= m_vox_spacing[2];

            vtkMath::Add(a, center, center);
//            std::vector<double> aa (a[0], a[1], a[2]);
            scaledPoints.push_back( std::vector<double> (a, a + sizeof(a) / sizeof(double)) );
        }
        length += edge_lengths->GetValue( m_graphs[cyclesGraph[i]]->GetEdgeId(cycles[i][cycles[i].size()-1], cycles[i][0]));
        lengths.push_back(length);
        num_isec_nodes.push_back(num_isec);

        center[0] /= cycles[i].size();
        center[1] /= cycles[i].size();
        center[2] /= cycles[i].size();

//        std::cout << " has length = " << length;
//        std::cout << " and nodes = " << cycles[i].size();
//        std::cout << " and num_isec_nodes = " << num_isec;
//        std::cout << " and has center of mass at (" << center[0] << ", " << center[1] << ", " << center[2] << ")";

        double minDistToCenter, maxDistToCenter;

        for(unsigned int j=0; j<cycles[i].size(); j++) {
            double vec[3];

            double poi[3] = {scaledPoints[j][0], scaledPoints[j][1], scaledPoints[j][2]};
            vtkMath::Subtract(center, poi, vec);
            double dist = vtkMath::Norm(vec);

            if(j==0) {
                minDistToCenter = dist;
                maxDistToCenter = dist;
            }
            else if(dist < minDistToCenter) minDistToCenter = dist;
            else if(dist > maxDistToCenter) maxDistToCenter = dist;
        }
        minDistsToCenter.push_back(minDistToCenter);
        maxDistsToCenter.push_back(maxDistToCenter);

//        std::cout << " with a minDistToCenter = " << minDistToCenter << ", and a maxDistToCenter = " << maxDistToCenter << std::endl;
    }

    //write cycle properties to output file 4
    std::fstream file1, file2, tempfile;
    bool hasHeader = false;

    file1.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_4.txt").c_str(), fstream::in);

    std::string line;
    getline(file1, line);

    if(line.find("dataSet")!=std::string::npos) {
        hasHeader = true;

        tempfile.open((m_networkDataSetPath + "../" + m_output_prefix + "tempfile_4.txt").c_str(), fstream::out);
        tempfile << line;

        while(!file1.eof()) {
            getline(file1, line);
            if(line.find(m_dataSetID)!=1) {
                tempfile << std::endl;
                tempfile << line;
            }
        }
        file1.close();
        tempfile.close();

        std::remove((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_4.txt").c_str());
        std::rename((m_networkDataSetPath + "../" + m_output_prefix + "tempfile_4.txt").c_str(), (m_networkDataSetPath + "../" + m_output_prefix + "Output_file_4.txt").c_str());
    }
    else
        file1.close();

    if(!hasHeader)
        file2.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_4.txt").c_str(), fstream::out);
    else
        file2.open((m_networkDataSetPath + "../" + m_output_prefix + "Output_file_4.txt").c_str(), fstream::out | std::ios::app);

    file2.flags(fstream::left | fstream::scientific);
    if(!hasHeader) {
        file2.width(40);
        file2 << "dataSet";
        file2.width(25);
        file2 << "length";
        file2.width(10);
        file2 << "#nodes";
        file2.width(10);
        file2 << "#isecs";
        file2.width(10);
        file2 << "#polyEdges";
        file2.width(20);
        file2 << "minRad";
        file2.width(20);
        file2 << "maxRad";
        file2.width(10);
        file2 << "avg_angle";
        file2.width(10);
        file2 << "avg_branch_length";
        file2.width(10);
        file2 << "avg_polyEdges_length"<< std::endl;
    }

    for(unsigned int i=0; i<cycles.size(); i++)
    {
        file2.width(40);
        file2 << m_dataSetID;
        file2.width(25);
        file2 << lengths[i];
        file2.width(10);
        file2 << cycles[i].size();
        file2.width(10);
        file2 << num_isec_nodes[i];
        file2.width(10);
        file2 << -1;                     //TODO
        file2.width(20);
        file2 << minDistsToCenter[i];
        file2.width(20);
        file2 << maxDistsToCenter[i] << std::endl;
    }
}


void AnalyzeBileNetworkFilter::Update()
{
    ParseParameterContext();
	
	CollectBasicImageInformation();

    BasePipeline<3>::ReadImage(m_networkDataSetPath + m_networkDataSetName + m_networkDataSetFileExtension, m_networkBin, m_vox_spacing);
	if(m_withNecReg) {
	    long double necRegVolume = BasePipeline<3>::ComputeVolumeOfRegions(m_necRegDataSetPath + m_necRegDataSetName + m_necRegDataSetFileExtension, m_vox_spacing, m_voxel_volume);
	    m_effective_dataset_volume -= necRegVolume;
	}

	if(m_withCell)
	    BasePipeline<3>::ReadImage(m_cellDataSetPath + m_cellDataSetName + m_cellDataSetFileExtension, m_cellBin, m_vox_spacing);

	if(m_withCV) {
	    BasePipeline<3>::ReadImage(m_cvDataSetPath + m_cvDataSetName + m_cvDataSetFileExtension, m_cvBin, m_vox_spacing);
	    m_cv_volume = BasePipeline<3>::ComputeVolumeOfRegions(m_cvBin, m_voxel_volume);
	    m_effective_dataset_volume -= m_cv_volume;
	}
	else {
        m_cvBin = CScalarVoImageType::New();
        m_cvDistMap = FScalarVoImageType::New();
        BasePipeline<3>::CreateImage(m_cvBin, m_networkBin->GetLargestPossibleRegion(), 0, m_vox_spacing);
        m_cv_volume = 0;
	}

	if(m_withPV) {
	    BasePipeline<3>::ReadImage(m_pvDataSetPath + m_pvDataSetName + m_pvDataSetFileExtension, m_pvBin, m_vox_spacing);
	    m_pv_volume = BasePipeline<3>::ComputeVolumeOfRegions(m_pvBin, m_voxel_volume);
	    m_effective_dataset_volume -= m_pv_volume;
	}
	else {
	    m_pvBin = CScalarVoImageType::New();
	    m_pvDistMap = FScalarVoImageType::New();
	    BasePipeline<3>::CreateImage(m_pvBin, m_networkBin->GetLargestPossibleRegion(), 0, m_vox_spacing);
	    m_pv_volume = 0;
	}

	m_network_vol_count = BasePipeline<3>::ComputeVolumeOfRegions(m_networkBin, m_voxel_volume);
	m_network_vol_count_per_vol = m_network_vol_count / m_effective_dataset_volume;

    if(m_withRadiusMeasurement) BasePipeline<3>::BuildDistanceMap(m_networkBin, m_networkDistMap, m_vox_spacing, false);
    if(m_withCell)              BasePipeline<3>::BuildLabelImage(m_cellBin, m_cellLabelImage, m_vox_spacing, false);
    if(m_withCV)                BasePipeline<3>::BuildDistanceMap(m_cvBin, m_cvDistMap, m_vox_spacing, true);
    else                        BasePipeline<3>::CreateImage(m_cvDistMap, m_networkBin->GetLargestPossibleRegion(), m_maxCVeinConnectednessDistance, m_vox_spacing);
    if(m_withPV)                BasePipeline<3>::BuildDistanceMap(m_pvBin, m_pvDistMap, m_vox_spacing, true);
    else                        BasePipeline<3>::CreateImage(m_pvDistMap, m_networkBin->GetLargestPossibleRegion(), m_maxPVeinConnectednessDistance, m_vox_spacing);

    mpAnno = new GraphAnnotationHelper();
    BuildBranchMapsAndPopulateGraphBranchArrays();

	CollectBasicNodeEdgeInformation();

	CollectFirstOrderBranchInformation();
	CollectSecondOrderBranchInformation();
//	CollectCycleInformation();  //TODO

	std::cout << "Save graphs?" << std::endl;
	if(m_saveGraphsWithArrays) {
	    std::cout << "yes!" << std::endl;
	    std::cout << "m_graphDataSetPath = " << m_graphDataSetPath << std::endl;

	    std::string path, filename, filenamePrefix, ext;
	    int startIndex = 0;
	    FilenameParser::ParseFilename(m_graphDataSetPath, path, filename, ext);
	    std::cout << "path = " << path << ", filename = " << filename << ", ext = " << ext << std::endl;
	    FilenameParser::IsSequence(filename, filenamePrefix, startIndex);
	    std::cout << "filenamePrefix = " << filenamePrefix << std::endl;
	    SaveGraphs(path, filenamePrefix, ext);
	}

	if(m_saveAnalysisFile)
	    WriteToFile();

	std::cout << "num_branches " << m_num_branches << "; num_branches_per_vol " << m_num_branches_per_vol << "; num_edges_per_branch " << m_num_edges_per_branch << "; mean_branch_length " << m_mean_branch_length << std::endl;
	std::cout << "mean radius (micron) " << m_mean_radius << "; m_num_neg_radius " << m_num_neg_radius << std::endl;
	std::cout << "num_firstOrderIsecBranches " << m_num_branches_between_intersections << "; mean_length_firstOrderIsecBranches " << m_mean_branch_between_intersections_length << std::endl;
	std::cout << "num_firstOrderDEndBranches " << m_num_branches_with_dead_end << "; mean_length_firstOrderDEndBranches "<< m_mean_branch_with_dead_end_length << "; num_firstOrderRestBranches (dEnd in contact with dataset border) " << m_num_branches_rest << std::endl;
	std::cout << "num_secOrderBranches " << m_num_secOrderBranches << "\n" << std::endl;

	m_networkBin->ReleaseData();
	m_networkBin = NULL;
	m_networkDistMap->ReleaseData();
	m_networkDistMap = NULL;

	for(unsigned int i=0; i<m_deadEndMasks.size(); i++) {
	    m_deadEndMasks[i]->ReleaseData();
	    m_deadEndMasks[i] = NULL;
	}

	if(m_withCell) {
	    m_cellBin->ReleaseData();
	    m_cellBin = NULL;
	    m_cellLabelImage->ReleaseData();
	    m_cellLabelImage = NULL;
	}
	if(m_withCV) {
	    m_cvBin->ReleaseData();
	    m_cvBin = NULL;
	    m_cvDistMap->ReleaseData();
	    m_cvDistMap = NULL;
	}
	if(m_withPV) {
	    m_pvBin->ReleaseData();
	    m_pvBin = NULL;
	    m_pvDistMap->ReleaseData();
	    m_pvDistMap = NULL;
	}
}
