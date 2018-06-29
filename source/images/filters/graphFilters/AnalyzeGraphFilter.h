///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeGraphFilter.h                                                 //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-01-17                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ANALYZEGRAPHFILTER_H_
#define ANALYZEGRAPHFILTER_H_

#include <string>
#include <vector>

#include <QFileInfo>
#include <QString>

#include "itkImage.h"

#include <vtkSmartPointer.h>
#include <vtkUndirectedGraph.h>

#include "GraphHandlingHelper.h"
#include "../../pipelines/BasePipeline.h"
#include "../../tools/GraphAnnotationHelper.h"
#include "../../../tools/input/FilenameParser.h"
#include "../../../tools/parameters/CSParameterChoice.h"

class CSParameterContext;


class AnalyzeGraphFilter
{
protected:
    typedef unsigned char                               CScalarPixelType;
    typedef float                                       FScalarPixelType;

    typedef itk::Image<CScalarPixelType, 3>             CScalarVoImageType;
    typedef itk::Image<FScalarPixelType, 3>             FScalarVoImageType;

public:
	AnalyzeGraphFilter();
	virtual ~AnalyzeGraphFilter();

	void SetParameterContext(CSParameterContext *paramContext) { m_paramContext = paramContext; };
	void SetGraphsToAnalyze(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs) { m_graphs = graphs; };
	void SaveGraphsWithArrays(bool save) { m_saveGraphsWithArrays = save; };
	void SetGraphFilename(std::string graphFilename) { m_graphDataSetPath = graphFilename; };
	void SetRadiusMeasurement(bool active) { m_withRadiusMeasurement = active; };

	void AddDeadEndMaskImage(std::string file);

	std::vector< vtkSmartPointer<vtkUndirectedGraph> > GetGraphs() { return m_graphs; };

	void Update();

protected:
	void ParseParameterContext();

	double MeasureDistance(vtkGraph *input, vtkIdType vertex, FScalarVoImageType::Pointer image);

	void CollectBasicImageInformation();
	void CollectBasicNodeEdgeInformation();

	void SaveGraphs(std::string path, std::string filename , std::string ext);

	CSParameterContext *m_paramContext;

	std::string m_dataSetID;

	QString m_networkDataSetFullFilename;
	QFileInfo m_infoFile;
    std::string m_networkDataSetPath;
    std::string m_networkDataSetName;
    std::string m_networkDataSetFileExtension;

    bool m_withRadiusMeasurement;

    bool m_saveGraphsWithArrays;
    std::string m_graphDataSetPath;

    bool m_saveAnalysisFile;

    CScalarVoImageType::Pointer m_networkBin;
    FScalarVoImageType::Pointer m_networkDistMap;
	std::vector< vtkSmartPointer<vtkUndirectedGraph> > m_graphs;

	GraphAnnotationHelper* mpAnno;

	long int m_num_dim;

	CScalarVoImageType::SizeType m_dim;
	CScalarVoImageType::SpacingType m_vox_spacing;
	long double m_voxel_volume;
	long double m_dataset_volume;
	long double m_effective_dataset_volume;

    long double m_rad;

    bool m_withDeadEndMasking;
    std::vector< CScalarVoImageType::Pointer > m_deadEndMasks;

	long int m_num_nodes;
	long int m_num_intersec_nodes;
	long int m_num_reg_nodes;
	long int m_num_total_dead_end_nodes;
	long int m_num_within_dead_end_nodes;
	long int m_num_zero_nodes;
	long int m_num_edges;

	long double m_num_nodes_per_vol;
	long double m_num_intersec_nodes_per_vol;
	long double m_num_within_dead_end_nodes_per_vol;
	long double m_mean_edges_per_intersec_node;
	long double m_vari_edges_per_intersec_node;
	long double m_total_edge_length;
	long double m_total_edge_length_per_vol;
	long double m_mean_radius;
	long double m_num_neg_radius;

	long double m_network_vol_calc;
    long double m_network_vol_calc_per_vol;
    long double m_network_vol_count;
    long double m_network_vol_count_per_vol;

    long double m_mean_isecBranch_minAngle;
    long double m_mean_isecBranch_minAvgAngle;
};

#endif /* ANALYZEGRAPHFILTER_H_ */
