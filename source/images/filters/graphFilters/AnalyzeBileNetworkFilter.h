///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeBileNetworkFilter.h                                           //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-01-17                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ANALYZEBILENETWORKFILTER_H_
#define ANALYZEBILENETWORKFILTER_H_

#include "AnalyzeGraphFilter.h"

#include <stack>

#include <itkConstantBoundaryCondition.h>
#include <itkConstShapedNeighborhoodIterator.h>
#include <itkConstNeighborhoodIterator.h>


class AnalyzeBileNetworkFilter : public AnalyzeGraphFilter
{
    typedef long                            LScalarPixelType;
    typedef itk::Image<LScalarPixelType, 3> LScalarVoImageType;

    typedef itk::ConstantBoundaryCondition<LScalarVoImageType>                              BoundaryConditionType;
    typedef itk::ConstNeighborhoodIterator<LScalarVoImageType, BoundaryConditionType>       IteratorType2;

public:
	AnalyzeBileNetworkFilter();
	~AnalyzeBileNetworkFilter();

	void SetVeinConnectednessMeasurement(bool active) { m_withVeinConnectednessMeasurement = active; };

	int GetNumberDeadEndBranches() { return m_innerDeadEndBranches.size(); };
	int GetNumberIntersectionBranches() { return m_isecBranches.size(); };
	int GetNumberSecondOrderBranches() { return m_secondOrderBranches.size(); };

	std::vector< std::vector<vtkIdType> > AssembleCycles(int subgraph, vtkIdType node, int cycleSize);

	void Update();

private:
	void ParseParameterContext();

	void BuildBranchMapsAndPopulateGraphBranchArrays();
	void IdentifyAnalysableDeadEnds(vtkGraph *input, std::vector<FirstOrderBranch> &deadEndBranches,
	        std::vector<FirstOrderBranch> &analysableDeadEndBranches, std::vector<FirstOrderBranch> &nonAnalysableDeadEndBranches);
	void AnalyzeBranchesAndPopulateGraphArrays(vtkGraph *input, std::vector<FirstOrderBranch> &isecBranches, std::vector<FirstOrderBranch> &deadEndBranches, std::map<vtkIdType, int> edgeIdToBranchId);
	void AnalyzeBranchesAndPopulateGraphArrays(vtkGraph *input, std::vector<SecondOrderBranch> &secOrderBranches, std::map<vtkIdType, int> &vertexIdToSecOState, std::map<vtkIdType, int> &edgeIdToSecOBranchId);
	void AddBranchStateGraphArrays(vtkGraph *input, std::vector<FirstOrderBranch> &isecBranches, std::vector<FirstOrderBranch> &innerDeadEndBranches);

	void MeasureBranchVeinDistances(vtkGraph *input, FirstOrderBranch &branch);
	void MeasureBranchVeinDistances(SecondOrderBranch &branch);
	std::set<int> FindAdjacentCells(vtkGraph *input, vtkIdType vertex, float radius);

	void CollectFirstOrderBranchInformation();
	void CollectSecondOrderBranchInformation();
	void CollectCycleInformation();

	void WriteToFile();
	void WriteGeneralInfoToFile();
	void WriteNodeInfoToFile();
	void WriteBranchInfoToFile();
	void WriteSecOrderBranchInfoToFile();
	void WriteCycleInfoToFile();


	std::string m_output_prefix;

	QString m_necRegDataSetFullFilename;
	QFileInfo m_infoNecRegDataSetFullFilename;
    std::string m_necRegDataSetPath;
    std::string m_necRegDataSetName;
    std::string m_necRegDataSetFileExtension;

    QString m_CellDataSetFullFilename;
    QFileInfo m_infoCellDataSetFullFilename;
    std::string m_cellDataSetPath;
    std::string m_cellDataSetName;
    std::string m_cellDataSetFileExtension;

	QString m_CVDataSetFullFilename;
	QFileInfo m_infoCVDataSetFullFilename;
    std::string m_cvDataSetPath;
    std::string m_cvDataSetName;
    std::string m_cvDataSetFileExtension;

	QString m_PVDataSetFullFilename;
	QFileInfo m_infoPVDataSetFullFilename;
    std::string m_pvDataSetPath;
    std::string m_pvDataSetName;
    std::string m_pvDataSetFileExtension;

    CScalarVoImageType::Pointer m_cellBin;
    LScalarVoImageType::Pointer m_cellLabelImage;

    CScalarVoImageType::Pointer m_cvBin;
    CScalarVoImageType::Pointer m_pvBin;

    FScalarVoImageType::Pointer m_cvDistMap;
    FScalarVoImageType::Pointer m_pvDistMap;

	bool m_withNecReg;
	bool m_withCell;
	bool m_withCV;
	bool m_withPV;
    bool m_withVeinConnectednessMeasurement;

    double m_maxCVeinConnectednessDistance;
    double m_maxPVeinConnectednessDistance;

    int m_maxAdjacentCells;
    int m_maxMinAngles;
    std::vector< std::multimap<vtkIdType,int> > m_adjacentCellsPerGraph;
    std::vector< std::multimap<vtkIdType,float> > m_minAnglesPerGraph;

    //First & second order branches------------------
    std::map< int, std::vector<FirstOrderBranch> > m_isecBranches;          //graphNr, isecBranches
    std::map< int, std::vector<FirstOrderBranch> > m_innerDeadEndBranches;   //graphNr, deBranches
    std::map< int, std::vector<FirstOrderBranch> > m_maskedDeadEndBranches;          //graphNr, restBranches

    std::map< int, std::vector<SecondOrderBranch> > m_secondOrderBranches;     //graphNr, secOBranches
	//-----------------------------------------------

	long double m_cv_volume;
    long double m_pv_volume;

	long int m_num_branches;
	long int m_num_branches_with_dead_end;
	long int m_num_branches_between_intersections;
	long int m_num_branches_rest;                               //dead end branches in contact with data set border and therefore discarded from quantification
	long int m_num_secOrderBranches;                            //branches between intersection nodes with at least three intersection branches attached to

	long double m_num_branches_per_vol;
	long double m_num_branches_with_dead_end_per_vol;
	long double m_num_branches_between_intersections_per_vol;
	long double m_num_secOrderBranches_per_vol;

	long double m_num_isecBranches_per_secOrderBranch;
	long double m_num_edges_per_branch;

	long double m_mean_branch_length;
	long double m_mean_branch_between_intersections_length;
	long double m_mean_branch_with_dead_end_length;
	long double m_mean_secOrderBranch_length;
	long double m_mean_minSecOrderBranch_length;
	long double m_mean_maxSecOrderBranch_length;

	long double m_mean_branch_between_intersections_radius;
	long double m_mean_branch_with_dead_end_radius;
	long double m_mean_secOrderBranch_radius;
};

#endif /* ANALYZEBILENETWORKFILTER_H_ */
