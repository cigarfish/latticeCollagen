///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CollapseIntersectionNodesFilter.h                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-01-31                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef COLLAPSEINTERSECTIONNODESFILTER_H_
#define COLLAPSEINTERSECTIONNODESFILTER_H_

#include <stack>

#include "itkImage.h"
#include "itkSmartPointer.h"

#include "vtkGraphAlgorithm.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"

#include "GraphHandlingHelper.h"


class CollapseIntersectionNodesFilter  : public vtkGraphAlgorithm
{
    typedef float                           FScalarPixelType;
    typedef itk::Image<FScalarPixelType, 3> FScalarVoImageType;


public:
    static CollapseIntersectionNodesFilter* New();
    vtkTypeMacro(CollapseIntersectionNodesFilter, vtkGraphAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    vtkGetMacro(IsecDistanceThreshold, int);
    vtkSetMacro(IsecDistanceThreshold, int);

    void SetNetworkDistanceMap(FScalarVoImageType::Pointer networkDistanceMap) { mpNetworkDistMap = networkDistanceMap; };
    void SetImageSpacing(FScalarVoImageType::SpacingType spacing) { mSpacing = spacing; };

    vtkGetMacro(TestOutput, bool);
    vtkSetMacro(TestOutput, bool);

    // Description:
    // Specify the first vtkGraph input and the second vtkSelection input.
    int FillInputPortInformation(int port, vtkInformation* info);
protected:
    CollapseIntersectionNodesFilter();
    ~CollapseIntersectionNodesFilter();

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    void BuildBranchMaps(vtkGraph *input);
    void IdentifyCollapsableIsecBranches(vtkGraph *input, std::vector<FirstOrderBranch> &isecBranches, std::vector<FirstOrderBranch> &collapseableBranches,
            std::vector<FirstOrderBranch> &unalteredBranches);
    double MeasureDistance(FScalarVoImageType::Pointer image, double *pos);

    CollapseIntersectionNodesFilter(const CollapseIntersectionNodesFilter&); // Not implemented
    void operator=(const CollapseIntersectionNodesFilter&);   // Not implemented

    double IsecDistanceThreshold;
    bool TestOutput;

    FScalarVoImageType::SpacingType mSpacing;
    FScalarVoImageType::Pointer mpNetworkDistMap;

    std::vector<FirstOrderBranch> mIsecCollapseBranch;
    std::vector<FirstOrderBranch> mIsecUnalteredBranch;
    std::vector<FirstOrderBranch> mDeadEndBranch;
};

#endif /* COLLAPSEINTERSECTIONNODESFILTER_H_ */
