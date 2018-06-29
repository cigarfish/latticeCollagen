///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  TopologicalPruneGraphFilter.h                                        //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-02-08                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef TOPOLOGICALPRUNEGRAPHFILTER_H_
#define TOPOLOGICALPRUNEGRAPHFILTER_H_

#include <vector>

#include "itkImage.h"
#include "itkSmartPointer.h"

#include "vtkGraphAlgorithm.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include "GraphHandlingHelper.h"


class TopologicalPruneGraphFilter  : public vtkGraphAlgorithm
{
    typedef unsigned char                   CScalarPixelType;
    typedef float                           FScalarPixelType;

    typedef itk::Image<CScalarPixelType, 3> CScalarVoImageType;
    typedef itk::Image<FScalarPixelType, 3> FScalarVoImageType;


public:
    static TopologicalPruneGraphFilter* New();
    vtkTypeMacro(TopologicalPruneGraphFilter, vtkGraphAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    vtkGetMacro(DeadEndLengthThreshold, double);
    vtkSetMacro(DeadEndLengthThreshold, double);

    vtkGetMacro(PruneAll, bool);
    vtkSetMacro(PruneAll, bool);

    vtkGetMacro(TestOutput, bool);
    vtkSetMacro(TestOutput, bool);

    void SetVeinDistanceMaps(FScalarVoImageType::Pointer cvDistMap, FScalarVoImageType::Pointer pvDistMap) { mpCVDistMap = cvDistMap; mpPVDistMap = pvDistMap; mSize = mpCVDistMap->GetLargestPossibleRegion().GetSize(); };
    void SetVeinConnectednessDistances(double maxCVDist, double maxPVDist) { m_maxCVeinConnectednessDistance = maxCVDist; m_maxPVeinConnectednessDistance = maxPVDist; };
    void SetImageSpacing(FScalarVoImageType::SpacingType spacing) { mSpacing = spacing; };
    void SetAssumedStandardRadius(double rad) { mAssumedStandardRadius = rad; };

    // Description:
    // Specify the first vtkGraph input and the second vtkSelection input.
    int FillInputPortInformation(int port, vtkInformation* info);
protected:
    TopologicalPruneGraphFilter();
    ~TopologicalPruneGraphFilter();

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    void BuildBranchMaps(vtkGraph *input, bool removeBorderDeadEnds);
    void IdentifyRemovableDeadEnds(vtkGraph *input, std::vector<FirstOrderBranch> &deadEndBranch, std::vector<FirstOrderBranch> &removableDeadEndBranch,
            std::vector<FirstOrderBranch> &nonRemovableDeadEndBranch, bool removeBorderDeadEnds);

    TopologicalPruneGraphFilter(const TopologicalPruneGraphFilter&); // Not implemented
    void operator=(const TopologicalPruneGraphFilter&);   // Not implemented

    double DeadEndLengthThreshold;
    bool PruneAll;
    bool TestOutput;

    double mAssumedStandardRadius;

    FScalarVoImageType::SizeType mSize;
    FScalarVoImageType::SpacingType mSpacing;
    FScalarVoImageType::Pointer mpCVDistMap;
    FScalarVoImageType::Pointer mpPVDistMap;

    double m_maxCVeinConnectednessDistance;
    double m_maxPVeinConnectednessDistance;

    std::vector<FirstOrderBranch> mRemovableDeadEndBranches;
    std::vector<FirstOrderBranch> mNonRemovableDeadEndBranches;
    std::vector<FirstOrderBranch> mIsecBranches;
};

#endif /* TOPOLOGICALPRUNEGRAPHFILTER_H_ */
