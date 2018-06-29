///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ResampleUndirectedGraphFilter.h                                      //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-11-21                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef RESAMPLEUNDIRECTEDGRAPHFILTER_H_
#define RESAMPLEUNDIRECTEDGRAPHFILTER_H_

#include <vector>

#include "vtkGraphAlgorithm.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include "GraphHandlingHelper.h"

/*!
  \brief A class for resampling vtkUndirectedGraph objects.
  A resampling factor can be specified, default is 2.
  \n The filter preserves vertices which have not exactly two edges (non-regular vertices). On paths composed of regular nodes the algorithm preserves every _resampling_factor_th vertex.
  \n Therefore the filter is preserving the topological structure, but not the geometric structure.
*/
class ResampleUndirectedGraphFilter : public vtkGraphAlgorithm
{
public:
    static ResampleUndirectedGraphFilter* New();
    vtkTypeMacro(ResampleUndirectedGraphFilter, vtkGraphAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // Get/Set lower threshold. This would be the value against which
    // edge or vertex data array value will be compared.
    vtkGetMacro(ResamplingFactor, int);
    vtkSetMacro(ResamplingFactor, int);

    vtkGetMacro(MaxResamplingDist, double);
    vtkSetMacro(MaxResamplingDist, double);

    vtkGetMacro(TestOutput, bool);
    vtkSetMacro(TestOutput, bool);

    // Description:
    // Specify the first vtkGraph input and the second vtkSelection input.
    int FillInputPortInformation(int port, vtkInformation* info);
protected:
    ResampleUndirectedGraphFilter();
    ~ResampleUndirectedGraphFilter();

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    void BuildBranchMaps(vtkGraph *input);

    ResampleUndirectedGraphFilter(const ResampleUndirectedGraphFilter&); // Not implemented
    void operator=(const ResampleUndirectedGraphFilter&);   // Not implemented

    int ResamplingFactor;
    double MaxResamplingDist;
    bool TestOutput;

    std::vector<FirstOrderBranch> branches;
};


#endif /* RESAMPLEUNDIRECTEDGRAPHFILTER_H_ */
