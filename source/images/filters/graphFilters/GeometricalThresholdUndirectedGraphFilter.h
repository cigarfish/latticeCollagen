///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GeometricalThresholdUndirectedGraphFilter.h                          //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-12-03                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GEOMETRICALTHRESHOLDUNDIRECTEDGRAPHFILTER_H_
#define GEOMETRICALTHRESHOLDUNDIRECTEDGRAPHFILTER_H_

#include "vtkGraphAlgorithm.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

class GeometricalThresholdUndirectedGraphFilter : public vtkGraphAlgorithm
{
public:
    static GeometricalThresholdUndirectedGraphFilter* New();
    vtkTypeMacro(GeometricalThresholdUndirectedGraphFilter, vtkGraphAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // Get/Set lower threshold. This would be the value against which
    // edge or vertex data array value will be compared.
    vtkGetMacro(AngleThreshold, float);
    vtkSetMacro(AngleThreshold, float);

    vtkGetMacro(TestOutput, bool);
    vtkSetMacro(TestOutput, bool);

    vtkGetMacro(WithAngleAnnotation, bool);
    vtkSetMacro(WithAngleAnnotation, bool);

    // Description:
    // Specify the first vtkGraph input and the second vtkSelection input.
    int FillInputPortInformation(int port, vtkInformation* info);
protected:
    GeometricalThresholdUndirectedGraphFilter();
    ~GeometricalThresholdUndirectedGraphFilter();

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    float AngleThreshold;
    bool TestOutput;
    bool WithAngleAnnotation;

    GeometricalThresholdUndirectedGraphFilter(const GeometricalThresholdUndirectedGraphFilter&);  // Not implemented
    void operator=(const GeometricalThresholdUndirectedGraphFilter&);                             // Not implemented
};


#endif /* GEOMETRICALTHRESHOLDUNDIRECTEDGRAPHFILTER_H_ */
