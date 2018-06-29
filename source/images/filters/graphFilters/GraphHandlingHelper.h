///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphHandlingHelper.h                                                //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-07-07                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHHANDLINGHELPER_H_
#define GRAPHHANDLINGHELPER_H_

#include <stack>
#include <string>
#include <vector>

#include <itkImage.h>
#include <itkSmartPointer.h>

#include <vtkSmartPointer.h>
#include <vtkGraph.h>

#include "../tools/GraphAnnotationHelper.h"


struct FirstOrderBranch
{
public:
    unsigned int id;
//    unsigned int graphNr;

    std::vector<vtkIdType> nodes;
    std::vector<unsigned long int> nodesPIDs;
    vtkIdType firstVert;
    vtkIdType lastVert;
    unsigned long int firstVertPID;
    unsigned long int lastVertPID;

    double middlePoint[3];
    long double length;
    long double radiusMean;
    long double distanceToCV;
    long double distanceToPV;
};

struct SecondOrderBranch
{
public:
    unsigned int id;
    unsigned int graphNr;

    std::vector<int> firstOrderBranches;     //position in m_isecBranches vector

    vtkIdType firstVert;
    vtkIdType lastVert;

    long double length;
    long double radiusMean;
    long double distanceToCV;
    long double distanceToPV;
};


class GraphHandlingHelper
{
    typedef unsigned char                   CScalarPixelType;
    typedef float                           FScalarPixelType;

    typedef itk::Image<CScalarPixelType, 3> CScalarVoImageType;
    typedef itk::Image<FScalarPixelType, 3> FScalarVoImageType;

public:
    GraphHandlingHelper();
    virtual ~GraphHandlingHelper();

    static std::map<vtkIdType, int> AssembleFirstOrderBranches(vtkSmartPointer<vtkGraph> input, std::vector<FirstOrderBranch> &deadEndBranches, std::vector<FirstOrderBranch> &isecBranches);
    static void AssembleSecondOrderBranches(vtkSmartPointer<vtkGraph> input, std::vector<SecondOrderBranch> &secOrderBranches, std::vector<FirstOrderBranch> &intersectionBranches,
            std::vector<FirstOrderBranch> &innerDeadEndBranches, std::vector<FirstOrderBranch> &maskedDeadEndBranches, std::map<vtkIdType, int> &vertexIdToSecOState, std::map<vtkIdType, int> &edgeIdToBranchId);

    static void GetMiddlepointOfVertices(vtkSmartPointer<vtkGraph> input, vtkIdType first, vtkIdType last, double &x, double &y, double &z);
    static bool IsVertexAtDatasetBorder(vtkSmartPointer<vtkGraph> input, vtkIdType vertexID, double minDistanceInUnit, FScalarVoImageType::SpacingType spacing, FScalarVoImageType::SizeType datasetSize);
    static bool IsVertexMasked(vtkSmartPointer<vtkGraph> input, vtkIdType vertexID, CScalarVoImageType::Pointer maskImage, CScalarPixelType maskingValue = 255);
    static int GetVeinConnectednessStatus(vtkSmartPointer<vtkGraph> input, vtkIdType vertexID, FScalarVoImageType::Pointer cvDistMap, FScalarVoImageType::Pointer pvDistMap,
            double maxCVDist, double maxPVDist, bool hasCV, bool hasPV);
    static void ComputeMinAngles(vtkSmartPointer<vtkGraph> input, vtkIdType vertexA, float lengthOfBranchVector, std::vector<float> &angles);

};

#endif /* GRAPHHANDLINGHELPER_H_ */
