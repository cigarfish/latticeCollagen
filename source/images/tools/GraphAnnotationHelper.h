///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphAnnotationHelper.h                                              //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-01-06                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHANNOTATIONHELPER_H_
#define GRAPHANNOTATIONHELPER_H_

#include <map>
#include <string>
#include <vector>

#include "vtkGraph.h"
#include "vtkUndirectedGraph.h"
#include "vtkSmartPointer.h"

class GraphAnnotationHelper
{
public:
    GraphAnnotationHelper();
    virtual ~GraphAnnotationHelper();

    void AddPredefinedAnnotations(vtkSmartPointer<vtkGraph> graph);

    //TODO: write the overload Custom function you need
    void AddCustomVertexAnnotation(vtkSmartPointer<vtkGraph> graph, std::string arrayName, std::map<vtkIdType,int> &data, int default_value = 0);
    void AddCustomVertexAnnotation(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, std::string arrayNameSuffix, std::vector< std::multimap<vtkIdType,int> > &data, int maxArrays, int default_value = 0);
    void AddCustomVertexAnnotation(vtkSmartPointer<vtkGraph> graph, std::string arrayName, std::map<vtkIdType,float> &data, float default_value = 0);
    void AddCustomVertexAnnotation(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, std::string arrayNameSuffix, std::vector< std::multimap<vtkIdType,float> > &data, int maxArrays, float default_value = 0);
    void AddCustomEdgeAnnotation(vtkSmartPointer<vtkGraph> graph, std::string arrayName, std::map<vtkIdType, int> &data, int default_value = 0);
    void AddCustomEdgeAnnotation(vtkSmartPointer<vtkGraph> graph, std::string arrayName, std::map<vtkIdType, float> &data, float default_value = 0);
    void AddCustomEdgeAnnotation(vtkSmartPointer<vtkGraph> graph, std::string arrayName, std::map<vtkIdType, double> &data, double default_value = 0);

    void EnableVertexIDAnnotation() { addVertexID = true; };
    void EnableVertexDegreeAnnotation() { addVertexDegree = true; };
    void EnableEdgeIDAnnotation() { addEdgeID = true; };
    void EnableEdgeLengthAnnotation(double x_size, double y_size, double z_size);

    void SetArrayNameVertexID(std::string name) { arrayName_vertexID = name; };
    void SetArrayNameVertexDegree(std::string name) { arrayName_vertexDegree = name; };
    void SetArrayNameEdgeID(std::string name) { arrayName_edgeID = name; };
    void SetArrayNameEdgeLength(std::string name) { arrayName_edgeLength = name; };

    void SetVertexIDOffset(unsigned long offset) { vertexID_offset = offset; };

    std::string GetArrayNameVertexID() { return arrayName_vertexID; };
    std::string GetArrayNameVertexDegree() { return arrayName_vertexDegree; };
    std::string GetArrayNameEdgeID() { return arrayName_edgeID; };
    std::string GetArrayNameEdgeLength() { return arrayName_edgeLength; };

private:
    bool addVertexID;
    bool addVertexDegree;

    bool addEdgeID;
    bool addEdgeLength;

    double m_vox_size_x;
    double m_vox_size_y;
    double m_vox_size_z;

    std::string arrayName_vertexID;
    std::string arrayName_vertexDegree;

    std::string arrayName_edgeID;
    std::string arrayName_edgeLength;

    unsigned long vertexID_offset;
};

#endif /* GRAPHANNOTATIONHELPER_H_ */
