///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphAnnotationHelper.cpp                                            //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-01-06                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "GraphAnnotationHelper.h"

#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkEdgeListIterator.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkUnsignedLongArray.h"
#include "vtkVertexListIterator.h"

#include <set>
#include <sstream>


GraphAnnotationHelper::GraphAnnotationHelper()
{
	addVertexID 	        = false;
	addVertexDegree	        = false;
	addEdgeID 		        = false;
	addEdgeLength	        = false;

	arrayName_vertexID 		        = "Pedigree IDs";
	arrayName_vertexDegree 	        = "degrees";
	arrayName_edgeID		        = "EdgeIDs";
	arrayName_edgeLength	        = "EdgeLengths";

	vertexID_offset = 0;
}


GraphAnnotationHelper::~GraphAnnotationHelper()
{
	//Auto-generated destructor stub
}


void GraphAnnotationHelper::EnableEdgeLengthAnnotation(double x_size, double y_size, double z_size)
{
    m_vox_size_x = x_size;
    m_vox_size_y = y_size;
    m_vox_size_z = z_size;

    addEdgeLength = true;
}


void GraphAnnotationHelper::AddPredefinedAnnotations(vtkSmartPointer<vtkGraph> graph)
{
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIds = vtkSmartPointer<vtkUnsignedLongArray>::New();               // Create an array for the vertex id labels
    pedigreeIds->SetName(arrayName_vertexID.c_str());

	vtkSmartPointer<vtkIntArray> vertexDegrees = vtkSmartPointer<vtkIntArray>::New();								// Create an array for the vertex degrees labels
	vertexDegrees->SetNumberOfComponents(1);
	vertexDegrees->SetName(arrayName_vertexDegree.c_str());

	vtkSmartPointer<vtkIntArray> edgeIDs = vtkSmartPointer<vtkIntArray>::New();										// Create an array for the edge id labels
	edgeIDs->SetNumberOfComponents(1);
	edgeIDs->SetName(arrayName_edgeID.c_str());

	vtkSmartPointer<vtkDoubleArray> edgeLengths = vtkSmartPointer<vtkDoubleArray>::New();							// Create an array for the edge length labels
	edgeLengths->SetNumberOfComponents(1);
	edgeLengths->SetName(arrayName_edgeLength.c_str());

	for(vtkIdType j=0; j<graph->GetNumberOfVertices(); j++) {
		if(addVertexID)			pedigreeIds->InsertNextValue(vertexID_offset+j);

		if(addVertexDegree)		vertexDegrees->InsertNextValue(graph->GetDegree(j));
	}

	for(vtkIdType j=0; j<graph->GetNumberOfEdges(); j++) {
		if(addEdgeID)			edgeIDs->InsertNextValue(j);
	}

	vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
	graph->GetEdges(it);
	while(it->HasNext()) {
		vtkEdgeType e = it->Next();

		if(addEdgeLength) {
			double a[3], b[3], vec[3];
			graph->GetPoint(e.Target, a);
			graph->GetPoint(e.Source, b);
			vtkMath::Subtract(b, a, vec);

			vec[0] *= m_vox_size_x;
			vec[1] *= m_vox_size_y;
			vec[2] *= m_vox_size_z;

			edgeLengths->InsertValue(e.Id, vtkMath::Norm(vec));
		}
	}

	if(addVertexID)			graph->GetVertexData()->SetPedigreeIds(pedigreeIds);;
	if(addVertexDegree)		graph->GetVertexData()->AddArray(vertexDegrees);
	if(addEdgeID)			graph->GetEdgeData()->AddArray(edgeIDs);
	if(addEdgeLength)		graph->GetEdgeData()->AddArray(edgeLengths);
}


void GraphAnnotationHelper::AddCustomVertexAnnotation(vtkSmartPointer<vtkGraph> graph, std::string arrayName, std::map<vtkIdType, int> &data, int default_value)
{
    vtkSmartPointer<vtkIntArray> vertexData = vtkSmartPointer<vtkIntArray>::New();                                  // Create an array for the vertex id labels
    vertexData->SetNumberOfComponents(1);
    vertexData->SetName(arrayName.c_str());

    vtkSmartPointer<vtkVertexListIterator> it = vtkSmartPointer<vtkVertexListIterator>::New();
    graph->GetVertices(it);

    if(data.size() != 0) {
        while(it->HasNext()) {
            vtkIdType e = it->Next();

            if(data.count(e)!=0)    vertexData->InsertValue(e, data[e]);
            else                    vertexData->InsertValue(e, default_value);
        }
    }
    else {
        while(it->HasNext()) {
            vtkIdType e = it->Next();
            vertexData->InsertValue(e, default_value);
        }
    }

    graph->GetVertexData()->AddArray(vertexData);
}


void GraphAnnotationHelper::AddCustomVertexAnnotation(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, std::string arrayNameSuffix, std::vector< std::multimap<vtkIdType,int> > &data, int maxArrays, int default_value)
{
    for(unsigned int i=0; i<graphs.size(); i++) {
        for(int j=1; j<=maxArrays; j++) {
            std::stringstream arrayName;
            arrayName << arrayNameSuffix << "_" << j;

            std::set<vtkIdType> uniqueIds;
            std::map<vtkIdType,int> dataArray;

            for(std::multimap<vtkIdType,int>::iterator it=data[i].begin(); it!=data[i].end(); ++it) {
                if(uniqueIds.count((*it).first)==0)
                    uniqueIds.insert((*it).first);
            }

            for(std::set<vtkIdType>::iterator it=uniqueIds.begin(); it!=uniqueIds.end(); ++it) {
                if(data[i].count(*it)<j)
                    dataArray.insert(std::pair<vtkIdType,int>(*it, default_value));
                else {
                    std::pair<std::multimap<vtkIdType,int>::iterator, std::multimap<vtkIdType,int>::iterator> ret;
                    ret = data[i].equal_range(*it);
                    int counter = 1;
                    for(std::multimap<vtkIdType,int>::iterator iter=ret.first; iter!=ret.second; ++iter) {
                        if(counter == j) {
                            dataArray.insert(std::pair<vtkIdType,int>(*it, iter->second));
                            break;
                        }
                        counter++;
                    }
                }
            }

            AddCustomVertexAnnotation(graphs[i], arrayName.str(), dataArray, default_value);
        }
    }
}


void GraphAnnotationHelper::AddCustomVertexAnnotation(vtkSmartPointer<vtkGraph> graph, std::string arrayName, std::map<vtkIdType, float> &data, float default_value)
{
	vtkSmartPointer<vtkFloatArray> vertexData = vtkSmartPointer<vtkFloatArray>::New();									// Create an array for the vertex id labels
	vertexData->SetNumberOfComponents(1);
	vertexData->SetName(arrayName.c_str());

	vtkSmartPointer<vtkVertexListIterator> it = vtkSmartPointer<vtkVertexListIterator>::New();
	graph->GetVertices(it);

	if(data.size() != 0) {
	    while(it->HasNext()) {
	        vtkIdType e = it->Next();

	        if(data.count(e)!=0)    vertexData->InsertValue(e, data[e]);
	        else                    vertexData->InsertValue(e, default_value);
	    }
	}
	else {
	    while(it->HasNext()) {
	        vtkIdType e = it->Next();
	        vertexData->InsertValue(e, default_value);
	    }
	}

	graph->GetVertexData()->AddArray(vertexData);
}


void GraphAnnotationHelper::AddCustomVertexAnnotation(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, std::string arrayNameSuffix, std::vector< std::multimap<vtkIdType,float> > &data, int maxArrays, float default_value)
{
    for(unsigned int i=0; i<graphs.size(); i++) {
        for(int j=1; j<=maxArrays; j++) {
            std::stringstream arrayName;
            arrayName << arrayNameSuffix << "_" << j;

            std::set<vtkIdType> uniqueIds;
            std::map<vtkIdType,float> dataArray;

            for(std::multimap<vtkIdType,float>::iterator it=data[i].begin(); it!=data[i].end(); ++it) {
                if(uniqueIds.count((*it).first)==0)
                    uniqueIds.insert((*it).first);
            }

            for(std::set<vtkIdType>::iterator it=uniqueIds.begin(); it!=uniqueIds.end(); ++it) {
                if(data[i].count(*it)<j)
                    dataArray.insert(std::pair<vtkIdType,float>(*it, default_value));
                else {
                    std::pair<std::multimap<vtkIdType,float>::iterator, std::multimap<vtkIdType,float>::iterator> ret;
                    ret = data[i].equal_range(*it);
                    int counter = 1;
                    for(std::multimap<vtkIdType,float>::iterator iter=ret.first; iter!=ret.second; ++iter) {
                        if(counter == j) {
                            dataArray.insert(std::pair<vtkIdType,float>(*it, iter->second));
                            break;
                        }
                        counter++;
                    }
                }
            }

            AddCustomVertexAnnotation(graphs[i], arrayName.str(), dataArray, default_value);
        }
    }
}


void GraphAnnotationHelper::AddCustomEdgeAnnotation(vtkSmartPointer<vtkGraph> graph, std::string arrayName, std::map<vtkIdType, int> &data, int default_value)
{
	vtkSmartPointer<vtkIntArray> edgeData = vtkSmartPointer<vtkIntArray>::New();									// Create an array for the vertex id labels
	edgeData->SetNumberOfComponents(1);
	edgeData->SetName(arrayName.c_str());

	vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
	graph->GetEdges(it);

	if(data.size() != 0) {
	    while(it->HasNext()) {
	        vtkEdgeType e = it->Next();

	        if(data.count(e.Id)!=0) edgeData->InsertValue(e.Id, data[e.Id]);
	        else                    edgeData->InsertValue(e.Id, default_value);
	    }
	}
	else {
	    while(it->HasNext()) {
	        vtkEdgeType e = it->Next();
	        edgeData->InsertValue(e.Id, default_value);
	    }
	}

	graph->GetEdgeData()->AddArray(edgeData);
}


void GraphAnnotationHelper::AddCustomEdgeAnnotation(vtkSmartPointer<vtkGraph> graph, std::string arrayName, std::map<vtkIdType, float> &data, float default_value)
{
    vtkSmartPointer<vtkFloatArray> edgeData = vtkSmartPointer<vtkFloatArray>::New();                                  // Create an array for the vertex id labels
    edgeData->SetNumberOfComponents(1);
    edgeData->SetName(arrayName.c_str());

    vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
    graph->GetEdges(it);

    if(data.size() != 0) {
        while(it->HasNext()) {
            vtkEdgeType e = it->Next();

            if(data.count(e.Id)!=0) edgeData->InsertValue(e.Id, data[e.Id]);
            else                    edgeData->InsertValue(e.Id, default_value);
        }
    }
    else {
        while(it->HasNext()) {
            vtkEdgeType e = it->Next();
            edgeData->InsertValue(e.Id, default_value);
        }
    }

    graph->GetEdgeData()->AddArray(edgeData);
}


void GraphAnnotationHelper::AddCustomEdgeAnnotation(vtkSmartPointer<vtkGraph> graph, std::string arrayName, std::map<vtkIdType, double> &data, double default_value)
{
	vtkSmartPointer<vtkDoubleArray> edgeData = vtkSmartPointer<vtkDoubleArray>::New();									// Create an array for the vertex id labels
	edgeData->SetNumberOfComponents(1);
	edgeData->SetName(arrayName.c_str());

	vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
	graph->GetEdges(it);

	if(data.size() != 0) {
	    while(it->HasNext()) {
	        vtkEdgeType e = it->Next();

	        if(data.count(e.Id)!=0) edgeData->InsertValue(e.Id, data[e.Id]);
	        else                    edgeData->InsertValue(e.Id, default_value);
	    }
	}
	else {
	    while(it->HasNext()) {
	        vtkEdgeType e = it->Next();
	        edgeData->InsertValue(e.Id, default_value);
	    }
	}

	graph->GetEdgeData()->AddArray(edgeData);
}
