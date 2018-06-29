///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GeometricalThresholdUndirectedGraphFilter.cpp                        //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-12-03                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "GeometricalThresholdUndirectedGraphFilter.h"
#include "../../tools/GraphAnnotationHelper.h"

#include <set>

#include "vtkAdjacentVertexIterator.h"
#include "vtkDataSetAttributes.h"
#include "vtkMath.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkVertexListIterator.h"

vtkStandardNewMacro(GeometricalThresholdUndirectedGraphFilter);


//-----------------------------------------------------------------------------
GeometricalThresholdUndirectedGraphFilter::GeometricalThresholdUndirectedGraphFilter() : vtkGraphAlgorithm()
{
	this->SetNumberOfInputPorts(1);
	this->AngleThreshold = 360;
	this->TestOutput = false;
	this->WithAngleAnnotation = false;
}


//-----------------------------------------------------------------------------
GeometricalThresholdUndirectedGraphFilter::~GeometricalThresholdUndirectedGraphFilter()
{
  // Do nothing.
}


int GeometricalThresholdUndirectedGraphFilter::FillInputPortInformation(int port, vtkInformation* info)
{
	if(port==0) {
		info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkGraph");
		return 1;
	}
	return 0;
}


//-----------------------------------------------------------------------------
void GeometricalThresholdUndirectedGraphFilter::PrintSelf(ostream &os, vtkIndent indent)
{
  // Base class print.
  vtkGraphAlgorithm::PrintSelf(os, indent);

  os << indent << "AngleThreshold: "  << this->AngleThreshold << endl;
}


//----------------------------------------------------------------------------
int GeometricalThresholdUndirectedGraphFilter::RequestDataObject(vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	if(!inInfo)
		return 0;

	vtkGraph *input = vtkGraph::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	if(input) {
		vtkInformation* info = outputVector->GetInformationObject(0);
		vtkGraph *output = vtkGraph::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));

		if(!output) {
			output = input->NewInstance();

			output->SetPipelineInformation(info);
			output->Delete();
		}
		return 1;
	}
	return 0;
}


//----------------------------------------------------------------------------
int GeometricalThresholdUndirectedGraphFilter::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	vtkGraph* input = vtkGraph::GetData(inputVector[0]);
	vtkGraph* output = vtkGraph::GetData(outputVector);

	std::map<vtkIdType, float> angles;

	if(input->GetNumberOfVertices()==0)
		return 1;							// If input graph is empty, return an empty graph

	if(AngleThreshold == 0) {
	    if(!output->CheckedShallowCopy(input)) {
	        vtkErrorMacro(<<"Invalid graph structure.");
	        return 0;
	    }
	    output->GetFieldData()->PassData(input->GetFieldData());
	    return 1;
	}


	vtkSmartPointer<vtkMutableUndirectedGraph> builder = vtkSmartPointer<vtkMutableUndirectedGraph>::New();

	vtkDataSetAttributes* vdIn = input->GetVertexData();
	vtkDataSetAttributes* edIn = input->GetEdgeData();
	vtkDataSetAttributes* vdOut = builder->GetVertexData();
	vtkDataSetAttributes* edOut = builder->GetEdgeData();
	vtkPoints* ptsIn = input->GetPoints();
	vtkPoints* ptsOut = builder->GetPoints();
	vdOut->CopyAllocate(vdIn);
	edOut->CopyAllocate(edIn);

	std::map<vtkIdType, vtkIdType> vertex_map;
	vtkIdType out_vertex = 0;
	std::set<vtkIdType>	discovered;
	std::multiset<vtkIdType> will_be_processed;

	//Find start vertex; (and detect pure circle graphs & handle them)
	vtkIdType start_vertex;
	bool real_start_vertex_found = false;

	vtkSmartPointer<vtkVertexListIterator> vertexIter = vtkSmartPointer<vtkVertexListIterator>::New();
	input->GetVertices(vertexIter);

	while(vertexIter->HasNext() && !real_start_vertex_found) {
		vtkIdType start_cand = vertexIter->Next();

		if(input->GetDegree(start_cand) != 2) {
			start_vertex = start_cand;
			real_start_vertex_found = true;
		}
	}
	if(!real_start_vertex_found) {
		input->GetVertices(vertexIter);

		if(vertexIter->HasNext())
			start_vertex = vertexIter->Next();
	}

	//Let's start the fun
	vtkSmartPointer<vtkAdjacentVertexIterator> neighIterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();

	will_be_processed.insert(start_vertex);
	discovered.insert(start_vertex);

	while(will_be_processed.size()>0) {
		if(TestOutput) std::cout << "will_be_processed.size() = " << will_be_processed.size() << std::endl;

		vtkIdType this_vertex = (*will_be_processed.begin());
		will_be_processed.erase(will_be_processed.begin());
		if(TestOutput) std::cout << "crawler: this_vertex = " << this_vertex << std::endl;

		if(vertex_map.count(this_vertex)==0) {
			builder->AddVertex();
			vdOut->CopyData(vdIn, this_vertex, out_vertex);
			ptsOut->InsertNextPoint(ptsIn->GetPoint(this_vertex));

			vertex_map[this_vertex] = out_vertex;
			out_vertex++;
		}

		input->GetAdjacentVertices(this_vertex, neighIterator);
		while(neighIterator->HasNext()) {
			vtkIdType direction_vertex = neighIterator->Next();
			if(TestOutput) std::cout << "crawler: direction_vertex candidates= " << direction_vertex << std::endl;

			if(discovered.count(direction_vertex) == 0) {
				if(TestOutput) std::cout << "crawler: direction_vertex = " << direction_vertex << std::endl;

				vtkIdType vertex1 = this_vertex;
				vtkIdType vertex2 = direction_vertex;
				vtkIdType vertex3 = vertex2;

				vtkIdType start_vertex = vertex1;
				vtkIdType end_vertex = vertex3;
				bool go_on = true;

				while(go_on) {
					if(input->GetDegree(vertex2) == 2) {
						if(discovered.count(vertex2)==0)
							discovered.insert(vertex2);

						vtkSmartPointer<vtkAdjacentVertexIterator> neighJterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
						input->GetAdjacentVertices(vertex2, neighJterator);
						vtkIdType v3a = neighJterator->Next();
						vtkIdType v3b = neighJterator->Next();

						if(v3a != vertex1)
							vertex3 = v3a;
						else
							vertex3 = v3b;

						double p1[3], p2[3], p3[3], vec1[3], vec2[3];
						ptsIn->GetPoint(vertex1, p1);
						ptsIn->GetPoint(vertex2, p2);
						ptsIn->GetPoint(vertex3, p3);
						vtkMath::Subtract(p2,p1,vec1);
						vtkMath::Subtract(p3,p2,vec2);

						double dot = vtkMath::Dot(vec1, vec2);
						double norms = vtkMath::Norm(vec1) * vtkMath::Norm(vec2);
						double dn = (double)((int)(dot/norms*1000))/1000.0;
						double angle = vtkMath::DegreesFromRadians(acos(dn));
//						if(TestOutput) std::cout << "Angle between " << vertex1 << "," << vertex2 << " and " << vertex2 << "," << vertex3 << " is " << angle << " (dot/norms=" << dn << ") versus angleThreshold of " << AngleThreshold << std::endl;

						if(WithAngleAnnotation)
							angles[vertex2] = angle;

						if(angle<=AngleThreshold) {
							end_vertex = vertex3;

							vertex1 = vertex2;
							vertex2 = vertex3;
						}
						else
							go_on = false;
					}
					else
						go_on = false;
				}

				if(vertex_map.count(end_vertex)==0) {
					will_be_processed.insert(end_vertex);
					if(discovered.count(end_vertex)==0)
						discovered.insert(end_vertex);
					builder->AddVertex();
					vdOut->CopyData(vdIn, end_vertex, out_vertex);
					ptsOut->InsertNextPoint(ptsIn->GetPoint(end_vertex));

					vertex_map[end_vertex] = out_vertex;
					out_vertex++;
				}

				vtkEdgeType f;
				f = builder->AddEdge(vertex_map[start_vertex], vertex_map[end_vertex]);
				vtkIdType edgeId = input->GetEdgeId(start_vertex, end_vertex);
				if(edgeId != -1) {
					edOut->CopyData(edIn, edgeId, f.Id);
					// Copy edge layout to the output.
					vtkIdType npts;
					double* pts;
					input->GetEdgePoints(edgeId, npts, pts);
					builder->SetEdgePoints(f.Id, npts, pts);
				}
			}
			else if(will_be_processed.count(direction_vertex)>0) {
				vtkEdgeType f;
				f = builder->AddEdge(vertex_map[this_vertex], vertex_map[direction_vertex]);
				vtkIdType edgeId = input->GetEdgeId(this_vertex, direction_vertex);
				if(edgeId != -1) {
					edOut->CopyData(edIn, edgeId, f.Id);
					// Copy edge layout to the output.
					vtkIdType npts;
					double* pts;
					input->GetEdgePoints(edgeId, npts, pts);
					builder->SetEdgePoints(f.Id, npts, pts);
				}
			}
		}
	}

	if(WithAngleAnnotation) {
		GraphAnnotationHelper anno;
		anno.AddCustomVertexAnnotation(input, "VertexAngles", angles, -1);
	}


	// Pass constructed graph to output.
	if(!output->CheckedShallowCopy(builder)) {
		vtkErrorMacro(<<"Invalid graph structure.");
		return 0;
	}
	output->GetFieldData()->PassData(input->GetFieldData());

	// Clean up
	output->Squeeze();

	return 1;
}

