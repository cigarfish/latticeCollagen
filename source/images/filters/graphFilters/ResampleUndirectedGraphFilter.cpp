///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ResampleUndirectedGraphFilter.cpp                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-11-21                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "ResampleUndirectedGraphFilter.h"

#include <map>
#include <set>

#include "vtkAdjacentVertexIterator.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkVertexListIterator.h"


vtkStandardNewMacro(ResampleUndirectedGraphFilter);

/*!
  \brief Constructor.
  Sets the resampling factor to 2.
*/
ResampleUndirectedGraphFilter::ResampleUndirectedGraphFilter() : vtkGraphAlgorithm()
{
	this->SetNumberOfInputPorts(1);
	this->ResamplingFactor = 2;
	this->TestOutput = false;
}


/*!
  \brief Destructor.
*/
ResampleUndirectedGraphFilter::~ResampleUndirectedGraphFilter()
{
  // Do nothing.
}


int ResampleUndirectedGraphFilter::FillInputPortInformation(int port, vtkInformation* info)
{
	if(port==0) {
		info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkGraph");
		return 1;
	}
	return 0;
}


/*!
  \brief Prints itself (resampling factor).
*/
void ResampleUndirectedGraphFilter::PrintSelf(ostream &os, vtkIndent indent)
{
  // Base class print.
  vtkGraphAlgorithm::PrintSelf(os, indent);

  os << indent << "ResamplingFactor: "  << this->ResamplingFactor << endl;
}



int ResampleUndirectedGraphFilter::RequestDataObject(vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
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


void ResampleUndirectedGraphFilter::BuildBranchMaps(vtkGraph *input)
{
    branches.clear();

    std::vector<FirstOrderBranch> deadEndBranches;

    GraphHandlingHelper::AssembleFirstOrderBranches(input, deadEndBranches, branches);

    for(unsigned int i=0; i<deadEndBranches.size(); i++)
        branches.push_back(deadEndBranches[i]);
}


/*!
  \brief Method that implements the resampling.
*/
int ResampleUndirectedGraphFilter::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
    vtkGraph* input = vtkGraph::GetData(inputVector[0]);
    vtkGraph* output = vtkGraph::GetData(outputVector);

    if(input->GetNumberOfVertices()==0)
        return 1;                           // If input graph is empty, return an empty graph


    //TODO invoke edge length calculate if input->GetEdgeData()->HasArray("EdgeLengths") == 0

    BuildBranchMaps(input);

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

    vtkSmartPointer<vtkAdjacentVertexIterator> neighIterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
    vtkDoubleArray* edge_lengths = vtkDoubleArray::SafeDownCast( input->GetEdgeData()->GetArray("EdgeLengths") );

    for(unsigned int i=0; i<branches.size(); i++) {
        FirstOrderBranch b = branches[i];

        vtkIdType w = b.nodes[0];

        if(TestOutput) std::cout << "new branch starts with w = " << w << std::endl;

        if(vertex_map.count(w)==0) {
            builder->AddVertex();
            vdOut->CopyData(vdIn, w, out_vertex);
            ptsOut->InsertNextPoint(ptsIn->GetPoint(w));

            vertex_map[w] = out_vertex;
            out_vertex++;
        }

        for(unsigned int j=1; j<b.nodes.size(); j++) {
            vtkIdType v = b.nodes[j];
            vtkIdType u = w;

            int c=0;
            double length = edge_lengths->GetValue(input->GetEdgeId(w, v));
            if(TestOutput) std::cout << "length = " << length << std::endl;

            while(c<ResamplingFactor-1 && input->GetDegree(v) == 2) {
                if(TestOutput) std::cout << "skip v = " << v << std::endl;

                u = v;
                v = b.nodes[j+1];

                length += edge_lengths->GetValue(input->GetEdgeId(u, v));
                if(TestOutput) std::cout << "length = " << length << std::endl;

                if(length>MaxResamplingDist) {
                    length -= edge_lengths->GetValue(input->GetEdgeId(u, v));
                    v = u;
                    break;
                }
                else
                    j++;
                c++;
            }

            if(TestOutput) std::cout << "stop node v = " << v << std::endl;

            if(vertex_map.count(v)==0) {
                builder->AddVertex();
                vdOut->CopyData(vdIn, v, out_vertex);
                ptsOut->InsertNextPoint(ptsIn->GetPoint(v));

                vertex_map[v] = out_vertex;
                out_vertex++;
            }

            vtkEdgeType f;
            f = builder->AddEdge(vertex_map[v], vertex_map[w]);
            vtkIdType edgeId = input->GetEdgeId(v, w);

            if(TestOutput) std::cout << "connect w = " << w << " and v = " << v << std::endl;

            if(edgeId != -1) {
                edOut->CopyData(edIn, edgeId, f.Id);
                // Copy edge layout to the output.
                vtkIdType npts;
                double* pts;
                input->GetEdgePoints(edgeId, npts, pts);
                builder->SetEdgePoints(f.Id, npts, pts);
            }
            w = v;
            if(TestOutput) std::cout << "new start is w = " << w << std::endl;
        }
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
