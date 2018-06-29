///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  TopologicalPruneGraphFilter.cpp                                      //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-02-08                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "TopologicalPruneGraphFilter.h"

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
#include "vtkUnsignedLongArray.h"
#include "vtkVertexListIterator.h"

#include "../../tools/GraphAnnotationHelper.h"



vtkStandardNewMacro(TopologicalPruneGraphFilter);


TopologicalPruneGraphFilter::TopologicalPruneGraphFilter() : vtkGraphAlgorithm()
{
	this->SetNumberOfInputPorts(1);
	this->DeadEndLengthThreshold = 2.0;
	this->TestOutput = false;
}


/*!
  \brief Destructor.
*/
TopologicalPruneGraphFilter::~TopologicalPruneGraphFilter()
{
  // Do nothing.
}


int TopologicalPruneGraphFilter::FillInputPortInformation(int port, vtkInformation* info)
{
	if(port==0) {
		info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkGraph");
		return 1;
	}
	return 0;
}


/*!
  \brief Prints itself
*/
void TopologicalPruneGraphFilter::PrintSelf(ostream &os, vtkIndent indent)
{
  // Base class print.
  vtkGraphAlgorithm::PrintSelf(os, indent);

  os << indent << "DeadEndLengthThreshold: "  << this->DeadEndLengthThreshold << endl;
}



int TopologicalPruneGraphFilter::RequestDataObject(vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
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


void TopologicalPruneGraphFilter::IdentifyRemovableDeadEnds(vtkGraph *input, std::vector<FirstOrderBranch> &deadEndBranches, std::vector<FirstOrderBranch> &removableDeadEndBranches,
        std::vector<FirstOrderBranch> &nonRemovableDeadEndBranches, bool removeBorderDeadEnds)
{
    for(unsigned int i=0; i<deadEndBranches.size(); i++) {
        bool isVertexAtDatasetBorder = false;
        bool isVertexConnectedToVein = false;

        if(input->GetDegree(deadEndBranches[i].firstVert) == 1) {
            isVertexAtDatasetBorder |= GraphHandlingHelper::IsVertexAtDatasetBorder(input, deadEndBranches[i].firstVert, (2.0*mAssumedStandardRadius), mSpacing, mSize);
            int veinStatus = GraphHandlingHelper::GetVeinConnectednessStatus(input, deadEndBranches[i].firstVert, mpCVDistMap, mpPVDistMap,
                    m_maxCVeinConnectednessDistance, m_maxPVeinConnectednessDistance, true, true);
            if(veinStatus!=0)  isVertexConnectedToVein = true;
        }
        if(input->GetDegree(deadEndBranches[i].lastVert) == 1) {
            isVertexAtDatasetBorder |= GraphHandlingHelper::IsVertexAtDatasetBorder(input, deadEndBranches[i].lastVert, (2.0*mAssumedStandardRadius), mSpacing, mSize);
            int veinStatus = GraphHandlingHelper::GetVeinConnectednessStatus(input, deadEndBranches[i].lastVert, mpCVDistMap, mpPVDistMap,
                    m_maxCVeinConnectednessDistance, m_maxPVeinConnectednessDistance, true, true);
            if(veinStatus!=0)  isVertexConnectedToVein = true;
        }

        if(removeBorderDeadEnds) {                                                              //Standard case for first iteration: for all non-border & non-masked dead ends decision is length dependent
            if(isVertexAtDatasetBorder || isVertexConnectedToVein || deadEndBranches[i].length > DeadEndLengthThreshold)
                nonRemovableDeadEndBranches.push_back(deadEndBranches[i]);
            else
                removableDeadEndBranches.push_back(deadEndBranches[i]);
        }
        else {                                                                                  //In further iterations (PruneAll==true) remove all non-border & non-masked dead ends
            if(isVertexAtDatasetBorder || isVertexConnectedToVein)
                nonRemovableDeadEndBranches.push_back(deadEndBranches[i]);
            else
                removableDeadEndBranches.push_back(deadEndBranches[i]);
        }
    }
}


void TopologicalPruneGraphFilter::BuildBranchMaps(vtkGraph *input, bool removeBorderDeadEnds)
{
    mIsecBranches.clear();
    mRemovableDeadEndBranches.clear();
    mNonRemovableDeadEndBranches.clear();

    std::vector<FirstOrderBranch> deadEndBranches;

    GraphHandlingHelper::AssembleFirstOrderBranches(input, deadEndBranches, mIsecBranches);
    IdentifyRemovableDeadEnds(input, deadEndBranches, mRemovableDeadEndBranches, mNonRemovableDeadEndBranches, removeBorderDeadEnds);
}


/*!
  \brief Method that implements the thresholding of dead end branches
*/
int TopologicalPruneGraphFilter::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	vtkGraph* input = vtkGraph::GetData(inputVector[0]);
	vtkGraph* output = vtkGraph::GetData(outputVector);

	vtkGraph* inter = vtkGraph::GetData(outputVector);
	vtkSmartPointer<vtkMutableUndirectedGraph> builder;
	bool first_iter = true;

	if(input->GetNumberOfVertices()==0)
	    return 1;                           // If input graph is empty, return an empty graph

	if(!inter->CheckedShallowCopy(input)) {
	    vtkErrorMacro(<<"Invalid graph structure.");
	    return 0;
	}

	GraphAnnotationHelper anno;
	anno.EnableVertexIDAnnotation();
	anno.EnableEdgeLengthAnnotation(mSpacing[0], mSpacing[1], mSpacing[2]);
	anno.AddPredefinedAnnotations(inter);

	BuildBranchMaps(inter, first_iter);

	while(mRemovableDeadEndBranches.size()>0 || first_iter)
	{
	    builder = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	    first_iter = false;

	    vtkDataSetAttributes* vdIn = inter->GetVertexData();
	    vtkDataSetAttributes* edIn = inter->GetEdgeData();
	    vtkDataSetAttributes* vdOut = builder->GetVertexData();
	    vtkDataSetAttributes* edOut = builder->GetEdgeData();
	    vtkPoints* ptsIn = inter->GetPoints();
	    vtkPoints* ptsOut = builder->GetPoints();
	    vdOut->CopyAllocate(vdIn);
	    edOut->CopyAllocate(edIn);

	    std::map<vtkIdType, vtkIdType> vertex_map;
	    vtkIdType out_vertex = 0;

	    vtkSmartPointer<vtkAdjacentVertexIterator> neighIterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();

	    for(unsigned int i=0; i<mIsecBranches.size(); i++) {
	        FirstOrderBranch b = mIsecBranches[i];

	        vtkIdType w = -1;
	        for(unsigned int j=0; j<b.nodes.size(); j++) {
	            vtkIdType v = b.nodes[j];

	            if(inter->GetDegree(v) < 3) {
	                builder->AddVertex();
	                vdOut->CopyData(vdIn, v, out_vertex);
	                ptsOut->InsertNextPoint(ptsIn->GetPoint(v));

	                vertex_map[v] = out_vertex;
	                out_vertex++;

	                if(w!=-1) {
	                    vtkEdgeType f;
	                    f = builder->AddEdge(vertex_map[v], vertex_map[w]);
	                    vtkIdType edgeId = inter->GetEdgeId(v, w);

	                    if(edgeId != -1) {
	                        edOut->CopyData(edIn, edgeId, f.Id);
	                        // Copy edge layout to the output.
	                        vtkIdType npts;
	                        double* pts;
	                        inter->GetEdgePoints(edgeId, npts, pts);
	                        builder->SetEdgePoints(f.Id, npts, pts);
	                    }
	                }
	                w = v;
	            }
	            else {
	                if(vertex_map.count(v)==0) {
	                    builder->AddVertex();
	                    vdOut->CopyData(vdIn, v, out_vertex);
	                    ptsOut->InsertNextPoint(ptsIn->GetPoint(v));

	                    vertex_map[v] = out_vertex;
	                    out_vertex++;
	                }

	                inter->GetAdjacentVertices(v, neighIterator);
	                while(neighIterator->HasNext()) {
	                    vtkIdType d = neighIterator->Next();

	                    if(vertex_map.count(d)!=0 && builder->GetEdgeId(vertex_map[v], vertex_map[d])==-1) {
	                        vtkEdgeType f;
	                        f = builder->AddEdge(vertex_map[v], vertex_map[d]);
	                        vtkIdType edgeId = inter->GetEdgeId(v, d);

	                        edOut->CopyData(edIn, edgeId, f.Id);
	                        // Copy edge layout to the output.
	                        vtkIdType npts;
	                        double* pts;
	                        inter->GetEdgePoints(edgeId, npts, pts);
	                        builder->SetEdgePoints(f.Id, npts, pts);
	                    }
	                }
	                w = v;
	            }
	        }
	    }

	    for(unsigned int i=0; i<mNonRemovableDeadEndBranches.size(); i++) {
	        FirstOrderBranch b = mNonRemovableDeadEndBranches[i];

	        vtkIdType w = -1;
	        for(unsigned int j=0; j<b.nodes.size(); j++) {
	            vtkIdType v = b.nodes[j];

	            if(vertex_map.count(v)==0) {
	                builder->AddVertex();
	                vdOut->CopyData(vdIn, v, out_vertex);
	                ptsOut->InsertNextPoint(ptsIn->GetPoint(v));

	                vertex_map[v] = out_vertex;
	                out_vertex++;
	            }

	            inter->GetAdjacentVertices(v, neighIterator);
	            while(neighIterator->HasNext()) {
	                vtkIdType d = neighIterator->Next();

	                if(vertex_map.count(d)!=0 && builder->GetEdgeId(vertex_map[v], vertex_map[d])==-1) {
	                    vtkEdgeType f;
	                    f = builder->AddEdge(vertex_map[v], vertex_map[d]);
	                    vtkIdType edgeId = inter->GetEdgeId(v, d);

	                    edOut->CopyData(edIn, edgeId, f.Id);
	                    // Copy edge layout to the output.
	                    vtkIdType npts;
	                    double* pts;
	                    inter->GetEdgePoints(edgeId, npts, pts);
	                    builder->SetEdgePoints(f.Id, npts, pts);
	                }
	            }
	        }
	    }

        if(!inter->CheckedShallowCopy(builder)) {
            vtkErrorMacro(<<"Invalid graph structure.");
            return 0;
        }

        if(!PruneAll)                                                               //normal mode PruneAll = False = only one dead end removal iteration, otherwise iterative removal !!
            break;

        GraphAnnotationHelper anno;
        anno.EnableVertexIDAnnotation();
        anno.EnableEdgeLengthAnnotation(mSpacing[0], mSpacing[1], mSpacing[2]);
        anno.AddPredefinedAnnotations(inter);

        BuildBranchMaps(inter, first_iter);
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
