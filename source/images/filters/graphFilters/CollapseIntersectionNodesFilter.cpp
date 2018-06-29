///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CollapseIntersectionNodesFilter.cpp                                  //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-01-31                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "CollapseIntersectionNodesFilter.h"

#include <map>
#include <set>

#include <vtkAdjacentVertexIterator.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkEdgeListIterator.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkObjectFactory.h>
#include <vtkOutEdgeIterator.h>
#include <vtkPoints.h>
#include <vtkUnsignedLongArray.h>
#include <vtkVertexListIterator.h>

#include "../../tools/GraphAnnotationHelper.h"


vtkStandardNewMacro(CollapseIntersectionNodesFilter);


CollapseIntersectionNodesFilter::CollapseIntersectionNodesFilter() : vtkGraphAlgorithm()
{
	this->SetNumberOfInputPorts(1);
	this->IsecDistanceThreshold = 5;
	this->TestOutput = false;
}


/*!
  \brief Destructor.
*/
CollapseIntersectionNodesFilter::~CollapseIntersectionNodesFilter()
{
  // Do nothing.
}


int CollapseIntersectionNodesFilter::FillInputPortInformation(int port, vtkInformation* info)
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
void CollapseIntersectionNodesFilter::PrintSelf(ostream &os, vtkIndent indent)
{
  // Base class print.
  vtkGraphAlgorithm::PrintSelf(os, indent);

  os << indent << "IsecDistanceThreshold: "  << this->IsecDistanceThreshold << endl;
}



int CollapseIntersectionNodesFilter::RequestDataObject(vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
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


double CollapseIntersectionNodesFilter::MeasureDistance(FScalarVoImageType::Pointer image, double *pos)
{
    FScalarVoImageType::IndexType idx;
    idx[0] = pos[0];
    idx[1] = pos[1];
    idx[2] = pos[2];

    return image->GetPixel(idx);
}


void CollapseIntersectionNodesFilter::IdentifyCollapsableIsecBranches(vtkGraph *input, std::vector<FirstOrderBranch> &isecBranches, std::vector<FirstOrderBranch> &collapseableBranches,
        std::vector<FirstOrderBranch> &unalteredBranches)
{
    for(unsigned int i=0; i<isecBranches.size(); i++) {
        if(isecBranches[i].length <= IsecDistanceThreshold && MeasureDistance(mpNetworkDistMap, isecBranches[i].middlePoint) > 0)
            collapseableBranches.push_back(isecBranches[i]);
        else
            unalteredBranches.push_back(isecBranches[i]);
    }
}


void CollapseIntersectionNodesFilter::BuildBranchMaps(vtkGraph *input)
{
    mIsecCollapseBranch.clear();
    mIsecUnalteredBranch.clear();
    mDeadEndBranch.clear();

    std::vector<FirstOrderBranch> isecBranches;

    GraphHandlingHelper::AssembleFirstOrderBranches(input, mDeadEndBranch, isecBranches);
    IdentifyCollapsableIsecBranches(input, isecBranches, mIsecCollapseBranch, mIsecUnalteredBranch);
}


/*!
  \brief Method that implements the collapsing of adjacent intersection nodes
*/
int CollapseIntersectionNodesFilter::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector, vtkInformationVector* outputVector)
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

    BuildBranchMaps(inter);

    while(mIsecCollapseBranch.size()>0 || first_iter)
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
        std::map<vtkIdType, vtkIdType> collapsedTo;
        vtkIdType out_vertex = 0;

        vtkSmartPointer<vtkAdjacentVertexIterator> neighIterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();

        for(unsigned int i=0; i<mIsecCollapseBranch.size(); i++) {
            FirstOrderBranch b = mIsecCollapseBranch[i];

            vtkIdType inId = b.firstVert;
            vtkIdType outId = b.lastVert;

            if(collapsedTo.count(inId)==0 && collapsedTo.count(outId)==0) {
                builder->AddVertex();
                double pos[3];
                pos[0] = b.middlePoint[0]; pos[1] = b.middlePoint[1]; pos[2] = b.middlePoint[2];
                ptsOut->InsertNextPoint(pos);

                vertex_map[inId] = out_vertex;
                vertex_map[outId] = out_vertex;
                collapsedTo[inId] = out_vertex;
                collapsedTo[outId] = out_vertex;
                out_vertex++;

                inter->GetAdjacentVertices(inId, neighIterator);
                while(neighIterator->HasNext()) {
                    vtkIdType d = neighIterator->Next();

                    if(vertex_map.count(d)!=0 && builder->GetEdgeId(vertex_map[inId], vertex_map[d])==-1 && d!=outId) {
                        vtkEdgeType f;
                        f = builder->AddEdge(vertex_map[inId], vertex_map[d]);
                        vtkIdType edgeId = inter->GetEdgeId(inId, d);

                        edOut->CopyData(edIn, edgeId, f.Id);
                        // Copy edge layout to the output.
                        vtkIdType npts;
                        double* pts;
                        inter->GetEdgePoints(edgeId, npts, pts);
                        builder->SetEdgePoints(f.Id, npts, pts);
                    }
                }

                inter->GetAdjacentVertices(outId, neighIterator);
                while(neighIterator->HasNext()) {
                    vtkIdType d = neighIterator->Next();

                    if(vertex_map.count(d)!=0 && builder->GetEdgeId(vertex_map[outId], vertex_map[d])==-1 && d!=inId) {
                        vtkEdgeType f;
                        f = builder->AddEdge(vertex_map[outId], vertex_map[d]);
                        vtkIdType edgeId = inter->GetEdgeId(outId, d);

                        edOut->CopyData(edIn, edgeId, f.Id);
                        // Copy edge layout to the output.
                        vtkIdType npts;
                        double* pts;
                        inter->GetEdgePoints(edgeId, npts, pts);
                        builder->SetEdgePoints(f.Id, npts, pts);
                    }
                }
            }
            else
                mIsecUnalteredBranch.push_back(mIsecCollapseBranch[i]);
        }

        for(unsigned int i=0; i<mDeadEndBranch.size(); i++) {
            FirstOrderBranch b = mDeadEndBranch[i];

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

                        if(vertex_map.count(d)!=0 && builder->GetEdgeId(vertex_map[v], vertex_map[d])==-1 && !(collapsedTo.count(v)!=0 && collapsedTo.count(d)!=0 && collapsedTo[v]==collapsedTo[d])) {
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

        for(unsigned int i=0; i<mIsecUnalteredBranch.size(); i++) {
            FirstOrderBranch b = mIsecUnalteredBranch[i];

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

                        if(vertex_map.count(d)!=0 && builder->GetEdgeId(vertex_map[v], vertex_map[d])==-1 && !(collapsedTo.count(v)!=0 && collapsedTo.count(d)!=0 && collapsedTo[v]==collapsedTo[d])) {
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

        if(!inter->CheckedShallowCopy(builder)) {
            vtkErrorMacro(<<"Invalid graph structure.");
            return 0;
        }

        GraphAnnotationHelper anno;
        anno.EnableVertexIDAnnotation();
        anno.EnableEdgeLengthAnnotation(mSpacing[0], mSpacing[1], mSpacing[2]);
        anno.AddPredefinedAnnotations(inter);

        BuildBranchMaps(inter);
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

