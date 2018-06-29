///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphHandlingHelper.cpp                                              //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-07-07                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "GraphHandlingHelper.h"

#include "vtkAdjacentVertexIterator.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkPoints.h"
#include "vtkUnsignedLongArray.h"
#include "vtkVertexListIterator.h"


GraphHandlingHelper::GraphHandlingHelper()
{
    // TODO Auto-generated constructor stub

}


GraphHandlingHelper::~GraphHandlingHelper()
{
    // TODO Auto-generated destructor stub
}


std::map<vtkIdType, int> GraphHandlingHelper::AssembleFirstOrderBranches(vtkSmartPointer<vtkGraph> input, std::vector<FirstOrderBranch> &deadEndBranches, std::vector<FirstOrderBranch> &isecBranches)
{
    int num_branches = 0;
    int branchState = 0;

    std::map<vtkIdType, int> edgeIdToBranchId;

    vtkSmartPointer<vtkVertexListIterator> vertexIter = vtkSmartPointer<vtkVertexListIterator>::New();
    input->GetVertices(vertexIter);

    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIds = vtkUnsignedLongArray::SafeDownCast( input->GetVertexData()->GetPedigreeIds() );
    vtkSmartPointer<vtkDoubleArray> edge_lengths = vtkDoubleArray::SafeDownCast( input->GetEdgeData()->GetArray("EdgeLengths") );

    while(vertexIter->HasNext()) {                                                                              //go through all vertices
        vtkIdType v = vertexIter->Next();
        unsigned long int vPID = pedigreeIds->GetValue(v);

        if(input->GetDegree(v) != 2) {                                                                          //found a non-regular-vertex?
            branchState = 0;

            for(int j=0; j<input->GetDegree(v); j++) {                                                          //then go over all its edges
                vtkOutEdgeType e = input->GetOutEdge(v, j);                                                     //take one of its edges

                if(edgeIdToBranchId.count(e.Id)==0) {                                                           //has the edge already a branch number?
                    FirstOrderBranch branch;
                    branch.nodes.push_back(v);
                    branch.nodesPIDs.push_back(vPID);
                    branch.firstVert = v;
                    branch.firstVertPID = vPID;

                    edgeIdToBranchId[e.Id] = num_branches;                                                      //if not give it to it

                    if(input->GetDegree(v) == 1)                                                                //if the journey starts at a dead end
                        branchState = 1;                                                                        //we have no 'intersection-branch'
                    else
                        branchState = 2;

                    vtkOutEdgeType a = e;                                                                       //save the actual edge in this variable

                    double branch_length = edge_lengths->GetValue(a.Id);                                        //and add its length to the branch length variable

                    while(input->GetDegree(a.Target)==2)    {                                                   //if the other node of the edge is a regular node the branch continues
                        branch.nodes.push_back(a.Target);
                        branch.nodesPIDs.push_back( pedigreeIds->GetValue(a.Target) );

                        for(int k=0; k<2; k++) {
                            vtkOutEdgeType ec = input->GetOutEdge(a.Target, k);                                 //look at both edges
                            if(ec.Id!=a.Id) {                                                                   //and chose the one which is not the actual edge
                                edgeIdToBranchId[ec.Id] = num_branches;                                         //give this edge its branch number
                                a = ec;                                                                         //and make this edge the actual edge

                                branch_length += edge_lengths->GetValue(a.Id);                                  //and add its length to the branch length variable

                                break;                                                                          //and leave the for loop
                            }
                        }
                    }                                                                                           //now we have reached a non regular node which make the branch end
                    if(input->GetDegree(a.Target) == 1)                                                         //if the journey ends at a dead end
                        branchState += 1;                                                                       //we have no 'intersection-branch'
                    else                                                                                        //otherwise
                        branchState += 2;

                    branch.nodes.push_back(a.Target);
                    branch.nodesPIDs.push_back( pedigreeIds->GetValue(a.Target) );
                    branch.lastVert = a.Target;
                    branch.lastVertPID = pedigreeIds->GetValue(a.Target);
                    branch.length = branch_length;

                    double outPos[3];
                    GetMiddlepointOfVertices(input, branch.firstVert, branch.lastVert, outPos[0], outPos[1], outPos[2]);
                    branch.middlePoint[0] = outPos[0];
                    branch.middlePoint[1] = outPos[1];
                    branch.middlePoint[2] = outPos[2];

                    if(branchState==4) {
                        branch.id = isecBranches.size();
                        isecBranches.push_back(branch);
                    }
                    else {
                        branch.id = deadEndBranches.size();
                        deadEndBranches.push_back(branch);
                    }

                    num_branches++;                                                                             //therefore add one to the branch counter
                }
            }
        }
    }
    return edgeIdToBranchId;
}


void GraphHandlingHelper::AssembleSecondOrderBranches(vtkSmartPointer<vtkGraph> input, std::vector<SecondOrderBranch> &secOrderBranches,
        std::vector<FirstOrderBranch> &intersectionBranches, std::vector<FirstOrderBranch> &innerDeadEndBranches, std::vector<FirstOrderBranch> &maskedDeadEndBranches,
        std::map<vtkIdType, int> &vertexIdToSecOState, std::map<vtkIdType, int> &edgeIdToBranchId)
{
    std::multimap<vtkIdType, FirstOrderBranch*> vertIdToIsecBranches;               //build lookup maps
    std::multimap<vtkIdType, FirstOrderBranch*> vertIdToDeadEndBranches;
    std::multimap<vtkIdType, FirstOrderBranch*> vertIdToRestBranches;

    for(unsigned int i=0; i<intersectionBranches.size(); i++)
        for(unsigned int j=0; j<intersectionBranches[i].nodes.size(); j++)
            vertIdToIsecBranches.insert( std::pair<vtkIdType, FirstOrderBranch*>(intersectionBranches[i].nodes[j], &(intersectionBranches[i])) );

    for(unsigned int i=0; i<innerDeadEndBranches.size(); i++)
        for(unsigned int j=0; j<innerDeadEndBranches[i].nodes.size(); j++)
            vertIdToDeadEndBranches.insert( std::pair<vtkIdType, FirstOrderBranch*>(innerDeadEndBranches[i].nodes[j], &(innerDeadEndBranches[i])) );

    for(unsigned int i=0; i<maskedDeadEndBranches.size(); i++)
        for(unsigned int j=0; j<maskedDeadEndBranches[i].nodes.size(); j++)
            vertIdToRestBranches.insert( std::pair<vtkIdType, FirstOrderBranch*>(maskedDeadEndBranches[i].nodes[j], &(maskedDeadEndBranches[i])) );
                                                                                    //identify second order vertices (= vertices which branch in >=3 intersection branches)
    std::set<vtkIdType> secondOrderVertices;
    for(unsigned int i=0; i<input->GetNumberOfVertices(); i++) {
        if(vertIdToIsecBranches.count(i)>=3) {                                      //sec order vertex -> state 2
            vertexIdToSecOState.insert( std::pair<vtkIdType,int>(i, 2) );
            secondOrderVertices.insert(i);                                          //and add it id to second order vertex set
        }
        else if(vertIdToIsecBranches.count(i)==2) {                                 //first order intersection vertex
            if(vertIdToRestBranches.count(i)>0)                                     //is a dead end branch touching the dataset border attached -> state -1
                vertexIdToSecOState.insert( std::pair<vtkIdType,int>(i, -1) );
            else if(vertIdToDeadEndBranches.count(i)>0)                             //is a dead end branch not touching the dataset border attached -> state 0
                vertexIdToSecOState.insert( std::pair<vtkIdType,int>(i, 0) );
        }
        else if(vertIdToIsecBranches.count(i)==1)                                   //vertex somewhere on intersection branch -> state 1
            vertexIdToSecOState.insert( std::pair<vtkIdType,int>(i, 1) );
        else if(vertIdToIsecBranches.count(i)==0) {                                 //vertex is not on an intersection branch at all
            if(vertIdToRestBranches.count(i)>0)                                     //is it on a dead end branch touching the dataset border -> state -1
                vertexIdToSecOState.insert( std::pair<vtkIdType,int>(i, -1) );
            else if(vertIdToDeadEndBranches.count(i)>0)                             //is it on a dead end branch not touching the dataset border -> state 0
                vertexIdToSecOState.insert( std::pair<vtkIdType,int>(i, 0) );
            else
                vertexIdToSecOState.insert( std::pair<vtkIdType,int>(i, -10) );     //should not happen
        }
    }

    for(std::set<vtkIdType>::iterator itStart=secondOrderVertices.begin(); itStart!=secondOrderVertices.end(); ++itStart) {             //iterate over the second order vertices
        std::pair< std::multimap<vtkIdType, FirstOrderBranch*>::iterator, std::multimap<vtkIdType, FirstOrderBranch*>::iterator > ret;
        ret = vertIdToIsecBranches.equal_range(*itStart);

        for(std::multimap<vtkIdType, FirstOrderBranch*>::iterator it=ret.first; it!=ret.second; ++it) {                                 //iterate over all first order branches that start at a second order vertex
            FirstOrderBranch *itBranch = it->second;                                                                                    //get the first order branch

            if(edgeIdToBranchId.count( input->GetEdgeId(itBranch->nodes[0], itBranch->nodes[1]) ) == 0) {                               //check if this branch is not already part of a second order branch
                SecondOrderBranch b;                                                                                                    //start second order branch and initialize it
                b.firstOrderBranches.push_back(itBranch->id);
                b.length = itBranch->length;
                if(secondOrderVertices.count(itBranch->firstVert)!=0)
                    b.firstVert = itBranch->firstVert;
                else
                    b.firstVert = itBranch->lastVert;

                int dir = 0;                                                                                                            //check whether the second order vertex is first or last vertex in the first order branch
                if(itBranch->nodes[0] != it->first)                                                                                     //based on that determine direction (technically speaking)
                    dir = 1;

                vtkIdType nextStop;                                                                                                     //determine next first order intersection vertex
                if(dir==0)
                    nextStop = itBranch->nodes[itBranch->nodes.size()-1];
                else
                    nextStop = itBranch->nodes[0];

                bool foundBorderDeadEnd = false;                                                                                        //variable to track if we encountered dead end in contact with dataset border
                bool endAtDeadEndBranchingPoint = false;                                                                                //variable to track if intersection vertex branches in dead ends ONLY

                while( secondOrderVertices.count(nextStop)==0 ) {                                                                       //follow first order intersection branches until arrival at next second order vertex
                    if(vertIdToIsecBranches.count(nextStop)==2) {                                                                       //check whether dead end in contact with dataset border is attached to this vertex
                        if(vertIdToRestBranches.count(nextStop)>0)
                            foundBorderDeadEnd = true;
                    }
                    else {                                                                                                              //in case we haven't stopped at a second order vertex, and intersection vertex has not exactly to attached intersection branches
                        endAtDeadEndBranchingPoint = true;                                                                              //we ended up at a vertex that has one intersection branch attached plus some dead ends only
                        break;                                                                                                          //-> its a blind alley, so quit loop
                    }
                                                                                                                                        //iterate over the first order intersection branches beginning here
                    for(std::multimap<vtkIdType, FirstOrderBranch*>::iterator itIsecs = vertIdToIsecBranches.equal_range(nextStop).first; itIsecs != vertIdToIsecBranches.equal_range(nextStop).second; ++itIsecs) {
                        if(itIsecs->second->id != itBranch->id) {
                            itBranch = itIsecs->second;                                                                                 //and catch the one we haven't visited already
                            break;
                        }
                    }
                    b.firstOrderBranches.push_back(itBranch->id);                                                                       //update our second order branch
                    b.length += itBranch->length;

                    dir = 0;
                    if(itBranch->nodes[0] != nextStop)                                                                                  //based on that determine direction (technically speaking)
                        dir = 1;
                    if(dir==0)                                                                                                          //and retrieve the next stop vertex, at other side of branch
                        nextStop = itBranch->nodes[itBranch->nodes.size()-1];
                    else
                        nextStop = itBranch->nodes[0];
                }

                if(!foundBorderDeadEnd) {                                                                                               //after arrival at next second order vertex check whether border dead end was on the way
                    if(!endAtDeadEndBranchingPoint) {
                        for(unsigned int i=0; i<b.firstOrderBranches.size(); i++)                                                       //if not, iterate over all composing first order intersection branches
                            for(unsigned int j=1; j<intersectionBranches[b.firstOrderBranches[i]].nodes.size(); j++)                 //and their composing vertices, to retrieve sec order branch composing edges, and attach second order branch id to them
                                edgeIdToBranchId.insert( std::pair<vtkIdType,int>(input->GetEdgeId(intersectionBranches[b.firstOrderBranches[i]].nodes[j-1], intersectionBranches[b.firstOrderBranches[i]].nodes[j]), secOrderBranches.size()) );

                        if(b.firstOrderBranches.size()!=1) {
                            if(secondOrderVertices.count(itBranch->firstVert)!=0)
                                b.lastVert = itBranch->firstVert;
                            else
                                b.lastVert = itBranch->lastVert;
                        }
                        else
                            b.lastVert = itBranch->lastVert;

                        b.id = secOrderBranches.size();
                        secOrderBranches.push_back(b);                                                                                  //finally, remember second order branch
                    }
                    else {
                        for(unsigned int i=0; i<b.firstOrderBranches.size(); i++)                                                       //so attach '-1' to edges, composing this formal, but rejected second order branch
                            for(unsigned int j=1; j<intersectionBranches[b.firstOrderBranches[i]].nodes.size(); j++)
                                edgeIdToBranchId.insert( std::pair<vtkIdType,int>(input->GetEdgeId(intersectionBranches[b.firstOrderBranches[i]].nodes[j-1], intersectionBranches[b.firstOrderBranches[i]].nodes[j]), -1) );
                    }
                }
                else {                                                                                                                  //if dataset border touching dead end was discovered, we can not tell for sure that this is really a second order branch
                    for(unsigned int i=0; i<b.firstOrderBranches.size(); i++)                                                           //so attach '-2' to edges, composing this formal, but rejected second order branch
                        for(unsigned int j=1; j<intersectionBranches[b.firstOrderBranches[i]].nodes.size(); j++)
                            edgeIdToBranchId.insert( std::pair<vtkIdType,int>(input->GetEdgeId(intersectionBranches[b.firstOrderBranches[i]].nodes[j-1], intersectionBranches[b.firstOrderBranches[i]].nodes[j]), -2) );
                }
            }
        }
    }
}


void GraphHandlingHelper::GetMiddlepointOfVertices(vtkSmartPointer<vtkGraph> input, vtkIdType first, vtkIdType last, double &x, double &y, double &z)
{
    double inPos1[3], inPos2[3];
    input->GetPoint( first, inPos1);
    input->GetPoint( last, inPos2);

    x = (inPos1[0] + inPos2[0]) / 2.;
    y = (inPos1[1] + inPos2[1]) / 2.;
    z = (inPos1[2] + inPos2[2]) / 2.;

    if(x>0) x = floor(x);
    else    x = ceil(x);

    if(y>0) y = floor(y);
    else    y = ceil(y);

    if(z>0) z = floor(z);
    else    z = ceil(z);
}


bool GraphHandlingHelper::IsVertexAtDatasetBorder(vtkSmartPointer<vtkGraph> input, vtkIdType vertexID, double minDistanceInUnit, FScalarVoImageType::SpacingType spacing, FScalarVoImageType::SizeType datasetSize)
{
    double pos[3];
    input->GetPoint(vertexID, pos);

    CScalarVoImageType::IndexType idx;
    idx[0] = pos[0];
    idx[1] = pos[1];
    idx[2] = pos[2];

    int pixelToDataSetBorder[3];

    pixelToDataSetBorder[0] = ceil(minDistanceInUnit / spacing[0]);
    pixelToDataSetBorder[1] = ceil(minDistanceInUnit / spacing[1]);
    pixelToDataSetBorder[2] = ceil(minDistanceInUnit / spacing[2]);

    for(int i=0; i<3; i++)
        if(idx[i] <= pixelToDataSetBorder[i] || idx[i] >= (datasetSize[i]-1-pixelToDataSetBorder[i]))
            return true;

    return false;
}


bool GraphHandlingHelper::IsVertexMasked(vtkSmartPointer<vtkGraph> input, vtkIdType vertexID, CScalarVoImageType::Pointer maskImage, CScalarPixelType maskingValue)
{
    double pos[3];
    input->GetPoint(vertexID, pos);

    CScalarVoImageType::IndexType idx;
    idx[0] = pos[0];
    idx[1] = pos[1];
    idx[2] = pos[2];

    if(maskImage->GetPixel(idx) == maskingValue)
        return true;

    return false;
}


int GraphHandlingHelper::GetVeinConnectednessStatus(vtkSmartPointer<vtkGraph> input, vtkIdType vertexID, FScalarVoImageType::Pointer cvDistMap, FScalarVoImageType::Pointer pvDistMap,
        double maxCVDist, double maxPVDist, bool hasCV, bool hasPV)
{
    double pos[3];
    input->GetPoint(vertexID, pos);

    CScalarVoImageType::IndexType idx;
    idx[0] = pos[0];
    idx[1] = pos[1];
    idx[2] = pos[2];

    if(hasCV && cvDistMap->GetPixel(idx) < maxCVDist)
        return 1;
    else if(hasPV && pvDistMap->GetPixel(idx) < maxPVDist)
        return 2;
    else
        return 0;
}

//Computes for each branch connecting to a node the minimal angle considering all other connected branches
void GraphHandlingHelper::ComputeMinAngles(vtkSmartPointer<vtkGraph> input, vtkIdType vertexA, float lengthOfBranchVector, std::vector<float> &angles)
{
    vtkPoints* ptsIn = input->GetPoints();
    vtkSmartPointer<vtkDoubleArray> edge_lengths = vtkDoubleArray::SafeDownCast( input->GetEdgeData()->GetArray("EdgeLengths") );

    std::vector<vtkIdType> endVertices;
    std::map<vtkIdType,std::vector<vtkIdType> > branchPairs;

    for(int j=0; j<input->GetDegree(vertexA); j++) {
        vtkOutEdgeType e = input->GetOutEdge(vertexA, j);
        float length = edge_lengths->GetValue(e.Id);

        while(lengthOfBranchVector > length && input->GetDegree(e.Target)==2) {
            for(int k=0; k<2; k++) {
                vtkOutEdgeType ec = input->GetOutEdge(e.Target, k);
                if(ec.Id!=e.Id) {
                    e = ec;
                    length += edge_lengths->GetValue(e.Id);
                    break;
                }
            }
        }
        endVertices.push_back(e.Target);
    }

    for(unsigned int i=0; i<endVertices.size(); i++) {
        std::vector<vtkIdType> ends;
        for(unsigned int j=0; j<endVertices.size(); j++)
            if(endVertices[i]!=endVertices[j])
                ends.push_back(endVertices[j]);
        branchPairs.insert( std::pair<vtkIdType,std::vector<vtkIdType> >(endVertices[i],ends) );
    }

    for(std::map<vtkIdType,std::vector<vtkIdType> >::iterator it=branchPairs.begin(); it!=branchPairs.end(); ++it) {
        vtkIdType vertexB = (*it).first;
        float minAngle = 360;
        for(unsigned int i=0; i<(*it).second.size(); i++) {
            vtkIdType vertexC = (*it).second[i];

            double pA[3], pB[3], pC[3], vec1[3], vec2[3];
            ptsIn->GetPoint(vertexA, pA);
            ptsIn->GetPoint(vertexB, pB);
            ptsIn->GetPoint(vertexC, pC);
            vtkMath::Subtract(pB,pA,vec1);
            vtkMath::Subtract(pC,pA,vec2);

            float dot = vtkMath::Dot(vec1, vec2);
            float norms = vtkMath::Norm(vec1) * vtkMath::Norm(vec2);
            float dn = (double)((int)(dot/norms*1000))/1000.0;
            float angle = vtkMath::DegreesFromRadians(acos(dn));

            if(minAngle > angle)
                minAngle = angle;
        }
        angles.push_back(minAngle);
    }
}
