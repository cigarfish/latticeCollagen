///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelImageToGraphFilter.tpp                                          //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-08-07                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "LabelImageToGraphFilter.h"

#include <map>
#include <set>

#include "itkConstNeighborhoodIterator.h"
#include "vtkDataSetAttributes.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkPoints.h"
#include "vtkUnsignedLongArray.h"



template< unsigned int VImageDimension > LabelImageToGraphFilter< VImageDimension >::LabelImageToGraphFilter()
{
    // TODO Auto-generated constructor stub

}


template< unsigned int VImageDimension > LabelImageToGraphFilter< VImageDimension >::~LabelImageToGraphFilter()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > vtkSmartPointer<vtkUndirectedGraph> LabelImageToGraphFilter< VImageDimension >::LabelImageToGraph(itk::SmartPointer< itk::Image<unsigned long, VImageDimension> > image, bool ignoreBackgroundLabel)
{
    typedef typename itk::Image<unsigned long, VImageDimension>::IndexType  IndexType;
    typedef typename itk::Image<unsigned long, VImageDimension>::OffsetType OffsetType;
    typedef typename itk::Image<unsigned long, VImageDimension>::SizeType   SizeType;

    vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIds = vtkSmartPointer<vtkUnsignedLongArray>::New();
    pedigreeIds->SetName("Pedigree IDs");

    // Define offsets
    int numOffInDimMode = 4;
    if(ImageDimension == 3) numOffInDimMode = 6;
    OffsetType *off;
    if(ImageDimension == 2) {
        off = new OffsetType[4];        // face connected
        off[0][0] = -1; off[0][1] =  0; // 0 - NM
        off[1][0] =  1, off[1][1] =  0; // 1 - SM
        off[2][0] =  0; off[2][1] =  1; // 2 - MS
        off[3][0] =  0, off[3][1] = -1; // 3 - MN
    }
    else if(ImageDimension == 3) {
        off = new OffsetType[6];            // face connected
        off[0][0] = -1; off[0][1] =  0; off[0][2] =  0;    // 0 - NMM
        off[1][0] =  1; off[1][1] =  0; off[1][2] =  0;    // 1 - SMM
        off[2][0] =  0; off[2][1] =  1; off[2][2] =  0;    // 2 - MSM
        off[3][0] =  0; off[3][1] = -1; off[3][2] =  0;    // 3 - MNM
        off[4][0] =  0; off[4][1] =  0; off[4][2] = -1;    // 4 - MMN
        off[5][0] =  0; off[5][1] =  0; off[5][2] =  1;    // 5 - MMS
    }

    std::map<unsigned long, vtkIdType>          labelIdToVertexId;
    std::multimap<unsigned long, unsigned long> labelNeighborhood;
    std::map<vtkIdType, unsigned long>          vertexIdToLabelId;
    std::map<vtkIdType, IndexType>              vertexIdToCenterOfMass;
    std::map<vtkIdType, unsigned int>           vertexIdToNumberPixel;

    SizeType radius;
    radius.Fill(1);

    itk::ConstNeighborhoodIterator< itk::Image<unsigned long, VImageDimension> > it(radius, image, image->GetLargestPossibleRegion());
    it.GoToBegin();

    while(!it.IsAtEnd())
    {
        if(ignoreBackgroundLabel && it.GetCenterPixel()==0)
            ;
        else {
            if(labelIdToVertexId.count(it.GetCenterPixel()) == 0) {
                std::cout << "LabelImageToGraphFilter: add graph label " << it.GetCenterPixel() << std::endl;
                vtkIdType v = graph->AddVertex();
                pedigreeIds->InsertValue(v, it.GetCenterPixel());

                labelIdToVertexId.insert(std::pair<unsigned long, vtkIdType>(it.GetCenterPixel(), v));
                vertexIdToLabelId.insert(std::pair<vtkIdType, unsigned long>(v, it.GetCenterPixel()));

                IndexType idx = it.GetIndex();

                vertexIdToCenterOfMass.insert(std::pair<vtkIdType, IndexType>(v, idx));
                vertexIdToNumberPixel.insert(std::pair<vtkIdType, unsigned int>(v, 1));
            }
            else {
                vertexIdToCenterOfMass.at(labelIdToVertexId.at(it.GetCenterPixel()))[0] += it.GetIndex()[0];
                vertexIdToCenterOfMass.at(labelIdToVertexId.at(it.GetCenterPixel()))[1] += it.GetIndex()[1];
                if(ImageDimension==3)
                    vertexIdToCenterOfMass.at(labelIdToVertexId.at(it.GetCenterPixel()))[2] += it.GetIndex()[2];
                vertexIdToNumberPixel.at(labelIdToVertexId.at(it.GetCenterPixel())) += 1;
            }

            for(unsigned int i=0; i<numOffInDimMode; i++) {
                if(it.GetCenterPixel() != it.GetPixel(off[i]) && it.GetPixel(off[i])!=0) {
                    bool discovered = false;

                    unsigned long c = it.GetCenterPixel();

                    std::multimap<unsigned long, unsigned long>::iterator iter;
                    for(iter = labelNeighborhood.equal_range(c).first; iter!=labelNeighborhood.equal_range(c).second; ++iter)
                        if(iter->second == it.GetPixel(off[i]))
                            discovered = true;

                    if(!discovered) {
                        c = it.GetPixel(off[i]);

                        for(iter = labelNeighborhood.equal_range(c).first; iter!=labelNeighborhood.equal_range(c).second; ++iter)
                            if(iter->second == it.GetCenterPixel())
                                discovered = true;
                    }

                    if(!discovered)
                        labelNeighborhood.insert(std::pair<unsigned long, unsigned long>(it.GetCenterPixel(), it.GetPixel(off[i])));
                }
            }
        }

        ++it;
    }
    for(std::multimap<unsigned long, unsigned long>::iterator iter = labelNeighborhood.begin(); iter != labelNeighborhood.end(); ++iter) {
        vtkVariant id1 = static_cast<unsigned long>(iter->first);
        vtkVariant id2 = static_cast<unsigned long>(iter->second);

        graph->AddEdge(labelIdToVertexId.at(iter->first), labelIdToVertexId.at(iter->second));
    }

    std::cout << "LabelImageToGraphFilter: number vertices in label map graph = " << graph->GetNumberOfVertices() << std::endl;

    typedef typename std::map<vtkIdType, IndexType>::iterator MapIteratorType;
    for(MapIteratorType iter = vertexIdToCenterOfMass.begin(); iter != vertexIdToCenterOfMass.end(); ++iter) {
        iter->second[0] = iter->second[0] / vertexIdToNumberPixel.at(iter->first);
        iter->second[1] = iter->second[1] / vertexIdToNumberPixel.at(iter->first);
        if(ImageDimension==3)
            iter->second[2] = iter->second[2] / vertexIdToNumberPixel.at(iter->first);

        if(ImageDimension==2)
            points->InsertPoint(iter->first, iter->second[0], iter->second[1], 0);
        else if(ImageDimension==3)
            points->InsertPoint(iter->first, iter->second[0], iter->second[1], iter->second[2]);
    }
    std::cout << "LabelImageToGraphFilter: number edges in label map graph = " << graph->GetNumberOfEdges() << std::endl;

    graph->GetVertexData()->SetPedigreeIds(pedigreeIds);
    graph->SetPoints(points);

    //TODO: handle label 0 (=background) should be part of the graph

    return graph;
}
