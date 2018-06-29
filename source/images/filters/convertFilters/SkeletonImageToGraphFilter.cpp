///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SkeletonImageToGraphFilter.cpp                                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-11-09                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "SkeletonImageToGraphFilter.h"

#include <map>
#include <set>

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMap.h"
#include "itkShapeLabelObject.h"

#include "vtkPoints.h"



typedef itk::LabelMap<itk::ShapeLabelObject<itk::SizeValueType, 2> > LabelMap2DType;
typedef itk::LabelMap<itk::ShapeLabelObject<itk::SizeValueType, 3> > LabelMap3DType;

typedef LabelMap2DType::LabelObjectType::IndexType				IndexType2D;
typedef LabelMap2DType::LabelObjectType::IndexType::OffsetType OffsetType2D;

typedef LabelMap3DType::LabelObjectType::IndexType				IndexType3D;
typedef LabelMap3DType::LabelObjectType::IndexType::OffsetType OffsetType3D;


//--Helper class: defines order on IndexType; needed for map<IndexType>
class IndexCompare2D
{
   public:
      bool operator()(IndexType2D x, IndexType2D y) const
      {
    	  for(unsigned int i=0; i<x.GetIndexDimension(); i++)
    	  {
    		  if(x.GetIndex()[i] > y.GetIndex()[i])
    			  return true;

    		  else if(x.GetIndex()[i] < y.GetIndex()[i])
    			  return false;
    	  }
    	  return false;
      }
};

//--Helper class: defines order on IndexType; needed for map<IndexType>
class IndexCompare3D
{
   public:
      bool operator()(IndexType3D x, IndexType3D y) const
      {
    	  for(unsigned int i=0; i<x.GetIndexDimension(); i++)
    	  {
    		  if(x.GetIndex()[i] > y.GetIndex()[i])
    			  return true;

    		  else if(x.GetIndex()[i] < y.GetIndex()[i])
    			  return false;
    	  }
    	  return false;
      }
};


/*!
  \brief Method for converting skeletons in 2D images into graphs.
  Conversion assumes skeletons as 8-connected structures. The skeleton image is assumed to be a binary image where pixel value 1 indicate a skeleton pixel.
  \n Each pixel is converted into a graph vertex, to avoid zero-loops at crossings the pixel with the highest order induces the crossing vertex (face connectedness > edge connectedness > vertex connectedness).
  \n Skeletons that are not connected result in different graphs.

  \param image A 2-dimensional itk::Image which contains one or more skeletons.
  \returns A vector of vtkMutableUndirectedGraph objects that represents the skeletons in the input image.
*/
std::vector< vtkSmartPointer<vtkUndirectedGraph> > SkeletonImageToGraphFilter::skeletonToGraph2D(const itk::SmartPointer< itk::Image<bool,2> > image)
{
	typedef itk::BinaryImageToShapeLabelMapFilter< itk::Image<bool,2> > ImageToShapeLabelMapFilterType;

	ImageToShapeLabelMapFilterType::Pointer skeletonLabelMapFilter = ImageToShapeLabelMapFilterType::New();
	skeletonLabelMapFilter->ReleaseDataFlagOn();
	skeletonLabelMapFilter->SetFullyConnected(true);
	skeletonLabelMapFilter->SetInput(image);
	skeletonLabelMapFilter->Update();
	LabelMap2DType * skeletonMap = dynamic_cast<LabelMap2DType *>(skeletonLabelMapFilter->GetOutput());
	//itk::SmartPointer<LabelMap2DType> skeletonMap = skeletonLabelMapFilter->GetOutput();


	std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs;

	// Define offsets
	OffsetType2D of[8] = {
			// edge connected
			{{-1, 0}},	// 0 - north
			{{ 1, 0}},	// 1 - south
			{{ 0, 1}},	// 2 - east
			{{ 0,-1}},	// 3 - west
			// corner connected
			{{-1, 1}},	// 4 - northeast
			{{-1,-1}},  // 5 - northwest
			{{ 1, 1}},	// 6 - southeast
			{{ 1,-1}}	// 7 - southwest
	};

	bool con[8];
	bool conClose[8];
	int con_edge, con_vertex;


	for(unsigned int i=0; i<skeletonMap->GetNumberOfLabelObjects(); ++i) {
		vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

		LabelMap2DType::LabelObjectType* label_obj = skeletonMap->GetNthLabelObject(i);
		std::set<IndexType2D, IndexCompare2D>	not_discovered;
		std::set<IndexType2D, IndexCompare2D> will_be_processed;
		std::map<IndexType2D, vtkIdType, IndexCompare2D> indexVertexLink;
#if (ITK_VERSION_MAJOR < 4)
		for(unsigned int j=0; j<label_obj->GetSize(); j++)
#else
		for(unsigned int j=0; j<label_obj->GetNumberOfPixels(); j++)
#endif
			not_discovered.insert(label_obj->GetIndex(j));


		IndexType2D c = label_obj->GetIndex(0);
		vtkIdType v = graph->AddVertex();
		points->InsertNextPoint(c[0], c[1], 0.0);

		indexVertexLink.insert( std::pair<IndexType2D,vtkIdType>(c,v) );
		not_discovered.erase(c);

		con_edge = 0;
		for(int j=0; j<4; j++)													//edge connected have highest priority
		{
			IndexType2D t = c + of[j];
			con[j] = not_discovered.count(t);

			if(con[j])
			{
				con_edge++;

				v = graph->AddVertex();
				points->InsertNextPoint(t[0], t[1], 0.0);
				indexVertexLink.insert( std::pair<IndexType2D,vtkIdType>(t,v) );

				graph->AddEdge(indexVertexLink[c], indexVertexLink[t]);

				not_discovered.erase(t);
				will_be_processed.insert(t);
			}
		}

		con_vertex = 0;
		for(int j=4; j<8; j++)													//edge connected have lowest priority
		{
			IndexType2D t = c + of[j];
			con[j] = not_discovered.count(t);

			if(con[j]) {
				bool is_real_neighbor = true;

				for(int k=0; k<4; k++)
				{
					if(con[k])
					{
						int length = 0;

						for(unsigned int l=0; l<of[j].GetOffsetDimension(); l++)
							length = length + vcl_abs(of[j][l] - of[k][l]);

						if(length==1)
						{
							is_real_neighbor = false;
							break;
						}
					}
				}

				if(is_real_neighbor)
				{
					con_vertex++;

					v = graph->AddVertex();
					points->InsertNextPoint(t[0], t[1], 0.0);
					indexVertexLink.insert( std::pair<IndexType2D,vtkIdType>(t,v) );

					graph->AddEdge(indexVertexLink[c], indexVertexLink[t]);

					not_discovered.erase(t);
					will_be_processed.insert(t);
				}
			}
		}

		while(will_be_processed.size()>0)
		{
			IndexType2D c = *(will_be_processed.begin());
			vtkIdType v = indexVertexLink[c];
			will_be_processed.erase(c);

			con_edge = 0;
			for(int j=0; j<4; j++)
			{
				IndexType2D t = c + of[j];
				con[j] = not_discovered.count(t);
				conClose[j] = will_be_processed.count(t);

				if(con[j])
				{
					con_edge++;

					v = graph->AddVertex();
					points->InsertNextPoint(t[0], t[1], 0.0);
					indexVertexLink.insert( std::pair<IndexType2D,vtkIdType>(t,v) );

					graph->AddEdge(indexVertexLink[c], indexVertexLink[t]);

					not_discovered.erase(t);
					will_be_processed.insert(t);
				}
				if(conClose[j])
					graph->AddEdge(indexVertexLink[c], indexVertexLink[t]);
			}

			con_vertex = 0;
			for(int j=4; j<8; j++)
			{
				IndexType2D t = c + of[j];
				con[j] = not_discovered.count(t);
				conClose[j] = will_be_processed.count(t);

				if(con[j])
				{
					bool is_real_neighbor = true;

					for(int k=0; k<4; k++)
					{
						if(con[k] || conClose[k])
						{
							int length = 0;

							for(unsigned int l=0; l<of[j].GetOffsetDimension(); l++)
								length = length + vcl_abs(of[j][l] - of[k][l]);

							if(length==1)
							{
								is_real_neighbor = false;
								break;
							}
						}
					}

					if(is_real_neighbor)
					{
						con_vertex++;

						v = graph->AddVertex();
						points->InsertNextPoint(t[0], t[1], 0.0);
						indexVertexLink.insert( std::pair<IndexType2D,vtkIdType>(t,v) );

						graph->AddEdge(indexVertexLink[c], indexVertexLink[t]);

						not_discovered.erase(t);
						will_be_processed.insert(t);
					}
				}

				if(conClose[j])
				{
					IndexType2D nei1, nei2;
					nei1[0] = c[0]; nei1[1] = t[1];
					nei2[0] = t[0]; nei2[1] = c[1];

					if(!indexVertexLink.count(nei1) && !indexVertexLink.count(nei2))
						graph->AddEdge(indexVertexLink[c], indexVertexLink[t]);
				}
			}
		}
		graph->SetPoints(points);
		graphs.push_back(graph);
	}


	return graphs;
}


/*!
  \brief Method for converting skeletons in 3D images into graphs.
  Conversion assumes skeletons as 26-connected structures. The skeleton image is assumed to be a binary image where pixel value 1 indicate a skeleton pixel.
  \n Each pixel is converted into a graph vertex, to avoid zero-loops at crossings the pixel with the highest order induces the crossing vertex (face connectedness > edge connectedness > vertex connectedness).
  \n Skeletons that are not connected result in different graphs.

  \param image A 3-dimensional itk::Image which contains one or more skeletons.
  \returns A vector of vtkMutableUndirectedGraph objects that represents the skeletons in the input image.
*/
std::vector< vtkSmartPointer<vtkUndirectedGraph> > SkeletonImageToGraphFilter::skeletonToGraph3D(const itk::SmartPointer< itk::Image<bool,3> > image, bool fullConnectednessOn)
{
	typedef itk::BinaryImageToShapeLabelMapFilter< itk::Image<bool,3> > ImageToShapeLabelMapFilterType;

	ImageToShapeLabelMapFilterType::Pointer skeletonLabelMapFilter = ImageToShapeLabelMapFilterType::New();
	skeletonLabelMapFilter->ReleaseDataFlagOn();
	skeletonLabelMapFilter->SetFullyConnected(true);
	skeletonLabelMapFilter->SetInput(image);
	skeletonLabelMapFilter->Update();
	LabelMap3DType * skeletonMap = dynamic_cast<LabelMap3DType *>(skeletonLabelMapFilter->GetOutput());
	//itk::SmartPointer<LabelMap2DType> skeletonMap = skeletonLabelMapFilter->GetOutput();

	std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs;

	// Define offsets
	OffsetType3D of[26] = {	//     xyz
			// face connected
			{{-1, 0, 0}},	// 0 - NMM
			{{ 1, 0, 0}},	// 1 - SMM
			{{ 0, 1, 0}},	// 2 - MSM
			{{ 0,-1, 0}},	// 3 - MNM
			{{ 0, 0,-1}},	// 4 - MMN
			{{ 0, 0, 1}},	// 5 - MMS
			// edge connected
			{{-1, 1, 0}},	// 6 - NSM
			{{-1,-1, 0}},  	// 7 - NNM
			{{ 1, 1, 0}},	// 8 - SSM
			{{ 1,-1, 0}},	// 9 - SNM
			{{-1, 0,-1}},	//10 - NMN
			{{ 1, 0,-1}},	//11 - SMN
			{{ 0, 1,-1}},	//12 - MSN
			{{ 0,-1,-1}},	//13 - MNN
			{{-1, 0, 1}},	//14 - NMS
			{{ 1, 0, 1}},	//15 - SMS
			{{ 0, 1, 1}},	//16 - MSS
			{{ 0,-1, 1}},	//17 - MNS
			// corner connected
			{{-1, 1,-1}},	//18 - NSN
			{{-1,-1,-1}},  	//19 - NNN
			{{ 1, 1,-1}},	//20 - SSN
			{{ 1,-1,-1}},	//21 - SNN
			{{-1, 1, 1}},	//22 - NSS
			{{-1,-1, 1}},  	//23 - NNS
			{{ 1, 1, 1}},	//24 - SSS
			{{ 1,-1, 1}}	//25 - SNS
	};

	for(unsigned int i=0; i<skeletonMap->GetNumberOfLabelObjects(); ++i) {
		vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

		LabelMap3DType::LabelObjectType* label_obj = skeletonMap->GetNthLabelObject(i);
		std::map<IndexType3D, vtkIdType, IndexCompare3D> indexVertexLink;
		std::map<IndexType3D, vtkIdType, IndexCompare3D>::iterator it;
#if (ITK_VERSION_MAJOR < 4)
		for(unsigned int j=0; j<label_obj->GetSize(); j++) {
#else
		for(unsigned int j=0; j<label_obj->GetNumberOfPixels(); j++) {
#endif
			IndexType3D zIdx = label_obj->GetIndex(j);

			vtkIdType zVtx = graph->AddVertex();
			points->InsertNextPoint(zIdx[0], zIdx[1], zIdx[2]);
			indexVertexLink[zIdx] = zVtx;
		}

		for(it = indexVertexLink.begin(); it != indexVertexLink.end(); it++) {
			IndexType3D zIdx = (*it).first;

			for(int k=0; k<6; k++) {
				IndexType3D tIdx = zIdx + of[k];

				if(indexVertexLink.count(tIdx)!=0 && graph->GetEdgeId(indexVertexLink[zIdx], indexVertexLink[tIdx])==-1)
					graph->AddEdge(indexVertexLink[zIdx], indexVertexLink[tIdx]);
			}

			if(fullConnectednessOn) {
			    for(int j=6; j<18; j++) {
			        IndexType3D tIdx = zIdx + of[j];

			        if(indexVertexLink.count(tIdx)!=0 && graph->GetEdgeId(indexVertexLink[zIdx], indexVertexLink[tIdx])==-1) {
			            IndexType3D nei1, nei2, nei3, nei4, nei5, nei6;				//neighbor1 & 2 are neighbor voxels-that are face connected to both zIdx and tIdx (-> they have higher rank)
			            bool first_disp = true;										//rest neighbors are neighbor voxels that are edge connected to both zIdx and tIdx (-> have same rank)

			            for(int k=0; k<3; k++) {
			                if(zIdx[k] != tIdx[k] && !first_disp) {
			                    nei1[k] = zIdx[k];
			                    nei2[k] = tIdx[k];
			                    nei3[k] = tIdx[k];
			                    nei4[k] = tIdx[k];
			                    nei5[k] = zIdx[k];
			                    nei6[k] = zIdx[k];
			                }
			                if(zIdx[k] != tIdx[k] && first_disp) {
			                    nei1[k] = tIdx[k];
			                    nei2[k] = zIdx[k];
			                    nei3[k] = zIdx[k];
			                    nei4[k] = zIdx[k];
			                    nei5[k] = tIdx[k];
			                    nei6[k] = tIdx[k];

			                    first_disp = false;
			                }
			                if(zIdx[k] == tIdx[k]) {
			                    nei1[k] = zIdx[k];
			                    nei2[k] = zIdx[k];
			                    nei3[k] = zIdx[k]+1;
			                    nei4[k] = zIdx[k]-1;
			                    nei5[k] = zIdx[k]+1;
			                    nei6[k] = zIdx[k]-1;
			                }
			            }
			            //					std::cout << "zIdx(" << zIdx[0] << "," << zIdx[1] << "," << zIdx[2] << "); tIdx(" << tIdx[0] << "," << tIdx[1] << "," << tIdx[2] << "); nei1(" << nei1[0] << "," << nei1[1] << "," << nei1[2] << "); nei2(" << nei2[0] << "," << nei2[1] << "," << nei2[2];
			            //					std::cout << "); nei3(" << nei3[0] << "," << nei3[1] << "," << nei3[2] << "); nei4(" << nei4[0] << "," << nei4[1] << "," << nei4[2] << "); nei5(" << nei5[0] << "," << nei5[1] << "," << nei5[2] << "); nei6(" << nei6[0] << "," << nei6[1] << "," << nei6[2] << ");" << std::endl;

			            if( indexVertexLink.count(nei1)==0 && indexVertexLink.count(nei2)==0 &&																															//if there are higher rank neighbor voxels no edge is needed
			                    (indexVertexLink.count(nei3)==0 || graph->GetEdgeId(indexVertexLink[zIdx], indexVertexLink[nei3]) == -1 || graph->GetEdgeId(indexVertexLink[tIdx], indexVertexLink[nei3]) == -1) &&		//if there are same rank neighbor voxels that are
			                    (indexVertexLink.count(nei4)==0 || graph->GetEdgeId(indexVertexLink[zIdx], indexVertexLink[nei4]) == -1 || graph->GetEdgeId(indexVertexLink[tIdx], indexVertexLink[nei4]) == -1) &&		//already connected to both zIdx and tIdx no edge is needed
			                    (indexVertexLink.count(nei5)==0 || graph->GetEdgeId(indexVertexLink[zIdx], indexVertexLink[nei5]) == -1 || graph->GetEdgeId(indexVertexLink[tIdx], indexVertexLink[nei5]) == -1) &&
			                    (indexVertexLink.count(nei6)==0 || graph->GetEdgeId(indexVertexLink[zIdx], indexVertexLink[nei6]) == -1 || graph->GetEdgeId(indexVertexLink[tIdx], indexVertexLink[nei6]) == -1) )
			                graph->AddEdge(indexVertexLink[zIdx], indexVertexLink[tIdx]);
			        }
			    }

			    for(int j=18; j<26; j++) {
			        IndexType3D tIdx = zIdx + of[j];

			        if(indexVertexLink.count(tIdx)!=0 && graph->GetEdgeId(indexVertexLink[zIdx], indexVertexLink[tIdx])==-1) {
			            IndexType3D nei1, nei2, nei3, nei4, nei5, nei6;			//all neighbor voxels have higher rank (-> if one of them is present no edge is needed)

			            nei1[0] = tIdx[0]; nei1[1] = zIdx[1]; nei1[2] = zIdx[2];
			            nei2[0] = zIdx[0]; nei2[1] = tIdx[1]; nei2[2] = zIdx[2];
			            nei3[0] = zIdx[0]; nei3[1] = zIdx[1]; nei3[2] = tIdx[2];
			            nei4[0] = zIdx[0]; nei4[1] = tIdx[1]; nei4[2] = tIdx[2];
			            nei5[0] = tIdx[0]; nei5[1] = zIdx[1]; nei5[2] = tIdx[2];
			            nei6[0] = tIdx[0]; nei6[1] = tIdx[1]; nei6[2] = zIdx[2];

			            //					std::cout << "zIdx(" << zIdx[0] << "," << zIdx[1] << "," << zIdx[2] << "); tIdx(" << tIdx[0] << "," << tIdx[1] << "," << tIdx[2] << "); nei1(" << nei1[0] << "," << nei1[1] << "," << nei1[2] << "); nei2(" << nei2[0] << "," << nei2[1] << "," << nei2[2];
			            //					std::cout << "); nei3(" << nei3[0] << "," << nei3[1] << "," << nei3[2] << "); nei4(" << nei4[0] << "," << nei4[1] << "," << nei4[2] << "); nei5(" << nei5[0] << "," << nei5[1] << "," << nei5[2] << "); nei6(" << nei6[0] << "," << nei6[1] << "," << nei6[2] << ");" << std::endl;

			            if(indexVertexLink.count(nei1)==0 && indexVertexLink.count(nei2)==0 && indexVertexLink.count(nei3)==0 && indexVertexLink.count(nei4)==0 && indexVertexLink.count(nei5)==0 && indexVertexLink.count(nei6)==0)
			                graph->AddEdge(indexVertexLink[zIdx], indexVertexLink[tIdx]);
			        }
			    }
			}
		}

		graph->SetPoints(points);
		graphs.push_back(graph);
	}

	return graphs;
}




