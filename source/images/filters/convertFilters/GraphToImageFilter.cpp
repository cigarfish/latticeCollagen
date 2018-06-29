///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphToImageFilter.cpp                                               //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-08-13                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "GraphToImageFilter.h"

#include "itkGroupSpatialObject.h"
#include "itkMaskedSpatialObjectToImageFilter.h"
#include "itkSphereSpatialFunction.h"
#include "itkTubeSpatialObject.h"
#include "itkVesselTubeSpatialObject.h"

#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkEdgeListIterator.h"
//#include "vtkPoints.h"


GraphToImageFilter::GraphToImageFilter()
{
    // TODO Auto-generated constructor stub

}


GraphToImageFilter::~GraphToImageFilter()
{
    // TODO Auto-generated destructor stub
}


void GraphToImageFilter::ComputeBounds(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, double frameFactor, double*& boundsAll, double*& offset, double *& rim)
{
    std::vector<double*> bounds;
 //   double rim[3];

    for(unsigned int i=0; i<graphs.size(); i++) {
        double b[6];

        graphs[i]->ComputeBounds();
        graphs[i]->GetBounds(b);

        bounds.push_back(b);
    }

    for(unsigned int i=0; i<6; i++)
        boundsAll[i] = bounds[0][i];
/*
    for(unsigned int i=0; i<6; i=i+2){
      if( boundsAll[i] > 5 )  
         boundsAll[i] -= 6;
      else 
        boundsAll[i] = 0;
    }
    for(unsigned int i=1; i<6; i=i+2)
        boundsAll[i] += 6;
*/

    for(unsigned int i=1; i<bounds.size(); i++) {
        boundsAll[0] = std::min(bounds[i][0], boundsAll[0]);
        boundsAll[2] = std::min(bounds[i][2], boundsAll[2]);
        boundsAll[4] = std::min(bounds[i][4], boundsAll[4]);
        boundsAll[1] = std::max(bounds[i][1], boundsAll[1]);
        boundsAll[3] = std::max(bounds[i][3], boundsAll[3]);
        boundsAll[5] = std::max(bounds[i][5], boundsAll[5]);
    }

    std::cout << "Original bounds: bounds[0] = " << boundsAll[0] << " bounds[1] = " << boundsAll[1] << " bounds[2] = " << boundsAll[2] << " bounds[3] = " << boundsAll[3] << " bounds[4] = " << boundsAll[4] <<
            " bounds[5] = " << boundsAll[5] << std::endl;

    ofstream F;
    F.open( "../../input/image.txt", ios::out|ios::app );

      F << "Original bounds: bounds[0] = " << boundsAll[0] << " bounds[1] = " << boundsAll[1] << " bounds[2] = " << boundsAll[2] << " bounds[3] = " << boundsAll[3] << " bounds[4] = " << boundsAll[4] <<
            " bounds[5] = " << boundsAll[5] << std::endl;

    rim[0] = (boundsAll[1]-boundsAll[0])/frameFactor;
    rim[1] = (boundsAll[3]-boundsAll[2])/frameFactor;
    rim[2] = (boundsAll[5]-boundsAll[4])/frameFactor;

//    if(boundsAll[0]-rim[0]<0)
//        offset[0] = 0;
//    else
        offset[0] = boundsAll[0]-rim[0];

//    if(boundsAll[2]-rim[1]<0)
//        offset[1] = 0;
//    else
        offset[1] = boundsAll[2]-rim[1];

//    if(boundsAll[4]-rim[2]<0)
//        offset[2] = 0;
//    else
        offset[2] = boundsAll[4]-rim[2];

    boundsAll[0] -= offset[0];
    boundsAll[1] -= offset[0];
    boundsAll[2] -= offset[1];
    boundsAll[3] -= offset[1];
    boundsAll[4] -= offset[2];
    boundsAll[5] -= offset[2];

    std::cout << "Corrected bounds: bounds[0] = " << boundsAll[0] << " bounds[1] = " << boundsAll[1] << " bounds[2] = " << boundsAll[2] << " bounds[3] = " << boundsAll[3] << " bounds[4] = " << boundsAll[4] <<
            " bounds[5] = " << boundsAll[5] << std::endl;


    
    F << "Corrected bounds: bounds[0] = " << boundsAll[0] << " bounds[1] = " << boundsAll[1] << " bounds[2] = " << boundsAll[2] << " bounds[3] = " << boundsAll[3] << " bounds[4] = " << boundsAll[4] <<
            " bounds[5] = " << boundsAll[5] << std::endl;

    F << "offset: x = " << offset[0] << " y = " << offset[1] << " z = " << offset[2] << std::endl;
     F << "rim: x = " << rim[0] << " y = " << rim[1] << " z = " << rim[2] << std::endl;
F.close();
}


void GraphToImageFilter::GraphToImage(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, itk::SmartPointer< itk::Image<unsigned char, 3> > &image, std::string radiusArrayName, double sampleFactor, double frameFactor,
        unsigned char foreground, unsigned char background)
{
    typedef unsigned char                           PixelType;

    typedef itk::Image<PixelType, 3>                ImageType;

    typedef itk::TubeSpatialObject<3>               TubeType;
    typedef itk::TubeSpatialObjectPoint<3>          TubePointType;

    typedef itk::SphereSpatialFunction<3>           SphereType;
    
    typedef TubePointType::CovariantVectorType      VectorType;
    typedef itk::GroupSpatialObject<3>              GroupType;
    typedef GroupType::TransformType                TransformType;

    typedef itk::MaskedSpatialObjectToImageFilter<GroupType, ImageType>     SpatialObjectToImageFilterType;

    double* bounds = new double[6];
    double* offset = new double[3];
    double* rim    = new double[3];

    ComputeBounds(graphs, frameFactor, bounds, offset, rim);

    ImageType::SizeType size;
    size[0] = (bounds[1]+2*rim[0])*sampleFactor;
    size[1] = (bounds[3]+2*rim[1])*sampleFactor;
    size[2] = (bounds[5]+2*rim[2])*sampleFactor;

    ImageType::SpacingType spacing;
    spacing[0] =  1. / sampleFactor;
    spacing[1] =  1. / sampleFactor;
    spacing[2] =  1. / sampleFactor;

    GroupType::Pointer group = GroupType::New();

        ofstream F;
    F.open( "../../input/image.txt", ios::out|ios::app );


    unsigned int numEdges = 0;

    for(unsigned int i=0; i<graphs.size(); i++) {
        vtkDoubleArray* edgeRadius = vtkDoubleArray::SafeDownCast( graphs[i]->GetEdgeData()->GetArray(radiusArrayName.c_str()) );

        vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
        graphs[i]->GetEdges(it);

        while(it->HasNext()) {
            vtkEdgeType e = it->Next();

            TubeType::PointListType list;

            TubePointType source, target;
            source.SetPosition(graphs[i]->GetPoint(e.Source)[0]-offset[0], graphs[i]->GetPoint(e.Source)[1]-offset[1], graphs[i]->GetPoint(e.Source)[2]-offset[2]);
            source.SetRadius(edgeRadius->GetValue(e.Id));
            list.push_back(source);
//            std::cout << "tube " << numEdges << " starting point = (" << (graphs[i]->GetPoint(e.Source)[0]-offset[0])*sampleFactor << ", " << (graphs[i]->GetPoint(e.Source)[1]-offset[1])*sampleFactor
//                    << ", " << (graphs[i]->GetPoint(e.Source)[2]-offset[2])*sampleFactor << ")" << std::endl;
            F << "tube " << numEdges << " starting point = (" << (graphs[i]->GetPoint(e.Source)[0]-offset[0])*sampleFactor << ", " << (graphs[i]->GetPoint(e.Source)[1]-offset[1])*sampleFactor
                    << ", " << (graphs[i]->GetPoint(e.Source)[2]-offset[2])*sampleFactor << ")" << std::endl;

            target.SetPosition(graphs[i]->GetPoint(e.Target)[0]-offset[0], graphs[i]->GetPoint(e.Target)[1]-offset[1], graphs[i]->GetPoint(e.Target)[2]-offset[2]);
            target.SetRadius(edgeRadius->GetValue(e.Id));
            list.push_back(target);
//            std::cout << "tube " << numEdges << " end point = (" << (graphs[i]->GetPoint(e.Target)[0]-offset[0])*sampleFactor << ", " << (graphs[i]->GetPoint(e.Target)[1]-offset[1])*sampleFactor
//                    << ", " << (graphs[i]->GetPoint(e.Target)[2]-offset[2])*sampleFactor << ")" << std::endl;
//            std::cout << "tube " << numEdges << " radius = " << edgeRadius->GetValue(e.Id) << "(*" << sampleFactor << ")" << std::endl;
            F << "tube " << numEdges << " end point = (" << (graphs[i]->GetPoint(e.Target)[0]-offset[0])*sampleFactor << ", " << (graphs[i]->GetPoint(e.Target)[1]-offset[1])*sampleFactor
                    << ", " << (graphs[i]->GetPoint(e.Target)[2]-offset[2])*sampleFactor << ")" << std::endl;
            F << "tube " << numEdges << " radius = " << edgeRadius->GetValue(e.Id) << "(*" << sampleFactor << ")" << std::endl;

            TubeType::Pointer tube1 = TubeType::New();
            TubeType::Pointer tube2 = TubeType::New();
            tube1->SetPoints(list);
            tube2->SetPoints(list);

//            std::cout << "Add tube number " << numEdges << " to group" << std::endl;
            group->AddSpatialObject(tube1);

            tube2->SetEndType(1);
            group->AddSpatialObject(tube2);

            SphereType::Pointer sphere1 = SphereType::New();

            SphereType::InputType center;
              center[0] = 0;
              center[1] = 0;
              center[2] = 0;
            sphere1->SetCenter( center );
            sphere1->SetRadius( 0.25 );


//            group->AddSpatialObject(sphere1);


            numEdges++;
        }
    }

    SpatialObjectToImageFilterType::Pointer converter = SpatialObjectToImageFilterType::New();
    converter->SetNumberOfThreads(4);
    converter->SetInput(group);
    converter->SetInsideValue(foreground);
    converter->SetOutsideValue(background);
    converter->SetSize(size);
    converter->SetSpacing(spacing);
    converter->SetMaskResampleFactor(4);
    converter->SetMaskDilationSize(2);
    converter->Update();

    image = converter->GetOutput();
    image->DisconnectPipeline();

    delete [] bounds;
    delete [] offset;
}
