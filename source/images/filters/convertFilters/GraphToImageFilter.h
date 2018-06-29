///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphToImageFilter.h                                                 //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-08-13                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPHTOIMAGEFILTER_H_
#define GRAPHTOIMAGEFILTER_H_

#include "itkImage.h"
#include "itkSmartPointer.h"

#include "vtkSmartPointer.h"
#include "vtkUndirectedGraph.h"


class GraphToImageFilter
{
public:
    GraphToImageFilter();
    virtual ~GraphToImageFilter();

    static void GraphToImage(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, itk::SmartPointer< itk::Image<unsigned char, 3> > &image, std::string radiusArrayName = "radius", double sampleFactor = 1., double frameFactor = 50.,
            unsigned char foreground = 255, unsigned char background = 0);

private:
    static void ComputeBounds(std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs, double frameFactor, double*& boundsAll, double*& offset, double *& rim);
};

#endif /* GRAPHTOIMAGEFILTER_H_ */
