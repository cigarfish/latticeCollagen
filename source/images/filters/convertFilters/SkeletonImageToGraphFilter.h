///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SkeletonImageToGraphFilter.h                                         //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-11-09                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SKELETONIMAGETOGRAPHFILTER_H_
#define SKELETONIMAGETOGRAPHFILTER_H_

#include <vector>

#include "itkImage.h"
#include "itkSmartPointer.h"

#include "vtkMutableUndirectedGraph.h"
#include "vtkSmartPointer.h"

class vtkImageData;
class vtkMutableUndirectedGraph;


/*!
  \brief A class for converting a skeleton to a vtkGraph structure.
  This class provides methods for converting skeleton structures in 2D or 3D itk::Image objects into vtkGraph objects.
  \n Make sure that the input image contains only skeleton structures (use this converter after thinning).
  \n Conversion assumes skeletons as 8-connected structures in the 2D case.
  \n Each pixel is converted into a graph vertex, to avoid zero-loops at crossings the pixel with the highest order induces the crossing vertex (face connectedness > edge connectedness > vertex connectedness).
  \n Skeletons that are not connected result in different graphs.
*/
class SkeletonImageToGraphFilter
{
public:
    SkeletonImageToGraphFilter();
    ~SkeletonImageToGraphFilter();

    static std::vector< vtkSmartPointer<vtkUndirectedGraph> > skeletonToGraph2D(itk::SmartPointer< itk::Image<bool,2> > image);
    static std::vector< vtkSmartPointer<vtkUndirectedGraph> > skeletonToGraph3D(itk::SmartPointer< itk::Image<bool,3> > image, bool fullConnectednessOn);
};



#endif /* SKELETONIMAGETOGRAPHFILTER_H_ */
