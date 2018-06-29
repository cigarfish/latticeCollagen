///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelImageToGraphFilter.h                                            //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-08-07                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LABELIMAGETOGRAPHFILTER_H_
#define LABELIMAGETOGRAPHFILTER_H_

#include "itkImage.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkSmartPointer.h"

#include "vtkSmartPointer.h"
#include "vtkUndirectedGraph.h"


/*!
  \brief Method for constructing a neighborhood graph based on label objects within an image.
*/
template< unsigned int VImageDimension > class LabelImageToGraphFilter
{
public:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

    LabelImageToGraphFilter();
    ~LabelImageToGraphFilter();

    static vtkSmartPointer<vtkUndirectedGraph> LabelImageToGraph(itk::SmartPointer< itk::Image<unsigned long, VImageDimension> > image, bool ignoreBackgroundLabel=false);
};

#include "LabelImageToGraphFilter.tpp"

#endif /* LABELIMAGETOGRAPHFILTER_H_ */
