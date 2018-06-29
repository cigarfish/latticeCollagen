///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Lobule2DHandlingFilter.h                                             //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-04-24                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LOBULE2DHANDLINGFILTER_H_
#define LOBULE2DHANDLINGFILTER_H_

#include <itkImage.h>
#include <itkEllipsoidInteriorExteriorSpatialFunction.h>
#include <itkLabelImageToShapeLabelMapFilter.h>
#include <itkLineIterator.h>


class Lobule2DHandlingFilter
{
public:
    typedef unsigned char                                               CScalarPixelType;
    typedef float                                                       FScalarPixelType;
    typedef long                                                        LScalarPixelType;

    typedef itk::Image<CScalarPixelType, 2>                             CScalar2DImageType;
    typedef itk::Image<FScalarPixelType, 2>                             FScalar2DImageType;
    typedef itk::Image<LScalarPixelType, 2>                             LScalar2DImageType;
    typedef itk::LabelMap<itk::ShapeLabelObject<CScalarPixelType, 2> >  LabelMap2DType;

    typedef typename CScalar2DImageType::Pointer                        CScalar2DImagePointerType;
    typedef typename FScalar2DImageType::Pointer                        FScalar2DImagePointerType;
    typedef typename LScalar2DImageType::Pointer                        LScalar2DImagePointerType;

    typedef typename CScalar2DImageType::RegionType                     CScalar2DRegionType;
    typedef typename CScalar2DImageType::SpacingType                    CScalar2DSpacingType;
    typedef typename CScalar2DImageType::SizeType                       CScalar2DSizeType;
    typedef typename CScalar2DImageType::IndexType                      CScalar2DIndexType;
    typedef typename CScalar2DIndexType::IndexValueType                 CScalar2DIndexValueType;
    typedef itk::Point<double, 2>                                       PointType;
    typedef itk::Vector<double, 2>                                      VectorType;

    typedef itk::EllipsoidInteriorExteriorSpatialFunction<2>            EllipsoidFunctionType;
    typedef EllipsoidFunctionType::InputType                            EllipsoidFunctionVectorType;
    typedef itk::LineIterator<CScalar2DImageType>                       LineIteratorType;


    typedef itk::LabelImageToShapeLabelMapFilter<CScalar2DImageType, LabelMap2DType>    LabelImageToShapeLabelMapFilterType;

    Lobule2DHandlingFilter();
    virtual ~Lobule2DHandlingFilter();

    void ComputeVeinCenterPointsIn2DImage(CScalar2DImagePointerType veinLabelImage, std::map<int, CScalar2DIndexType> &veinIdToCenterPoints);
    void DrawLine(CScalar2DImagePointerType &image, CScalar2DIndexType p1, CScalar2DIndexType p2);
    void DrawEllipse(CScalar2DImagePointerType &image, CScalar2DIndexType p1, CScalar2DIndexType p2, double radius2);
};

#endif /* LOBULE2DHANDLINGFILTER_H_ */
