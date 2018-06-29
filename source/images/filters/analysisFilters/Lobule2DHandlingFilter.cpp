///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Lobule2DHandlingFilter.cpp                                           //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-04-24                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "Lobule2DHandlingFilter.h"

#include <itkFloodFilledSpatialFunctionConditionalIterator.h>


Lobule2DHandlingFilter::Lobule2DHandlingFilter()
{
    // TODO Auto-generated constructor stub

}


Lobule2DHandlingFilter::~Lobule2DHandlingFilter()
{
    // TODO Auto-generated destructor stub
}


void Lobule2DHandlingFilter::ComputeVeinCenterPointsIn2DImage(CScalar2DImagePointerType veinLabelImage, std::map<int, CScalar2DIndexType> &veinIdToCenterPoints)
{
    typename LabelImageToShapeLabelMapFilterType::Pointer labelImageToLabelMap = LabelImageToShapeLabelMapFilterType::New();
    labelImageToLabelMap->SetInput(veinLabelImage);
    labelImageToLabelMap->Update();

    for(unsigned int j=0; j<labelImageToLabelMap->GetOutput()->GetNumberOfLabelObjects(); j++) {
        CScalar2DIndexType idx;
        veinLabelImage->TransformPhysicalPointToIndex(labelImageToLabelMap->GetOutput()->GetNthLabelObject(j)->GetCentroid(), idx);

        veinIdToCenterPoints.insert(std::pair<int, CScalar2DIndexType>(labelImageToLabelMap->GetOutput()->GetNthLabelObject(j)->GetLabel(), idx));
    }
}


void Lobule2DHandlingFilter::DrawLine(CScalar2DImagePointerType &image, CScalar2DIndexType p1, CScalar2DIndexType p2)
{
    LineIteratorType it1(image, p1, p2);
    while(!it1.IsAtEnd()) {
        it1.Set(255);
        ++it1;
    }
}


void Lobule2DHandlingFilter::DrawEllipse(CScalar2DImagePointerType &image, CScalar2DIndexType p1, CScalar2DIndexType p2, double radius2)
{
    PointType start, end;
    start[0] = p1[0];
    start[1] = p1[1];
    end[0] = p2[0];
    end[1] = p2[1];

    EllipsoidFunctionVectorType axes;
    axes[0] = start.EuclideanDistanceTo(end);                                   //axis length in pxl
    axes[1] = 2*radius2;

    EllipsoidFunctionVectorType center;
    center[0] = (start[0] + end[0]) / (2.);                                     //center coord in pxl
    center[1] = (start[1] + end[1]) / (2.);

    VectorType vec0, vec1;
    vec0 = end - start;
    vec0.Normalize();

    vec1[0] = -vec0[1];                                                         //vec1 orthogonal to vec0
    vec1[1] = vec0[0];

    double data[] = {vec0[0], vec0[1], vec1[0], vec1[1]};
    vnl_matrix<double> orientations (data, 2, 2);

//    std::cout << "lineSegStart = " << start << " lineSegEnd = " << end << std::endl;
//    std::cout << "axes = " << axes << std::endl;
//    std::cout << "center = " << center << std::endl;
//    std::cout << "vec0 = " << vec0 << "vec1 = " << vec1 << std::endl;

    typename EllipsoidFunctionType::Pointer spatialFunc = EllipsoidFunctionType::New();
    spatialFunc->SetAxes(axes);
    spatialFunc->SetCenter(center);
    spatialFunc->SetOrientations(orientations);

    CScalar2DIndexType seedPos;
    seedPos[0] = static_cast<CScalar2DIndexValueType>(center[0]);
    seedPos[1] = static_cast<CScalar2DIndexValueType>(center[1]);

    int numInteriorPixels1 = 0;
    unsigned char interiorPixelValue = 255;

    itk::FloodFilledSpatialFunctionConditionalIterator<CScalar2DImageType, EllipsoidFunctionType> sfi =
            itk::FloodFilledSpatialFunctionConditionalIterator<CScalar2DImageType, EllipsoidFunctionType>(image, spatialFunc, seedPos);

    for(; !sfi.IsAtEnd(); ++sfi) {
        sfi.Set(interiorPixelValue);
        ++numInteriorPixels1;
    }
}
