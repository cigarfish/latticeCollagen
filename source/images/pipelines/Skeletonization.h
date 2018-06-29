///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Skeletonization.h                                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-12-14                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SKELETONIZATION_H_
#define SKELETONIZATION_H_

#include "BasePipeline.h"

#include "../filters/imageFilters/itkBinaryThinningImageFilter3D.h"

//#include "../filters/imageFilters/itkAverageOutwardFluxImageFilter.h"
//#include "../filters/imageFilters/itkMedialCurveImageFilter.h"


class SkeletonizationPipeline : public BasePipeline<3>
{
    typedef itk::BinaryThinningImageFilter3D<CScalarVoImageType,CScalarVoImageType>     ThinningImageFilterType;

public:
    SkeletonizationPipeline();
    virtual ~SkeletonizationPipeline();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);


    std::string m_fullFilename;
    std::string m_path;
    std::string m_filename;
    std::string m_fileExtension;
};

#endif /* SKELETONIZATION_H_ */
