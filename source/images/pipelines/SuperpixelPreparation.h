///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SuperpixelPreparation.h                                              //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-03-31                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SUPERPIXELPREPARATION_H_
#define SUPERPIXELPREPARATION_H_

#include "LabelMapGraphBasePipeline.h"

#include "../filters/imageFilters/SLICImageFilter.h"


template< unsigned int VImageDimension > class SuperpixelPreparation : public LabelMapGraphBasePipeline<VImageDimension>
{
private:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

    typedef LabelMapGraphBasePipeline<VImageDimension>                  Superclass;

    typedef typename Superclass::LScalarImageType                       LScalarImageType;
    typedef typename Superclass::CRGBImageType                          CRGBImageType;

    typedef typename Superclass::SpacingType                            SpacingType;

    typedef typename Superclass::LabelMapType                           LabelMapType;
    typedef typename Superclass::UnaryFeatureType                       UnaryFeatureType;
    typedef typename Superclass::BinaryFeatureType                      BinaryFeatureType;

    typedef typename Superclass::LabelImageToGraphFilterType            LabelImageToGraphFilterType;
    typedef itk::SLICImageFilter<CRGBImageType, LScalarImageType>       SLICFilterType;

public:
    SuperpixelPreparation();
    virtual ~SuperpixelPreparation();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);

    void AddFeaturesToLabelMap();

    int mEntrypoint;

    bool mSuperpixelBySize;
    SpacingType mSuperpixelSpacing;
    unsigned int mSuperpixelNumber;

    bool* mWithFeature;

    std::string mLogFilenameSave;
    std::string mFilenameSave;
};

#include "SuperpixelPreparation.tpp"

#endif /* SUPERPIXELPREPARATION_H_ */
