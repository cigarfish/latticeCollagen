///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ObjectsToMatrixFilter.h                                              //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-04-14                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef OBJECTSTOMATRIXFILTER_H_
#define OBJECTSTOMATRIXFILTER_H_

#include "BasePipeline.h"

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkLabelMapToLabelImageFilter.h"


class ObjectsToMatrixFilter : public BasePipeline<3>
{
    typedef long                            LScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType> CRGBPixelType;

    typedef itk::Image<LScalarPixelType, 3>             LScalarVoImageType;
    typedef itk::Image<CRGBPixelType, 3>                CRGBVoImageType;

    typedef itk::ImageRegion<3>                         ImageRegionType;
    typedef typename LScalarVoImageType::SpacingType    SpacingType;

    typedef itk::ImageFileReader<CRGBVoImageType>       CRGBVoReaderType;

    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>               ImageToShapeLabelMapFilterType;
    typedef ImageToShapeLabelMapFilterType::OutputImageType                         ShapeLabelMapType;
    typedef ShapeLabelMapType::LabelObjectType                                      ShapeLabelObjectType;
    typedef itk::LabelMapToLabelImageFilter<ShapeLabelMapType, LScalarVoImageType>  LabelMapToLabelImageFilterType;
    typedef itk::BinaryThresholdImageFilter<LScalarVoImageType, CScalarVoImageType> ThresholdFilterType;

public:
    ObjectsToMatrixFilter();
    virtual ~ObjectsToMatrixFilter();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteMetaFile();

    void SetupMatrixImage(itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> labelMap, itk::SmartPointer<CRGBVoImageType> rgbImage);
    void SetupBorderRegions(CScalarVoImageType::Pointer &image);
    void RemoveLabelObjectsAtDatasetBorders(LScalarVoImageType::Pointer cellLabelImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap);
    void FilterLabelObjectsBySlide(LScalarVoImageType::Pointer cellLabelImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap);


    std::string mFullFilename;
    std::string mPath;
    std::string mFilename;
    std::string mExtension;

    std::string mOverlayFullFilename;

    std::string mLogFilenameSave;

    SpacingType mSpacing;
    double mVoxelVolume;

    int mRow;
    int mCol;

    int mDimX;
    int mDimY;
    int mDimZ;

    ImageRegionType mDataSetFaces[6];
    ImageRegionType mSlide;

    bool mFullyConnected;

    std::map<unsigned long, CRGBPixelType> mLabelObjectToColor;
    std::map<unsigned long, int*> mLabelObjectToMatrixPosition;
    std::map<unsigned long, itk::Index<3> > mLabelObjectToCentroid;
    std::map<unsigned long, double> mLabelObjectToVolume;

    bool mWithFiltering;
    unsigned int mFilterBySlide;
    bool mWriteMetaInfoFile;

    itk::SmartPointer<ShapeLabelMapType> mpMatrixMap;
};

#endif /* OBJECTSTOMATRIXFILTER_H_ */
