///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelMapGraphBasePipeline.h                                          //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-04-01                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef LABELMAPGRAPHBASEPIPELINE_H_
#define LABELMAPGRAPHBASEPIPELINE_H_

#include "BasePipeline.h"

#include "itkCastImageFilter.h"
#include "itkImageAdaptor.h"
#include "itkLabelContourImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkMaskNegatedImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

#include <vtkUndirectedGraph.h>
#include <vtkSmartPointer.h>

#include "../filters/convertFilters/LabelImageToLabelMapGraphFilter.h"
#include "../filters/convertFilters/LabelImageToGraphFilter.h"
#include "../tools/LabelMapGraph.h"


template< unsigned int VImageDimension > class LabelMapGraphBasePipeline : public BasePipeline<VImageDimension>
{
protected:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

    typedef typename BasePipeline<VImageDimension>::CScalarPixelType    CScalarPixelType;
    typedef unsigned long int                                           LScalarPixelType;
    typedef unsigned int                                                IScalarPixelType;
    typedef itk::RGBPixel<CScalarPixelType>                             CRGBPixelType;

    class BlueChannelPixelAccessor
    {
    public:
        typedef itk::RGBPixel<unsigned char>    InternalType;
        typedef unsigned char                   ExternalType;
        static ExternalType Get(const InternalType &input) { return static_cast<ExternalType>(input.GetBlue()); }
    };

    typedef typename BasePipeline<VImageDimension>::CScalarVoImageType  CScalarImageType;
    typedef itk::Image<LScalarPixelType, ImageDimension>                LScalarImageType;
    typedef itk::Image<IScalarPixelType, ImageDimension>                IScalarImageType;
    typedef itk::Image<CRGBPixelType, ImageDimension>                   CRGBImageType;

    typedef typename CScalarImageType::SpacingType                      SpacingType;

    typedef itk::ImageFileReader<CRGBImageType>                         CRGBReaderType;
    typedef itk::ImageFileReader<IScalarImageType>                      IScalarReaderType;
    typedef itk::ImageFileWriter<IScalarImageType>                      IScalarWriterType;
    typedef itk::ImageFileWriter<CRGBImageType>                         CRGBWriterType;


    typedef itk::CastImageFilter<LScalarImageType, IScalarImageType>        CastLIImageFilterType;
    typedef itk::CastImageFilter<IScalarImageType, LScalarImageType>        CastILImageFilterType;
    typedef itk::MaskNegatedImageFilter<CRGBImageType, LScalarImageType>    MaskNegatedFilterType;
    typedef itk::MinimumMaximumImageCalculator<LScalarImageType>            LMinMaxCalculatorType;
    typedef itk::MinimumMaximumImageCalculator<IScalarImageType>            IMinMaxCalculatorType;

    typedef LabelImageToGraphFilter<ImageDimension>                                         LabelImageToGraphFilterType;
    typedef itk::LabelImageToLabelMapGraphFilter<CRGBImageType, LScalarImageType>           LabelImageToLabelMapGraphFilterType;
    typedef typename LabelImageToLabelMapGraphFilterType::OutputImageType                   LabelMapType;
    typedef typename LabelMapType::Pointer                                                  LabelMapTypePointer;

    typedef typename LabelMapType::UnaryFeatures                                            UnaryFeatureType;
    typedef typename LabelMapType::BinaryFeatures                                           BinaryFeatureType;

    typedef itk::ImageAdaptor<CRGBImageType, BlueChannelPixelAccessor>                      BlueAdaptorType;
    typedef itk::LabelContourImageFilter<LScalarImageType, LScalarImageType>                LabelContourImageFilterType;
    typedef itk::LabelOverlayImageFilter<BlueAdaptorType, LScalarImageType, CRGBImageType>  LabelOverlayImageFilterType;

public:
    LabelMapGraphBasePipeline();
    virtual ~LabelMapGraphBasePipeline();

    vtkSmartPointer<vtkUndirectedGraph> GetObjectGraph() { return mObjectGraph; };

    virtual void Update() = 0;

protected:
    virtual void ParseParameterContext() = 0;
    virtual void WriteLogFile(std::string timeStamp) = 0;

    void ReadOriginalImage(std::string filename);
    void ReadObjectImage(std::string filename);
    void ReadObjectGraph(std::string filename);

    void InitLabelMap();

    void SaveObjectImage(std::string filename);
    void SaveObjectContourOverlayImage(std::string filename);
    void SaveObjectOverlayImage(std::string filename);
    void SaveObjectGraph(std::string filename);

    //Paramters used by the pipeline
    std::string mRawImageFullFilename;
    std::string mRawImagePath;
    std::string mRawImageFilename;
    std::string mRawImageFileExtension;

    std::string mObjectImageFullFilename;
    std::string mObjectImagePath;
    std::string mObjectImageFilename;
    std::string mObjectImageFileExtension;

    std::string mObjectGraphFullFilename;

    SpacingType mVoxelSpacing;

    typename LabelMapType::Pointer mLabelMap;
    itk::SmartPointer<CRGBImageType> mOriginalImage;
    itk::SmartPointer<LScalarImageType> mObjectImage;
    itk::SmartPointer<LScalarImageType> mObjectImage2;
    vtkSmartPointer<vtkUndirectedGraph> mObjectGraph;
};

#include "LabelMapGraphBasePipeline.tpp"

#endif /* LABELMAPGRAPHBASEPIPELINE_H_ */
