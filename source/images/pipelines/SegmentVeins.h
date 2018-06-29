///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentVeins.h                                                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-03-16                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SEGMENTVEINS_H_
#define SEGMENTVEINS_H_

#include "BasePipeline.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkIndex.h"

#include "itkAddImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#if (ITK_VERSION_MAJOR < 4 || ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR < 4 || ITK_VERSION_MAJOR > 4)
#include "itkCompose2DVectorImageFilter.h"
#include "itkCompose3DVectorImageFilter.h"
#else
#include "itkComposeImageFilter.h"
#endif
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkLabelContourImageFilter.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkNearestNeighborExtrapolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"

#include "../filters/imageFilters/MultiChannelConnectedThresholdImageFilter.h"


class SegmentVeins : public BasePipeline<3>
{
    typedef double                                      CoordRepType;
    typedef unsigned long                               LScalarPixelType;
    typedef float                                       FScalarPixelType;
    typedef itk::Vector<CScalarPixelType, 2>            CVector2DPixelType;
    typedef itk::Vector<CScalarPixelType, 3>            CVector3DPixelType;
    typedef itk::RGBPixel<CScalarPixelType>             RGBPixelType;

    typedef itk::Image<LScalarPixelType, 3>             LScalarVoImageType;
    typedef itk::Image<FScalarPixelType, 3>             FScalarVoImageType;
    typedef itk::Image<CVector2DPixelType, 3>           Vector2DVoImageType;
    typedef itk::Image<CVector3DPixelType, 3>           Vector3DVoImageType;
    typedef itk::Image<RGBPixelType, 3>                 RGBVoImageType;

    typedef itk::BinaryBallStructuringElement<CScalarVoImageType::PixelType, 3>     StructuringElementType;

    typedef itk::ImageFileWriter<RGBVoImageType>        RGBVoWriterType;

    typedef itk::AddImageFilter<FScalarVoImageType, FScalarVoImageType, FScalarVoImageType>                             AddImageFilterType;
    typedef itk::BinaryDilateImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>                DilateImageFilterType;
    typedef itk::BinaryMorphologicalClosingImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>  ClosingImageFilterType;
    typedef itk::BinaryMorphologicalOpeningImageFilter<CScalarVoImageType, CScalarVoImageType, StructuringElementType>  OpeningImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                                                   ImageToShapeLabelMapFilterType;
    typedef itk::BinaryThresholdImageFilter<CScalarVoImageType, CScalarVoImageType>                                     ThresholdFilterType;
#if (ITK_VERSION_MAJOR < 4 || ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR < 4 || ITK_VERSION_MAJOR > 4)
    typedef itk::Compose2DVectorImageFilter<CScalarVoImageType, Vector2DVoImageType>                                    Compose2DVectorImageFilterType;
    typedef itk::Compose3DVectorImageFilter<CScalarVoImageType, Vector3DVoImageType>                                    Compose3DVectorImageFilterType;
#else
    typedef itk::ComposeImageFilter<CScalarVoImageType, Vector2DVoImageType>                                            Compose2DVectorImageFilterType;
    typedef itk::ComposeImageFilter<CScalarVoImageType, Vector3DVoImageType>                                            Compose3DVectorImageFilterType;
#endif
    typedef itk::DanielssonDistanceMapImageFilter<CScalarVoImageType, FScalarVoImageType>                               DanielssonDistanceMapImageFilterType;
    typedef itk::DivideImageFilter<FScalarVoImageType, FScalarVoImageType, FScalarVoImageType>                          DivideImageFilterType;
    typedef itk::MaximumImageFilter<CScalarVoImageType, CScalarVoImageType>                                             MaximumImageFilter;
    typedef itk::MorphologicalWatershedImageFilter<FScalarVoImageType, LScalarVoImageType>                              MorphoWatershedImageFilterType;
    typedef itk::MultiChannelConnectedThresholdImageFilter<Vector2DVoImageType, CScalarVoImageType>                     Connected2DThresholdFilterType;
    typedef itk::MultiChannelConnectedThresholdImageFilter<Vector3DVoImageType, CScalarVoImageType>                     Connected3DThresholdFilterType;
    typedef itk::NearestNeighborExtrapolateImageFunction<CScalarVoImageType, CoordRepType>                              NearestNeighborExtrapolatorImageFunctionType;
    typedef itk::NearestNeighborInterpolateImageFunction<CScalarVoImageType, CoordRepType>                              NearestNeighborInterpolatorImageFunctionType;
    typedef itk::LabelContourImageFilter<LScalarVoImageType, CScalarVoImageType>                                        LabelContourImageFilterType;
    typedef itk::LabelMapToLabelImageFilter<ImageToShapeLabelMapFilterType::OutputImageType, CScalarVoImageType>        LabelMapToLabelImageFilterType;
    typedef itk::LabelOverlayImageFilter<CScalarVoImageType, CScalarVoImageType, RGBVoImageType>                        LabelOverlayImageFilterType;
    typedef itk::ResampleImageFilter<CScalarVoImageType, CScalarVoImageType>                                            ResampleImageFilterType;
    typedef itk::RescaleIntensityImageFilter<FScalarVoImageType, CScalarVoImageType>                                    RescaleImageFilterType;
    typedef itk::ShapeOpeningLabelMapFilter<ImageToShapeLabelMapFilterType::OutputImageType>                            ShapeOpeningLabelMapFilterType;

    typedef itk::LabelImageToLabelMapFilter<MorphoWatershedImageFilterType::OutputImageType>                            LabelImageToLabelMapFilterType;
    typedef itk::LabelMapToLabelImageFilter<LabelImageToLabelMapFilterType::OutputImageType, CScalarVoImageType>        LabelMapToLabelImageFilterTypeW;

public:
    SegmentVeins();
    virtual ~SegmentVeins();

    void SetCVSeedPoints(std::vector<itk::Index<3> > seedPoints)
    {
        m_CVSeedPoints.clear();
        m_CVSeedPoints = seedPoints;
    };

    void SetPVSeedPoints(std::vector<itk::Index<3> > seedPoints)
    {
        m_PVSeedPoints.clear();
        m_PVSeedPoints = seedPoints;
    };

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteDataSetSummary();

	QString m_fullFilenameChannel0;
    QFileInfo m_infoFullFilenameChannel0;
    std::string m_pathChannel0;
    std::string m_filenameChannel0;
    std::string m_fileExtensionChannel0;

    QString m_fullFilenameChannel1;
    QFileInfo m_infoFullFilenameChannel1;
    std::string m_pathChannel1;
    std::string m_filenameChannel1;
    std::string m_fileExtensionChannel1;

    QString m_fullFilenameChannel2;
    QFileInfo m_infoFullFilenameChannel2;
    std::string m_pathChannel2;
    std::string m_filenameChannel2;
    std::string m_fileExtensionChannel2;

	bool m_hasCh0;

    std::vector<itk::Index<3> > m_CVSeedPoints;
    std::vector<itk::Index<3> > m_PVSeedPoints;

    int m_thresholdDAPI;
    int m_thresholdDPPIV;
    int m_thresholdDM;

    bool m_withResampling;
    bool m_withSegmentation;
    bool m_withLobulePrediction;

    StructuringElementType m_openingStrucElem;
    StructuringElementType m_closingStrucElem;

//    StructuringElementType m_maskDilateStrucElem;

    unsigned int m_minimalVeinSize;

    float m_overlayOpacity;
	
    std::string m_cvSaveBin;
    std::string m_cvSaveOverlay;
    std::string m_pvSaveBin;
    std::string m_pvSaveOverlay;
//    std::string m_cvMaskSaveBin;
//    std::string m_cvMaskSaveOverlay;
//    std::string m_pvMaskSaveBin;
//    std::string m_pvMaskSaveOverlay;
};

#endif /* SEGMENTVEINS_H_ */
