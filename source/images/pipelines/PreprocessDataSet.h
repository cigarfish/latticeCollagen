///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  PreprocessDataSet.h                                                  //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-15                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef PREPROCESSDATASET_H_
#define PREPROCESSDATASET_H_

#include "BasePipeline.h"

#include "itkAddImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkHConvexImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkWhiteTopHatImageFilter.h"

#include "../filters/imageFilters/FastCLAHEImageFilter.h"


template< unsigned int VImageDimension > class PreprocessDataSet : public BasePipeline<VImageDimension>
{
public:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

private:
    typedef typename BasePipeline<VImageDimension>::CScalarPixelType                                CScalarPixelType;

    typedef typename BasePipeline<VImageDimension>::CScalarVoImageType                              CScalarImageType;

    typedef typename BasePipeline<VImageDimension>::ScalarVoReaderType                              CScalarImageReaderType;
    typedef typename BasePipeline<VImageDimension>::ScalarVoWriterType                              CScalarImageWriterType;

    typedef itk::BinaryBallStructuringElement<CScalarPixelType, VImageDimension>                    StructuringElementType;

    typedef itk::AddImageFilter<CScalarImageType, CScalarImageType>                                 AddFilterType;
    typedef itk::BinaryThresholdImageFilter<CScalarImageType, CScalarImageType>                     ThresholdFilterType;
    typedef itk::FastCLAHEImageFilter<CScalarImageType>                                             CLAHEImageFilterType;
    typedef itk::HConvexImageFilter<CScalarImageType, CScalarImageType>                             ConvexImageFilterType;
    typedef itk::MedianImageFilter<CScalarImageType, CScalarImageType >                             MedianImageFilterType;
    typedef itk::WhiteTopHatImageFilter<CScalarImageType, CScalarImageType, StructuringElementType> WhiteTopHatFilterType;

public:
    enum AlgorithmType
    {
        CLAHE,
        BackgroundElimination,
        HConvexImageFilter,
        AddImages
    };

    PreprocessDataSet();
    virtual ~PreprocessDataSet();

    void SetAlgorithm(AlgorithmType a) { algorithm = a; };

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);


    AlgorithmType algorithm;

	QString m_fullFilename;
	QString m_fullFilename1;
	QString m_fullFilename2;
    QFileInfo m_infoFullFilename;
    std::string m_path;
    std::string m_filename;
    std::string m_fileExtension;

    unsigned int m_CLAHE_histWinSize;
    unsigned int m_CLAHE_stepSize;
    double m_CLAHE_clipLevel;

    int m_medianRadius;
    unsigned int m_backgroundElimination_KernelRadius;

    int m_convexFilter_heightLevel;
    bool m_convexFilter_fullyConnected;

    int m_addImageFilter_intensity1;
    int m_addImageFilter_intensity2;
};

#include "PreprocessDataSet.tpp"

#endif /* PREPROCESSDATASET_H_ */
