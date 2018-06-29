///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CropDataSet.h                                                        //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-06                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef CROPDATASET_H_
#define CROPDATASET_H_

#include "BasePipeline.h"

#include "itkExtractImageFilter.h"


template< unsigned int VImageDimension > class CropDataSet : public BasePipeline<VImageDimension>
{
public:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

private:
    typedef typename BasePipeline<VImageDimension>::CScalarVoImageType      CScalarImageType;
    typedef typename BasePipeline<VImageDimension>::CScalarRegionType       RegionType;
    typedef typename BasePipeline<VImageDimension>::CScalarIndexType        IndexType;
    typedef typename BasePipeline<VImageDimension>::CScalarSizeType         SizeType;

    typedef typename BasePipeline<VImageDimension>::ScalarVoReaderType      ScalarReaderType;
    typedef typename BasePipeline<VImageDimension>::ScalarVoWriterType      ScalarWriterType;

    typedef itk::ExtractImageFilter<CScalarImageType, CScalarImageType>     ExtractImageFilterType;

public:
    CropDataSet();
    virtual ~CropDataSet();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);

	QString m_fullFilename;
    QFileInfo m_infoFullFilename;
    std::string m_path;
    std::string m_filename;
    std::string m_fileExtension;

    IndexType m_start;
    SizeType m_size;
};

#include "CropDataSet.tpp"

#endif /* CROPDATASET_H_ */
