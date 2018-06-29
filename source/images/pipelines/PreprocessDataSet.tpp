///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  PreprocessDataSet.tpp                                                //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-06                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "PreprocessDataSet.h"

#include <QDateTime>

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"


template< unsigned int VImageDimension > PreprocessDataSet< VImageDimension >::PreprocessDataSet()
{
    // TODO Auto-generated constructor stub
    algorithm = CLAHE;
}


template< unsigned int VImageDimension > PreprocessDataSet< VImageDimension >::~PreprocessDataSet()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void PreprocessDataSet< VImageDimension >::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    this->m_paramContext->dump(parameter);

    file.open((m_path + "log_preprocessing.txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-Preprocess data set------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-preprocess data set--------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


template< unsigned int VImageDimension > void PreprocessDataSet< VImageDimension >::ParseParameterContext()
{
    if(algorithm == 0) {
        if(this->m_paramContext->findContext("CLAHE",0)==NULL) {
            std::cout << "Error: PreprocessDataSet: Invalid parameter context: CLAHE not found" << std::endl;
            return;
        }
        m_fullFilename = QString::fromStdString( *(std::string*)(this->m_paramContext->findParameter("Image filename", 0)->dataPointer()) );
		m_infoFullFilename.setFile(m_fullFilename);

		if(!m_infoFullFilename.exists())
		throw std::string("Please specify a file");

        m_CLAHE_histWinSize = *(unsigned int*)(this->m_paramContext->findParameter("Histogram Window Size", 0)->dataPointer());
        m_CLAHE_stepSize  = *(unsigned int*)(this->m_paramContext->findParameter("Step Size", 0)->dataPointer());
        m_CLAHE_clipLevel = *(double*)(this->m_paramContext->findParameter("Clip Level", 0)->dataPointer());
    }
    else if(algorithm == 1) {
        if(this->m_paramContext->findContext("Background Elimination",0)==NULL) {
            std::cout << "Error: PreprocessDataSet: Invalid parameter context: Background Elimination not found" << std::endl;
            return;
        }
        m_fullFilename = QString::fromStdString( *(std::string*)(this->m_paramContext->findParameter("Image filename", 0)->dataPointer()) );
		m_infoFullFilename.setFile(m_fullFilename);

		if(!m_infoFullFilename.exists())
		throw std::string("Please specify a file");

        m_medianRadius = ( *(int*)(this->m_paramContext->findParameter("Median kernel radius", 0)->dataPointer()) );
        m_backgroundElimination_KernelRadius = *(unsigned int*)(this->m_paramContext->findParameter("Background Elimination kernel radius", 0)->dataPointer());
    }
    else if(algorithm == 2) {
        if(this->m_paramContext->findContext("HConvex Image Filter",0)==NULL) {
            std::cout << "Error: PreprocessDataSet: Invalid parameter context: HConvex Image Filter not found" << std::endl;
            return;
        }
        m_fullFilename = QString::fromStdString( *(std::string*)(this->m_paramContext->findParameter("Image filename", 0)->dataPointer()) );
		m_infoFullFilename.setFile(m_fullFilename);

		if(!m_infoFullFilename.exists())
		throw std::string("Please specify a file");

        m_convexFilter_heightLevel = ( *(int*)(this->m_paramContext->findParameter("Height Level for identification of local maxima", 0)->dataPointer()) );
        m_convexFilter_fullyConnected = *(bool*)(this->m_paramContext->findParameter("Fully Connected", 0)->dataPointer());
    }
    else if(algorithm == 3) {
        if(this->m_paramContext->findContext("Add Images",0)==NULL) {
            std::cout << "Error: PreprocessDataSet: Invalid parameter context: Add Image Filter not found" << std::endl;
            return;
        }
        m_fullFilename1 = QString::fromStdString( *(std::string*)(this->m_paramContext->findParameter("Image filename 1", 0)->dataPointer()) );
        m_fullFilename2 = QString::fromStdString( *(std::string*)(this->m_paramContext->findParameter("Image filename 2", 0)->dataPointer()) );
        m_infoFullFilename.setFile(m_fullFilename2);
        if(!m_infoFullFilename.exists())
            throw std::string("Please specify a file");

        m_infoFullFilename.setFile(m_fullFilename1);
        if(!m_infoFullFilename.exists())
            throw std::string("Please specify a file");

        m_addImageFilter_intensity1 = ( *(int*)(this->m_paramContext->findParameter("Intensity 1", 0)->dataPointer()) );
        m_addImageFilter_intensity2 = ( *(int*)(this->m_paramContext->findParameter("Intensity 2", 0)->dataPointer()) );
    }

	m_path = (m_infoFullFilename.path() + QString("/")).toStdString();
    m_filename = m_infoFullFilename.baseName().toStdString();
    m_fileExtension = (QString(".") + m_infoFullFilename.suffix()).toStdString();

    if(BasePipeline<ImageDimension>::GetNumberOfDimensions(m_path + m_filename + m_fileExtension) != ImageDimension)
        throw std::string("Data has different dimensionality than expected by pipeline.");
}


template< unsigned int VImageDimension > void PreprocessDataSet< VImageDimension >::Update()
{
    ParseParameterContext();

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();


    std::cout << "Preprocess data set: " << std::endl;
    std::cout << " dir: " << m_path << std::endl;
    std::cout << " file: " << m_filename << std::endl;
    std::cout << " ext: " << m_fileExtension << std::endl;


    typename CScalarImageReaderType::Pointer reader = CScalarImageReaderType::New();
    reader->SetFileName(m_path + m_filename + m_fileExtension);
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif

    if(algorithm == 0)
    {
        typename CLAHEImageFilterType::Pointer CLAHEFilter = CLAHEImageFilterType::New();
        CLAHEFilter->SetInput(reader->GetOutput());
        CLAHEFilter->SetRadius(m_CLAHE_histWinSize);
        CLAHEFilter->SetStepSizeRadius(m_CLAHE_stepSize);
        CLAHEFilter->SetClipLimit(m_CLAHE_clipLevel);

        typename CScalarImageWriterType::Pointer writer = CScalarImageWriterType::New();
        writer->SetFileName(m_path + m_filename + "_clahe" + m_fileExtension);
        writer->SetInput(CLAHEFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer->Update();
    }
    else if(algorithm == 1)
    {
        typename MedianImageFilterType::InputSizeType medianRadius;
        medianRadius.Fill(m_medianRadius);

        typename MedianImageFilterType::Pointer medianFilter = MedianImageFilterType::New();
        medianFilter->SetRadius(medianRadius);
        medianFilter->SetInput(reader->GetOutput());

        StructuringElementType topHatKernel;
        topHatKernel.SetRadius(m_backgroundElimination_KernelRadius);
        topHatKernel.CreateStructuringElement();

        typename WhiteTopHatFilterType::Pointer whiteTopHatFilter = WhiteTopHatFilterType::New();
        whiteTopHatFilter->SetInput(medianFilter->GetOutput());
        whiteTopHatFilter->SetKernel(topHatKernel);
        whiteTopHatFilter->SetSafeBorder(true);

        typename CScalarImageWriterType::Pointer writer = CScalarImageWriterType::New();
        writer->ReleaseDataFlagOn();
        writer->SetFileName(m_path + m_filename + "_bckEli" + m_fileExtension);
        writer->SetInput(whiteTopHatFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer->Update();
    }
    else if(algorithm == 2)
    {
        typename ConvexImageFilterType::Pointer convexFilter = ConvexImageFilterType::New();
        convexFilter->SetInput(reader->GetOutput());
        convexFilter->SetHeight(m_convexFilter_heightLevel);
        convexFilter->SetFullyConnected(m_convexFilter_fullyConnected);

        typename CScalarImageWriterType::Pointer writer = CScalarImageWriterType::New();
        writer->ReleaseDataFlagOn();
        writer->SetFileName(m_path + m_filename + "_cnvx" + m_fileExtension);
        writer->SetInput(convexFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer->Update();
    }
    else if(algorithm == 3)
    {
        typename CScalarImageReaderType::Pointer reader2 = CScalarImageReaderType::New();
        reader2->SetFileName(m_fullFilename2.toStdString());
#if (ITK_VERSION_MAJOR >= 4)
        reader2->SetImageIO( itk::TIFFImageIO::New() );
#endif

        typename ThresholdFilterType::Pointer threshold1 = ThresholdFilterType::New();
        threshold1->ReleaseDataFlagOn();
        threshold1->SetOutsideValue(0);
        threshold1->SetInsideValue(m_addImageFilter_intensity1);
        threshold1->SetLowerThreshold(255);
        threshold1->SetUpperThreshold(255);
        threshold1->SetInput(reader->GetOutput());

        typename ThresholdFilterType::Pointer threshold2 = ThresholdFilterType::New();
        threshold2->ReleaseDataFlagOn();
        threshold2->SetOutsideValue(0);
        threshold2->SetInsideValue(m_addImageFilter_intensity2);
        threshold2->SetLowerThreshold(255);
        threshold2->SetUpperThreshold(255);
        threshold2->SetInput(reader2->GetOutput());

        typename AddFilterType::Pointer addFilter = AddFilterType::New();
        addFilter->SetInput1(threshold1->GetOutput());
        addFilter->SetInput2(threshold2->GetOutput());

        typename CScalarImageWriterType::Pointer writer = CScalarImageWriterType::New();
        writer->ReleaseDataFlagOn();
        writer->SetFileName(m_path + m_filename + "_comb" + m_fileExtension);
        writer->SetInput(addFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer->Update();
    }

    WriteLogFile(timeStamp);
}

