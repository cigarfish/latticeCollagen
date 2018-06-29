///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CropDataSet.cpp                                                      //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-06                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "CropDataSet.h"

#include <QDateTime>

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"


template< unsigned int VImageDimension > CropDataSet< VImageDimension >::CropDataSet()
{
    // TODO Auto-generated constructor stub

}


template< unsigned int VImageDimension > CropDataSet< VImageDimension >::~CropDataSet()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void CropDataSet< VImageDimension >::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    this->m_paramContext->dump(parameter);

    file.open((m_path + "log_preprocessing.txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-crop data set-------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-crop data set---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


template< unsigned int VImageDimension > void CropDataSet< VImageDimension >::ParseParameterContext()
{
    if(this->m_paramContext->findContext("Crop",0)==NULL) {
        std::cout << "Error: CropDataSet: Invalid parameter context: " << std::endl;
        return;
    }
    m_fullFilename = QString::fromStdString( *(std::string*)(this->m_paramContext->findParameter("Image filename", 0)->dataPointer()) );
	m_infoFullFilename.setFile(m_fullFilename);

	if(!m_infoFullFilename.exists())
		throw std::string("Please specify a file");

	m_path = (m_infoFullFilename.path() + QString("/")).toStdString();
    m_filename = m_infoFullFilename.baseName().toStdString();
    m_fileExtension = (QString(".") + m_infoFullFilename.suffix()).toStdString();

    if(BasePipeline<ImageDimension>::GetNumberOfDimensions(m_path + m_filename + m_fileExtension) != ImageDimension)
        throw std::string("Data has different dimensionality than expected by pipeline.");

    m_start[0] = *(int*)(this->m_paramContext->findParameter("x start", 0)->dataPointer());
    m_start[1] = *(int*)(this->m_paramContext->findParameter("y start", 0)->dataPointer());
    if(ImageDimension==3)   m_start[2] = *(int*)(this->m_paramContext->findParameter("z start", 0)->dataPointer());

    m_size[0] = *(int*)(this->m_paramContext->findParameter("x end", 0)->dataPointer()) - *(int*)(this->m_paramContext->findParameter("x start", 0)->dataPointer());
    m_size[1] = *(int*)(this->m_paramContext->findParameter("y end", 0)->dataPointer()) - *(int*)(this->m_paramContext->findParameter("y start", 0)->dataPointer());
    if(ImageDimension==3)   m_size[2] = *(int*)(this->m_paramContext->findParameter("z end", 0)->dataPointer()) - *(int*)(this->m_paramContext->findParameter("z start", 0)->dataPointer());
}


template< unsigned int VImageDimension > void CropDataSet< VImageDimension >::Update()
{
    ParseParameterContext();

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();


    std::cout << "Crop data set: " << std::endl;
    std::cout << " dir: " << m_path << std::endl;
    std::cout << " file: " << m_filename << std::endl;
    std::cout << " ext: " << m_fileExtension << std::endl;


    typename ScalarReaderType::Pointer          reader;
    typename ExtractImageFilterType::Pointer    extractFilter;
    typename ScalarWriterType::Pointer          writer;

    reader = ScalarReaderType::New();
    reader->SetFileName(m_path + m_filename + m_fileExtension);
#if (ITK_VERSION_MAJOR >= 4)
    reader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    reader->Update();

    RegionType desiredRegion(m_start, m_size);

    extractFilter = ExtractImageFilterType::New();
    extractFilter->SetExtractionRegion(desiredRegion);
    extractFilter->SetInput(reader->GetOutput());

    writer = ScalarWriterType::New();
    writer->ReleaseDataFlagOn();
    writer->SetFileName(m_path + m_filename + "_cut" + m_fileExtension);
    writer->SetInput(extractFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();

    WriteLogFile(timeStamp);
}
