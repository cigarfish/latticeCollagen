///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeStellateCellsFilter.cpp                                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-01-30                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "AnalyzeStellateCellsFilter.h"

#include <math.h>
#include <fstream>
#include <string>
#include <sstream>

#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>

#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIOFactory.h>
#include <itkTIFFImageIO.h>
#endif // (ITK_VERSION_MAJOR >= 4)

#include <vtkDelimitedTextReader.h>
#include <vtkDelimitedTextWriter.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkVariantArray.h>

#include "../../../tools/parameters/CSParameter.h"
#include "../../../tools/parameters/CSParameterContext.h"



AnalyzeStellateCellsFilter::AnalyzeStellateCellsFilter()
{
    // TODO Auto-generated constructor stub
    m_cv_volume = 0;
    m_pv_volume = 0;

    m_normVeinDistRingStepSize = 50;
    m_numRings = 300 / m_normVeinDistRingStepSize;

    for(unsigned int i=m_normVeinDistRingStepSize; i<=300; i+=m_normVeinDistRingStepSize) {
        m_normCVDistRings[i]  = 0;
        m_normPVDistRings[i]  = 0;
    }

#if (ITK_VERSION_MAJOR >= 4)
    itk::TIFFImageIOFactory::RegisterOneFactory();
#endif // (ITK_VERSION_MAJOR >= 4)
}


AnalyzeStellateCellsFilter::~AnalyzeStellateCellsFilter()
{
}


void AnalyzeStellateCellsFilter::CollectBasicImageInformation()
{
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO((m_dataSetPathHepNuc + m_dataSetNameHepNuc + m_dataSetFileExtensionHepNuc).c_str(), itk::ImageIOFactory::ReadMode);
    if(imageIO->CanReadFile((m_dataSetPathHepNuc + m_dataSetNameHepNuc + m_dataSetFileExtensionHepNuc).c_str())) {
        imageIO->ReadImageInformation();

        m_num_dim = imageIO->GetNumberOfDimensions();
        m_dim = new long int[m_num_dim];
        m_voxel_volume = 1;
        m_dataset_volume = 1;

        for(unsigned int i=0; i<m_num_dim; i++) {
            m_dim[i] = imageIO->GetDimensions(i);
            //TODO: m_spacing[i] = imageIO->GetSpacing(i);
            m_voxel_volume *= m_spacing[i];
            m_dataset_volume *= m_dim[i]*m_spacing[i];
        }
        m_voxel_volume /= 1000000000;
        m_dataset_volume /= 1000000000;                                 //in mm³ (=1000³)

        m_effective_dataset_volume = m_dataset_volume;
    }
}


void AnalyzeStellateCellsFilter::CollectBasicNucleiInformation(bool withCV, bool withPV)
{
    //----------READER---------------------------------------------------------------------------------------------------------
    ScalarVoReaderType::Pointer hepNucReader = ScalarVoReaderType::New();
    hepNucReader->SetFileName(m_dataSetPathHepNuc + m_dataSetNameHepNuc + m_dataSetFileExtensionHepNuc);
    hepNucReader->ReleaseDataBeforeUpdateFlagOn();
    hepNucReader->Update();

    CScalarVoImageType::Pointer hepNucImage = hepNucReader->GetOutput();
    hepNucImage->DisconnectPipeline();
    hepNucImage->SetSpacing(m_spacing);
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER---CREATE-LABEL-MAP--------------------------------------------------------------------------------------
    ImageToShapeLabelMapFilterType::Pointer hepNucImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    hepNucImageToShaLabMapFilter->SetInput(hepNucImage);
    hepNucImageToShaLabMapFilter->Update();

    /*for(unsigned int i=0; i<hepNucImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
        std::cout << "roundness of label object " << i << ": " << hepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetRoundness() << std::endl;
        std::cout << "spherical radius of label object " << i << ": " << hepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetEquivalentSphericalRadius() << std::endl;
        std::cout << "physical size of label object " << i << ": " << hepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetPhysicalSize() << std::endl;
        std::cout << "number pixel of label object " << i << ": " << hepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixels() << std::endl;
    }*/

    m_num_hepaticNuclei = hepNucImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
    m_num_hepaticNucleiPerEffVol = (double)m_num_hepaticNuclei / m_effective_dataset_volume;


    //----------READER---------------------------------------------------------------------------------------------------------
    ScalarVoReaderType::Pointer nonHepNucReader = ScalarVoReaderType::New();
    nonHepNucReader->SetFileName(m_dataSetPathNonHepNuc + m_dataSetNameNonHepNuc + m_dataSetFileExtensionNonHepNuc);
    nonHepNucReader->ReleaseDataBeforeUpdateFlagOn();
    nonHepNucReader->Update();

    CScalarVoImageType::Pointer nonHepNucImage = nonHepNucReader->GetOutput();
    nonHepNucImage->DisconnectPipeline();
    nonHepNucImage->SetSpacing(m_spacing);
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER---CREATE-LABEL-MAP--------------------------------------------------------------------------------------
    ImageToShapeLabelMapFilterType::Pointer nonHepNucImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    nonHepNucImageToShaLabMapFilter->SetInput(nonHepNucImage);
    nonHepNucImageToShaLabMapFilter->Update();

    /*for(unsigned int i=0; i<nonHepNucImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
        std::cout << "roundness of label object " << i << ": " << nonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetRoundness() << std::endl;
        std::cout << "spherical radius of label object " << i << ": " << nonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetEquivalentSphericalRadius() << std::endl;
        std::cout << "physical size of label object " << i << ": " << nonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetPhysicalSize() << std::endl;
        std::cout << "number pixel of label object " << i << ": " << nonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixels() << std::endl;
    }*/

    m_num_nonHepaticNuclei = nonHepNucImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
    m_num_nonHepaticNucleiPerEffVol = (double)m_num_nonHepaticNuclei / m_effective_dataset_volume;


    //----------READER---------------------------------------------------------------------------------------------------------
    ScalarVoReaderType::Pointer stcHepNucReader = ScalarVoReaderType::New();
    stcHepNucReader->SetFileName(m_dataSetPathStCNuc + m_dataSetNameStCNuc + m_dataSetFileExtensionStCNuc);
    stcHepNucReader->ReleaseDataBeforeUpdateFlagOn();
    stcHepNucReader->Update();

    CScalarVoImageType::Pointer stcHepNucImage = stcHepNucReader->GetOutput();
    stcHepNucImage->DisconnectPipeline();
    stcHepNucImage->SetSpacing(m_spacing);
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER---CREATE-LABEL-MAP--------------------------------------------------------------------------------------
    ImageToShapeLabelMapFilterType::Pointer stcHepNucImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    stcHepNucImageToShaLabMapFilter->SetInput(stcHepNucImage);
    stcHepNucImageToShaLabMapFilter->Update();

    m_num_stellateCellNuclei = stcHepNucImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
    m_num_stellateCellNucleiPerEffVol = (double)m_num_stellateCellNuclei / m_effective_dataset_volume;
    //-------------------------------------------------------------------------------------------------------------------------

    //----------READER---------------------------------------------------------------------------------------------------------
    ScalarVoReaderType::Pointer CVReader = ScalarVoReaderType::New();
    CScalarVoImageType::Pointer CVImage = CVReader->GetOutput();
    SignedMaurerDistanceMapImageFilterType::Pointer distanceMapCV = SignedMaurerDistanceMapImageFilterType::New();

    if(withCV && m_cv_volume!=0) {
        CVReader->SetFileName(m_dataSetPathCV + m_dataSetNameCV + m_dataSetFileExtensionCV);
        CVReader->ReleaseDataBeforeUpdateFlagOn();
        CVReader->Update();

        CVImage->DisconnectPipeline();
        CVImage->SetSpacing(m_spacing);

        distanceMapCV->ReleaseDataFlagOn();
        distanceMapCV->UseImageSpacingOn();
        distanceMapCV->SquaredDistanceOff();
        distanceMapCV->SetInput(CVImage);
        distanceMapCV->Update();

        ImageIteratorType it(distanceMapCV->GetOutput(), distanceMapCV->GetOutput()->GetLargestPossibleRegion());

        int i=0;
        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
            unsigned int bin = m_normVeinDistRingStepSize*((int)(it.Value() / m_normVeinDistRingStepSize) + 1);
            if(i<100) { std::cout << "Pixel with Dist = " <<  it.Value() << " sorted to bin = " << bin << std::endl; i++; }

            m_normCVDistRings[bin] += 1;
        }
        std::map<unsigned int, double>::iterator iter;
        for(iter = m_normCVDistRings.begin(); iter!=m_normCVDistRings.end(); iter++) {
            iter->second = iter->second * m_voxel_volume;
            std::cout << "Ring to " << iter->first << " has volume = " << iter->second << std::endl;
        }
    }

    ScalarVoReaderType::Pointer PVReader = ScalarVoReaderType::New();
    CScalarVoImageType::Pointer PVImage = PVReader->GetOutput();
    SignedMaurerDistanceMapImageFilterType::Pointer distanceMapPV = SignedMaurerDistanceMapImageFilterType::New();

    if(withPV && m_pv_volume!=0) {
        PVReader->SetFileName(m_dataSetPathPV + m_dataSetNamePV + m_dataSetFileExtensionPV);
        PVReader->ReleaseDataBeforeUpdateFlagOn();
        PVReader->Update();

        PVImage->DisconnectPipeline();
        PVImage->SetSpacing(m_spacing);

        distanceMapPV->ReleaseDataFlagOn();
        distanceMapPV->UseImageSpacingOn();
        distanceMapPV->SquaredDistanceOff();
        distanceMapPV->SetInput(PVImage);
        distanceMapPV->Update();

        ImageIteratorType it(distanceMapPV->GetOutput(), distanceMapPV->GetOutput()->GetLargestPossibleRegion());

        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
            unsigned int bin = m_normVeinDistRingStepSize*((int)(it.Value() / m_normVeinDistRingStepSize) + 1);

            m_normPVDistRings[bin] += 1;
        }
        std::map<unsigned int, double>::iterator iter;
        for(iter = m_normPVDistRings.begin(); iter!=m_normPVDistRings.end(); iter++) {
            iter->second = iter->second * m_voxel_volume;
        }
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------COMPUTE-STELLATE-CELL-TO-VEIN-DISTANCES------------------------------------------------------------------------
    for(unsigned int i=0; i<stcHepNucImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
        StellateCell cell;

        cell.mLabel = stcHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetLabel();
        if(withCV && m_cv_volume!=0) cell.distToCV = distanceMapCV->GetOutput()->GetPixel( stcHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetIndex(0) );
        if(withPV && m_pv_volume!=0) cell.distToPV = distanceMapPV->GetOutput()->GetPixel( stcHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetIndex(0) );

        m_stellateCells.insert(std::pair<unsigned long int, StellateCell>(cell.mLabel, cell));
    }
    //-------------------------------------------------------------------------------------------------------------------------
}


void AnalyzeStellateCellsFilter::CollectBasicStellateCellInformation()
{
    //----------READER---------------------------------------------------------------------------------------------------------
    ScalarVoReaderType::Pointer stcBodiesReader = ScalarVoReaderType::New();
    stcBodiesReader->SetFileName(m_dataSetPathStCBody + m_dataSetNameStCBody + m_dataSetFileExtensionStCBody);
    stcBodiesReader->ReleaseDataBeforeUpdateFlagOn();
    stcBodiesReader->Update();

    CScalarVoImageType::Pointer stcBodiesImage = stcBodiesReader->GetOutput();
    stcBodiesImage->DisconnectPipeline();
    stcBodiesImage->SetSpacing(m_spacing);
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER---CREATE-LABEL-MAP--------------------------------------------------------------------------------------
    ImageToShapeLabelMapFilterType::Pointer stcBodiesImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    stcBodiesImageToShaLabMapFilter->SetInput(stcBodiesImage);
    stcBodiesImageToShaLabMapFilter->Update();

    m_volume_stellateCellBodies = 0;
    for(unsigned int i=0; i<stcBodiesImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++)
        m_volume_stellateCellBodies += stcBodiesImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetPhysicalSize();
    m_volume_stellateCellBodies /= 1000000000;

    m_num_stellateCellBodies = stcBodiesImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
    m_num_stellateCellBodiesPerEffVol = (double)m_num_stellateCellBodies / m_effective_dataset_volume;

    m_volume_stellateCellBodiesPerEffVol = m_volume_stellateCellBodies / m_effective_dataset_volume;
}


void AnalyzeStellateCellsFilter::CollectBasicHepatocyteInformation()
{
    ScalarVoReaderType::Pointer cellReader = ScalarVoReaderType::New();
    cellReader->SetFileName(m_dataSetPathHepCell + m_dataSetNameHepCell + m_dataSetFileExtensionHepCell);
#if (ITK_VERSION_MAJOR >= 4)
    cellReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    cellReader->Update();

    CScalarVoImageType::Pointer cellImage = cellReader->GetOutput();
    cellImage->DisconnectPipeline();
    cellImage->SetSpacing(m_spacing);

    ImageToShapeLabelMapFilterType::Pointer cellImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    cellImageToShaLabMapFilter->SetInput(cellImage);
    cellImageToShaLabMapFilter->SetFullyConnected(false);
    cellImageToShaLabMapFilter->Update();

    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap = cellImageToShaLabMapFilter->GetOutput();
    cellLabelMap->DisconnectPipeline();

    m_num_hepatocytes = cellLabelMap->GetNumberOfLabelObjects();
    m_num_hepatocytesPerEffVol = (double)m_num_hepatocytes / m_effective_dataset_volume;

    for(unsigned int i=0; i<cellLabelMap->GetNumberOfLabelObjects(); ++i)
        m_hepatocyteVolume.push_back( cellLabelMap->GetNthLabelObject(i)->GetNumberOfPixels() * m_voxel_volume );
}


long double AnalyzeStellateCellsFilter::ComputeVolumeOfRegions(std::string path, std::string filename, std::string ext)
{
    long double volume = 0;

    ifstream file((path + filename + ext).c_str());
    if(file.good())
    {
        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO((path + filename + ext).c_str(), itk::ImageIOFactory::ReadMode);
        imageIO->ReadImageInformation();

        ScalarVoReaderType::Pointer                     reader;
        ImageToShapeLabelMapFilterType::Pointer         imageToShaLabMapFilter;

        //----------READER---------------------------------------------------------------------------------------------------------
        reader = ScalarVoReaderType::New();
        reader->SetFileName(path + filename + ext);
        reader->ReleaseDataBeforeUpdateFlagOn();
        //-------------------------------------------------------------------------------------------------------------------------

        //----------FILTER---CREATE-LABEL-MAP--------------------------------------------------------------------------------------
        imageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
        imageToShaLabMapFilter->SetInput(reader->GetOutput());
        imageToShaLabMapFilter->Update();

        int numPixel = 0;

        int numLabObjects = imageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
        for(unsigned int i=0; i<numLabObjects; ++i) {
            imageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->Optimize();
            numPixel += imageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixels();
        }

        volume = numPixel * m_voxel_volume;
        //-------------------------------------------------------------------------------------------------------------------------
    }
    return volume;
}


void AnalyzeStellateCellsFilter::WriteBasicInformationFile()
{
    std::fstream file1, file2, tempfile;
    bool hasHeader = false;

    file1.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_1.txt").c_str(), fstream::in);

    std::string line;
    getline(file1, line);

    if(line.find("dataSet")!=std::string::npos) {
        hasHeader = true;

        tempfile.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_1_tempfile.txt").c_str(), fstream::out);
        tempfile << line;

        while(!file1.eof()) {
            getline(file1, line);
            if(line.find(m_dataSetID)!=0) {
                tempfile << std::endl;
                tempfile << line;
            }
        }
        file1.close();
        tempfile.close();

        std::remove((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_1.txt").c_str());
        std::rename((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_1_tempfile.txt").c_str(), (m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_1.txt").c_str());
    }
    else
        file1.close();

    if(!hasHeader)
        file2.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_1.txt").c_str(), fstream::out);
    else
        file2.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_1.txt").c_str(), fstream::out | fstream::app);

    file2.flags(fstream::left | fstream::scientific);
    if(!hasHeader) {
        file2.width(40);
        file2 << "dataSet";
        file2.width(15);
        file2 << "dim";
        file2.width(15);
        file2 << "dimX";
        file2.width(15);
        file2 << "dimY";
        file2.width(15);
        file2 << "dimZ";
        file2.width(20);
        file2 << "dataSetVolume";
        file2.width(20);
        file2 << "veinVolume";
        file2.width(20);
        file2 << "effectiveVolume";
        file2.width(20);
        file2 << "numberHepaticNuclei";
        file2.width(35);
        file2 << "numberNonHepaticNuclei";
        file2.width(35);
        file2 << "numberStellateCellNuclei";
        file2.width(35);
        file2 << "numberStellateCellBodies";
        file2.width(35);
        file2 << "numberHepatocytes";
        file2.width(35);
        file2 << "volumeStellateCellBodies";
        file2.width(35);
        file2 << "volumeStellateCellBodies/effVol";
        file2.width(35);
        file2 << "numberHepaticNuclei/effVol";
        file2.width(35);
        file2 << "numberNonHepaticNuclei/effVol";
        file2.width(35);
        file2 << "numberStellateCellNuclei/effVol";
        file2.width(35);
        file2 << "numberStellateCellBodies/effVol";
        file2.width(35);
        file2 << "numberHepatocytes/effVol";
        for(unsigned int i=0; i<m_numRings; i++) {
            std::stringstream s;
            s <<  "cvRingVol_" << i;

            file2.width(20);
            file2 << s.str();
        }
        for(unsigned int i=0; i<m_numRings; i++) {
            std::stringstream s;
            s <<  "pvRingVol_" << i;

            file2.width(20);
            file2 << s.str();
        }
        file2 << std::endl;
    }

    file2.width(40);
    file2 << m_dataSetID;
    file2.width(15);
    file2 << m_num_dim;
    file2.width(15);
    file2 << m_dim[0];
    file2.width(15);
    file2 << m_dim[1];
    file2.width(15);
    file2 << m_dim[2];
    file2.width(20);
    file2 << (double)m_dataset_volume;
    file2.width(20);
    file2 << (double)(m_cv_volume + m_pv_volume);
    file2.width(20);
    file2 << (double)m_effective_dataset_volume;
    file2.width(20);
    file2 << m_num_hepaticNuclei;
    file2.width(35);
    file2 << m_num_nonHepaticNuclei;
    file2.width(35);
    file2 << m_num_stellateCellNuclei;
    file2.width(35);
    file2 << m_num_stellateCellBodies;
    file2.width(35);
    file2 << m_num_hepatocytes;
    file2.width(35);
    file2 << m_volume_stellateCellBodies;
    file2.width(35);
    file2 << m_volume_stellateCellBodiesPerEffVol;
    file2.width(35);
    file2 << m_num_hepaticNucleiPerEffVol;
    file2.width(35);
    file2 << m_num_nonHepaticNucleiPerEffVol;
    file2.width(35);
    file2 << m_num_stellateCellNucleiPerEffVol;
    file2.width(35);
    file2 << m_num_stellateCellBodiesPerEffVol;
    file2.width(35);
    file2 << m_num_hepatocytesPerEffVol;
    std::map<unsigned int, double>::iterator iter;
    for(iter = m_normCVDistRings.begin(); iter!=m_normCVDistRings.end(); iter++) {
        file2.width(20);
        file2 << iter->second;
    }
    for(iter = m_normPVDistRings.begin(); iter!=m_normPVDistRings.end(); iter++) {
        file2.width(20);
        file2 << iter->second;
    }
    file2 << std::endl;

    file2.close();
}


void AnalyzeStellateCellsFilter::WriteBasicNucleiInformationFile()
{
    std::fstream file1, file2, tempfile;
    bool hasHeader = false;

    file1.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_2.txt").c_str(), fstream::in);

    std::string line;
    getline(file1, line);

    if(line.find("dataSet")!=std::string::npos) {
        hasHeader = true;

        tempfile.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_2_tempfile.txt").c_str(), fstream::out);
        tempfile << line;

        while(!file1.eof()) {
            getline(file1, line);
            if(line.find(m_dataSetID)!=0) {
                tempfile << std::endl;
                tempfile << line;
            }
        }
        file1.close();
        tempfile.close();

        std::remove((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_2.txt").c_str());
        std::rename((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_2_tempfile.txt").c_str(), (m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_2.txt").c_str());
    }
    else
        file1.close();

    if(!hasHeader)
        file2.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_2.txt").c_str(), fstream::out);
    else
        file2.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_2.txt").c_str(), fstream::out | fstream::app);

    file2.flags(fstream::left | fstream::scientific);
    if(!hasHeader) {
        file2.width(40);
        file2 << "dataSet";
        file2.width(15);
        file2 << "label";
        file2.width(15);
        file2 << "distToCV";
        file2.width(15);
        file2 << "distToPV" << std::endl;
    }

    for(std::map<unsigned long int, StellateCell>::iterator cellIt = m_stellateCells.begin(); cellIt != m_stellateCells.end(); ++cellIt) {
        file2.width(40);
        file2 << m_dataSetID;
        file2.width(15);
        file2 << cellIt->second.mLabel;
        file2.width(15);
        file2 << cellIt->second.distToCV;
        file2.width(15);
        file2 << cellIt->second.distToPV << std::endl;
    }

    file2.close();
}


void AnalyzeStellateCellsFilter::WriteBasicHepatocyteInformationFile()
{
    std::fstream file1, file2, tempfile;
    bool hasHeader = false;

    file1.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_3.txt").c_str(), fstream::in);

    std::string line;
    getline(file1, line);

    if(line.find("dataSet")!=std::string::npos) {
        hasHeader = true;

        tempfile.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_3_tempfile.txt").c_str(), fstream::out);
        tempfile << line;

        while(!file1.eof()) {
            getline(file1, line);
            if(line.find(m_dataSetID)!=0) {
                tempfile << std::endl;
                tempfile << line;
            }
        }
        file1.close();
        tempfile.close();

        std::remove((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_3.txt").c_str());
        std::rename((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_3_tempfile.txt").c_str(), (m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_3.txt").c_str());
    }
    else
        file1.close();

    if(!hasHeader)
        file2.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_3.txt").c_str(), fstream::out);
    else
        file2.open((m_dataSetPathHepNuc + "../" + "StellateCell_Output_file_3.txt").c_str(), fstream::out | fstream::app);

    file2.flags(fstream::left | fstream::scientific);
    if(!hasHeader) {
        file2.width(40);
        file2 << "dataSet";
        file2.width(30);
        file2 << "volume" << std::endl;
    }

    for(unsigned int i = 0; i<m_hepatocyteVolume.size(); ++i) {
        file2.width(40);
        file2 << m_dataSetID;
        file2.width(30);
        file2 << m_hepatocyteVolume[i] << std::endl;
    }

    file2.close();
}


void AnalyzeStellateCellsFilter::ParseParameterContext()
{
    if(m_paramContext->findContext("Analyze Stellate Cells",0)==NULL) {
        std::cout << "Error: AnalyzeStellateCells: Invalid parameter context" << std::endl;
        return;
    }

    m_dataSetID = *(std::string*)(m_paramContext->findParameter("Dataset ID", 0)->dataPointer());

    m_dataSetFullFilenameHepNuc = *(std::string*)(m_paramContext->findParameter("Hepatic nuclei segmentation", 0)->dataPointer());
    m_dataSetFullFilenameNonHepNuc = *(std::string*)(m_paramContext->findParameter("Non-Hepatic nuclei segmentation", 0)->dataPointer());
    m_dataSetFullFilenameStCNuc = *(std::string*)(m_paramContext->findParameter("Stellate nuclei segmentation", 0)->dataPointer());
    m_dataSetFullFilenameStCBody = *(std::string*)(m_paramContext->findParameter("Stellate cell body segmentation", 0)->dataPointer());
    m_dataSetFullFilenameHepCell = *(std::string*)(m_paramContext->findParameter("Hepatocyte segmentation", 0)->dataPointer());
    m_dataSetFullFilenameCV = *(std::string*)(m_paramContext->findParameter("Central Vein segmentation", 0)->dataPointer());
    m_dataSetFullFilenamePV = *(std::string*)(m_paramContext->findParameter("Portal Vein segmentation", 0)->dataPointer());

    m_spacing[0] = *(double*)(m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    m_spacing[1] = *(double*)(m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    m_spacing[2] = *(double*)(m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());
}


void AnalyzeStellateCellsFilter::Update()
{
    ParseParameterContext();

    bool f0 = FilenameParser::ParseFilename(m_dataSetFullFilenameHepNuc, m_dataSetPathHepNuc, m_dataSetNameHepNuc, m_dataSetFileExtensionHepNuc);
    bool f1 = FilenameParser::ParseFilename(m_dataSetFullFilenameNonHepNuc, m_dataSetPathNonHepNuc, m_dataSetNameNonHepNuc, m_dataSetFileExtensionNonHepNuc);
    bool f2 = FilenameParser::ParseFilename(m_dataSetFullFilenameStCNuc, m_dataSetPathStCNuc, m_dataSetNameStCNuc, m_dataSetFileExtensionStCNuc);
    bool f3 = FilenameParser::ParseFilename(m_dataSetFullFilenameStCBody, m_dataSetPathStCBody, m_dataSetNameStCBody, m_dataSetFileExtensionStCBody);
    bool f4 = FilenameParser::ParseFilename(m_dataSetFullFilenameHepCell, m_dataSetPathHepCell, m_dataSetNameHepCell, m_dataSetFileExtensionHepCell);
    bool f5 = FilenameParser::ParseFilename(m_dataSetFullFilenameCV, m_dataSetPathCV, m_dataSetNameCV, m_dataSetFileExtensionCV);
    bool f6 = FilenameParser::ParseFilename(m_dataSetFullFilenamePV, m_dataSetPathPV, m_dataSetNamePV, m_dataSetFileExtensionPV);

    if(!f0 || !f1 || !f2 || !f3 || !f4)
        return;

    CollectBasicImageInformation();

    if(f5) {
        m_cv_volume = ComputeVolumeOfRegions(m_dataSetPathCV, m_dataSetNameCV, m_dataSetFileExtensionCV);
        m_effective_dataset_volume -= m_cv_volume;
    }
    if(f6) {
        m_pv_volume = ComputeVolumeOfRegions(m_dataSetPathPV, m_dataSetNamePV, m_dataSetFileExtensionPV);
        m_effective_dataset_volume -= m_pv_volume;
    }

    CollectBasicNucleiInformation(f5, f6);
    CollectBasicStellateCellInformation();
    CollectBasicHepatocyteInformation();


    WriteBasicInformationFile();
    WriteBasicNucleiInformationFile();
    WriteBasicHepatocyteInformationFile();
}
