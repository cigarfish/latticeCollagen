///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  EstimateCellShape.cpp                                                //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-06-04                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "EstimateCellShape.h"

#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIO.h>
#include <itkNrrdImageIO.h>
#endif

#include <QDateTime>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"



EstimateCellShape::EstimateCellShape()
{
    m_overlayOpacity = 0.5;

    m_saveSuffixesForFinals[0] = "_step3_bin";
    m_saveSuffixesForFinals[1] = "_step3_gray";
    m_saveSuffixesForFinals[2] = "_step3_dppiv_cellOverlay";
    m_saveSuffixesForFinals[3] = "_step3_dppiv_borderOverlay";
    m_saveSuffixesForFinals[4] = "_step3_dppiv_monoNucOverlay";
    m_saveSuffixesForFinals[5] = "_step3_dppiv_biNucOverlay";
    m_saveSuffixesForFinals[6] = "_step3_dppiv_multiNucOverlay";
}


EstimateCellShape::~EstimateCellShape()
{
    // TODO Auto-generated destructor stub
}


void EstimateCellShape::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

    file.open((m_path + m_logFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-cell-shape-estimation-------------------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-cell-shape-estimation---------------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


void EstimateCellShape::WriteDataSetSummary()
{
    ImageAnalysisSummaryFileIO::AddEntry(CellShapeBin, m_path, m_path + m_filenameSave + m_saveSuffixesForFinals[0] + m_filenameExtension);
    ImageAnalysisSummaryFileIO::AddEntry(CellShapeOverlay, m_path, m_path + m_filenameSave + m_saveSuffixesForFinals[2] + m_filenameExtension);
}


void EstimateCellShape::ParseParameterContext()
{
    if(m_paramContext->findContext("Approximate Cell Shape",0)==NULL) {
        std::cout << "Error: EstimateCellShape: Invalid parameter context" << std::endl;
        return;
    }

    m_filenameDPPIVChannel = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("DPPIV channel", 0)->dataPointer()) );
    m_filenameNucleiBin = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Hepatic nuclei segmentation", 0)->dataPointer()) );
    m_filenameBileBin = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Bile segmentation", 0)->dataPointer()) );
    m_filenameSinusoidBin = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Sinusoid segmentation", 0)->dataPointer()) );
    m_filenameCVBin = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Central vein segmentation", 0)->dataPointer()) );
    m_filenamePVBin = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Portal vein segmentation", 0)->dataPointer()) );

	m_infoFile.setFile(m_filenameDPPIVChannel);
	if(!m_infoFile.exists())
		throw std::string("Please specify DPPIV channel");

	m_infoFile.setFile(m_filenameNucleiBin);
	if(!m_infoFile.exists())
		throw std::string("Please specify nuclei segmentation file");

	m_path = (m_infoFile.path() + QString("/")).toStdString();
    m_filenameExtension = (QString(".") + m_infoFile.suffix()).toStdString();

	m_infoFile.setFile(m_filenameBileBin);
	if(!m_infoFile.exists())
		throw std::string("Please specify bile segmentation file");

	m_infoFile.setFile(m_filenameSinusoidBin);
	if(!m_infoFile.exists())
		throw std::string("Please specify sinusoid segmentation file");

	m_infoFile.setFile(m_filenameCVBin);
	m_hasCV = m_infoFile.exists();

	m_infoFile.setFile(m_filenamePVBin);
	m_hasPV = m_infoFile.exists();

    m_hasNR = *(bool*)(m_paramContext->findParameter("Is there a necrotic region", 0)->dataPointer());
	m_fullFilenameNecroticRegion = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Necrotic region", 0)->dataPointer()) );

	m_infoFile.setFile(m_fullFilenameNecroticRegion);
	if(m_hasNR && !m_infoFile.exists())
		throw std::string("Please specify Necrotic Region binary mask");

    m_spacing[0] = *(double*)(m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    m_spacing[1] = *(double*)(m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    m_spacing[2] = *(double*)(m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    m_bileSinWeight = *(double*)(m_paramContext->findParameter("Bile weight [0,1]", 0)->dataPointer());
    m_nucleiWeight = 1. - m_bileSinWeight;

    m_minimalCellDiameter = *(double*)(m_paramContext->findParameter("Minimal cell diameter", 0)->dataPointer());
    m_minimalCellVolume = 1./6. * itk::Math::pi * pow(m_minimalCellDiameter, 3);

    m_watershedFloodLevel = *(double*)(m_paramContext->findParameter("Alpha", 0)->dataPointer());

    m_withPositionFilter = *(bool*)(m_paramContext->findParameter("With position filter", 0)->dataPointer());
    m_numberCells = *(int*)(m_paramContext->findParameter("Number of cells to keep", 0)->dataPointer());
    m_pos[0] = (*(double*)(m_paramContext->findParameter("Coordinate x", 0)->dataPointer()))*m_spacing[0];
    m_pos[1] = (*(double*)(m_paramContext->findParameter("Coordinate y", 0)->dataPointer()))*m_spacing[1];
    m_pos[2] = (*(double*)(m_paramContext->findParameter("Coordinate z", 0)->dataPointer()))*m_spacing[2];

    std::string saveMode = ( (CSParameterChoice*)(m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        m_saveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        m_saveEverything = 0;
    m_logFilenameSave = *(std::string*)(m_paramContext->findParameter("Log file name", 0)->dataPointer());
    m_filenameSave = *(std::string*)(m_paramContext->findParameter("Save prefix", 0)->dataPointer());
    m_writeCellOutlineToFile = *(bool*)(m_paramContext->findParameter("Write cell outline to file", 0)->dataPointer());
    m_writeSinusoidOutlineToFile = *(bool*)(m_paramContext->findParameter("Write sinusoid outline to file", 0)->dataPointer());
#ifndef CS_TI_QUANT_ONLY
    m_writeSinusoidGraphToFile = *(bool*)(m_paramContext->findParameter("Write sinusoid graph to file", 0)->dataPointer());
    m_filenameSinusoidSkeleton = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Sinusoid skeleton image", 0)->dataPointer()) );
#else
	m_writeSinusoidGraphToFile = false;
	m_filenameSinusoidSkeleton = QString("");
#endif
}


void EstimateCellShape::SetupBorderRegions(CScalarVoImageType::Pointer &image)
{
    itk::Index<3> indices[6];
    indices[0][0] = 0;  indices[1][0] = 0;
    indices[0][1] = 0;  indices[1][1] = 0;
    indices[0][2] = 0;  indices[1][2] = image->GetLargestPossibleRegion().GetSize()[2]-1;

    indices[2][0] = 0;  indices[3][0] = 0;
    indices[2][1] = 0;  indices[3][1] = image->GetLargestPossibleRegion().GetSize()[1]-1;
    indices[2][2] = 1;  indices[3][2] = 1;

    indices[4][0] = 0;  indices[5][0] = image->GetLargestPossibleRegion().GetSize()[0]-1;
    indices[4][1] = 1;  indices[5][1] = 1;
    indices[4][2] = 1;  indices[5][2] = 1;

    itk::Size<3> sizes[6];
    sizes[0][0] = image->GetLargestPossibleRegion().GetSize()[0];   sizes[1][0] = image->GetLargestPossibleRegion().GetSize()[0];
    sizes[0][1] = image->GetLargestPossibleRegion().GetSize()[1];   sizes[1][1] = image->GetLargestPossibleRegion().GetSize()[1];
    sizes[0][2] = 1;                                                sizes[1][2] = 1;

    sizes[2][0] = image->GetLargestPossibleRegion().GetSize()[0];   sizes[3][0] = image->GetLargestPossibleRegion().GetSize()[0];
    sizes[2][1] = 1;                                                sizes[3][1] = 1;
    sizes[2][2] = image->GetLargestPossibleRegion().GetSize()[2]-2; sizes[3][2] = image->GetLargestPossibleRegion().GetSize()[2]-2;

    sizes[4][0] = 1;                                                sizes[5][0] = 1;
    sizes[4][1] = image->GetLargestPossibleRegion().GetSize()[1]-2; sizes[5][1] = image->GetLargestPossibleRegion().GetSize()[1]-2;
    sizes[4][2] = image->GetLargestPossibleRegion().GetSize()[2]-2; sizes[5][2] = image->GetLargestPossibleRegion().GetSize()[2]-2;


    mDataSetFaces[0].SetIndex(indices[0]);   mDataSetFaces[0].SetSize(sizes[0]);
    mDataSetFaces[1].SetIndex(indices[1]);   mDataSetFaces[1].SetSize(sizes[1]);
    mDataSetFaces[2].SetIndex(indices[2]);   mDataSetFaces[2].SetSize(sizes[2]);
    mDataSetFaces[3].SetIndex(indices[3]);   mDataSetFaces[3].SetSize(sizes[3]);
    mDataSetFaces[4].SetIndex(indices[4]);   mDataSetFaces[4].SetSize(sizes[4]);
    mDataSetFaces[5].SetIndex(indices[5]);   mDataSetFaces[5].SetSize(sizes[5]);
}


void EstimateCellShape::SetBorderPixel(CScalarVoImageType::Pointer &image, CScalarPixelType value)
{
    for(unsigned int i=0; i<6; i++) {
        itk::ImageRegionIterator<CScalarVoImageType> iter(image, mDataSetFaces[i]);

        while(!iter.IsAtEnd()) {
            iter.Set(value);
            ++iter;
        }
    }
}


void EstimateCellShape::SetBorderPixelAtLocalMaxima(CScalarVoImageType::Pointer &image, FScalarVoImageType::Pointer distMap, CScalarPixelType value)
{
    RegionalMaximaImageFilterType::Pointer maximaFilter = RegionalMaximaImageFilterType::New();
    maximaFilter->SetInput(distMap);

    CScalarVoImageType::Pointer distMaxImage = CScalarVoImageType::New();
    distMaxImage->SetRegions(image->GetLargestPossibleRegion());
    distMaxImage->Allocate();
    distMaxImage->FillBuffer(0);
    distMaxImage->SetSpacing(m_spacing);

    PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New();
    for(unsigned int i=0; i<6; i++) {
        maximaFilter->GetOutput()->SetRequestedRegion(mDataSetFaces[i]);

        pasteFilter->SetSourceImage(maximaFilter->GetOutput());
        pasteFilter->SetSourceRegion(maximaFilter->GetOutput()->GetRequestedRegion());
        pasteFilter->SetDestinationImage(distMaxImage);
        pasteFilter->SetDestinationIndex(mDataSetFaces[i].GetIndex());
        pasteFilter->Update();

        distMaxImage = pasteFilter->GetOutput();
        distMaxImage->DisconnectPipeline();
    }

    AddCImageFilterType::Pointer addFilter = AddCImageFilterType::New();
    addFilter->SetInput1(distMaxImage);
    addFilter->SetInput2(image);
    addFilter->Update();

    image = addFilter->GetOutput();
    image->DisconnectPipeline();
}


void EstimateCellShape::SetBorderPixelAtLocalMaximaExceedingThreshold(CScalarVoImageType::Pointer &image, FScalarVoImageType::Pointer distMap, FScalarPixelType threshold, CScalarPixelType value)
{
    for(unsigned int i=0; i<6; i++) {
        itk::ImageRegionIterator<FScalarVoImageType> iter(distMap, mDataSetFaces[i]);

        while(!iter.IsAtEnd()) {
            if(iter.Get()>threshold)
                image->SetPixel(iter.GetIndex(), value);
            ++iter;
        }
    }
}


//TODO: move this into nuclei class as optional hep/non-hep classification helper
void EstimateCellShape::RemoveNucleiInContactWithSinusoids(CScalarVoImageType::Pointer &nucImage, CScalarVoImageType::Pointer sinImage)
{
    ImageToShapeLabelMapFilterType::Pointer nucLabelMapFilter = ImageToShapeLabelMapFilterType::New();
    nucLabelMapFilter->SetInput(nucImage);
    nucLabelMapFilter->SetFullyConnected(false);
    nucLabelMapFilter->Update();

    ImageToShapeLabelMapFilterType::OutputImageType::Pointer nucLabelMap = nucLabelMapFilter->GetOutput();
    nucLabelMap->DisconnectPipeline();

    LabelMapToLabelImageFilterType::Pointer nucLabelMapToLabelImage = LabelMapToLabelImageFilterType::New();
    nucLabelMapToLabelImage->SetInput(nucLabelMap);
    nucLabelMapToLabelImage->Update();

    LScalarVoImageType::Pointer nucLabelImage = nucLabelMapToLabelImage->GetOutput();
    nucLabelImage->DisconnectPipeline();

    MaskImageFilterType::Pointer maskFilter = MaskImageFilterType::New();
    maskFilter->SetMaskImage(sinImage);
    maskFilter->SetInput(nucLabelImage);
    maskFilter->SetOutsideValue(0);
    maskFilter->Update();

    LScalarVoImageType::Pointer maskedNucLabelImage = maskFilter->GetOutput();
    maskedNucLabelImage->DisconnectPipeline();

    RescaleLIImageFilterType::Pointer rescaler = RescaleLIImageFilterType::New();
    rescaler->SetInput(maskedNucLabelImage);

    IScalarVoWriterType::Pointer maskWriter = IScalarVoWriterType::New();
    maskWriter->SetFileName(m_path + m_filenameSave + "_step0_maskSinNucOverlap.nrrd");
    maskWriter->SetInput(rescaler->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    maskWriter->SetImageIO( itk::NrrdImageIO::New() );
#endif
    maskWriter->Update();

    std::map<unsigned long, float> nucleiOverlap;
    itk::ImageRegionConstIterator<LScalarVoImageType> iter(maskedNucLabelImage, maskedNucLabelImage->GetLargestPossibleRegion());

    while(!iter.IsAtEnd()) {
        if(iter.Value() != 0) {
            if(nucleiOverlap.count(iter.Value()) == 0)
                nucleiOverlap.insert( std::pair<unsigned long, int>(iter.Value(), 1) );
            else
                nucleiOverlap[iter.Value()]+=1;
        }
        ++iter;
    }

    for(std::map<unsigned long, float>::iterator it = nucleiOverlap.begin(); it != nucleiOverlap.end(); ++it)
        nucleiOverlap[it->first] = (float)it->second / (float)nucLabelMap->GetLabelObject(it->first)->Size();

    for(std::map<unsigned long, float>::iterator it = nucleiOverlap.begin(); it != nucleiOverlap.end(); ++it)
        if(it->second>0.1)
            nucLabelMap->RemoveLabel(it->first);

    nucLabelMapToLabelImage->SetInput(nucLabelMap);
    nucLabelMapToLabelImage->Update();

    ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
    thresholdFilter->SetInput(nucLabelMapToLabelImage->GetOutput());
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->SetInsideValue(itk::NumericTraits<CScalarPixelType>::max());
    thresholdFilter->SetLowerThreshold(1);
    thresholdFilter->SetUpperThreshold(itk::NumericTraits<LabelMapToLabelImageFilterType::OutputImagePixelType>::max());
    thresholdFilter->Update();

    nucImage = thresholdFilter->GetOutput();
    nucImage->DisconnectPipeline();

    ScalarVoWriterType::Pointer writerTest2 = ScalarVoWriterType::New();
    writerTest2->SetFileName(m_path + m_filenameSave + "_step0_hepsOnly.tif");
    writerTest2->SetInput(nucImage);
#if (ITK_VERSION_MAJOR >= 4)
    writerTest2->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerTest2->Update();
}


void EstimateCellShape::RemoveCellsAtDatasetBorder(LScalarVoImageType::Pointer cellLabelImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithOneNucleusLabelMap,
        itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithTwoNucleiLabelMap, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithMoreNucleiLabelMap,
        double allowedBorderPixel)
{
    std::map<unsigned long, int> cellsToNumberBorderPixel;

    for(unsigned int i=0; i<6; i++) {
        itk::ImageRegionConstIterator<LScalarVoImageType> iter(cellLabelImage, mDataSetFaces[i]);

        while(!iter.IsAtEnd()) {
            if(iter.Value() != 0) {
                if(cellsToNumberBorderPixel.count(iter.Value()) == 0)
                    cellsToNumberBorderPixel.insert( std::pair<unsigned long, int>(iter.Value(), 1) );
                else
                    cellsToNumberBorderPixel[iter.Value()]++;
            }
            ++iter;
        }
    }

    for(std::map<unsigned long, int>::iterator it = cellsToNumberBorderPixel.begin(); it != cellsToNumberBorderPixel.end(); ++it) {
        if(it->second>allowedBorderPixel) {
            cellsWithOneNucleusLabelMap->RemoveLabel(it->first);
            cellsWithTwoNucleiLabelMap->RemoveLabel(it->first);
            cellsWithMoreNucleiLabelMap->RemoveLabel(it->first);
        }
    }
}


void EstimateCellShape::BuildCellNucleiAlignmentMaps(LScalarVoImageType::Pointer nucleiLabelImage, LScalarVoImageType::Pointer cellLabelImage,
        itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> nucleiLabelMap, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap)
{
    std::cout << "nucleiLabelImage size " << nucleiLabelImage->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "cellLabelImage size " << cellLabelImage->GetLargestPossibleRegion().GetSize() << std::endl;

    itk::ImageRegionConstIterator<LScalarVoImageType> iterNucl(nucleiLabelImage, nucleiLabelImage->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<LScalarVoImageType> iterCell(cellLabelImage, cellLabelImage->GetLargestPossibleRegion());
    for(iterNucl = iterNucl.Begin(), iterCell = iterCell.Begin(); !iterNucl.IsAtEnd(); ++iterNucl, ++iterCell) {
        if(iterNucl.Value() != 0 && iterCell.Value() != 0) {
            if(m_cellToNuclei.count(iterCell.Value()) == 0) {
                std::pair<unsigned long int, double> nucleus(iterNucl.Value(), 1.);
                m_cellToNuclei.insert(std::pair<unsigned long int, std::pair<unsigned long int, double> >(iterCell.Value(), nucleus));
            }
            else {
                std::multimap<unsigned long int, std::pair<unsigned long int, double> >::iterator alignIter;
                alignIter = m_cellToNuclei.find(iterCell.Value());

                bool alreadyFound = false;
                while(alignIter != m_cellToNuclei.end() && alignIter->first == iterCell.Value()) {
                    if(alignIter->second.first == iterNucl.Value()) {
                        alignIter->second.second += 1.;
                        alreadyFound = true;
                        break;
                    }
                    ++alignIter;
                }
                if(!alreadyFound) {
                    std::pair<unsigned long int, double> nucleus(iterNucl.Value(), 1.);
                    m_cellToNuclei.insert(std::pair<unsigned long int, std::pair<unsigned long int, double> >(iterCell.Value(), nucleus));
                }
            }
        }
    }

    std::cout << "Start cell evaluation " << std::endl;
    std::multimap<unsigned long int, std::pair<unsigned long int, double> >::iterator alignIter = m_cellToNuclei.begin();
    while(alignIter != m_cellToNuclei.end()) {
        alignIter->second.second = alignIter->second.second / (double)(nucleiLabelMap->GetLabelObject(alignIter->second.first)->GetNumberOfPixels());

        std::cout << "Cell " << alignIter->first << " which starts at " << cellLabelMap->GetLabelObject(alignIter->first)->GetIndex(0) << " contains nucleus " << alignIter->second.first <<
                " which starts at " << nucleiLabelMap->GetLabelObject(alignIter->second.first)->GetIndex(0) << " to factor " << alignIter->second.second << std::endl;

        if(alignIter->second.second <= 0.5) {
            std::cout << "This cell-nuclei pair is removed, due to minor overlapping." << std::endl;
            m_cellToNuclei.erase(alignIter++);
        }
        else
            ++alignIter;
    }

    for(unsigned int i=0; i<cellLabelMap->GetNumberOfLabelObjects(); i++) {
        if( m_cellToNuclei.count( cellLabelMap->GetNthLabelObject(i)->GetLabel() ) == 0 ) {
            std::cout << "cell " << cellLabelMap->GetNthLabelObject(i)->GetLabel() << " has no nucleus!" << std::endl;
            m_labelsWithZeroNuclei.push_back(cellLabelMap->GetNthLabelObject(i)->GetLabel());
        }
        else if( m_cellToNuclei.count( cellLabelMap->GetNthLabelObject(i)->GetLabel() ) == 1 ) {
            std::cout << "cell " << cellLabelMap->GetNthLabelObject(i)->GetLabel() << " has one nucleus!" << std::endl;
            m_labelsWithOneNucleus.push_back(cellLabelMap->GetNthLabelObject(i)->GetLabel());
        }
        else if( m_cellToNuclei.count( cellLabelMap->GetNthLabelObject(i)->GetLabel() ) == 2 ) {
            std::cout << "cell " << cellLabelMap->GetNthLabelObject(i)->GetLabel() << " has two nuclei!" << std::endl;
            m_labelsWithTwoNuclei.push_back(cellLabelMap->GetNthLabelObject(i)->GetLabel());
        }
        else if( m_cellToNuclei.count( cellLabelMap->GetNthLabelObject(i)->GetLabel() ) > 2 ) {
            std::cout << "cell " << cellLabelMap->GetNthLabelObject(i)->GetLabel() << " has more than two nuclei!" << std::endl;
            m_labelsWithMoreNuclei.push_back(cellLabelMap->GetNthLabelObject(i)->GetLabel());
        }
    }
}


void EstimateCellShape::WriteNumNucleiOutlineFiles(itk::SmartPointer<CScalarVoImageType> dppivImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithOneNucleusLabelMap,
        itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithTwoNucleiLabelMap, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithMoreNucleiLabelMap)
{
    LabelOverlayImageFilterType::Pointer    overlayImage3, overlayImage4, overlayImage5;
    RGBVoWriterType::Pointer                writerRGB3, writerRGB4, writerRGB5;

    for(unsigned int i=0; i<m_labelsWithZeroNuclei.size(); i++)
        cellsWithOneNucleusLabelMap->RemoveLabel(m_labelsWithZeroNuclei[i]);
    for(unsigned int i=0; i<m_labelsWithTwoNuclei.size(); i++)
        cellsWithOneNucleusLabelMap->RemoveLabel(m_labelsWithTwoNuclei[i]);
    for(unsigned int i=0; i<m_labelsWithMoreNuclei.size(); i++)
        cellsWithOneNucleusLabelMap->RemoveLabel(m_labelsWithMoreNuclei[i]);
    std::cout << "after removal of cells with other than one nucleus " << cellsWithOneNucleusLabelMap->GetNumberOfLabelObjects() << " cells left" << std::endl;

    for(unsigned int i=0; i<m_labelsWithZeroNuclei.size(); i++)
        cellsWithTwoNucleiLabelMap->RemoveLabel(m_labelsWithZeroNuclei[i]);
    for(unsigned int i=0; i<m_labelsWithOneNucleus.size(); i++)
        cellsWithTwoNucleiLabelMap->RemoveLabel(m_labelsWithOneNucleus[i]);
    for(unsigned int i=0; i<m_labelsWithMoreNuclei.size(); i++)
        cellsWithTwoNucleiLabelMap->RemoveLabel(m_labelsWithMoreNuclei[i]);
    std::cout << "after removal of cells with other than two nuclei " << cellsWithTwoNucleiLabelMap->GetNumberOfLabelObjects() << " cells left" << std::endl;

    for(unsigned int i=0; i<m_labelsWithZeroNuclei.size(); i++)
        cellsWithMoreNucleiLabelMap->RemoveLabel(m_labelsWithZeroNuclei[i]);
    for(unsigned int i=0; i<m_labelsWithOneNucleus.size(); i++)
        cellsWithMoreNucleiLabelMap->RemoveLabel(m_labelsWithOneNucleus[i]);
    for(unsigned int i=0; i<m_labelsWithTwoNuclei.size(); i++)
        cellsWithMoreNucleiLabelMap->RemoveLabel(m_labelsWithTwoNuclei[i]);
    std::cout << "after removal of cells with less than two nuclei " << cellsWithMoreNucleiLabelMap->GetNumberOfLabelObjects() << " cells left" << std::endl;

    LabelMapToLabelImageFilterType::Pointer labelMapToImageA = LabelMapToLabelImageFilterType::New();
    labelMapToImageA->SetInput(cellsWithOneNucleusLabelMap);

    LabelMapToLabelImageFilterType::Pointer labelMapToImageB = LabelMapToLabelImageFilterType::New();
    labelMapToImageB->SetInput(cellsWithTwoNucleiLabelMap);

    LabelMapToLabelImageFilterType::Pointer labelMapToImageC = LabelMapToLabelImageFilterType::New();
    labelMapToImageC->SetInput(cellsWithMoreNucleiLabelMap);

    overlayImage3 = LabelOverlayImageFilterType::New();
    overlayImage3->ReleaseDataFlagOn();
    overlayImage3->SetOpacity(m_overlayOpacity);
    overlayImage3->SetInput(dppivImage);
    overlayImage3->SetLabelImage(labelMapToImageA->GetOutput());

    overlayImage4 = LabelOverlayImageFilterType::New();
    overlayImage4->ReleaseDataFlagOn();
    overlayImage4->SetOpacity(m_overlayOpacity);
    overlayImage4->SetInput(dppivImage);
    overlayImage4->SetLabelImage(labelMapToImageB->GetOutput());

    overlayImage5 = LabelOverlayImageFilterType::New();
    overlayImage5->ReleaseDataFlagOn();
    overlayImage5->SetOpacity(m_overlayOpacity);
    overlayImage5->SetInput(dppivImage);
    overlayImage5->SetLabelImage(labelMapToImageC->GetOutput());

    writerRGB3 = RGBVoWriterType::New();
    writerRGB3->SetFileName(m_path + m_filenameSave + m_saveSuffixesForFinals[4] + m_filenameExtension);
    writerRGB3->SetInput(overlayImage3->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writerRGB3->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerRGB3->Update();

    writerRGB4 = RGBVoWriterType::New();
    writerRGB4->SetFileName(m_path + m_filenameSave + m_saveSuffixesForFinals[5] + m_filenameExtension);
    writerRGB4->SetInput(overlayImage4->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writerRGB4->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerRGB4->Update();

    writerRGB5 = RGBVoWriterType::New();
    writerRGB5->SetFileName(m_path + m_filenameSave + m_saveSuffixesForFinals[6] + m_filenameExtension);
    writerRGB5->SetInput(overlayImage5->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writerRGB5->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerRGB5->Update();
}


void EstimateCellShape::WriteSinusoidGraphToFile()
{
    BScalarVoReaderType::Pointer                skeletonReader;
    ScalarVoReaderType::Pointer                 readerSin;
    ImageToShapeLabelMapFilterType::Pointer     sinImageToShaLabMapFilter;
    ShapeOpeningLabelMapFilterType2::Pointer    sinShapeOpeningLabMapFilter;

    std::string dirSkl, fileSkl, extSkl;
    m_infoFile.setFile(m_filenameSinusoidSkeleton);

	if(!m_infoFile.exists())
		throw std::string("Please specify sinusoid skeleton");

	dirSkl = (m_infoFile.path() + QString("/")).toStdString();
    fileSkl = m_infoFile.baseName().toStdString();
    extSkl = (QString(".") + m_infoFile.suffix()).toStdString();

    //----------READER---------------------------------------------------------------------------------------------------------
    skeletonReader = BScalarVoReaderType::New();
    skeletonReader->SetFileName(m_path + fileSkl + extSkl);
#if (ITK_VERSION_MAJOR >= 4)
    skeletonReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    skeletonReader->Update();

    readerSin = ScalarVoReaderType::New();
	readerSin->SetFileName(m_filenameSinusoidBin.toStdString());
    readerSin->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerSin->SetImageIO( itk::TIFFImageIO::New() );
#endif
    readerSin->Update();

    //----------FILTER---CREATE-LABEL-MAP--------------------------------------------------------------------------------------
    sinImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    sinImageToShaLabMapFilter->SetInput(readerSin->GetOutput());
    sinImageToShaLabMapFilter->SetFullyConnected(true);
    sinImageToShaLabMapFilter->Update();

    std::cout << "before removal " << sinImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects() << " sinusoid objects" << std::endl;

    sinShapeOpeningLabMapFilter = ShapeOpeningLabelMapFilterType2::New();
    sinShapeOpeningLabMapFilter->SetLambda(100);                                           //attribute value
    sinShapeOpeningLabMapFilter->ReverseOrderingOff();                                     //removes objects with attribute smaller than lambda
    sinShapeOpeningLabMapFilter->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::NUMBER_OF_PIXELS);
    sinShapeOpeningLabMapFilter->SetInput(sinImageToShaLabMapFilter->GetOutput());
    sinShapeOpeningLabMapFilter->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    std::vector< vtkSmartPointer<vtkUndirectedGraph> > graphs = SkeletonImageToGraphFilter::skeletonToGraph3D(skeletonReader->GetOutput(), true);

    GraphAnnotationHelper anno;
    anno.EnableEdgeLengthAnnotation(1, 1, 1);

    for(unsigned int i=0; i<graphs.size(); i++)
        anno.AddPredefinedAnnotations(graphs[i]);

    std::vector< vtkSmartPointer<vtkUndirectedGraph> > resampledGraphs;

    for(unsigned int i=0; i<graphs.size(); i++) {
        vtkSmartPointer<ResampleUndirectedGraphFilter> resampleGraph = vtkSmartPointer<ResampleUndirectedGraphFilter>::New();
        resampleGraph->SetInput(graphs[i]);
        resampleGraph->SetResamplingFactor(9);
        resampleGraph->SetMaxResamplingDist(40);
        resampleGraph->SetTestOutput(false);
        resampleGraph->Update();

        vtkSmartPointer<vtkUndirectedGraph> resampledUndirectedGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
        resampleGraph->GetOutput()->ToUndirectedGraph(resampledUndirectedGraph);
        resampledGraphs.push_back(resampledUndirectedGraph);
    }

    std::stringstream nameOfFirstOne;
    nameOfFirstOne << m_path << m_filenameSave << "_supportGraph" << 0 << ".txt";
    std::string nameOfFirstGraph = nameOfFirstOne.str();

    for(unsigned int i=0; i<resampledGraphs.size(); i++) {
        vtkSmartPointer<vtkGraphWriter> writer = vtkSmartPointer<vtkGraphWriter>::New();

        std::stringstream name;
        name << m_path << m_filenameSave << "_supportGraph" << i << ".txt";
        std::cout << "save final graph to " << name.str() << std::endl;

        writer->SetFileName(name.str().c_str());
        writer->SetInput(resampledGraphs[i]);
        writer->Update();
    }

    std::fstream file2;
    file2.open((m_path + m_filenameSave + "_skeletonSinusoids.txt").c_str(), std::fstream::out);

    for(unsigned int i=0; i<resampledGraphs.size(); i++) {
        vtkSmartPointer<vtkVertexListIterator> vertexListIterator = vtkSmartPointer<vtkVertexListIterator>::New();
        resampledGraphs[i]->GetVertices(vertexListIterator);

        if(vertexListIterator->HasNext()) {
            double *a;
            a = resampledGraphs[i]->GetPoint(vertexListIterator->Next());
            itk::Index<3> idx;
            idx[0] = a[0]; idx[1] = a[1]; idx[2] = a[2];

            bool found = false;
            for(unsigned int j=0; j<sinShapeOpeningLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); j++){
                if( sinShapeOpeningLabMapFilter->GetOutput()->GetNthLabelObject(j)->HasIndex(idx) ) {
                    file2 << j << ": " << std::endl;
                    found = true;
                    break;
                }
            }
            if(!found)
                file2 << "not found: " << std::endl;
            file2 << a[0] << ", " << a[1] << ", " << a[2] << "; ";
        }

        while(vertexListIterator->HasNext())
        {
            double *a;
            a = resampledGraphs[i]->GetPoint(vertexListIterator->Next());
            file2 << a[0] << ", " << a[1] << ", " << a[2] << "; ";
        }
        file2 << std::endl;
    }
    file2.close();
}


void EstimateCellShape::WriteSiunsoidOutlineToFile()
{
    ScalarVoReaderType::Pointer                     readerSin;
    ImageToShapeLabelMapFilterType::Pointer         sinImageToShaLabMapFilter;
    ShapeOpeningLabelMapFilterType2::Pointer        sinShapeOpeningLabMapFilter;
    LabelMapToLabelImageFilterType2::Pointer        sinLabelMapToImage;
    LabelContourImageFilterType2::Pointer           sinLabelContour;
    LabelImageToShapeLabelMapFilterType2::Pointer   sinLabelImageToShapeLabelMap;

    readerSin = ScalarVoReaderType::New();
	readerSin->SetFileName(m_filenameSinusoidBin.toStdString());
    readerSin->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerSin->SetImageIO( itk::TIFFImageIO::New() );
#endif
    readerSin->Update();

    //----------FILTER---CREATE-LABEL-MAP--------------------------------------------------------------------------------------
    sinImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    sinImageToShaLabMapFilter->SetInput(readerSin->GetOutput());
    sinImageToShaLabMapFilter->SetFullyConnected(true);
    sinImageToShaLabMapFilter->Update();

    std::cout << "before removal " << sinImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects() << " sinusoid objects" << std::endl;

    sinShapeOpeningLabMapFilter = ShapeOpeningLabelMapFilterType2::New();
    sinShapeOpeningLabMapFilter->SetLambda(100);                                           //attribute value
    sinShapeOpeningLabMapFilter->ReverseOrderingOff();                                     //removes objects with attribute smaller than lambda
    sinShapeOpeningLabMapFilter->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::NUMBER_OF_PIXELS);
    sinShapeOpeningLabMapFilter->SetInput(sinImageToShaLabMapFilter->GetOutput());
    sinShapeOpeningLabMapFilter->Update();

    std::cout << "after removal " << sinShapeOpeningLabMapFilter->GetOutput()->GetNumberOfLabelObjects() << " sinusoid objects" << std::endl;

    std::vector< std::vector<ShapeOpeningLabelMapFilterType::IndexType> > sinBorderVoxel;
    ShapeOpeningLabelMapFilterType::IndexType idx;

    for(unsigned int i=0; i<sinShapeOpeningLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
        int numVoxel = sinShapeOpeningLabMapFilter->GetOutput()->GetNthLabelObject(i)->Size();
        std::vector<ShapeOpeningLabelMapFilterType::IndexType> bV;

        for(unsigned int j=0; j<numVoxel; j++) {
            idx = sinShapeOpeningLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetIndex(j);

            if(idx[0] == 0 || idx[1] == 0 || idx[2] == 0 || idx[0] == m_XMax || idx[1] == m_YMax || idx[2] == m_ZMax)
                bV.push_back(idx);
        }
        sinBorderVoxel.push_back(bV);
    }

    sinLabelMapToImage = LabelMapToLabelImageFilterType2::New();
    sinLabelMapToImage->SetInput(sinShapeOpeningLabMapFilter->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------CONTOUR-OF-SINUSOIDS------------------------------------------------------------------------------------------
    sinLabelContour = LabelContourImageFilterType2::New();
    sinLabelContour->SetFullyConnected(true);
    sinLabelContour->SetInput(sinLabelMapToImage->GetOutput());

    sinLabelImageToShapeLabelMap = LabelImageToShapeLabelMapFilterType2::New();
    sinLabelImageToShapeLabelMap->SetBackgroundValue(0);
    sinLabelImageToShapeLabelMap->SetInput(sinLabelContour->GetOutput());
    sinLabelImageToShapeLabelMap->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    std::fstream file1;

    file1.open((m_path + m_filenameSave + "_outlineSinusoids.txt").c_str(), std::fstream::out);

    for(unsigned int i=0; i<sinLabelImageToShapeLabelMap->GetOutput()->GetNumberOfLabelObjects(); i++) {
        file1 << i << ": " << std::endl;

        int numVoxel = sinLabelImageToShapeLabelMap->GetOutput()->GetNthLabelObject(i)->Size();

        for(unsigned int j=0; j<numVoxel; j++) {
            idx = sinLabelImageToShapeLabelMap->GetOutput()->GetNthLabelObject(i)->GetIndex(j);

            file1 << idx[0] << ", " << idx[1] << ", " << idx[2] << "; ";
        }
        for(unsigned int j=0; j<sinBorderVoxel[i].size(); j++)
            file1 << sinBorderVoxel[i][j][0] << ", " << sinBorderVoxel[i][j][1] << ", " << sinBorderVoxel[i][j][2] << "; ";

        file1 << std::endl;
    }
    file1.close();
}


void EstimateCellShape::Update()
{
    ParseParameterContext();

    std::string dirDPPIV, fileDPPIV, extDPPIV, dirNuc, fileNuc, extNuc, dirBile, fileBile, extBile, dirSin, fileSin, extSin, dirCV, fileCV, extCV, dirPV, filePV, extPV, dirNecReg, fileNecReg, extNecReg;

	m_infoFile.setFile(m_filenameDPPIVChannel);
	dirDPPIV = (m_infoFile.path() + QString("/")).toStdString();
    fileDPPIV = m_infoFile.baseName().toStdString();
    extDPPIV = (QString(".") + m_infoFile.suffix()).toStdString();

	m_infoFile.setFile(m_filenameNucleiBin);
	dirNuc = (m_infoFile.path() + QString("/")).toStdString();
    fileNuc = m_infoFile.baseName().toStdString();
    extNuc = (QString(".") + m_infoFile.suffix()).toStdString();

	m_infoFile.setFile(m_filenameBileBin);
	dirBile = (m_infoFile.path() + QString("/")).toStdString();
    fileBile = m_infoFile.baseName().toStdString();
    extBile = (QString(".") + m_infoFile.suffix()).toStdString();

	m_infoFile.setFile(m_filenameSinusoidBin);
	dirSin = (m_infoFile.path() + QString("/")).toStdString();
    fileSin = m_infoFile.baseName().toStdString();
    extSin = (QString(".") + m_infoFile.suffix()).toStdString();

	if(m_hasCV) {
		m_infoFile.setFile(m_filenameCVBin);
		dirCV = (m_infoFile.path() + QString("/")).toStdString();
		fileCV = m_infoFile.baseName().toStdString();
		extCV = (QString(".") + m_infoFile.suffix()).toStdString();
	}

	if(m_hasPV) {
		m_infoFile.setFile(m_filenamePVBin);
		dirPV = (m_infoFile.path() + QString("/")).toStdString();
		filePV = m_infoFile.baseName().toStdString();
		extPV = (QString(".") + m_infoFile.suffix()).toStdString();
	}

	if(m_hasNR) {
		m_infoFile.setFile(m_fullFilenameNecroticRegion);
		dirNecReg = (m_infoFile.path() + QString("/")).toStdString();
		fileNecReg = m_infoFile.baseName().toStdString();
		extNecReg = (QString(".") + m_infoFile.suffix()).toStdString();
	}

    m_path = dirNuc;
    m_filenameExtension = extNuc;

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    std::cout << "Segment estimate cell shape: " << std::endl;
    std::cout << " dirDPPIV: " << dirDPPIV << std::endl;
    std::cout << " fileDPPIV: " << fileDPPIV << std::endl;
    std::cout << " extDPPIV: " << extDPPIV << std::endl;

    std::cout << " dirNuc: " << dirNuc << std::endl;
    std::cout << " fileNuc: " << fileNuc << std::endl;
    std::cout << " extNuc: " << extNuc << std::endl;

    std::cout << " dirBile: " << dirBile << std::endl;
    std::cout << " fileBile: " << fileBile << std::endl;
    std::cout << " extBile: " << extBile << std::endl;

    std::cout << " dirSin: " << dirSin << std::endl;
    std::cout << " fileSin: " << fileSin << std::endl;
    std::cout << " extSin: " << extSin << std::endl;

    std::cout << " dirCV: " << dirCV << std::endl;
    std::cout << " fileCV: " << fileCV << std::endl;
    std::cout << " extCV: " << extCV << std::endl;

    std::cout << " dirPV: " << dirPV << std::endl;
    std::cout << " filePV: " << filePV << std::endl;
    std::cout << " extPV: " << extPV << std::endl;

    ScalarVoReaderType::Pointer                     readerNuc;
    ScalarVoReaderType::Pointer                     readerBile;
    ScalarVoReaderType::Pointer                     readerSin;
    ScalarVoReaderType::Pointer                     readerCV;
    ScalarVoReaderType::Pointer                     readerPV;
    ScalarVoReaderType::Pointer                     reader1, reader2;
    ScalarVoReaderType::Pointer                     readerDPPIV;
    ScalarVoReaderType::Pointer                     readerNecReg;
    SignedMaurerDistanceMapImageFilterType::Pointer distanceMapNuc, distanceMapBileSin;
    ImageCalculatorFilterType::Pointer              distanceMapNucCalc, distanceMapBileSinCalc;
    IntensityWindowingImageFilter::Pointer          distanceMapNucCapped, distanceMapBileSinCapped;
    AddCImageFilterType::Pointer                    addVeins;
    AddCImageFilterType::Pointer                    addVeinsSin;
    AddCImageFilterType::Pointer                    addVeinsSinBile;
    AddCImageFilterType::Pointer                    addVeinsSinNecReg;
    MultiplyImageFilterType::Pointer                multNuc;
    MultiplyImageFilterType::Pointer                multBS;
    AddFImageFilterType::Pointer                    add2;
    DivideImageFilterType::Pointer                  divide;
    RescaleImageFilterType::Pointer                 rescaler1, rescaler2, rescaler3, rescaler4;
    MorphoWatershedImageFilterType::Pointer         morphWatershed;
    RelabelImageFilterType::Pointer                 relabelImage;
    MaskNegatedImageFilterType::Pointer             maskSinNegatedImageFilter;
    ThresholdFilterType::Pointer                    cellImageToBin;
    ImageToShapeLabelMapFilterType::Pointer         cellBinImageToShapeLabelMap;
    LabelImageToShapeLabelMapFilterType::Pointer    labelImageToShapeLabelMap;
    LabelImageToShapeLabelMapFilterType2::Pointer   labelImageToShapeLabelMap2;
    ShapeOpeningLabelMapFilterType::Pointer         shapeOpeningLabMapFilter;
    ImageToShapeLabelMapFilterType::Pointer         nucleiBinToLabMapFilter;
    LabelMapToLabelImageFilterType::Pointer         nucleiLabelMapToImageFilter;
    LabelMapToLabelImageFilterType::Pointer         cellLabelMapToImageFilter;
    LabelMapToLabelImageFilterType::Pointer         labelMapToImage, labelMapToImage2;
    ThresholdFilterType::Pointer                    thresholdImageFilter;
    LabelContourImageFilterType::Pointer            labelContour, labelContour2;
    ScalarVoWriterType::Pointer                     writer1, writer2, writer3, writer4, writer5, writer6;
    LabelOverlayImageFilterType::Pointer            overlayImage1, overlayImage2, overlayImage3, overlayImage4, overlayImage5;
    RGBVoWriterType::Pointer                        writerRGB1, writerRGB2, writerRGB3, writerRGB4, writerRGB5;
    KeepObjectsNearMiddleFilterType::Pointer        keepNearMiddle;

    //----------READER--------------------------------------------------------------------------------------------------------
    if(GetNumberOfDimensions(dirNuc + fileNuc + extNuc) != 3)
        throw std::string("Please specify a nuclei segmentation file with three dimensions.");

    readerNuc = ScalarVoReaderType::New();
    readerNuc->SetFileName(dirNuc + fileNuc + extNuc);
    readerNuc->ReleaseDataFlagOn();
    readerNuc->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerNuc->SetImageIO( itk::TIFFImageIO::New() );
#endif
    readerNuc->Update();

    itk::SmartPointer<CScalarVoImageType> nucImage = readerNuc->GetOutput();
    nucImage->DisconnectPipeline();

    nucImage->SetSpacing(m_spacing);

    SetupBorderRegions(nucImage);

    if(GetNumberOfDimensions(dirBile + fileBile + extBile) != 3)
        throw std::string("Please specify a bile segmentation file with three dimensions.");

    readerBile = ScalarVoReaderType::New();
    readerBile->SetFileName(dirBile + fileBile + extBile);
    readerBile->ReleaseDataFlagOn();
    readerBile->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerBile->SetImageIO( itk::TIFFImageIO::New() );
#endif

    if(GetNumberOfDimensions(dirSin + fileSin + extSin) != 3)
        throw std::string("Please specify a sinusoid segmentation file with three dimensions.");

    readerSin = ScalarVoReaderType::New();
    readerSin->SetFileName(dirSin + fileSin + extSin);
    readerSin->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerSin->SetImageIO( itk::TIFFImageIO::New() );
#endif
    readerSin->Update();

    if(m_hasCV) {
        if(GetNumberOfDimensions(dirCV + fileCV + extCV) != 3)
            throw std::string("Please specify a central vein segmentation file with three dimensions.");

        readerCV = ScalarVoReaderType::New();
        readerCV->SetFileName(dirCV + fileCV + extCV);
        readerCV->ReleaseDataFlagOn();
        readerCV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
        readerCV->SetImageIO( itk::TIFFImageIO::New() );
#endif
    }

    if(m_hasPV) {
        if(GetNumberOfDimensions(dirPV + filePV + extPV) != 3)
            throw std::string("Please specify a portal vein segmentation file with three dimensions.");

        readerPV = ScalarVoReaderType::New();
        readerPV->SetFileName(dirPV + filePV + extPV);
        readerPV->ReleaseDataFlagOn();
        readerPV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
        readerPV->SetImageIO( itk::TIFFImageIO::New() );
#endif
    }

    if(GetNumberOfDimensions(dirDPPIV + fileDPPIV + extDPPIV) != 3)
        throw std::string("Please specify a DPPIV channel file with three dimensions.");

    readerDPPIV = ScalarVoReaderType::New();
    readerDPPIV->SetFileName(dirDPPIV + fileDPPIV + extDPPIV);
    readerDPPIV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerDPPIV->SetImageIO( itk::TIFFImageIO::New() );
#endif
    readerDPPIV->Update();

    itk::SmartPointer<CScalarVoImageType> dppivImage = readerDPPIV->GetOutput();
    dppivImage->DisconnectPipeline();

    dppivImage->SetSpacing(m_spacing);

    m_XMax = nucImage->GetLargestPossibleRegion().GetSize()[0] - 1;
    m_YMax = nucImage->GetLargestPossibleRegion().GetSize()[1] - 1;
    m_ZMax = nucImage->GetLargestPossibleRegion().GetSize()[2] - 1;

    itk::SmartPointer<CScalarVoImageType> necRegImage;
    if(m_hasNR) {
        if(GetNumberOfDimensions(dirNecReg + fileNecReg + extNecReg) != 3)
            throw std::string("Please specify a necrotic region file with three dimensions.");

        readerNecReg = ScalarVoReaderType::New();
        readerNecReg->SetFileName(dirNecReg + fileNecReg + extNecReg);
        readerNecReg->ReleaseDataBeforeUpdateFlagOn();
        readerNecReg->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
        readerNecReg->SetImageIO( itk::TIFFImageIO::New() );
#endif
        readerNecReg->Update();

        necRegImage = readerNecReg->GetOutput();
        necRegImage->DisconnectPipeline();

        necRegImage->SetSpacing(m_spacing);
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------ADD-SEGMENTATIONS----------------------------------------------------------------------------------------------
    if(m_hasCV && m_hasPV) {
        addVeins = AddCImageFilterType::New();
        addVeins->ReleaseDataFlagOn();
        addVeins->SetInput1(readerCV->GetOutput());
        addVeins->SetInput2(readerPV->GetOutput());
    }

    if(m_hasCV || m_hasPV) {
        addVeinsSin = AddCImageFilterType::New();
        addVeinsSin->ReleaseDataFlagOn();
        if(m_hasCV && m_hasPV)  addVeinsSin->SetInput1(addVeins->GetOutput());
        else if(m_hasCV)		addVeinsSin->SetInput1(readerCV->GetOutput());
        else if(m_hasPV)		addVeinsSin->SetInput1(readerPV->GetOutput());
        addVeinsSin->SetInput2(readerSin->GetOutput());
        addVeinsSin->Update();
    }

    itk::SmartPointer<CScalarVoImageType> addVeinSinImage;
    if(m_hasCV || m_hasPV)	addVeinSinImage  = addVeinsSin->GetOutput();
    else					addVeinSinImage  = readerSin->GetOutput();
    addVeinSinImage->DisconnectPipeline();

    addVeinsSinBile = AddCImageFilterType::New();
    addVeinsSinBile->ReleaseDataFlagOn();
    addVeinsSinBile->SetInput1(addVeinSinImage);
    addVeinsSinBile->SetInput2(readerBile->GetOutput());
    addVeinsSinBile->Update();

    itk::SmartPointer<CScalarVoImageType> addVeinSinBileImage = addVeinsSinBile->GetOutput();
    addVeinSinBileImage->DisconnectPipeline();
    addVeinSinBileImage->SetSpacing(m_spacing);
    addVeinSinImage->SetSpacing(m_spacing);

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writerVeinSinBile = ScalarVoWriterType::New();
        writerVeinSinBile->SetFileName(m_path + m_filenameSave + "_step0_addVeinSinBile.tif");
        writerVeinSinBile->SetInput(addVeinSinBileImage);
#if (ITK_VERSION_MAJOR >= 4)
        writerVeinSinBile->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerVeinSinBile->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------DISTANCE-MAP---------------------------------------------------------------------------------------------------
    distanceMapBileSin = SignedMaurerDistanceMapImageFilterType::New();
    distanceMapBileSin->ReleaseDataFlagOn();
    distanceMapBileSin->UseImageSpacingOn();
    distanceMapBileSin->SquaredDistanceOff();
    distanceMapBileSin->SetInput(addVeinSinBileImage);
    distanceMapBileSin->Update();

    distanceMapBileSinCalc = ImageCalculatorFilterType::New();
    distanceMapBileSinCalc->SetImage(distanceMapBileSin->GetOutput());
    distanceMapBileSinCalc->Compute();

    distanceMapBileSinCapped = IntensityWindowingImageFilter::New();
    distanceMapBileSinCapped->ReleaseDataFlagOn();
    distanceMapBileSinCapped->SetWindowMinimum(0);
    distanceMapBileSinCapped->SetWindowMaximum(distanceMapBileSinCalc->GetMaximum());
    distanceMapBileSinCapped->SetOutputMinimum(0);
    distanceMapBileSinCapped->SetOutputMaximum(distanceMapBileSinCalc->GetMaximum());
    distanceMapBileSinCapped->SetInput(distanceMapBileSin->GetOutput());
    distanceMapBileSinCapped->Update();

    itk::SmartPointer<FScalarVoImageType> distMapBileSinImage = distanceMapBileSinCapped->GetOutput();
    distMapBileSinImage->DisconnectPipeline();

//    RemoveNucleiInContactWithSinusoids(nucImage, addVeinSinBileImage);

    distanceMapNuc = SignedMaurerDistanceMapImageFilterType::New();
    distanceMapNuc->ReleaseDataFlagOn();
    distanceMapNuc->UseImageSpacingOn();
    distanceMapNuc->SquaredDistanceOff();
    distanceMapNuc->SetInput(nucImage);
    distanceMapNuc->Update();

//    itk::SmartPointer<FScalarVoImageType> distMapNucImage = distanceMapNuc->GetOutput();
//    distMapNucImage->DisconnectPipeline();
//
//    //----------ADD-NUCLEUS-RIM-TO-IMAGE---------------------------------------------------------------------------------------
//    SetBorderPixelAtLocalMaximaExceedingThreshold(nucImage, distMapNucImage, 10., 255);             //TODO: distance 10 has to be subject to parametrization
//
//    if(m_saveEverything) {
//        ScalarVoWriterType::Pointer writerNucBorder = ScalarVoWriterType::New();
//        writerNucBorder->SetFileName(m_path + m_filenameSave + "_step0_addBorderNuclei.tif");
//        writerNucBorder->SetInput(nucImage);
//#if (ITK_VERSION_MAJOR >= 4)
//        writerNucBorder->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//        writerNucBorder->Update();
//    }
//
//    distanceMapNuc->ReleaseDataFlagOn();
//    distanceMapNuc->UseImageSpacingOn();
//    distanceMapNuc->SquaredDistanceOff();
//    distanceMapNuc->SetInput(nucImage);
//    distanceMapNuc->Update();

    distanceMapNucCalc = ImageCalculatorFilterType::New();
    distanceMapNucCalc->SetImage(distanceMapNuc->GetOutput());
    distanceMapNucCalc->Compute();

    distanceMapNucCapped = IntensityWindowingImageFilter::New();
    distanceMapNucCapped->ReleaseDataFlagOn();
    distanceMapNucCapped->SetWindowMinimum(0);
    distanceMapNucCapped->SetWindowMaximum(distanceMapNucCalc->GetMaximum());
    distanceMapNucCapped->SetOutputMinimum(0);
    distanceMapNucCapped->SetOutputMaximum(distanceMapNucCalc->GetMaximum());
    distanceMapNucCapped->SetInput(distanceMapNuc->GetOutput());
    distanceMapNucCapped->Update();

    itk::SmartPointer<FScalarVoImageType> distMapNucImage = distanceMapNucCapped->GetOutput();
    distMapNucImage->DisconnectPipeline();
    //-------------------------------------------------------------------------------------------------------------------------

    std::cout << "Write distMap intermediate results " << m_saveEverything << std::endl;
    if(m_saveEverything) {
        rescaler1 = RescaleImageFilterType::New();
        rescaler1->ReleaseDataFlagOn();
        rescaler1->SetInput(distMapNucImage);

        rescaler2 = RescaleImageFilterType::New();
        rescaler2->ReleaseDataFlagOn();
        rescaler2->SetInput(distMapBileSinImage);

        writer2 = ScalarVoWriterType::New();
        writer2->ReleaseDataFlagOn();
        writer2->SetFileName(m_path + m_filenameSave + "_step1_distMapNuc.tif");
        writer2->SetInput(rescaler1->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer2->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer2->Update();

        writer3 = ScalarVoWriterType::New();
        writer3->ReleaseDataFlagOn();
        writer3->SetFileName(m_path + m_filenameSave + "_step1_distMapBileSin.tif");
        writer3->SetInput(rescaler2->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer3->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer3->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------MULTIPLY-WITH-WEIGHTS------------------------------------------------------------------------------------------
    multNuc = MultiplyImageFilterType::New();
    multNuc->SetInput(distMapNucImage);
    multNuc->SetConstant(m_nucleiWeight);

    multBS = MultiplyImageFilterType::New();
    multBS->ReleaseDataFlagOn();
    multBS->SetInput(distMapBileSinImage);
    multBS->SetConstant(m_bileSinWeight);
    //------------------------------------------------------------------------------------------------------------------------

    //----------CALCULATE-REFERENCE-IMAGE-------------------------------------------------------------------------------------
    add2 = AddFImageFilterType::New();
    add2->SetInput1(multNuc->GetOutput());
    add2->SetInput2(multBS->GetOutput());
    add2->Update();

    ImageCalculatorFilterType::Pointer add2Calc = ImageCalculatorFilterType::New();
    add2Calc->SetImage(add2->GetOutput());
    add2Calc->Compute();

    IntensityWindowingImageFilter::Pointer add2Capped = IntensityWindowingImageFilter::New();
    add2Capped->SetWindowMinimum(0.0000001);
    add2Capped->SetWindowMaximum(add2Calc->GetMaximum());
    add2Capped->SetOutputMinimum(0.0000001);
    add2Capped->SetOutputMaximum(add2Calc->GetMaximum());
    add2Capped->SetInput(add2->GetOutput());

    divide = DivideImageFilterType::New();
    divide->SetInput1(multNuc->GetOutput());
    divide->SetInput2(add2Capped->GetOutput());
    divide->Update();

    ImageCalculatorFilterType::Pointer divideCalc = ImageCalculatorFilterType::New();
    divideCalc->SetImage(divide->GetOutput());
    divideCalc->Compute();

    IntensityWindowingImageFilter::Pointer divideCapped = IntensityWindowingImageFilter::New();
    divideCapped->SetWindowMinimum(divideCalc->GetMinimum());
    divideCapped->SetWindowMaximum(divideCalc->GetMaximum());
    divideCapped->SetOutputMinimum(0.);
    divideCapped->SetOutputMaximum(100.);
    divideCapped->SetInput(divide->GetOutput());
    divideCapped->Update();

    std::cout << "Write formula intermediate results " << m_saveEverything << std::endl;
    if(m_saveEverything) {
        rescaler3 = RescaleImageFilterType::New();
        rescaler3->ReleaseDataFlagOn();
        rescaler3->SetInput(add2Capped->GetOutput());

        rescaler4 = RescaleImageFilterType::New();
        rescaler4->ReleaseDataFlagOn();
        rescaler4->SetInput(divideCapped->GetOutput());

        writer4 = ScalarVoWriterType::New();
        writer4->SetFileName(m_path + m_filenameSave + "_step2_add2.tif");
        writer4->SetInput(rescaler3->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer4->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer4->Update();

        writer5 = ScalarVoWriterType::New();
        writer5->SetFileName(m_path + m_filenameSave + "_step2_divide.tif");
        writer5->SetInput(rescaler4->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer5->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer5->Update();
    }

    distMapNucImage->ReleaseData();
    distMapBileSinImage->ReleaseData();
    //------------------------------------------------------------------------------------------------------------------------

    //----------WATERSHED-----------------------------------------------------------------------------------------------------
    morphWatershed = MorphoWatershedImageFilterType::New();
    morphWatershed->ReleaseDataFlagOn();
    morphWatershed->SetLevel(m_watershedFloodLevel);
    morphWatershed->FullyConnectedOff();
    morphWatershed->SetInput(divideCapped->GetOutput());
    //------------------------------------------------------------------------------------------------------------------------

    //----------ADD-SEGMENTATIONS----------------------------------------------------------------------------------------------
    //Begin: This code segment only valid for special case of dataset 250; has to be deleted for other datasets
//    itk::Size<3> rad;
//    rad[0] = 8;
//    rad[1] = 8;
//    rad[2] = 4;
//
//    StructuringElementType openingStructuringElement;
//    openingStructuringElement.SetRadius(rad);
//    openingStructuringElement.CreateStructuringElement();
//
//    OpeningImageFilterType::Pointer openingFilter = OpeningImageFilterType::New();
//    openingFilter->SetKernel(openingStructuringElement);
//    openingFilter->SetBackgroundValue(0);
//    openingFilter->SetForegroundValue(255);
//    openingFilter->SetInput(readerSin->GetOutput());
//
//    AddCImageFilterType::Pointer addVeinsSin2 = AddCImageFilterType::New();
//    addVeinsSin2->ReleaseDataFlagOn();
//    if(f4 && f5)    addVeinsSin2->SetInput1(addVeins->GetOutput());
//    else if(f4)     addVeinsSin2->SetInput1(readerCV->GetOutput());
//    else if(f5)     addVeinsSin2->SetInput1(readerPV->GetOutput());
//    addVeinsSin2->SetInput2(openingFilter->GetOutput());
//    addVeinsSin2->Update();
//
//    itk::SmartPointer<CScalarVoImageType> addVeinSinOpenImage = addVeinsSin2->GetOutput();
//    addVeinSinOpenImage->DisconnectPipeline();
//    addVeinSinOpenImage->SetSpacing(m_spacing);
    //END: This code segment only valid for special case of dataset 250; has to be deleted for other datasets

    itk::SmartPointer<CScalarVoImageType> addVeinSinNecRegImage;

    if(m_hasNR) {
        addVeinsSinNecReg = AddCImageFilterType::New();
        addVeinsSinNecReg->ReleaseDataFlagOn();
        addVeinsSinNecReg->SetInput1(addVeinSinImage);
        addVeinsSinNecReg->SetInput2(necRegImage);
        addVeinsSinNecReg->Update();

        addVeinSinNecRegImage = addVeinsSinNecReg->GetOutput();
        addVeinSinNecRegImage->DisconnectPipeline();

        addVeinSinNecRegImage->SetSpacing(m_spacing);
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------REMOVE-SINUSOID-AREA------------------------------------------------------------------------------------------
    maskSinNegatedImageFilter = MaskNegatedImageFilterType::New();
    maskSinNegatedImageFilter->ReleaseDataFlagOn();
    maskSinNegatedImageFilter->SetInput(morphWatershed->GetOutput());
    if(m_hasNR)     maskSinNegatedImageFilter->SetMaskImage(addVeinSinNecRegImage);
//    else            maskSinNegatedImageFilter->SetMaskImage(addVeinSinOpenImage);
    else            maskSinNegatedImageFilter->SetMaskImage(addVeinSinImage);
    //------------------------------------------------------------------------------------------------------------------------

    //----------RELABEL-WATERSHED-LABEL-OBJECTS-------------------------------------------------------------------------------
    cellImageToBin = ThresholdFilterType::New();
    cellImageToBin->SetOutsideValue(0);
    cellImageToBin->SetInsideValue(itk::NumericTraits<LScalarPixelType>::max());
    cellImageToBin->SetLowerThreshold(1);
    cellImageToBin->SetUpperThreshold(itk::NumericTraits<LabelMapToLabelImageFilterType::OutputImagePixelType>::max());
    cellImageToBin->SetInput(maskSinNegatedImageFilter->GetOutput());
    cellImageToBin->Update();

    itk::SmartPointer<CScalarVoImageType> cellBinImage = cellImageToBin->GetOutput();
    cellBinImage->DisconnectPipeline();

    cellBinImage->SetSpacing(m_spacing);

    cellBinImageToShapeLabelMap = ImageToShapeLabelMapFilterType::New();
    cellBinImageToShapeLabelMap->SetFullyConnected(false);
    cellBinImageToShapeLabelMap->SetInput(cellBinImage);
    cellBinImageToShapeLabelMap->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------------
    std::cout << "remove objects with volume smaller than " << m_minimalCellVolume << std::endl;

    std::cout << "before volume based removal " << cellBinImageToShapeLabelMap->GetOutput()->GetNumberOfLabelObjects() << " objects" << std::endl;

    shapeOpeningLabMapFilter = ShapeOpeningLabelMapFilterType::New();
    shapeOpeningLabMapFilter->SetLambda(m_minimalCellVolume);                                           //attribute value
    shapeOpeningLabMapFilter->ReverseOrderingOff();                                                       //removes objects with attribute smaller than lambda
    shapeOpeningLabMapFilter->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::PHYSICAL_SIZE);
    shapeOpeningLabMapFilter->SetInput(cellBinImageToShapeLabelMap->GetOutput());
    shapeOpeningLabMapFilter->Update();

    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap = shapeOpeningLabMapFilter->GetOutput();
    cellLabelMap->DisconnectPipeline();
    std::cout << "after volume based removal " << cellLabelMap->GetNumberOfLabelObjects() << " objects" << std::endl;

    shapeOpeningLabMapFilter->Update();
    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithOneNucleusLabelMap = shapeOpeningLabMapFilter->GetOutput();
    cellsWithOneNucleusLabelMap->DisconnectPipeline();

    shapeOpeningLabMapFilter->Update();
    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithTwoNucleiLabelMap = shapeOpeningLabMapFilter->GetOutput();
    cellsWithTwoNucleiLabelMap->DisconnectPipeline();

    shapeOpeningLabMapFilter->Update();
    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithMoreNucleiLabelMap = shapeOpeningLabMapFilter->GetOutput();
    cellsWithMoreNucleiLabelMap->DisconnectPipeline();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------REMOVE-OBECTS-WITHOUT-NUCLEI-----------------------------------------------------------------------------------
    nucleiBinToLabMapFilter = ImageToShapeLabelMapFilterType::New();
    nucleiBinToLabMapFilter->SetInput(nucImage);
    nucleiBinToLabMapFilter->Update();

    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> nucleiLabelMap = nucleiBinToLabMapFilter->GetOutput();
    nucleiLabelMap->DisconnectPipeline();

    nucleiLabelMapToImageFilter = LabelMapToLabelImageFilterType::New();
    nucleiLabelMapToImageFilter->SetInput(nucleiLabelMap);
    nucleiLabelMapToImageFilter->Update();

    cellLabelMapToImageFilter = LabelMapToLabelImageFilterType::New();
    cellLabelMapToImageFilter->SetInput(cellLabelMap);
    cellLabelMapToImageFilter->Update();

//    LScalarVoImageType::Pointer nucleiLabelImage = nucleiLabelMapToImageFilter->GetOutput();
//    LScalarVoImageType::Pointer cellLabelImage = cellLabelMapToImageFilter->GetOutput();

    BuildCellNucleiAlignmentMaps(nucleiLabelMapToImageFilter->GetOutput(), cellLabelMapToImageFilter->GetOutput(), nucleiLabelMap, cellLabelMap);

    std::cout << "before removal " << cellLabelMap->GetNumberOfLabelObjects() << " cells" << std::endl;
    for(unsigned int i=0; i<m_labelsWithZeroNuclei.size(); i++)
        cellLabelMap->RemoveLabel(m_labelsWithZeroNuclei[i]);
    std::cout << "after removal of cells without a nucleus " << cellLabelMap->GetNumberOfLabelObjects() << " cells left" << std::endl;

    RemoveCellsAtDatasetBorder(cellLabelMapToImageFilter->GetOutput(), cellsWithOneNucleusLabelMap, cellsWithTwoNucleiLabelMap, cellsWithMoreNucleiLabelMap, 0);
#ifndef CS_TI_QUANT_ONLY
    WriteNumNucleiOutlineFiles(dppivImage, cellsWithOneNucleusLabelMap, cellsWithTwoNucleiLabelMap, cellsWithMoreNucleiLabelMap); //keep near position filter not applicable to this atm
#endif
    //-------------------------------------------------------------------------------------------------------------------------

    //----------KEEP-OBJECTS-NEAR-POSITION------------------------------------------------------------------------------------
    std::cout << "Use position filter " << m_withPositionFilter << std::endl;

    if(m_withPositionFilter) {
        keepNearMiddle = KeepObjectsNearMiddleFilterType::New();
        keepNearMiddle->ReverseOrderingOff();
        keepNearMiddle->SetNumberOfObjects(m_numberCells);
        keepNearMiddle->SetPosition(m_pos);
        keepNearMiddle->SetInput(cellLabelMap);
        keepNearMiddle->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------REGISTER-BORDER-VOXEL-POSITIONS--------------------------------------------------------------------------------
    std::vector< std::vector<ShapeOpeningLabelMapFilterType::IndexType> > borderVoxel;

    std::cout << "Write cell outline to file " << m_writeCellOutlineToFile << std::endl;

    if(m_writeCellOutlineToFile) {
        ShapeOpeningLabelMapFilterType::IndexType idx;
        int numLabelObjects;
        if(m_withPositionFilter)    numLabelObjects = keepNearMiddle->GetOutput()->GetNumberOfLabelObjects();
        else                        numLabelObjects = cellLabelMap->GetNumberOfLabelObjects();

        for(int i=0; i<numLabelObjects; i++) {
            int numVoxel;
            if(m_withPositionFilter)    numVoxel = keepNearMiddle->GetOutput()->GetNthLabelObject(i)->Size();
            else                        numVoxel = cellLabelMap->GetNthLabelObject(i)->Size();
            std::vector<ShapeOpeningLabelMapFilterType::IndexType> bV;

            for(int j=0; j<numVoxel; j++) {
                if(m_withPositionFilter)    idx = keepNearMiddle->GetOutput()->GetNthLabelObject(i)->GetIndex(j);
                else                        idx = cellLabelMap->GetNthLabelObject(i)->GetIndex(j);

                if(idx[0] == 0 || idx[1] == 0 || idx[2] == 0 || idx[0] == m_XMax || idx[1] == m_YMax || idx[2] == m_ZMax)
                    bV.push_back(idx);
            }
            borderVoxel.push_back(bV);
        }
    }

    labelMapToImage = LabelMapToLabelImageFilterType::New();
    if(m_withPositionFilter)    keepNearMiddle->SetInput(cellLabelMap);
    else                        labelMapToImage->SetInput(cellLabelMap);
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
    thresholdImageFilter = ThresholdFilterType::New();
    thresholdImageFilter->SetOutsideValue(0);
    thresholdImageFilter->SetInsideValue(itk::NumericTraits<CScalarPixelType>::max());
    thresholdImageFilter->SetLowerThreshold(1);
    thresholdImageFilter->SetUpperThreshold(itk::NumericTraits<LabelMapToLabelImageFilterType::OutputImagePixelType>::max());
    thresholdImageFilter->SetInput(labelMapToImage->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------CONTOUR-WATERSHED-REGIONS-------------------------------------------------------------------------------------
    labelContour = LabelContourImageFilterType::New();
    labelContour->SetFullyConnected(true);
    labelContour->SetInput(labelMapToImage->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    labelImageToShapeLabelMap2 = LabelImageToShapeLabelMapFilterType2::New();
    labelImageToShapeLabelMap2->SetBackgroundValue(0);
    labelImageToShapeLabelMap2->SetInput(labelContour->GetOutput());
    labelImageToShapeLabelMap2->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------WRITE-LABEL-OBJECTS-TO-FILE-----------------------------------------------------------------------------------
    if(m_writeCellOutlineToFile) {
        std::fstream file;

        //TODO: write also center of mass to file
        file.open((m_path + m_filenameSave + "_outlineCells.txt").c_str(), std::fstream::out);

        LabelImageToShapeLabelMapFilterType2::IndexType idx2;

        for(unsigned int i=0; i<labelImageToShapeLabelMap2->GetOutput()->GetNumberOfLabelObjects(); i++) {
            file << i << ": " << std::endl;

            int numVoxel = labelImageToShapeLabelMap2->GetOutput()->GetNthLabelObject(i)->Size();

            for(int j=0; j<numVoxel; j++) {
                idx2 = labelImageToShapeLabelMap2->GetOutput()->GetNthLabelObject(i)->GetIndex(j);

                file << idx2[0] << ", " << idx2[1] << ", " << idx2[2] << "; ";
            }
            for(unsigned int j=0; j<borderVoxel[i].size(); j++)
                file << borderVoxel[i][j][0] << ", " << borderVoxel[i][j][1] << ", " << borderVoxel[i][j][2] << "; ";

            file << std::endl;
        }
        file.close();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------VISULALIZATION-ONLY-FILTER---NOT-NEEDED-BY-OTHER-PIPELINE-FILTERS----------------------------------------------
    overlayImage1 = LabelOverlayImageFilterType::New();
    overlayImage1->ReleaseDataFlagOn();
    overlayImage1->SetOpacity(m_overlayOpacity);
    overlayImage1->SetInput(dppivImage);
    overlayImage1->SetLabelImage(labelMapToImage->GetOutput());

    overlayImage2 = LabelOverlayImageFilterType::New();
    overlayImage2->ReleaseDataFlagOn();
    overlayImage2->SetOpacity(m_overlayOpacity);
    overlayImage2->SetInput(dppivImage);
    overlayImage2->SetLabelImage(labelContour->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------RESCALE-IMAGES------------------------------------------------------------------------------------------------
    RescaleLImageFilterType::Pointer rescaler5 = RescaleLImageFilterType::New();
    rescaler5->ReleaseDataFlagOn();
    rescaler5->SetInput(labelMapToImage->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------WRITER--------------------------------------------------------------------------------------------------------
    writer1 = ScalarVoWriterType::New();
    writer1->SetFileName(m_path + m_filenameSave + m_saveSuffixesForFinals[0] + m_filenameExtension);
    writer1->SetInput(thresholdImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer1->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer1->Update();

    writer6 = ScalarVoWriterType::New();
    writer6->SetFileName(m_path + m_filenameSave + m_saveSuffixesForFinals[1] + m_filenameExtension);
    writer6->SetInput(rescaler5->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer6->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer6->Update();

    writerRGB1 = RGBVoWriterType::New();
    writerRGB1->SetFileName(m_path + m_filenameSave + m_saveSuffixesForFinals[2] + m_filenameExtension);
    writerRGB1->SetInput(overlayImage1->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writerRGB1->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerRGB1->Update();

    writerRGB2 = RGBVoWriterType::New();
    writerRGB2->SetFileName(m_path + m_filenameSave + m_saveSuffixesForFinals[3] + m_filenameExtension);
    writerRGB2->SetInput(overlayImage2->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writerRGB2->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerRGB2->Update();
    //-------------------------------------------------------------------------------------------------------------------------


    std::cout << "Write sinusoid outline to file " << m_writeSinusoidOutlineToFile << std::endl;
    if(m_writeSinusoidOutlineToFile)
        WriteSiunsoidOutlineToFile();

    std::cout << "Write sinusoid graph to file " << m_writeSinusoidGraphToFile << std::endl;
    if(m_writeSinusoidGraphToFile)
        WriteSinusoidGraphToFile();

    WriteLogFile(timeStamp);
    WriteDataSetSummary();


    nucImage->ReleaseData();
    nucImage = NULL;

    dppivImage->ReleaseData();
    dppivImage = NULL;

    if(m_hasNR) {
        necRegImage->ReleaseData();
        necRegImage = NULL;
    }

    addVeinSinImage->ReleaseData();
    addVeinSinImage = NULL;

    addVeinSinBileImage->ReleaseData();
    addVeinSinBileImage = NULL;

    distMapBileSinImage = NULL;

    distMapNucImage = NULL;

    if(m_hasNR) {
        addVeinSinNecRegImage->ReleaseData();
        addVeinSinNecRegImage = NULL;
    }

    cellBinImage->ReleaseData();
    cellBinImage = NULL;
}
