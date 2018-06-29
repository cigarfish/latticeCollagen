///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SegmentCellMembraneOnBCat60x.cpp                                     //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-09-02                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "SegmentCellMembraneOnBCat60x.h"

#if (ITK_VERSION_MAJOR < 4)
#define NUMBER_OF_PIXELS SIZE
#endif

#include <QDateTime>

#include "../tools/ImageAnalysisSummaryFileIO.h"

#include "../../tools/input/FilenameParser.h"
#include "../../tools/parameters/CSParameter.h"
#include "../../tools/parameters/CSParameterContext.h"

#include "itkPointSet.h"


SegmentCellMembraneOnBCat60x::SegmentCellMembraneOnBCat60x()
{
    m_overlayOpacity = 0.5;

    m_saveSuffixesForFinals[0] = "_step7_cell_bin";
    m_saveSuffixesForFinals[1] = "_step7_bcat_cell_overlay";
    m_saveSuffixesForFinals[2] = "_step7_nuc_cell_overlay";
    m_saveSuffixesForFinals[3] = "_step7_bcat_monoNucOverlay";
    m_saveSuffixesForFinals[4] = "_step7_bcat_biNucOverlay";
    m_saveSuffixesForFinals[5] = "_step7_bcat_multNucOverlay";
}


SegmentCellMembraneOnBCat60x::~SegmentCellMembraneOnBCat60x()
{
    // TODO Auto-generated destructor stub
}


void SegmentCellMembraneOnBCat60x::WriteLogFile(std::string timeStamp)
{
    std::fstream file;

    std::stringstream parameter;
    m_paramContext->dump(parameter);

    file.open((m_pathChannelBCat + m_logFilenameSave + ".txt").c_str(), std::ios::out | std::ios::app);

    file << "Begin-segment-cell-membrane-60x-----------------------------------------------------------------------------------------------------------\n";
    file << "-------------------------time stamp start: " << timeStamp << "\n";
    file << parameter.str();
    file << "End-segment-cell-membrane-60x-------------------------------------------------------------------------------------------------------------\n\n";

    file.close();
}


void SegmentCellMembraneOnBCat60x::WriteDataSetSummary()
{
//    ImageAnalysisSummaryFileIO::AddEntry(NonHepNucleiSegmentationBin, m_pathChannelDAPI, m_pathChannelDAPI + m_filenameSave + m_saveSuffixesForFinals[0] + m_fileExtensionChannelDAPI);
}


void SegmentCellMembraneOnBCat60x::ParseParameterContext()
{
    if(m_paramContext->findContext("Segment Cell Membrane",0)==NULL) {
        std::cout << "Error: SegmentCellmembraneOnBCat60xContext: Invalid parameter context" << std::endl;
        return;
    }

    m_fullFilenameChannelBCat = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("BCat channel", 0)->dataPointer()) );
    m_infoFullFilenameChannelDAPI.setFile(m_fullFilenameChannelBCat);

    if(!m_infoFullFilenameChannelDAPI.exists())
        throw std::string("Please specify BCat channel");

    m_pathChannelBCat = (m_infoFullFilenameChannelDAPI.path() + QString("/")).toStdString();
    m_filenameChannelBCat = m_infoFullFilenameChannelDAPI.baseName().toStdString();
    m_fileExtensionChannelBCat = (QString(".") + m_infoFullFilenameChannelDAPI.suffix()).toStdString();

    m_fullFilenameSegmentationSin = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Sinusoid segmentation", 0)->dataPointer()) );
    m_fullFilenameSegmentationNuc = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Hepatic nuclei segmentation", 0)->dataPointer()) );
    m_fullFilenameSegmentationCV = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Central vein segmentation", 0)->dataPointer()) );
    m_fullFilenameSegmentationPV = QString::fromStdString( *(std::string*)(m_paramContext->findParameter("Portal vein segmentation", 0)->dataPointer()) );

    m_spacing[0] = *(double*)(m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    m_spacing[1] = *(double*)(m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    m_spacing[2] = *(double*)(m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    std::string thresMode = ( (CSParameterChoice*)(m_paramContext->findParameter("Threshold mode", 0)->dataPointer()) )->currentString();
    if( thresMode.compare("Adaptive Otsu-Threshold")==0 )
        m_thresholdingMode = 0;
    else if( thresMode.compare("Otsu-Threshold")==0 )
        m_thresholdingMode = 1;
    else if( thresMode.compare("Manual Threshold")==0 )
        m_thresholdingMode = 2;

    m_adapOtsuRadius[0] = *(int*)(m_paramContext->findParameter("Sample region size x", 0)->dataPointer());
    m_adapOtsuRadius[1] = *(int*)(m_paramContext->findParameter("Sample region size y", 0)->dataPointer());
    m_adapOtsuRadius[2] = *(int*)(m_paramContext->findParameter("Sample region size z", 0)->dataPointer());
    m_adapOtsuSamplePoints = *(int*)(m_paramContext->findParameter("Number samples", 0)->dataPointer());
    m_lowerThreshold = *(int*)(m_paramContext->findParameter("BCat manual min threshold", 0)->dataPointer());
    m_upperThreshold = *(int*)(m_paramContext->findParameter("BCat manual max threshold", 0)->dataPointer());

    m_holeFillingNeighborhoodRadius[0] = *(int*)(m_paramContext->findContext("1.2) Hole Filling on 1.1",0)->findParameter("Radius x", 0)->dataPointer());
    m_holeFillingNeighborhoodRadius[1] = *(int*)(m_paramContext->findContext("1.2) Hole Filling on 1.1",0)->findParameter("Radius y", 0)->dataPointer());
    m_holeFillingNeighborhoodRadius[2] = *(int*)(m_paramContext->findContext("1.2) Hole Filling on 1.1",0)->findParameter("Radius z", 0)->dataPointer());
    m_holeFillingMajThreshold = *(int*)(m_paramContext->findContext("1.2) Hole Filling on 1.1",0)->findParameter("Majority threshold", 0)->dataPointer());

    m_invHoleFillingNeighborhoodRadius[0] = *(int*)(m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2",0)->findParameter("Radius x", 0)->dataPointer());
    m_invHoleFillingNeighborhoodRadius[1] = *(int*)(m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2",0)->findParameter("Radius y", 0)->dataPointer());
    m_invHoleFillingNeighborhoodRadius[2] = *(int*)(m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2",0)->findParameter("Radius z", 0)->dataPointer());
    m_invHoleFillingMajThreshold = *(int*)(m_paramContext->findContext("1.3) Inverse Hole Filling on 1.2",0)->findParameter("Majority threshold", 0)->dataPointer());

    itk::Size<3> rad;
    rad[0] = ( *(int*)(m_paramContext->findContext("1.4) Closing on 1.3", 0)->findParameter("Kernel radius x", 0)->dataPointer()) );
    rad[1] = ( *(int*)(m_paramContext->findContext("1.4) Closing on 1.3", 0)->findParameter("Kernel radius y", 0)->dataPointer()) );
    rad[2] = ( *(int*)(m_paramContext->findContext("1.4) Closing on 1.3", 0)->findParameter("Kernel radius z", 0)->dataPointer()) );
    m_closingStructuringElement.SetRadius(rad);
    m_closingStructuringElement.CreateStructuringElement();

    m_openingStructuringElement1.SetRadius(rad);
    m_openingStructuringElement1.CreateStructuringElement();

    rad[0] = ( *(int*)(m_paramContext->findContext("1.5) Opening on 1.4", 0)->findParameter("Kernel radius x", 0)->dataPointer()) );
    rad[1] = ( *(int*)(m_paramContext->findContext("1.5) Opening on 1.4", 0)->findParameter("Kernel radius y", 0)->dataPointer()) );
    rad[2] = ( *(int*)(m_paramContext->findContext("1.5) Opening on 1.4", 0)->findParameter("Kernel radius z", 0)->dataPointer()) );
    m_openingStructuringElement2.SetRadius(rad);
    m_openingStructuringElement2.CreateStructuringElement();

    m_floodLevel = *(double*)(m_paramContext->findParameter("Watershed flood level", 0)->dataPointer());
    m_minimalCellDiameter = *(double*)(m_paramContext->findParameter("Smallest cell diameter", 0)->dataPointer());
    m_minimalCellVolume = 1./6. * itk::Math::pi * pow(m_minimalCellDiameter, 3);

    std::string saveMode = ( (CSParameterChoice*)(m_paramContext->findParameter("Save mode", 0)->dataPointer()) )->currentString();
    if( saveMode.compare("Save everything")==0 )
        m_saveEverything = 1;
    else if( saveMode.compare("Save only essentials")==0 )
        m_saveEverything = 0;
    m_logFilenameSave = *(std::string*)(m_paramContext->findParameter("Log file name", 0)->dataPointer());
    m_filenameSave = *(std::string*)(m_paramContext->findParameter("Save prefix", 0)->dataPointer());
}


void SegmentCellMembraneOnBCat60x::BuildCellNucleiAlignmentMaps(LScalarVoImageType::Pointer nucleiLabelImage, LScalarVoImageType::Pointer cellLabelImage,
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


void SegmentCellMembraneOnBCat60x::RemoveCellsAtDatasetBorder(LScalarVoImageType::Pointer cellLabelImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithOneNucleusLabelMap,
        itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithTwoNucleiLabelMap, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithMoreNucleiLabelMap)
{
    std::set<unsigned long> cellsWithDatasetBorderContacts;

    itk::Index<3> topI, downI, frontI, rearI, rightI, leftI;
    topI[0] = 0;    downI[0] = 0;
    topI[1] = 0;    downI[1] = 0;
    topI[2] = 0;    downI[2] = cellLabelImage->GetLargestPossibleRegion().GetSize()[2]-1;

    frontI[0] = 0;  rearI[0] = 0;
    frontI[1] = 0;  rearI[1] = cellLabelImage->GetLargestPossibleRegion().GetSize()[1]-1;
    frontI[2] = 1;  rearI[2] = 1;

    rightI[0] = 0;   leftI[0] = cellLabelImage->GetLargestPossibleRegion().GetSize()[0]-1;
    rightI[1] = 1;   leftI[1] = 1;
    rightI[2] = 1;   leftI[2] = 1;

    itk::Size<3> topS, downS, frontS, rearS, rightS, leftS;
    topS[0] = cellLabelImage->GetLargestPossibleRegion().GetSize()[0]-1;        downS[0] = cellLabelImage->GetLargestPossibleRegion().GetSize()[0]-1;
    topS[1] = cellLabelImage->GetLargestPossibleRegion().GetSize()[1]-1;        downS[1] = cellLabelImage->GetLargestPossibleRegion().GetSize()[1]-1;
    topS[2] = 1;                                                                downS[2] = 1;

    frontS[0] = cellLabelImage->GetLargestPossibleRegion().GetSize()[0]-1;      rearS[0] = cellLabelImage->GetLargestPossibleRegion().GetSize()[0]-1;
    frontS[1] = 1;                                                              rearS[1] = 1;
    frontS[2] = cellLabelImage->GetLargestPossibleRegion().GetSize()[2]-1-2;    rearS[2] = cellLabelImage->GetLargestPossibleRegion().GetSize()[2]-1-2;

    rightS[0] = 1;                                                              leftS[0] = 1;
    rightS[1] = cellLabelImage->GetLargestPossibleRegion().GetSize()[1]-1-2;    leftS[1] = cellLabelImage->GetLargestPossibleRegion().GetSize()[1]-1-2;
    rightS[2] = cellLabelImage->GetLargestPossibleRegion().GetSize()[2]-1-2;    leftS[2] = cellLabelImage->GetLargestPossibleRegion().GetSize()[2]-1-2;     //-1-2 = starts with 0 -> -1; prevent double counting -> ignore first & last row of face -> -2

    typedef itk::ImageRegion<3> ImageRegionType;
    ImageRegionType dataSetFaces[6];
    dataSetFaces[0].SetIndex(topI);     dataSetFaces[0].SetSize(topS);
    dataSetFaces[1].SetIndex(downI);    dataSetFaces[1].SetSize(downS);
    dataSetFaces[2].SetIndex(frontI);   dataSetFaces[2].SetSize(frontS);
    dataSetFaces[3].SetIndex(rearI);    dataSetFaces[3].SetSize(rearS);
    dataSetFaces[4].SetIndex(rightI);   dataSetFaces[4].SetSize(rightS);
    dataSetFaces[5].SetIndex(leftI);    dataSetFaces[5].SetSize(leftS);

    for(unsigned int i=0; i<6; i++) {
        itk::ImageRegionConstIterator<LScalarVoImageType> iter(cellLabelImage, dataSetFaces[i]);

        for(iter = iter.Begin(); !iter.IsAtEnd(); ++iter) {
            if(iter.Value() != 0 && cellsWithDatasetBorderContacts.count(iter.Value()) == 0)
                cellsWithDatasetBorderContacts.insert(iter.Value());
        }
    }

    for(std::set<unsigned long>::iterator it = cellsWithDatasetBorderContacts.begin(); it != cellsWithDatasetBorderContacts.end(); ++it) {
        cellsWithOneNucleusLabelMap->RemoveLabel(*it);
        cellsWithTwoNucleiLabelMap->RemoveLabel(*it);
        cellsWithMoreNucleiLabelMap->RemoveLabel(*it);
    }
}


void SegmentCellMembraneOnBCat60x::WriteNumNucleiOutlineFiles(itk::SmartPointer<CScalarVoImageType> origImage, itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithOneNucleusLabelMap,
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
    overlayImage3->SetInput(origImage);
    overlayImage3->SetLabelImage(labelMapToImageA->GetOutput());

    overlayImage4 = LabelOverlayImageFilterType::New();
    overlayImage4->ReleaseDataFlagOn();
    overlayImage4->SetOpacity(m_overlayOpacity);
    overlayImage4->SetInput(origImage);
    overlayImage4->SetLabelImage(labelMapToImageB->GetOutput());

    overlayImage5 = LabelOverlayImageFilterType::New();
    overlayImage5->ReleaseDataFlagOn();
    overlayImage5->SetOpacity(m_overlayOpacity);
    overlayImage5->SetInput(origImage);
    overlayImage5->SetLabelImage(labelMapToImageC->GetOutput());

    writerRGB3 = RGBVoWriterType::New();
    writerRGB3->SetFileName(m_pathChannelBCat + m_filenameSave + m_saveSuffixesForFinals[3] + m_fileExtensionChannelBCat);
    writerRGB3->SetInput(overlayImage3->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writerRGB3->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerRGB3->Update();

    writerRGB4 = RGBVoWriterType::New();
    writerRGB4->SetFileName(m_pathChannelBCat + m_filenameSave + m_saveSuffixesForFinals[4] + m_fileExtensionChannelBCat);
    writerRGB4->SetInput(overlayImage4->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writerRGB4->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerRGB4->Update();

    writerRGB5 = RGBVoWriterType::New();
    writerRGB5->SetFileName(m_pathChannelBCat + m_filenameSave + m_saveSuffixesForFinals[5] + m_fileExtensionChannelBCat);
    writerRGB5->SetInput(overlayImage5->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writerRGB5->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writerRGB5->Update();
}


void SegmentCellMembraneOnBCat60x::Update()
{
    ParseParameterContext();

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    std::cout << "Segment cell membrane: " << std::endl;
    std::cout << " dir: " << m_pathChannelBCat << std::endl;
    std::cout << " file: " << m_filenameChannelBCat << std::endl;
    std::cout << " ext: " << m_fileExtensionChannelBCat << std::endl;


    ScalarVoReaderType::Pointer                     readerBCat;
    ScalarVoReaderType::Pointer                     readerSin;
    ScalarVoReaderType::Pointer                     readerNuc;
    ScalarVoReaderType::Pointer                     readerCV;
    ScalarVoReaderType::Pointer                     readerPV;
    AdaptiveOtsuThresholdImageFilterType::Pointer   adapOtsuFilter;
    OtsuThresholdImageFilterType::Pointer           otsuFilter;
    InvertIntensityImageFilterType::Pointer         invertFilter;
    ThresholdFilterType::Pointer                    thresBCatFilter;
    HoleFillingImageFilterType::Pointer             holeFillingFilter;
    HoleFillingImageFilterType::Pointer             invHoleFillingFilter;
    ClosingImageFilterType::Pointer                 closingFilter;
    OpeningImageFilterType::Pointer                 openingFilter1;
    InvertIntensityImageFilterType::Pointer         invertForWatershedFilter;
    ImageToShapeLabelMapFilterType::Pointer         imageToShapeLabelMap;
    ShapeOpeningLabelMapFilterType::Pointer         shapeOpeningLabMapFilter1;
    LabelMapToLabelImageFilterType::Pointer         labelMapToImage1;
    ThresholdLScalarFilterType::Pointer             thresholdLabelImageFilter1;
    SignedMaurerDistanceMapImageFilterType::Pointer distanceMap;
    RescaleImageFilterType::Pointer                 rescaler;
    MorphoWatershedImageFilterType::Pointer         morphWatershed;
    OpeningImageFilterType::Pointer                 openingFilter2;
    AddImageFilterType::Pointer                     addCVPVFilter;
    AddImageFilterType::Pointer                     addCVPVSinFilter;
    MaskImageFilterType::Pointer                    maskCellImageFilter;
    ThresholdLScalarFilterType::Pointer             cellImageToBin;
    ImageToShapeLabelMapFilterType::Pointer         cellImageToShapeLabelMap;
    ShapeOpeningLabelMapFilterType::Pointer         shapeOpeningLabMapFilter2;
    ImageToShapeLabelMapFilterType::Pointer         nucleiImageToShapeLabelMap;
    LabelMapToLabelImageFilterType::Pointer         nucleiLabelMapToImageFilter;
    LabelMapToLabelImageFilterType::Pointer         labelMapToImage2;
    ThresholdLScalarFilterType::Pointer             thresholdLabelImageFilter2;
    LabelOverlayImageFilterType::Pointer            overlayImageFilter1;
    LabelOverlayImageFilterType::Pointer            overlayImageFilter2;

    //----------READER---------------------------------------------------------------------------------------------------------
    itk::SmartPointer<CScalarVoImageType> readerBCatImage = CScalarVoImageType::New();
    this->ReadImage(m_pathChannelBCat + m_filenameChannelBCat + m_fileExtensionChannelBCat, readerBCatImage, m_spacing);

    readerSin = ScalarVoReaderType::New();
    readerSin->SetFileName(m_fullFilenameSegmentationSin.toStdString());
    readerSin->ReleaseDataBeforeUpdateFlagOn();
    readerSin->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerSin->SetImageIO( itk::TIFFImageIO::New() );
#endif

    itk::SmartPointer<CScalarVoImageType> nucImage = CScalarVoImageType::New();
    this->ReadImage(m_fullFilenameSegmentationNuc.toStdString(), nucImage, m_spacing);

    readerCV = ScalarVoReaderType::New();
    readerCV->SetFileName(m_fullFilenameSegmentationCV.toStdString());
    readerCV->ReleaseDataBeforeUpdateFlagOn();
    readerCV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerCV->SetImageIO( itk::TIFFImageIO::New() );
#endif

    readerPV = ScalarVoReaderType::New();
    readerPV->SetFileName(m_fullFilenameSegmentationCV.toStdString());
    readerPV->ReleaseDataBeforeUpdateFlagOn();
    readerPV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerPV->SetImageIO( itk::TIFFImageIO::New() );
#endif
    //-------------------------------------------------------------------------------------------------------------------------
//    StructuringElementType greyOpeningKernel;
//    greyOpeningKernel.SetRadius(m_greyscaleOpeningRadius);
//    greyOpeningKernel.CreateStructuringElement();
//
//    GrayscaleErodeImageFilterType::Pointer greyErodeFilter = GrayscaleErodeImageFilterType::New();
//    greyErodeFilter->SetInput(medianFilter->GetOutput());
//    greyErodeFilter->SetKernel(greyOpeningKernel);
//
//    if(m_saveEverything) {
//        ScalarVoWriterType::Pointer writerErode = ScalarVoWriterType::New();
//        writerErode->ReleaseDataFlagOn();
//        writerErode->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step0b_erode" + m_fileExtensionChannelDAPI);
//        writerErode->SetInput(greyErodeFilter->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//        writerErode->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//        writerErode->Update();
//    }
//
//    GrayscaleDilateImageFilterType::Pointer dilateFilter = GrayscaleDilateImageFilterType::New();
//    dilateFilter->SetInput(greyErodeFilter->GetOutput());
//    dilateFilter->SetKernel(greyOpeningKernel);
//
//    if(m_saveEverything) {
//        ScalarVoWriterType::Pointer writerDilate = ScalarVoWriterType::New();
//        writerDilate->ReleaseDataFlagOn();
//        writerDilate->SetFileName(m_pathChannelDAPI + m_filenameSave + "_step0c_dilate" + m_fileExtensionChannelDAPI);
//        writerDilate->SetInput(dilateFilter->GetOutput());
//#if (ITK_VERSION_MAJOR >= 4)
//        writerDilate->SetImageIO( itk::TIFFImageIO::New() );
//#endif
//        writerDilate->Update();
//    }

    //----------FILTER-THRESHOLD-ON-DAPI-CHANNEL-------------------------------------------------------------------------------
    switch(m_thresholdingMode)
    {
    case 0:
    {
        adapOtsuFilter = AdaptiveOtsuThresholdImageFilterType::New();
        adapOtsuFilter->SetInput(readerBCatImage);
        adapOtsuFilter->SetInsideValue(255);
        adapOtsuFilter->SetOutsideValue(0);
        adapOtsuFilter->SetNumberOfHistogramBins(256);
        adapOtsuFilter->SetSplineOrder(3);
        adapOtsuFilter->SetNumberOfControlPoints(5);
        adapOtsuFilter->SetNumberOfLevels(3);
        adapOtsuFilter->SetNumberOfSamples(m_adapOtsuSamplePoints);
        adapOtsuFilter->SetRadius(m_adapOtsuRadius);
        break;
    }
    case 1:
    {
        otsuFilter = OtsuThresholdImageFilterType::New();
        otsuFilter->SetInput(readerBCatImage);
        otsuFilter->Update();
        m_otsuThreshold = (CScalarPixelType)(otsuFilter->GetThreshold());

        invertFilter = InvertIntensityImageFilterType::New();
        invertFilter->SetInput(otsuFilter->GetOutput());
        invertFilter->SetMaximum(255);
        break;
    }
    default:
    {
        thresBCatFilter = ThresholdFilterType::New();
        thresBCatFilter->SetOutsideValue(0);
        thresBCatFilter->SetInsideValue(255);
        thresBCatFilter->SetLowerThreshold(m_lowerThreshold);
        thresBCatFilter->SetUpperThreshold(m_upperThreshold);
        thresBCatFilter->SetInput(readerBCatImage);
        break;
    }
    }

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer1 = ScalarVoWriterType::New();
        writer1->ReleaseDataFlagOn();
        writer1->SetFileName(m_pathChannelBCat + m_filenameSave + "_step1_bin" + m_fileExtensionChannelBCat);
        switch(m_thresholdingMode)
        {
        case 0:
            writer1->SetInput(adapOtsuFilter->GetOutput());
            break;
        case 1:
            writer1->SetInput(invertFilter->GetOutput());
            break;
        default:
            writer1->SetInput(thresBCatFilter->GetOutput());
            break;
        }
#if (ITK_VERSION_MAJOR >= 4)
        writer1->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer1->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    holeFillingFilter = HoleFillingImageFilterType::New();
    holeFillingFilter->SetRadius(m_holeFillingNeighborhoodRadius);
    holeFillingFilter->SetBackgroundValue(0);
    holeFillingFilter->SetForegroundValue(255);
    holeFillingFilter->SetMajorityThreshold(m_holeFillingMajThreshold);
    holeFillingFilter->SetMaximumNumberOfIterations(20);
    switch(m_thresholdingMode)
    {
    case 0:
        holeFillingFilter->SetInput(adapOtsuFilter->GetOutput());
        break;
    case 1:
        holeFillingFilter->SetInput(invertFilter->GetOutput());
        break;
    default:
        holeFillingFilter->SetInput(thresBCatFilter->GetOutput());
        break;
    }
    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writerHF = ScalarVoWriterType::New();
        writerHF->ReleaseDataFlagOn();
        writerHF->SetFileName(m_pathChannelBCat + m_filenameSave + "_step2_hole" + m_fileExtensionChannelBCat);
        writerHF->SetInput(holeFillingFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerHF->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerHF->Update();
    }

    invHoleFillingFilter = HoleFillingImageFilterType::New();
    invHoleFillingFilter->SetRadius(m_invHoleFillingNeighborhoodRadius);
    invHoleFillingFilter->SetBackgroundValue(255);
    invHoleFillingFilter->SetForegroundValue(0);
    invHoleFillingFilter->SetMajorityThreshold(m_invHoleFillingMajThreshold);
    invHoleFillingFilter->SetMaximumNumberOfIterations(20);
    invHoleFillingFilter->SetInput(holeFillingFilter->GetOutput());

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writerIHF = ScalarVoWriterType::New();
        writerIHF->ReleaseDataFlagOn();
        writerIHF->SetFileName(m_pathChannelBCat + m_filenameSave + "_step3_invHole" + m_fileExtensionChannelBCat);
        writerIHF->SetInput(invHoleFillingFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writerIHF->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writerIHF->Update();
    }

    //----------FILTER-CLOSING-ON-DAPI-CHANNEL---------------------------------------------------------------------------------
    closingFilter = ClosingImageFilterType::New();
    closingFilter->SetKernel(m_closingStructuringElement);
    closingFilter->SetForegroundValue(255);
    closingFilter->SetInput(invHoleFillingFilter->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-OPENING-ON-DAPI-CHANNEL---------------------------------------------------------------------------------
    openingFilter1 = OpeningImageFilterType::New();
    openingFilter1->SetKernel(m_openingStructuringElement1);
    openingFilter1->SetBackgroundValue(0);
    openingFilter1->SetForegroundValue(255);
    openingFilter1->SetInput(closingFilter->GetOutput());

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer3 = ScalarVoWriterType::New();
        writer3->ReleaseDataFlagOn();
        writer3->SetFileName(m_pathChannelBCat + m_filenameSave + "_step4_clopening" + m_fileExtensionChannelBCat);
        writer3->SetInput(openingFilter1->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer3->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer3->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------RELABEL-WATERSHED-LABEL-OBJECTS-------------------------------------------------------------------------------
    imageToShapeLabelMap = ImageToShapeLabelMapFilterType::New();
    imageToShapeLabelMap->SetFullyConnected(false);
    imageToShapeLabelMap->SetInput(openingFilter1->GetOutput());
    //-------------------------------------------------------------------------------------------------------------------------

    //----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------------
    shapeOpeningLabMapFilter1 = ShapeOpeningLabelMapFilterType::New();
    shapeOpeningLabMapFilter1->SetLambda(100);                                           //attribute value
    shapeOpeningLabMapFilter1->ReverseOrderingOff();                                    //removes objects with attribute smaller than lambda
    shapeOpeningLabMapFilter1->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::NUMBER_OF_PIXELS);
    shapeOpeningLabMapFilter1->SetInput(imageToShapeLabelMap->GetOutput());

    labelMapToImage1 = LabelMapToLabelImageFilterType::New();
    labelMapToImage1->SetInput(shapeOpeningLabMapFilter1->GetOutput());

    thresholdLabelImageFilter1 = ThresholdLScalarFilterType::New();
    thresholdLabelImageFilter1->SetOutsideValue(255);
    thresholdLabelImageFilter1->SetInsideValue(0);
    thresholdLabelImageFilter1->SetLowerThreshold(1);
    thresholdLabelImageFilter1->SetUpperThreshold(255);
    thresholdLabelImageFilter1->SetInput(labelMapToImage1->GetOutput());
    thresholdLabelImageFilter1->Update();

    itk::SmartPointer<CScalarVoImageType> borderImage = thresholdLabelImageFilter1->GetOutput();
    borderImage->DisconnectPipeline();

    borderImage->SetSpacing(m_spacing);

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer4 = ScalarVoWriterType::New();
        writer4->ReleaseDataFlagOn();
        writer4->SetFileName(m_pathChannelBCat + m_filenameSave + "_step5_remIsolPxl" + m_fileExtensionChannelBCat);
        writer4->SetInput(borderImage);
#if (ITK_VERSION_MAJOR >= 4)
        writer4->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer4->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------
//
//    //----------INVERT--------------------------------------------------------------------------------------------------------
//    invertForWatershedFilter = InvertIntensityImageFilterType::New();
//    invertForWatershedFilter->SetInput(openingFilter->GetOutput());
//    invertForWatershedFilter->SetMaximum(255);
//    //-------------------------------------------------------------------------------------------------------------------------

    //----------DISTANCE-MAP---------------------------------------------------------------------------------------------------
    distanceMap = SignedMaurerDistanceMapImageFilterType::New();
    distanceMap->ReleaseDataFlagOn();
    distanceMap->SquaredDistanceOff();
    distanceMap->UseImageSpacingOn();
    distanceMap->SetInput(borderImage);
    distanceMap->Update();

    itk::SmartPointer<FScalarVoImageType> distMapImage = distanceMap->GetOutput();
    distMapImage->DisconnectPipeline();

    std::cout << "Write distMap intermediate results " << m_saveEverything << std::endl;
    if(m_saveEverything) {
        rescaler = RescaleImageFilterType::New();
        rescaler->ReleaseDataFlagOn();
        rescaler->SetInput(distMapImage);

        ScalarVoWriterType::Pointer writer2 = ScalarVoWriterType::New();
        writer2->ReleaseDataFlagOn();
        writer2->SetFileName(m_pathChannelBCat + m_filenameSave + "_step6_distMap" + m_fileExtensionChannelBCat);
        writer2->SetInput(rescaler->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
        writer2->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer2->Update();
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------WATERSHED-FILTER-----------------------------------------------------------------------------------------------
    morphWatershed = MorphoWatershedImageFilterType::New();
    morphWatershed->SetLevel(m_floodLevel);
    morphWatershed->FullyConnectedOff();
    morphWatershed->ReleaseDataFlagOn();
    morphWatershed->SetInput(distMapImage);
    //-------------------------------------------------------------------------------------------------------------------------

    //----------MASK-VEINS-AND-SINUSOIDS---------------------------------------------------------------------------------------
    openingFilter2 = OpeningImageFilterType::New();
    openingFilter2->SetKernel(m_openingStructuringElement2);
    openingFilter2->SetBackgroundValue(0);
    openingFilter2->SetForegroundValue(255);
    openingFilter2->SetInput(readerSin->GetOutput());

    addCVPVFilter = AddImageFilterType::New();
    addCVPVFilter->SetInput1(readerCV->GetOutput());
    addCVPVFilter->SetInput2(readerPV->GetOutput());

    addCVPVSinFilter = AddImageFilterType::New();
    addCVPVSinFilter->SetInput1(addCVPVFilter->GetOutput());
    addCVPVSinFilter->SetInput2(readerSin->GetOutput());
    addCVPVSinFilter->Update();

    itk::SmartPointer<CScalarVoImageType> cvpvsinImage = addCVPVSinFilter->GetOutput();
    cvpvsinImage->DisconnectPipeline();

    cvpvsinImage->SetSpacing(m_spacing);

    if(m_saveEverything) {
        ScalarVoWriterType::Pointer writer5 = ScalarVoWriterType::New();
        writer5->ReleaseDataFlagOn();
        writer5->SetFileName(m_pathChannelBCat + m_filenameSave + "_step6_cvpvsinMask" + m_fileExtensionChannelBCat);
        writer5->SetInput(cvpvsinImage);
#if (ITK_VERSION_MAJOR >= 4)
        writer5->SetImageIO( itk::TIFFImageIO::New() );
#endif
        writer5->Update();
    }

    maskCellImageFilter = MaskImageFilterType::New();
    maskCellImageFilter->SetInput(morphWatershed->GetOutput());
    maskCellImageFilter->SetMaskImage(cvpvsinImage);
    //------------------------------------------------------------------------------------------------------------------------

    //----------RELABEL-WATERSHED-LABEL-OBJECTS-------------------------------------------------------------------------------
    cellImageToBin = ThresholdLScalarFilterType::New();
    cellImageToBin->SetOutsideValue(0);
    cellImageToBin->SetInsideValue(itk::NumericTraits<LScalarPixelType>::max());
    cellImageToBin->SetLowerThreshold(1);
    cellImageToBin->SetUpperThreshold(itk::NumericTraits<LabelMapToLabelImageFilterType::OutputImagePixelType>::max());
    cellImageToBin->SetInput(maskCellImageFilter->GetOutput());
    cellImageToBin->Update();

    itk::SmartPointer<CScalarVoImageType> cellBinImage = cellImageToBin->GetOutput();
    cellBinImage->DisconnectPipeline();

    cellBinImage->SetSpacing(m_spacing);

    cellImageToShapeLabelMap = ImageToShapeLabelMapFilterType::New();
    cellImageToShapeLabelMap->SetFullyConnected(false);
    cellImageToShapeLabelMap->SetInput(cellBinImage);
    cellImageToShapeLabelMap->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------REMOVE-SMALL-OBJECTS-------------------------------------------------------------------------------------------
    std::cout << "remove objects with volume smaller than " << m_minimalCellVolume << std::endl;

    std::cout << "before volume based removal " << cellImageToShapeLabelMap->GetOutput()->GetNumberOfLabelObjects() << " objects" << std::endl;

    shapeOpeningLabMapFilter2 = ShapeOpeningLabelMapFilterType::New();
    shapeOpeningLabMapFilter2->SetLambda(m_minimalCellVolume);                                           //attribute value
    shapeOpeningLabMapFilter2->ReverseOrderingOff();                                                       //removes objects with attribute smaller than lambda
    shapeOpeningLabMapFilter2->SetAttribute(ShapeOpeningLabelMapFilterType::LabelObjectType::PHYSICAL_SIZE);
    shapeOpeningLabMapFilter2->SetInput(cellImageToShapeLabelMap->GetOutput());
    shapeOpeningLabMapFilter2->Update();

    std::cout << "after volume based removal " << shapeOpeningLabMapFilter2->GetOutput()->GetNumberOfLabelObjects() << " objects" << std::endl;

    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap = shapeOpeningLabMapFilter2->GetOutput();
    cellLabelMap->DisconnectPipeline();
    std::cout << "after volume based removal " << cellLabelMap->GetNumberOfLabelObjects() << " objects" << std::endl;

    shapeOpeningLabMapFilter2->Update();
    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithOneNucleusLabelMap = shapeOpeningLabMapFilter2->GetOutput();
    cellsWithOneNucleusLabelMap->DisconnectPipeline();

    shapeOpeningLabMapFilter2->Update();
    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithTwoNucleiLabelMap = shapeOpeningLabMapFilter2->GetOutput();
    cellsWithTwoNucleiLabelMap->DisconnectPipeline();

    shapeOpeningLabMapFilter2->Update();
    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellsWithMoreNucleiLabelMap = shapeOpeningLabMapFilter2->GetOutput();
    cellsWithMoreNucleiLabelMap->DisconnectPipeline();

    labelMapToImage2 = LabelMapToLabelImageFilterType::New();
    labelMapToImage2->SetInput(shapeOpeningLabMapFilter2->GetOutput());
    labelMapToImage2->Update();

    nucleiImageToShapeLabelMap = ImageToShapeLabelMapFilterType::New();
    nucleiImageToShapeLabelMap->SetFullyConnected(false);
    nucleiImageToShapeLabelMap->SetInput(nucImage);
    nucleiImageToShapeLabelMap->Update();

    nucleiLabelMapToImageFilter = LabelMapToLabelImageFilterType::New();
    nucleiLabelMapToImageFilter->SetInput(nucleiImageToShapeLabelMap->GetOutput());
    nucleiLabelMapToImageFilter->Update();

    BuildCellNucleiAlignmentMaps(nucleiLabelMapToImageFilter->GetOutput(), labelMapToImage2->GetOutput(), nucleiImageToShapeLabelMap->GetOutput(), cellLabelMap);

    std::cout << "before removal " << cellLabelMap->GetNumberOfLabelObjects() << " cells" << std::endl;
    for(unsigned int i=0; i<m_labelsWithZeroNuclei.size(); i++)
        cellLabelMap->RemoveLabel(m_labelsWithZeroNuclei[i]);
    std::cout << "after removal of cells without a nucleus " << cellLabelMap->GetNumberOfLabelObjects() << " cells left" << std::endl;

    RemoveCellsAtDatasetBorder(labelMapToImage2->GetOutput(), cellsWithOneNucleusLabelMap, cellsWithTwoNucleiLabelMap, cellsWithMoreNucleiLabelMap);
    WriteNumNucleiOutlineFiles(readerBCatImage, cellsWithOneNucleusLabelMap, cellsWithTwoNucleiLabelMap, cellsWithMoreNucleiLabelMap); //keep near position filter not applicable to this atm
    //-------------------------------------------------------------------------------------------------------------------------

    //----------FILTER-THRESHOLD-ON-SEGMENTED-NETWORK--------------------------------------------------------------------------
    thresholdLabelImageFilter2 = ThresholdLScalarFilterType::New();
    thresholdLabelImageFilter2->SetOutsideValue(0);
    thresholdLabelImageFilter2->SetInsideValue(itk::NumericTraits<CScalarPixelType>::max());
    thresholdLabelImageFilter2->SetLowerThreshold(1);
    thresholdLabelImageFilter2->SetUpperThreshold(itk::NumericTraits<LabelMapToLabelImageFilterType::OutputImagePixelType>::max());
    thresholdLabelImageFilter2->SetInput(labelMapToImage2->GetOutput());

    ScalarVoWriterType::Pointer writer6 = ScalarVoWriterType::New();
    writer6->ReleaseDataFlagOn();
    writer6->SetFileName(m_pathChannelBCat + m_filenameSave + m_saveSuffixesForFinals[0] + m_fileExtensionChannelBCat);
    writer6->SetInput(thresholdLabelImageFilter2->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer6->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer6->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------VISULALIZATION-ONLY-FILTER---NOT-NEEDED-BY-OTHER-PIPELINE-FILTERS----------------------------------------------
    overlayImageFilter1 = LabelOverlayImageFilterType::New();
    overlayImageFilter1->SetOpacity(m_overlayOpacity);
    overlayImageFilter1->ReleaseDataFlagOn();
    overlayImageFilter1->SetInput(readerBCatImage);
    overlayImageFilter1->SetLabelImage(labelMapToImage2->GetOutput());

    RGBVoWriterType::Pointer writer7 = RGBVoWriterType::New();
    writer7->ReleaseDataFlagOn();
    writer7->SetFileName(m_pathChannelBCat + m_filenameSave + m_saveSuffixesForFinals[1] + m_fileExtensionChannelBCat);
    writer7->SetInput(overlayImageFilter1->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer7->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer7->Update();

    overlayImageFilter2 = LabelOverlayImageFilterType::New();
    overlayImageFilter2->SetOpacity(m_overlayOpacity);
    overlayImageFilter2->ReleaseDataFlagOn();
    overlayImageFilter2->SetInput(nucImage);
    overlayImageFilter2->SetLabelImage(labelMapToImage2->GetOutput());

    RGBVoWriterType::Pointer writer8 = RGBVoWriterType::New();
    writer8->ReleaseDataFlagOn();
    writer8->SetFileName(m_pathChannelBCat + m_filenameSave + m_saveSuffixesForFinals[2] + m_fileExtensionChannelBCat);
    writer8->SetInput(overlayImageFilter2->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer8->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer8->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    WriteLogFile(timeStamp);
    WriteDataSetSummary();
}
