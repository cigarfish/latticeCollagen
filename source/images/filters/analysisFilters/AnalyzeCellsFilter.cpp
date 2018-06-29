///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeCellsFilter.cpp                                               //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-05-03                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "AnalyzeCellsFilter.h"

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

#include <vtkDataSetAttributes.h>
#include <vtkDelimitedTextReader.h>
#include <vtkDelimitedTextWriter.h>
#include <vtkGraphWriter.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnsignedLongArray.h>
#include <vtkVariantArray.h>

#include "../../pipelines/BasePipeline.h"
#include "../../tools/GraphAnnotationHelper.h"
#include "../../../tools/parameters/CSParameter.h"
#include "../../../tools/parameters/CSParameterContext.h"


AnalyzeCellsFilter::AnalyzeCellsFilter()
{

}


AnalyzeCellsFilter::~AnalyzeCellsFilter()
{
    // TODO Auto-generated destructor stub
}


void AnalyzeCellsFilter::ParseParameterContext()
{
    if(mpParamContext->findContext("Analyze Cells",0)==NULL) {
        std::cout << "Error: AnalyzeCellsFilter: Invalid parameter context" << std::endl;
        return;
    }

    mDataSetID = *(std::string*)(mpParamContext->findParameter("Dataset ID", 0)->dataPointer());

    mDataSetFullFilenameCell = QString::fromStdString( *(std::string*)(mpParamContext->findParameter("Cell shape binary image", 0)->dataPointer()) );
	mInfoDataSetFullFilenameCell.setFile(mDataSetFullFilenameCell);
	if(!mInfoDataSetFullFilenameCell.exists())
		throw std::string("Please specify cell segmentation file");
	mDataSetPathCell = (mInfoDataSetFullFilenameCell.path() + QString("/")).toStdString();
    mDataSetNameCell = mInfoDataSetFullFilenameCell.baseName().toStdString();
    mDataSetFileExtensionCell = (QString(".") + mInfoDataSetFullFilenameCell.suffix()).toStdString();

    mDataSetFullFilenameNuclei = QString::fromStdString( *(std::string*)(mpParamContext->findParameter("Nuclei segmentation", 0)->dataPointer()) );
	mInfoDataSetFullFilenameNuclei.setFile(mDataSetFullFilenameNuclei);
	if(!mInfoDataSetFullFilenameNuclei.exists())
		throw std::string("Please specify nuclei segmentation file");
	mDataSetPathNuclei = (mInfoDataSetFullFilenameNuclei.path() + QString("/")).toStdString();
    mDataSetNameNuclei = mInfoDataSetFullFilenameNuclei.baseName().toStdString();
    mDataSetFileExtensionNuclei = (QString(".") + mInfoDataSetFullFilenameNuclei.suffix()).toStdString();

    mDataSetFullFilenameSinusoid = QString::fromStdString( *(std::string*)(mpParamContext->findParameter("Sinusoid segmentation", 0)->dataPointer()) );
	mInfoDataSetFullFilenameSinusoid.setFile(mDataSetFullFilenameSinusoid);
	if(!mInfoDataSetFullFilenameSinusoid.exists())
		throw std::string("Please specify sinusoid segmentation file");
	mDataSetPathSinusoid = (mInfoDataSetFullFilenameSinusoid.path() + QString("/")).toStdString();
    mDataSetNameSinusoid = mInfoDataSetFullFilenameSinusoid.baseName().toStdString();
    mDataSetFileExtensionSinusoid = (QString(".") + mInfoDataSetFullFilenameSinusoid.suffix()).toStdString();

    mDataSetFullFilenameBile = QString::fromStdString( *(std::string*)(mpParamContext->findParameter("Bile segmentation", 0)->dataPointer()) );
	mInfoDataSetFullFilenameBile.setFile(mDataSetFullFilenameBile);
	if(!mInfoDataSetFullFilenameBile.exists())
		throw std::string("Please specify bile segmentation file");
	mDataSetPathBile = (mInfoDataSetFullFilenameBile.path() + QString("/")).toStdString();
    mDataSetNameBile = mInfoDataSetFullFilenameBile.baseName().toStdString();
    mDataSetFileExtensionBile = (QString(".") + mInfoDataSetFullFilenameBile.suffix()).toStdString();

    mDataSetFullFilenameCV = QString::fromStdString( *(std::string*)(mpParamContext->findParameter("CV segmentation", 0)->dataPointer()) );
	mInfoDataSetFullFilenameCV.setFile(mDataSetFullFilenameCV);
	mHasCV = mInfoDataSetFullFilenameCV.exists();
	if(mHasCV) {
		mDataSetPathCV = (mInfoDataSetFullFilenameCV.path() + QString("/")).toStdString();
		mDataSetNameCV = mInfoDataSetFullFilenameCV.baseName().toStdString();
		mDataSetFileExtensionCV = (QString(".") + mInfoDataSetFullFilenameCV.suffix()).toStdString();
	}

    mDataSetFullFilenamePV = QString::fromStdString( *(std::string*)(mpParamContext->findParameter("PV segmentation", 0)->dataPointer()) );
	mInfoDataSetFullFilenamePV.setFile(mDataSetFullFilenamePV);
	mHasPV = mInfoDataSetFullFilenamePV.exists();
	if(mHasPV) {
		mDataSetPathPV = (mInfoDataSetFullFilenamePV.path() + QString("/")).toStdString();
		mDataSetNamePV = mInfoDataSetFullFilenamePV.baseName().toStdString();
		mDataSetFileExtensionPV = (QString(".") + mInfoDataSetFullFilenamePV.suffix()).toStdString();
	}

    mSpacing[0] = *(double*)(mpParamContext->findParameter("Voxel spacing x", 0)->dataPointer());
    mSpacing[1] = *(double*)(mpParamContext->findParameter("Voxel spacing y", 0)->dataPointer());
    mSpacing[2] = *(double*)(mpParamContext->findParameter("Voxel spacing z", 0)->dataPointer());

    mSaveAsGraph = *(bool*)(this->mpParamContext->findParameter("Save analysis result as graph", 0)->dataPointer());
}


void AnalyzeCellsFilter::PopulateCellAnalysisContainer()
{
    std::string timeStamp;
	
    std::cout << "Start contact analysis: " << std::endl;

    //----------READER---------------------------------------------------------------------------------------------------------
    if(BasePipeline<3>::GetNumberOfDimensions(mDataSetPathCell+ mDataSetNameCell + mDataSetFileExtensionCell) != 3)
        throw std::string("Please specify a cell segmentation file with three dimensions.");

    ScalarVoReaderType::Pointer cellReader = ScalarVoReaderType::New();
    cellReader->SetFileName(mDataSetPathCell+ mDataSetNameCell + mDataSetFileExtensionCell);
#if (ITK_VERSION_MAJOR >= 4)
    cellReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    cellReader->Update();

    if(BasePipeline<3>::GetNumberOfDimensions(mDataSetPathNuclei+ mDataSetNameNuclei + mDataSetFileExtensionNuclei) != 3)
        throw std::string("Please specify a nuclei segmentation file with three dimensions.");

    ScalarVoReaderType::Pointer nucleiReader = ScalarVoReaderType::New();
    nucleiReader->SetFileName(mDataSetPathNuclei+ mDataSetNameNuclei + mDataSetFileExtensionNuclei);
#if (ITK_VERSION_MAJOR >= 4)
    nucleiReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    nucleiReader->Update();

    if(BasePipeline<3>::GetNumberOfDimensions(mDataSetPathSinusoid + mDataSetNameSinusoid + mDataSetFileExtensionSinusoid) != 3)
        throw std::string("Please specify a sinusoid segmentation file with three dimensions.");

    ScalarVoReaderType::Pointer sinusoidReader = ScalarVoReaderType::New();
    sinusoidReader->SetFileName(mDataSetPathSinusoid + mDataSetNameSinusoid + mDataSetFileExtensionSinusoid);
#if (ITK_VERSION_MAJOR >= 4)
    sinusoidReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    sinusoidReader->Update();

    if(BasePipeline<3>::GetNumberOfDimensions(mDataSetPathBile + mDataSetNameBile + mDataSetFileExtensionBile) != 3)
        throw std::string("Please specify a bile segmentation file with three dimensions.");

    ScalarVoReaderType::Pointer bileReader = ScalarVoReaderType::New();
    bileReader->SetFileName(mDataSetPathBile + mDataSetNameBile + mDataSetFileExtensionBile);
#if (ITK_VERSION_MAJOR >= 4)
    bileReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
    bileReader->Update();

    ScalarVoReaderType::Pointer CVReader = ScalarVoReaderType::New();
    if(mHasCV) {
        if(BasePipeline<3>::GetNumberOfDimensions(mDataSetPathCV + mDataSetNameCV + mDataSetFileExtensionCV) != 3)
            throw std::string("Please specify a central vein segmentation file with three dimensions.");

        CVReader->SetFileName(mDataSetPathCV + mDataSetNameCV + mDataSetFileExtensionCV);
#if (ITK_VERSION_MAJOR >= 4)
        CVReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
        CVReader->Update();
    }

    ScalarVoReaderType::Pointer PVReader = ScalarVoReaderType::New();
    if(mHasPV) {
        if(BasePipeline<3>::GetNumberOfDimensions(mDataSetPathPV + mDataSetNamePV + mDataSetFileExtensionPV) != 3)
            throw std::string("Please specify a portal vein segmentation file with three dimensions.");

        PVReader->SetFileName(mDataSetPathPV + mDataSetNamePV + mDataSetFileExtensionPV);
#if (ITK_VERSION_MAJOR >= 4)
        PVReader->SetImageIO( itk::TIFFImageIO::New() );
#endif
        PVReader->Update();
    }

    CScalarVoImageType::Pointer cellImage = cellReader->GetOutput();
    cellImage->DisconnectPipeline();
    cellImage->SetSpacing(mSpacing);

    CScalarVoImageType::Pointer nucleiImage = nucleiReader->GetOutput();
    nucleiImage->DisconnectPipeline();
    nucleiImage->SetSpacing(mSpacing);

    CScalarVoImageType::Pointer sinusoidImage = sinusoidReader->GetOutput();
    sinusoidImage->DisconnectPipeline();
    sinusoidImage->SetSpacing(mSpacing);

    CScalarVoImageType::Pointer bileImage = bileReader->GetOutput();
    bileImage->DisconnectPipeline();
    bileImage->SetSpacing(mSpacing);

    CScalarVoImageType::Pointer CVImage = CScalarVoImageType::New();
    if(mHasCV) {
        CVImage = CVReader->GetOutput();
        CVImage->DisconnectPipeline();
    }
    else {
        CScalarVoImageType::RegionType region(cellImage->GetLargestPossibleRegion().GetIndex(), cellImage->GetLargestPossibleRegion().GetSize());

        CVImage->SetRegions(region);
        CVImage->Allocate();
        CVImage->FillBuffer(0);
    }
    CVImage->SetSpacing(mSpacing);

    CScalarVoImageType::Pointer PVImage = CScalarVoImageType::New();
    if(mHasPV) {
        PVImage = PVReader->GetOutput();
        PVImage->DisconnectPipeline();
    }
    else {
        CScalarVoImageType::RegionType region(cellImage->GetLargestPossibleRegion().GetIndex(), cellImage->GetLargestPossibleRegion().GetSize());

        PVImage->SetRegions(region);
        PVImage->Allocate();
        PVImage->FillBuffer(0);
    }
    PVImage->SetSpacing(mSpacing);
    //-------------------------------------------------------------------------------------------------------------------------

    //----------COMPUTE-CELL-VOLUME-AND-BOUNDARY-CONTACT-AREA------------------------------------------------------------------
    mVoxelVolume = mSpacing[0]*mSpacing[1]*mSpacing[2];

    ImageToShapeLabelMapFilterType::Pointer cellImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    cellImageToShaLabMapFilter->SetInput(cellImage);
    cellImageToShaLabMapFilter->SetFullyConnected(false);
    cellImageToShaLabMapFilter->Update();

    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> cellLabelMap = cellImageToShaLabMapFilter->GetOutput();
    cellLabelMap->DisconnectPipeline();

    std::cout << "cellLabelMap->GetNumberOfLabelObjects() = " << cellLabelMap->GetNumberOfLabelObjects() << std::endl;

    for(unsigned int i=0; i<cellLabelMap->GetNumberOfLabelObjects(); ++i) {
        CellAnalysisContainer cell;
        cell.mLabel = cellLabelMap->GetNthLabelObject(i)->GetLabel();
        PointType poi = cellLabelMap->GetNthLabelObject(i)->GetCentroid();
        cell.mX = poi[0];
        cell.mY = poi[1];
        cell.mZ = poi[2];

        cell.mVolume = cellLabelMap->GetNthLabelObject(i)->GetNumberOfPixels() * mVoxelVolume;

        if(cell.mLabel != 0)                                                                                //Label #0 is background label
            mCells.insert(std::pair<unsigned long int, CellAnalysisContainer>(cell.mLabel, cell) );
    }
    std::cout << "mCells.size = " << mCells.size() << std::endl;

    LabelMapToLabelImageFilterType::Pointer cellLabelMapToImage = LabelMapToLabelImageFilterType::New();
    cellLabelMapToImage->SetInput(cellLabelMap);
    cellLabelMapToImage->Update();

    LScalarVoImageType::Pointer cellLabelImage = cellLabelMapToImage->GetOutput();

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

    double voxelSurface[3];
    voxelSurface[0] = mSpacing[0]*mSpacing[1];        //z-plane voxel surface
    voxelSurface[1] = mSpacing[0]*mSpacing[2];        //y-plane voxel surface
    voxelSurface[2] = mSpacing[1]*mSpacing[2];        //x-plane voxel surface

    for(unsigned int i=0; i<6; i++) {
        itk::ImageRegionConstIterator<LScalarVoImageType> iter(cellLabelImage, dataSetFaces[i]);

        for(iter = iter.Begin(); !iter.IsAtEnd(); ++iter) {
            if(iter.Value() != 0)
                mCells[iter.Value()].mContactAreaWithDataSetBorder += voxelSurface[i/2];
//                mCells[iter.Value()].mContactAreaWithDataSetBorder++;
        }
    }

    //-------------------------------------------------------------------------------------------------------------------------

    //----------COMPUTE-CELL-CELL-CELL-BILE-AND-CELL-SINUSOID-CONTACT-AREA-----------------------------------------------------
    CScalarVoImageType::Pointer contactImage = CScalarVoImageType::New();
    contactImage->SetRegions(cellImage->GetLargestPossibleRegion());
    contactImage->Allocate();
    contactImage->FillBuffer(0);

    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);

    BoundaryConditionType boundaryCondition;
    NeighborhoodIteratorType it(radius, cellLabelImage, cellLabelImage->GetLargestPossibleRegion());
    it.SetBoundaryCondition(boundaryCondition);

	const unsigned int neighborhoodSize = 6;

    NeighborhoodIteratorType::OffsetType neighbor[neighborhoodSize];
    neighbor[0][0] = 1;     neighbor[1][0] = -1;
    neighbor[0][1] = 0;     neighbor[1][1] = 0;
    neighbor[0][2] = 0;     neighbor[1][2] = 0;

    neighbor[2][0] = 0;     neighbor[3][0] = 0;
    neighbor[2][1] = 1;     neighbor[3][1] = -1;
    neighbor[2][2] = 0;     neighbor[3][2] = 0;

    neighbor[4][0] = 0;     neighbor[5][0] = 0;
    neighbor[4][1] = 0;     neighbor[5][1] = 0;
    neighbor[4][2] = 1;     neighbor[5][2] = -1;

    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        if(it.GetCenterPixel() == 0) {
            std::map<long, double> cellHits;        //cell label, contact area

            for(unsigned int i=0; i<neighborhoodSize; i++) {
                if( cellLabelImage->GetLargestPossibleRegion().IsInside(it.GetIndex()+neighbor[i]) ) {
                    long p = cellLabelImage->GetPixel(it.GetIndex()+neighbor[i]);

                    if(p!=0) {
                        if(cellHits.count(p)==0) {
                            if(i==0 || i==1)
                                cellHits.insert(std::pair<long, double>(p,voxelSurface[2]));
                            else if(i==2 || i==3)
                                cellHits.insert(std::pair<long, double>(p,voxelSurface[1]));
                            else if(i==4 || i==5)
                                cellHits.insert(std::pair<long, double>(p,voxelSurface[0]));
                        }
                        else {
                            if(i==0 || i==1)
                                cellHits[p] += voxelSurface[2];
                            else if(i==2 || i==3)
                                cellHits[p] += voxelSurface[1];
                            else if(i==4 || i==5)
                                cellHits[p] += voxelSurface[0];
                        }
                    }
                }
            }
            if(cellHits.size()>1)
                contactImage->SetPixel(it.GetIndex(), 95);

            bool isSinusoid, isBile;

            if(sinusoidImage->GetPixel(it.GetIndex()) != 0 || CVImage->GetPixel(it.GetIndex()) != 0 || PVImage->GetPixel(it.GetIndex()) != 0) {
                isSinusoid = true;
                if(cellHits.size()>=1)
                    contactImage->SetPixel(it.GetIndex(), 175);
            }
            else
                isSinusoid = false;

            if(bileImage->GetPixel(it.GetIndex()) != 0) {
                isBile = true;
                if(cellHits.size()>=1)
                    contactImage->SetPixel(it.GetIndex(), 255);
            }
            else
                isBile = false;

            std::map<long, double>::iterator cHIter;
            for(cHIter = cellHits.begin(); cHIter != cellHits.end(); ++cHIter)
            {
                if(isBile)
                    mCells[cHIter->first].mContactAreaWithBile += cHIter->second;
                else if(isSinusoid)
                    mCells[cHIter->first].mContactAreaWithSinusoids += cHIter->second;
                else if(cellHits.size()>1)
                    mCells[cHIter->first].mContactAreaWithCells += cHIter->second;
                else
                    mCells[cHIter->first].mContactAreaNil += cHIter->second;
            }
        }
    }

    std::map<unsigned long int, CellAnalysisContainer>::iterator cellIter;
    for(cellIter = mCells.begin(); cellIter != mCells.end(); ++cellIter) {
//        double overallSurface = cellIter->second.mContactAreaWithCells + cellIter->second.mContactAreaWithSinusoids + cellIter->second.mContactAreaWithBile + cellIter->second.mContactAreaWithDataSetBorder
//                + cellIter->second.mContactAreaNil;
        double overallSurface = cellIter->second.mContactAreaWithCells + cellIter->second.mContactAreaWithSinusoids + cellIter->second.mContactAreaWithBile + cellIter->second.mContactAreaNil;

        cellIter->second.mContactAreaWithCellsRatio = cellIter->second.mContactAreaWithCells / overallSurface;
        cellIter->second.mContactAreaWithDataSetBorderRatio = 0;
        cellIter->second.mContactAreaWithSinusoidsRatio = cellIter->second.mContactAreaWithSinusoids / overallSurface;
        cellIter->second.mContactAreaWithBileRatio = cellIter->second.mContactAreaWithBile / overallSurface;
        cellIter->second.mContactAreaNilRatio = cellIter->second.mContactAreaNil / overallSurface;
    }

    CScalarVoWriterType::Pointer writer1 = CScalarVoWriterType::New();
    writer1->SetFileName(mDataSetPathCell+ "cellShape_step4_contactAnalysis" + mDataSetFileExtensionCell);
    writer1->SetInput(contactImage);
#if (ITK_VERSION_MAJOR >= 4)
    writer1->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer1->Update();
    //-------------------------------------------------------------------------------------------------------------------------

    //----------CELL-NUCLEI-ASSIGNMENT-----------------------------------------------------------------------------------------
    std::multimap<unsigned long int, std::pair<unsigned long int, double> > cellToNuclei;

    ImageToShapeLabelMapFilterType::Pointer nucleiBinToLabMapFilter = ImageToShapeLabelMapFilterType::New();
    nucleiBinToLabMapFilter->SetFullyConnected(true);
    nucleiBinToLabMapFilter->SetInput(nucleiImage);
    nucleiBinToLabMapFilter->Update();

    itk::SmartPointer<ImageToShapeLabelMapFilterType::OutputImageType> nucleiLabelMap = nucleiBinToLabMapFilter->GetOutput();
    nucleiLabelMap->DisconnectPipeline();

    LabelMapToLabelImageFilterType::Pointer nucleiLabelMapToImageFilter = LabelMapToLabelImageFilterType::New();
    nucleiLabelMapToImageFilter->SetInput(nucleiLabelMap);
    nucleiLabelMapToImageFilter->Update();

    LScalarVoImageType::Pointer nucleiLabelImage = nucleiLabelMapToImageFilter->GetOutput();
    std::cout << "nucleiLabelImage size " << nucleiLabelImage->GetLargestPossibleRegion().GetSize() << std::endl;

    itk::ImageRegionConstIterator<LScalarVoImageType> iterNucl(nucleiLabelImage, nucleiLabelImage->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<LScalarVoImageType> iterCell(cellLabelImage, cellLabelImage->GetLargestPossibleRegion());
    for(iterNucl = iterNucl.Begin(), iterCell = iterCell.Begin(); !iterNucl.IsAtEnd(); ++iterNucl, ++iterCell) {
        if(iterNucl.Value() != 0 && iterCell.Value() != 0) {
            if(cellToNuclei.count(iterCell.Value()) == 0) {
                std::pair<unsigned long int, double> nucleus(iterNucl.Value(), 1.);
                cellToNuclei.insert(std::pair<unsigned long int, std::pair<unsigned long int, double> >(iterCell.Value(), nucleus));
            }
            else {
                std::multimap<unsigned long int, std::pair<unsigned long int, double> >::iterator alignIter;
                alignIter = cellToNuclei.find(iterCell.Value());

                bool alreadyFound = false;
                while(alignIter != cellToNuclei.end() && alignIter->first == iterCell.Value()) {
                    if(alignIter->second.first == iterNucl.Value()) {
                        alignIter->second.second += 1.;
                        alreadyFound = true;
                        break;
                    }
                    ++alignIter;
                }
                if(!alreadyFound) {
                    std::pair<unsigned long int, double> nucleus(iterNucl.Value(), 1.);
                    cellToNuclei.insert(std::pair<unsigned long int, std::pair<unsigned long int, double> >(iterCell.Value(), nucleus));
                }
            }
        }
    }

    std::cout << "Start cell evaluation " << std::endl;
    std::multimap<unsigned long int, std::pair<unsigned long int, double> >::iterator alignIter = cellToNuclei.begin();

    std::cout << "before removal: cellToNuclei.size() = " << cellToNuclei.size() << std::endl;

    while(alignIter != cellToNuclei.end()) {
        alignIter->second.second = alignIter->second.second / (double)(nucleiLabelMap->GetLabelObject(alignIter->second.first)->GetNumberOfPixels());

//        std::cout << "Cell " << alignIter->first << " which starts at " << cellLabelMap->GetLabelObject(alignIter->first)->GetIndex(0) << " contains nucleus " << alignIter->second.first <<
//                " which starts at " << nucleiLabelMap->GetLabelObject(alignIter->second.first)->GetIndex(0) << " to factor " << alignIter->second.second << std::endl;

        if(alignIter->second.second <= 0.5) {
            cellToNuclei.erase(alignIter++);
//            std::cout << "remove the crappy bastard!" << std::endl;
        }
        else
            ++alignIter;
    }
    std::cout << "after removal: cellToNuclei.size() = " << cellToNuclei.size() << std::endl;

    for(cellIter = mCells.begin(); cellIter != mCells.end(); ++cellIter) {
        cellIter->second.mNumberNuclei = cellToNuclei.count(cellIter->first);

        alignIter = cellToNuclei.find(cellIter->first);
        for(int i=0; i<cellToNuclei.count(cellIter->first); i++) {
            cellIter->second.mNucleiVolume.push_back( (double)(nucleiLabelMap->GetLabelObject(alignIter->second.first)->GetNumberOfPixels())*mVoxelVolume );
            alignIter++;
        }
    }
    //-------------------------------------------------------------------------------------------------------------------------

    //----------COMPUTE-CELL-DISTANCE-TO-CV-PV---------------------------------------------------------------------------------
    CMinMaxCalculatorType::Pointer minMaxCalc = CMinMaxCalculatorType::New();
    minMaxCalc->SetImage(CVImage);
    minMaxCalc->Compute();

    if(minMaxCalc->GetMaximum() != 0) {
        SignedMaurerDistanceMapImageFilterType::Pointer distanceMapCV = SignedMaurerDistanceMapImageFilterType::New();
        distanceMapCV->ReleaseDataFlagOn();
        distanceMapCV->UseImageSpacingOn();
        distanceMapCV->SquaredDistanceOff();
        distanceMapCV->SetInput(CVImage);
        distanceMapCV->Update();

        std::map<unsigned long int, CellAnalysisContainer>::iterator cellIter;

        for(unsigned int i=0; i<cellLabelMap->GetNumberOfLabelObjects(); ++i) {
            cellIter = mCells.find(cellLabelMap->GetNthLabelObject(i)->GetLabel());

            itk::Point<double, 3> poi = cellLabelMap->GetNthLabelObject(i)->GetCentroid();
            itk::Index<3> idx;
            cellLabelMap->TransformPhysicalPointToIndex(poi, idx);
            cellIter->second.mDistToCV = distanceMapCV->GetOutput()->GetPixel(idx);
        }
    }
    else {
        std::map<unsigned long int, CellAnalysisContainer>::iterator cellIter;
        for(cellIter = mCells.begin(); cellIter != mCells.end(); ++cellIter)
            cellIter->second.mDistToCV = -1;
    }

    minMaxCalc->SetImage(PVImage);
    minMaxCalc->Compute();

    if(minMaxCalc->GetMaximum() != 0) {
        SignedMaurerDistanceMapImageFilterType::Pointer distanceMapPV = SignedMaurerDistanceMapImageFilterType::New();
        distanceMapPV->ReleaseDataFlagOn();
        distanceMapPV->UseImageSpacingOn();
        distanceMapPV->SquaredDistanceOff();
        distanceMapPV->SetInput(PVImage);
        distanceMapPV->Update();

        std::map<unsigned long int, CellAnalysisContainer>::iterator cellIter;

        for(unsigned int i=0; i<cellLabelMap->GetNumberOfLabelObjects(); ++i) {
            cellIter = mCells.find(cellLabelMap->GetNthLabelObject(i)->GetLabel());

            itk::Point<double, 3> poi = cellLabelMap->GetNthLabelObject(i)->GetCentroid();
            itk::Index<3> idx;
            cellLabelMap->TransformPhysicalPointToIndex(poi, idx);
            cellIter->second.mDistToPV = distanceMapPV->GetOutput()->GetPixel(idx);
        }
    }
    else {
        std::map<unsigned long int, CellAnalysisContainer>::iterator cellIter;
        for(cellIter = mCells.begin(); cellIter != mCells.end(); ++cellIter)
            cellIter->second.mDistToPV = -1;
    }
    //TODO: count CV/PV label objects -> if 0 then NA in table
    //-------------------------------------------------------------------------------------------------------------------------
}


void AnalyzeCellsFilter::WriteAsGraph(std::string filename)
{
    vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIds = vtkSmartPointer<vtkUnsignedLongArray>::New();
    pedigreeIds->SetName("Pedigree IDs");

    std::map<vtkIdType,float> vertexIdToVolume;

    GraphAnnotationHelper* anno = new GraphAnnotationHelper();

    for(std::map<unsigned long int, CellAnalysisContainer>::iterator cellIt = mCells.begin(); cellIt != mCells.end(); ++cellIt) {
        vtkIdType v = graph->AddVertex();
        pedigreeIds->InsertValue(v, cellIt->second.mLabel);
        points->InsertPoint(v, cellIt->second.mX, cellIt->second.mY, cellIt->second.mZ);

        vertexIdToVolume.insert(std::pair<vtkIdType,float>(v, cellIt->second.mVolume));
    }
    graph->GetVertexData()->SetPedigreeIds(pedigreeIds);
    graph->SetPoints(points);

    anno->AddCustomVertexAnnotation(graph, "volume", vertexIdToVolume, -1);

    vtkSmartPointer<vtkGraphWriter> writer = vtkSmartPointer<vtkGraphWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(graph);
    writer->Update();
}


void AnalyzeCellsFilter::WriteDataFile()
{
    std::fstream file1, file2, tempfile;
    bool hasHeader = false;

    file1.open((mDataSetPathCell + "../" + "Cell_Output_file.txt").c_str(), fstream::in);

    std::string line;
    getline(file1, line);

    if(line.find("dataSet")!=std::string::npos) {
        hasHeader = true;

        tempfile.open((mDataSetPathCell + "../" + "Cell_tempfile.txt").c_str(), fstream::out);
        tempfile << line;

        while(!file1.eof()) {
            getline(file1, line);
            if(line.find(mDataSetID)!=0) {
                tempfile << std::endl;
                tempfile << line;
            }
        }
        file1.close();
        tempfile.close();

        std::remove((mDataSetPathCell + "../" + "Cell_Output_file.txt").c_str());
        std::rename((mDataSetPathCell + "../" + "Cell_tempfile.txt").c_str(), (mDataSetPathCell + "../" + "Cell_Output_file.txt").c_str());
    }
    else
        file1.close();

    if(!hasHeader)
        file2.open((mDataSetPathCell + "../" + "Cell_Output_file.txt").c_str(), fstream::out);
    else
        file2.open((mDataSetPathCell + "../" + "Cell_Output_file.txt").c_str(), fstream::out | fstream::app);

    file2.flags(fstream::left | fstream::scientific);
    if(!hasHeader) {
        file2.width(40);
        file2 << "dataSet";
        file2.width(15);
        file2 << "label";
        file2.width(15);
        file2 << "centroidX";
        file2.width(15);
        file2 << "centroidY";
        file2.width(15);
        file2 << "centroidZ";
        file2.width(15);
        file2 << "distToCV";
        file2.width(15);
        file2 << "distToPV";
        file2.width(15);
        file2 << "volume";
        file2.width(15);
        file2 << "borderCATotal";
        file2.width(15);
        file2 << "borderCARatio";
        file2.width(15);
        file2 << "cellCATotal";
        file2.width(15);
        file2 << "cellCARatio";
        file2.width(20);
        file2 << "sinusoidCATotal";
        file2.width(20);
        file2 << "sinusoidCARatio";
        file2.width(15);
        file2 << "bileCATotal";
        file2.width(15);
        file2 << "bileCARatio";
        file2.width(15);
        file2 << "nilCATotal";
        file2.width(15);
        file2 << "nilCARatio";
        file2.width(15);
        file2 << "numNuclei";
        file2.width(15);
        file2 << "nuc1Vol";
        file2.width(15);
        file2 << "nuc2Vol";
        file2.width(15);
        file2 << "nuc3Vol";
        file2.width(15);
        file2 << "nuc4Vol" << std::endl;
    }

    for(std::map<unsigned long int, CellAnalysisContainer>::iterator cellIt = mCells.begin(); cellIt != mCells.end(); ++cellIt) {
        file2.width(40);
        file2 << mDataSetID;
        file2.width(15);
        file2 << cellIt->second.mLabel;
        file2.width(15);
        file2 << cellIt->second.mX;
        file2.width(15);
        file2 << cellIt->second.mY;
        file2.width(15);
        file2 << cellIt->second.mZ;
        file2.width(15);
        if(cellIt->second.mDistToCV != -1)
            file2 << cellIt->second.mDistToCV;
        else
            file2 << "NA";
        file2.width(15);
        if(cellIt->second.mDistToPV != -1)
            file2 << cellIt->second.mDistToPV;
        else
            file2 << "NA";
        file2.width(15);
        file2 << cellIt->second.mVolume;
        file2.width(15);
        file2 << cellIt->second.mContactAreaWithDataSetBorder;
        file2.width(15);
        file2 << cellIt->second.mContactAreaWithDataSetBorderRatio;
        file2.width(15);
        file2 << cellIt->second.mContactAreaWithCells;
        file2.width(15);
        file2 << cellIt->second.mContactAreaWithCellsRatio;
        file2.width(20);
        file2 << cellIt->second.mContactAreaWithSinusoids;
        file2.width(20);
        file2 << cellIt->second.mContactAreaWithSinusoidsRatio;
        file2.width(15);
        file2 << cellIt->second.mContactAreaWithBile;
        file2.width(15);
        file2 << cellIt->second.mContactAreaWithBileRatio;
        file2.width(15);
        file2 << cellIt->second.mContactAreaNil;
        file2.width(15);
        file2 << cellIt->second.mContactAreaNilRatio;
        file2.width(15);
        file2 << cellIt->second.mNumberNuclei;
        file2.width(15);
        if(cellIt->second.mNumberNuclei>=1)
            file2 << cellIt->second.mNucleiVolume[0];
        else
            file2 << "NA";
        file2.width(15);
        if(cellIt->second.mNumberNuclei>=2)
            file2 << cellIt->second.mNucleiVolume[1];
        else
            file2 << "NA";
        file2.width(15);
        if(cellIt->second.mNumberNuclei>=3)
            file2 << cellIt->second.mNucleiVolume[2];
        else
            file2 << "NA";
        file2.width(15);
        if(cellIt->second.mNumberNuclei>=4)
            file2 << cellIt->second.mNucleiVolume[3];
        else
            file2 << "NA";
        file2 << std::endl;
    }

    file2.close();
}


void AnalyzeCellsFilter::Update()
{
    ParseParameterContext();
    PopulateCellAnalysisContainer();
    WriteDataFile();

    if(mSaveAsGraph)
        WriteAsGraph(mDataSetPathCell + "cellAnalysisGraph.txt");
}

