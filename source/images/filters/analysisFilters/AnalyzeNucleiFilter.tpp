///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeNucleiFilter.cpp                                              //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-08                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "AnalyzeNucleiFilter.h"

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

#include "../../tools/GraphAnnotationHelper.h"

#define MMToMicron				1000



template< unsigned int VImageDimension > AnalyzeNucleiFilter< VImageDimension >::AnalyzeNucleiFilter()
{
    mDataSetID = "testDataSet";

    mNumLobuleSlices = 10;

	mVoxelVolume = 1;
	mDatasetVolume = 1;
	mCVVolume = 0;					//TODO input method for cv/pv similar to analyzeBileNetworkFilter
    mPVVolume = 0;
    mEffectiveDatasetVolume = 1;
    mTissueVolume = 0;
    mVoidVolume = 0;
    mNecRegVolume = 0;
    mAnalysedDatasetVolume = 0;
    mNumberLobules = 0;

	mSpacing[0] = 0.146;
	mSpacing[1] = 0.146;
	if(ImageDimension == 3) mSpacing[2] = 0.491;

	mSize[0] = 1024;
	mSize[1] = 1024;
	if(ImageDimension == 3) mSize[2] = 100;

#if (ITK_VERSION_MAJOR >= 4)
    itk::TIFFImageIOFactory::RegisterOneFactory();
#endif // (ITK_VERSION_MAJOR >= 4)

    mNumOverallNonProlifHepaticNuclei = 0;
    mNumOverallProlifHepaticNuclei = 0;
    mNumOverallNonProlifNonHepaticNuclei = 0;
    mNumOverallProlifNonHepaticNuclei = 0;
    mNumAnalysedNonProlifHepaticNuclei = 0;
    mNumAnalysedProlifHepaticNuclei = 0;
    mNumAnalysedNonProlifNonHepaticNuclei = 0;
    mNumAnalysedProlifNonHepaticNuclei = 0;
}


template< unsigned int VImageDimension > AnalyzeNucleiFilter< VImageDimension >::~AnalyzeNucleiFilter()
{
	// TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::WriteLogFile(std::string timeStamp)
{
    //Log not necessary
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::ParseParameterContext()
{
    if(this->m_paramContext->findContext("Analyze Nuclei",0)==NULL) {
        std::cout << "Error: AnalyzeNucleiFilter: Invalid parameter context" << std::endl;
        return;
    }

    mDataSetID = *(std::string*)(this->m_paramContext->findParameter("Dataset ID", 0)->dataPointer());

    mDataSetFullFilenameHepNuc = *(std::string*)(this->m_paramContext->findParameter("Hepatic nuclei segmentation", 0)->dataPointer());
    mDataSetFullFilenameNonHepNuc = *(std::string*)(this->m_paramContext->findParameter("Non-hepatic nuclei segmentation", 0)->dataPointer());
    mNonProlifNucleiThreshold = *(int*)(this->m_paramContext->findParameter("Non-proliferating nuclei intensity", 0)->dataPointer());
    mProlifNucleiThreshold = *(int*)(this->m_paramContext->findParameter("Proliferating nuclei intensity", 0)->dataPointer());
    mDataSetFullFilenameCV = *(std::string*)(this->m_paramContext->findParameter("CV segmentation", 0)->dataPointer());
    mCentralVeinThreshold = *(int*)(this->m_paramContext->findParameter("CV intensity", 0)->dataPointer());
    mDataSetFullFilenamePV = *(std::string*)(this->m_paramContext->findParameter("PV segmentation", 0)->dataPointer());
    mPortalVeinThreshold = *(int*)(this->m_paramContext->findParameter("PV intensity", 0)->dataPointer());
    mDataSetFullFilenameTissue = *(std::string*)(this->m_paramContext->findParameter("Tissue segmentation", 0)->dataPointer());
    mTissueThreshold = *(int*)(this->m_paramContext->findParameter("Tissue intensity", 0)->dataPointer());
    mDataSetFullFilenameNecReg = *(std::string*)(this->m_paramContext->findParameter("Necrotic Region segmentation", 0)->dataPointer());
    mNecRegThreshold = *(int*)(this->m_paramContext->findParameter("Necrotic Region intensity", 0)->dataPointer());
    mVoidThreshold = *(int*)(this->m_paramContext->findParameter("Void intensity", 0)->dataPointer());
    mDataSetFullFilenameLobules = *(std::string*)(this->m_paramContext->findParameter("Lobule segmentation", 0)->dataPointer());
    mLobulesThreshold = *(int*)(this->m_paramContext->findParameter("Lobule intensity", 0)->dataPointer());

    mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenameNonHepNuc) );
    mNonHepNucFileExists = mInfoFullFilename.exists();
    mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenameHepNuc) );
    mHepNucFileExists = mInfoFullFilename.exists();

    if(!mNonHepNucFileExists && !mHepNucFileExists)
        throw std::string("Please specify at least one nuclei segmentation file.");

    mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenameCV) );
    mCVFileExists = mInfoFullFilename.exists();

    mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenamePV) );
    mPVFileExists = mInfoFullFilename.exists();

    mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenameTissue) );
    mTissueFileExists = mInfoFullFilename.exists();

    mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenameNecReg) );
    mNecRegFileExists = mInfoFullFilename.exists();

    mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenameLobules) );
    mLobulesFileExists = mInfoFullFilename.exists();

    if(mNonHepNucFileExists) {
        mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenameNonHepNuc) );
        mPath = (mInfoFullFilename.path() + QString("/")).toStdString();
        mFilenameExtension = (QString(".") + mInfoFullFilename.suffix()).toStdString();
    }
    else {
        mInfoFullFilename.setFile(QString::fromStdString( mDataSetFullFilenameHepNuc) );
        mPath = (mInfoFullFilename.path() + QString("/")).toStdString();
        mFilenameExtension = (QString(".") + mInfoFullFilename.suffix()).toStdString();
    }

    mSpacing[0] = *(double*)(this->m_paramContext->findParameter("Voxel spacing x", 0)->dataPointer());
    mSpacing[1] = *(double*)(this->m_paramContext->findParameter("Voxel spacing y", 0)->dataPointer());
    if(ImageDimension == 3) mSpacing[2] = *(double*)(this->m_paramContext->findParameter("Voxel spacing z", 0)->dataPointer());

    mSaveAsGraph = *(bool*)(this->m_paramContext->findParameter("Save analysis result as graph", 0)->dataPointer());

    if(mNonHepNucFileExists && BasePipeline<ImageDimension>::GetNumberOfDimensions(mDataSetFullFilenameNonHepNuc) != ImageDimension)
        throw std::string("Data " + mDataSetFullFilenameNonHepNuc + " has different dimensionality than expected by pipeline.");
    if(mHepNucFileExists && BasePipeline<ImageDimension>::GetNumberOfDimensions(mDataSetFullFilenameHepNuc) != ImageDimension)
        throw std::string("Data " + mDataSetFullFilenameHepNuc + " has different dimensionality than expected by pipeline.");
    if(mCVFileExists && BasePipeline<ImageDimension>::GetNumberOfDimensions(mDataSetFullFilenameCV) != ImageDimension)
        throw std::string("Data " + mDataSetFullFilenameCV + " has different dimensionality than expected by pipeline.");
    if(mPVFileExists && BasePipeline<ImageDimension>::GetNumberOfDimensions(mDataSetFullFilenamePV) != ImageDimension)
        throw std::string("Data " + mDataSetFullFilenamePV + " has different dimensionality than expected by pipeline.");
    if(mTissueFileExists && BasePipeline<ImageDimension>::GetNumberOfDimensions(mDataSetFullFilenameTissue) != ImageDimension)
        throw std::string("Data " + mDataSetFullFilenameTissue + " has different dimensionality than expected by pipeline.");
    if(mNecRegFileExists && BasePipeline<ImageDimension>::GetNumberOfDimensions(mDataSetFullFilenameNecReg) != ImageDimension)
        throw std::string("Data " + mDataSetFullFilenameNecReg + " has different dimensionality than expected by pipeline.");
    if(mLobulesFileExists && BasePipeline<ImageDimension>::GetNumberOfDimensions(mDataSetFullFilenameLobules) != ImageDimension)
        throw std::string("Data " + mDataSetFullFilenameLobules + " has different dimensionality than expected by pipeline.");
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::ComputeRelativeLobularPositions(FScalarImagePointerType cvDistMap, FScalarImagePointerType lobuleDistMap)
{
    typename AddFImageFilterType::Pointer addDist;
    typename DivideImageFilterType::Pointer divide;

    addDist = AddFImageFilterType::New();
    addDist->ReleaseDataFlagOn();
    addDist->SetInput1(cvDistMap);
    addDist->SetInput2(lobuleDistMap);

    divide = DivideImageFilterType::New();
    divide->ReleaseDataFlagOn();
    divide->SetInput1(cvDistMap);
    divide->SetInput2(addDist->GetOutput());
    divide->Update();

    mRelativeLobularPositionMap = divide->GetOutput();
    mRelativeLobularPositionMap->DisconnectPipeline();
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::CollectBasicImageInformation()
{
    std::string pathToNucleiImage = mDataSetFullFilenameNonHepNuc;
    if(!mNonHepNucFileExists)
        pathToNucleiImage = mDataSetFullFilenameHepNuc;

	itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO((pathToNucleiImage).c_str(), itk::ImageIOFactory::ReadMode);
	if(imageIO->CanReadFile((pathToNucleiImage).c_str())) {
		imageIO->ReadImageInformation();

		mVoxelVolume = 1;
		mDatasetVolume = 1;

		for(unsigned int i=0; i<ImageDimension; i++) {
			mSize[i] = imageIO->GetDimensions(i);
			mVoxelVolume *= mSpacing[i];
			mDatasetVolume *= mSize[i]*mSpacing[i];
		}
		int dim = ImageDimension;
		mVoxelVolume /= std::pow(MMToMicron, dim);
		mDatasetVolume /= std::pow(MMToMicron, dim);
	}
	mEffectiveDatasetVolume = mDatasetVolume;
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::CollectBasicNucleiInformation()
{
    typename CScalarImageReaderType::Pointer hepNucReader, nonHepNucReader;
    typename ThresholdFilterType::Pointer thresholdHepNuc, thresholdNonHepNuc;
    typename ImageToShapeLabelMapFilterType::Pointer lobuleImageToShaLabMapFilter, nonProlifHepNucImageToShaLabMapFilter, prolifHepNucImageToShaLabMapFilter, nonProlifNonHepNucImageToShaLabMapFilter, prolifNonHepNucImageToShaLabMapFilter;
    typename LabelMapToLabelImageFilterType::Pointer lobuleShaLabMapToLabelImageFilter;


    typename LScalarImageType::Pointer lobuleImage;
    if(mLobulesFileExists) {
        lobuleImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
        lobuleImageToShaLabMapFilter->FullyConnectedOff();
        lobuleImageToShaLabMapFilter->SetInput(mLobuleBin);
        lobuleImageToShaLabMapFilter->Update();

        mNumberLobules = lobuleImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
        mAnalysedDatasetVolume = 0.;
        for(unsigned int i=0; i<lobuleImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++)
            if(lobuleImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetLabel()!=0)
                mAnalysedDatasetVolume += lobuleImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetPhysicalSize();

        lobuleShaLabMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
        lobuleShaLabMapToLabelImageFilter->SetInput(lobuleImageToShaLabMapFilter->GetOutput());
        lobuleShaLabMapToLabelImageFilter->Update();

        lobuleImage = lobuleShaLabMapToLabelImageFilter->GetOutput();
        lobuleImage->DisconnectPipeline();
    }
    else
        mAnalysedDatasetVolume = mEffectiveDatasetVolume;
    int dim = ImageDimension;
    mAnalysedDatasetVolume /= std::pow(MMToMicron, dim);

    if(mHepNucFileExists) {
        hepNucReader = CScalarImageReaderType::New();
        hepNucReader->SetFileName(mDataSetFullFilenameHepNuc);
        hepNucReader->ReleaseDataBeforeUpdateFlagOn();

        thresholdHepNuc = ThresholdFilterType::New();
        thresholdHepNuc->ReleaseDataFlagOn();
        thresholdHepNuc->SetOutsideValue(0);
        thresholdHepNuc->SetInsideValue(255);
        thresholdHepNuc->SetLowerThreshold(mNonProlifNucleiThreshold);
        thresholdHepNuc->SetUpperThreshold(mNonProlifNucleiThreshold);
        thresholdHepNuc->SetInput(hepNucReader->GetOutput());
        thresholdHepNuc->Update();

        itk::SmartPointer<CScalarImageType> nonProlifHepNucImage = thresholdHepNuc->GetOutput();
        nonProlifHepNucImage->DisconnectPipeline();
        nonProlifHepNucImage->SetSpacing(mSpacing);

        thresholdHepNuc->SetLowerThreshold(mProlifNucleiThreshold);
        thresholdHepNuc->SetUpperThreshold(mProlifNucleiThreshold);
        thresholdHepNuc->Update();

        itk::SmartPointer<CScalarImageType> prolifHepNucImage = thresholdHepNuc->GetOutput();
        prolifHepNucImage->DisconnectPipeline();
        prolifHepNucImage->SetSpacing(mSpacing);

        nonProlifHepNucImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
        nonProlifHepNucImageToShaLabMapFilter->FullyConnectedOn();
        nonProlifHepNucImageToShaLabMapFilter->SetInput(nonProlifHepNucImage);
        nonProlifHepNucImageToShaLabMapFilter->Update();

        prolifHepNucImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
        prolifHepNucImageToShaLabMapFilter->FullyConnectedOn();
        prolifHepNucImageToShaLabMapFilter->SetInput(prolifHepNucImage);
        prolifHepNucImageToShaLabMapFilter->Update();

        mNumOverallNonProlifHepaticNuclei = nonProlifHepNucImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
        mNumOverallProlifHepaticNuclei = prolifHepNucImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
        mNumAnalysedNonProlifHepaticNuclei = 0;
        mNumAnalysedProlifHepaticNuclei = 0;

        std::cout << "mNumOverallNonProlifHepaticNuclei = " << mNumOverallNonProlifHepaticNuclei << ", mNumOverallProlifHepaticNuclei = " << mNumOverallProlifHepaticNuclei << std::endl;

        for(unsigned int i=0; i<mNumOverallNonProlifHepaticNuclei; i++) {
            NucleiAnalysisContainer n;
            n.label = i;
            n.volume = nonProlifHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetPhysicalSize();
            n.roundness = nonProlifHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetRoundness();
            n.elongation = nonProlifHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetElongation();
            n.numPxlOnBorder = nonProlifHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixelsOnBorder();
            n.isProliferating = false;

            PointType poi = nonProlifHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetCentroid();
            n.x = poi[0];
            n.y = poi[1];
            n.z = poi[2];

            typename FScalarImageType::IndexType idx;
            nonProlifHepNucImageToShaLabMapFilter->GetOutput()->TransformPhysicalPointToIndex(poi, idx);
            if(mCVFileExists)                       n.distToCV = mCVDistMap->GetPixel( idx );
            if(mPVFileExists)                       n.distToPV = mPVDistMap->GetPixel( idx );
            if(mNecRegFileExists)                   n.distToNecReg = mNecRegDistMap->GetPixel( idx );
            if(mLobulesFileExists)                  n.lobuleLabel = lobuleImage->GetPixel( idx );
            if(mLobulesFileExists)                  n.distToLobularBoundary = mLobuleDistMap->GetPixel( idx );
            if(mLobulesFileExists && mCVFileExists) n.relativeLobularPos = mRelativeLobularPositionMap->GetPixel( idx );

            if((mLobulesFileExists && n.lobuleLabel!=0) || !mLobulesFileExists) {
                if((mNecRegFileExists && n.distToNecReg>0) || !mNecRegFileExists) {
                    mNumAnalysedNonProlifHepaticNuclei++;
                    mHepaticNuclei.push_back(n);
                }
            }
        }

        for(unsigned int i=0; i<mNumOverallProlifHepaticNuclei; i++) {
            NucleiAnalysisContainer n;
            n.label = i+mNumOverallNonProlifHepaticNuclei;
            n.volume = prolifHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetPhysicalSize();
            n.roundness = prolifHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetRoundness();
            n.elongation = prolifHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetElongation();
            n.numPxlOnBorder = prolifHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixelsOnBorder();
            n.isProliferating = true;

            PointType poi = prolifHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetCentroid();
            n.x = poi[0];
            n.y = poi[1];
            n.z = poi[2];

            typename FScalarImageType::IndexType idx;
            prolifHepNucImageToShaLabMapFilter->GetOutput()->TransformPhysicalPointToIndex(poi, idx);
            if(mCVFileExists)                       n.distToCV = mCVDistMap->GetPixel( idx );
            if(mPVFileExists)                       n.distToPV = mPVDistMap->GetPixel( idx );
            if(mNecRegFileExists)                   n.distToNecReg = mNecRegDistMap->GetPixel( idx );
            if(mLobulesFileExists)                  n.lobuleLabel = lobuleImage->GetPixel( idx );
            if(mLobulesFileExists)                  n.distToLobularBoundary = mLobuleDistMap->GetPixel( idx );
            if(mLobulesFileExists && mCVFileExists) n.relativeLobularPos = mRelativeLobularPositionMap->GetPixel( idx );

            if((mLobulesFileExists && n.lobuleLabel!=0) || !mLobulesFileExists) {
                if((mNecRegFileExists && n.distToNecReg>0) || !mNecRegFileExists) {
                    mNumAnalysedProlifHepaticNuclei++;
                    mHepaticNuclei.push_back(n);
                }
            }
        }
    }

    if(mNonHepNucFileExists) {
        nonHepNucReader = CScalarImageReaderType::New();
        nonHepNucReader->SetFileName(mDataSetFullFilenameNonHepNuc);
        nonHepNucReader->ReleaseDataBeforeUpdateFlagOn();

        thresholdNonHepNuc = ThresholdFilterType::New();
        thresholdNonHepNuc->ReleaseDataFlagOn();
        thresholdNonHepNuc->SetOutsideValue(0);
        thresholdNonHepNuc->SetInsideValue(255);
        thresholdNonHepNuc->SetLowerThreshold(mNonProlifNucleiThreshold);
        thresholdNonHepNuc->SetUpperThreshold(mNonProlifNucleiThreshold);
        thresholdNonHepNuc->SetInput(nonHepNucReader->GetOutput());
        thresholdNonHepNuc->Update();

        itk::SmartPointer<CScalarImageType> nonProlifNonHepNucImage = thresholdNonHepNuc->GetOutput();
        nonProlifNonHepNucImage->DisconnectPipeline();
        nonProlifNonHepNucImage->SetSpacing(mSpacing);

        thresholdNonHepNuc->SetLowerThreshold(mProlifNucleiThreshold);
        thresholdNonHepNuc->SetUpperThreshold(mProlifNucleiThreshold);
        thresholdNonHepNuc->Update();

        itk::SmartPointer<CScalarImageType> prolifNonHepNucImage = thresholdNonHepNuc->GetOutput();
        prolifNonHepNucImage->DisconnectPipeline();
        prolifNonHepNucImage->SetSpacing(mSpacing);
        //-------------------------------------------------------------------------------------------------------------------------

        //----------FILTER---CREATE-LABEL-MAP--------------------------------------------------------------------------------------
        nonProlifNonHepNucImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
        nonProlifNonHepNucImageToShaLabMapFilter->FullyConnectedOn();
        nonProlifNonHepNucImageToShaLabMapFilter->SetInput(nonProlifNonHepNucImage);
        nonProlifNonHepNucImageToShaLabMapFilter->Update();

        prolifNonHepNucImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
        prolifNonHepNucImageToShaLabMapFilter->FullyConnectedOn();
        prolifNonHepNucImageToShaLabMapFilter->SetInput(prolifNonHepNucImage);
        prolifNonHepNucImageToShaLabMapFilter->Update();

        mNumOverallNonProlifNonHepaticNuclei = nonProlifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
        mNumOverallProlifNonHepaticNuclei = prolifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
        mNumAnalysedNonProlifNonHepaticNuclei = 0;
        mNumAnalysedProlifNonHepaticNuclei = 0;

        std::cout << "mNumOverallNonProlifNonHepaticNuclei = " << mNumOverallNonProlifNonHepaticNuclei << ", mNumOverallProlifNonHepaticNuclei = " << mNumOverallProlifNonHepaticNuclei << std::endl;

        for(unsigned int i=0; i<mNumOverallNonProlifNonHepaticNuclei; i++) {
            NucleiAnalysisContainer n;
            n.label = i;
            n.volume = nonProlifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetPhysicalSize();
            n.roundness = nonProlifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetRoundness();
            n.elongation = nonProlifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetElongation();
            n.numPxlOnBorder = nonProlifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixelsOnBorder();
            n.isProliferating = false;

            n.x = nonProlifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetCentroid()[0];
            n.y = nonProlifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetCentroid()[1];
            if(ImageDimension == 3) n.z = nonProlifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetCentroid()[2];

            if((mLobulesFileExists && n.lobuleLabel!=0) || !mLobulesFileExists) {
                    mNumAnalysedNonProlifNonHepaticNuclei++;
                    mNonHepaticNuclei.push_back(n);
            }
        }

        for(unsigned int i=0; i<mNumOverallProlifNonHepaticNuclei; i++) {
            NucleiAnalysisContainer n;
            n.label = i+mNumOverallNonProlifNonHepaticNuclei;
            n.volume = prolifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetPhysicalSize();
            n.roundness = prolifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetRoundness();
            n.elongation = prolifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetElongation();
            n.numPxlOnBorder = prolifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixelsOnBorder();
            n.isProliferating = true;

            n.x = prolifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetCentroid()[0];
            n.y = prolifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetCentroid()[1];
            if(ImageDimension == 3) n.z = prolifNonHepNucImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetCentroid()[2];

            if((mLobulesFileExists && n.lobuleLabel!=0) || !mLobulesFileExists) {
                mNumAnalysedProlifNonHepaticNuclei++;
                mNonHepaticNuclei.push_back(n);
            }
        }
    }
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::CollectBasicLobuleInformation()
{
    typename ImageToShapeLabelMapFilterType::Pointer lobuleImageToShaLabMapFilter, necRegImageToShaLabMapFilter;
    typename LabelMapToLabelImageFilterType::Pointer lobuleShaLabMapToLabelImageFilter, necRegShaLabMapToLabelImageFilter;


    typename LScalarImageType::Pointer lobuleImage, necRegImage;
    if(mLobulesFileExists) {
        lobuleImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
        lobuleImageToShaLabMapFilter->FullyConnectedOff();
        lobuleImageToShaLabMapFilter->SetInput(mLobuleBin);
        lobuleImageToShaLabMapFilter->Update();

        lobuleShaLabMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
        lobuleShaLabMapToLabelImageFilter->SetInput(lobuleImageToShaLabMapFilter->GetOutput());
        lobuleShaLabMapToLabelImageFilter->Update();

        lobuleImage = lobuleShaLabMapToLabelImageFilter->GetOutput();
        lobuleImage->DisconnectPipeline();

        if(mNecRegFileExists) {
            necRegImageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
            necRegImageToShaLabMapFilter->FullyConnectedOff();
            necRegImageToShaLabMapFilter->SetInput(mNecRegBin);
            necRegImageToShaLabMapFilter->Update();

            necRegShaLabMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
            necRegShaLabMapToLabelImageFilter->SetInput(necRegImageToShaLabMapFilter->GetOutput());
            necRegShaLabMapToLabelImageFilter->Update();

            necRegImage = necRegShaLabMapToLabelImageFilter->GetOutput();
            necRegImage->DisconnectPipeline();
        }


        double sliceIntervals = 1. / (double)mNumLobuleSlices;

        std::map< double, std::vector<double> > lobuleToSliceVolume;

        itk::ImageRegionConstIterator<LScalarImageType> iter(lobuleImage, lobuleImage->GetLargestPossibleRegion());
        for(iter = iter.Begin(); !iter.IsAtEnd(); ++iter) {
            if(iter.Value() != 0) {
                if(lobuleToSliceVolume.count(iter.Value()) == 0) {
                    std::vector<double> sliceVolume(mNumLobuleSlices);
                    std::pair<double, std::vector<double> > n(iter.Value(), sliceVolume);
                    lobuleToSliceVolume.insert(n);
                }
                typename LScalarImageType::IndexType idx = iter.GetIndex();
                int sliceNr = mRelativeLobularPositionMap->GetPixel( idx ) / sliceIntervals;
                lobuleToSliceVolume[iter.Value()][sliceNr]++;
            }
        }

        for(unsigned int i=0; i<mNumberLobules; i++) {
            std::vector<double> sliceVolume = lobuleToSliceVolume[lobuleImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetLabel()];

            LobulesAnalysisContainer n;
            n.label = i;
            n.volume = lobuleImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetPhysicalSize();

            if(mNecRegFileExists) {
                for(unsigned int j=0; j<necRegImageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects(); j++) {
                    if(n.label == lobuleImage->GetPixel(necRegImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(j)->GetIndex(0))) {
                        n.volumeNecReg = necRegImageToShaLabMapFilter->GetOutput()->GetNthLabelObject(j)->GetPhysicalSize();
                        break;
                    }
                }
            }

            for(unsigned int j=0; j<sliceVolume.size(); j++)
                sliceVolume[j] = sliceVolume[j] * mVoxelVolume;

            n.volumePerRelativePositionSlice = sliceVolume;

            mLobules.push_back(n);
        }

    }
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::WriteAsGraph(std::string filename, std::vector<NucleiAnalysisContainer> nuclei)
{
    vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIds = vtkSmartPointer<vtkUnsignedLongArray>::New();
    pedigreeIds->SetName("Pedigree IDs");

    std::map<vtkIdType,int> vertexIdToLobuleLabel;
    std::map<vtkIdType,float> vertexIdToVolume;
    std::map<vtkIdType,float> vertexIdToRoundness;
    std::map<vtkIdType,float> vertexIdToElongation;
    std::map<vtkIdType,float> vertexIdToIsProliferating;
    std::map<vtkIdType,float> vertexIdToRelativeLobularPosition;
    std::map<vtkIdType,float> vertexIdToNumPxlOnBorder;
    std::map<vtkIdType,float> vertexIdToDistToCV;
    std::map<vtkIdType,float> vertexIdToDistToPV;
    std::map<vtkIdType,float> vertexIdToDistToNecReg;
    std::map<vtkIdType,float> vertexIdToDistToLobularBoundary;

    GraphAnnotationHelper* anno = new GraphAnnotationHelper();

    for(unsigned int i=0; i<nuclei.size(); i++) {
        vtkIdType v = graph->AddVertex();
        pedigreeIds->InsertValue(v, nuclei[i].label);
        points->InsertPoint(v, nuclei[i].x, nuclei[i].y, nuclei[i].z);

        vertexIdToLobuleLabel.insert(std::pair<vtkIdType,int>(v, nuclei[i].lobuleLabel));
        vertexIdToVolume.insert(std::pair<vtkIdType,float>(v, nuclei[i].volume));
        vertexIdToRoundness.insert(std::pair<vtkIdType,float>(v, nuclei[i].roundness));
        vertexIdToElongation.insert(std::pair<vtkIdType,float>(v, nuclei[i].elongation));
        vertexIdToIsProliferating.insert(std::pair<vtkIdType,float>(v, nuclei[i].isProliferating));
        vertexIdToRelativeLobularPosition.insert(std::pair<vtkIdType,float>(v, nuclei[i].relativeLobularPos));
        vertexIdToNumPxlOnBorder.insert(std::pair<vtkIdType,float>(v, nuclei[i].numPxlOnBorder));
        vertexIdToDistToCV.insert(std::pair<vtkIdType,float>(v, nuclei[i].distToCV));
        vertexIdToDistToPV.insert(std::pair<vtkIdType,float>(v, nuclei[i].distToPV));
        vertexIdToDistToNecReg.insert(std::pair<vtkIdType,float>(v, nuclei[i].distToNecReg));
        vertexIdToDistToLobularBoundary.insert(std::pair<vtkIdType,float>(v, nuclei[i].distToLobularBoundary));
    }
    graph->GetVertexData()->SetPedigreeIds(pedigreeIds);
    graph->SetPoints(points);

    anno->AddCustomVertexAnnotation(graph, "lobule label", vertexIdToLobuleLabel, -1);
    anno->AddCustomVertexAnnotation(graph, "volume", vertexIdToVolume, -1);
    anno->AddCustomVertexAnnotation(graph, "roundness", vertexIdToRoundness, -1);
    anno->AddCustomVertexAnnotation(graph, "elongation", vertexIdToElongation, -1);
    anno->AddCustomVertexAnnotation(graph, "is proliferating", vertexIdToIsProliferating, -0.1);
    anno->AddCustomVertexAnnotation(graph, "relative lobular position", vertexIdToRelativeLobularPosition, -0.1);
    anno->AddCustomVertexAnnotation(graph, "number pixel on border", vertexIdToNumPxlOnBorder, -1);
    anno->AddCustomVertexAnnotation(graph, "distance to CV", vertexIdToDistToCV, -1);
    anno->AddCustomVertexAnnotation(graph, "distance to PV", vertexIdToDistToPV, -1);
    anno->AddCustomVertexAnnotation(graph, "distance to Necrotic Region", vertexIdToDistToNecReg, -1);
    anno->AddCustomVertexAnnotation(graph, "distance to lobular boundary", vertexIdToDistToLobularBoundary, -1);

    vtkSmartPointer<vtkGraphWriter> writer = vtkSmartPointer<vtkGraphWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(graph);
    writer->Update();
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::WriteBasicInformationFile()
{
    std::string suffixVolORArea = "Area";
    if(ImageDimension == 3) suffixVolORArea = "Volume";

	vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
	reader->SetFileName((mPath + "../" + "Nuclei_Output_file_1.txt").c_str());
	reader->SetHaveHeaders(true);
	reader->DetectNumericColumnsOn();
	reader->SetFieldDelimiterCharacters(" ");
	reader->Update();

	vtkSmartPointer<vtkTable> readTable = reader->GetOutput();

	vtkSmartPointer<vtkVariantArray> head =  vtkVariantArray::SafeDownCast(readTable->GetColumnByName("header"));
	if(head==NULL) {
		vtkSmartPointer<vtkVariantArray> header = vtkSmartPointer<vtkVariantArray>::New();
		header->InsertNextValue( vtkVariant("numDim") );
		header->InsertNextValue( vtkVariant("dimX") );
		header->InsertNextValue( vtkVariant("dimY") );
		if(ImageDimension == 3) header->InsertNextValue( vtkVariant("dimZ") );
		header->InsertNextValue( vtkVariant("dataSet"+suffixVolORArea) );
		header->InsertNextValue( vtkVariant("effective"+suffixVolORArea) );
		header->InsertNextValue( vtkVariant("tissue"+suffixVolORArea) );
		header->InsertNextValue( vtkVariant("void"+suffixVolORArea) );
		header->InsertNextValue( vtkVariant("centralVein"+suffixVolORArea) );
		header->InsertNextValue( vtkVariant("portalVein"+suffixVolORArea) );
		header->InsertNextValue( vtkVariant("necReg"+suffixVolORArea) );
		header->InsertNextValue( vtkVariant("analysed"+suffixVolORArea) );
		header->InsertNextValue( vtkVariant("analysedLobules") );
		header->InsertNextValue( vtkVariant("numAnalysedNonProlifHepaticNuclei") );
		header->InsertNextValue( vtkVariant("numAnalysedProlifHepaticNuclei") );
		header->InsertNextValue( vtkVariant("numAnalysedNonProlifNonHepaticNuclei") );
		header->InsertNextValue( vtkVariant("numAnalysedProlifNonHepaticNuclei") );
		header->SetName("header");

		readTable->AddColumn(header);
	}

	vtkSmartPointer<vtkVariantArray> dataSet = vtkSmartPointer<vtkVariantArray>::New();
	dataSet->InsertNextValue( vtkVariant(ImageDimension) );
	dataSet->InsertNextValue( vtkVariant(mSize[0]) );
	dataSet->InsertNextValue( vtkVariant(mSize[1]) );
	if(ImageDimension == 3) dataSet->InsertNextValue( vtkVariant(mSize[2]) );
	dataSet->InsertNextValue( vtkVariant((double)mDatasetVolume) );
	dataSet->InsertNextValue( vtkVariant((double)mEffectiveDatasetVolume) );
	dataSet->InsertNextValue( vtkVariant((double)(mTissueVolume)) );
	dataSet->InsertNextValue( vtkVariant((double)(mVoidVolume)) );
	dataSet->InsertNextValue( vtkVariant((double)(mCVVolume)) );
	dataSet->InsertNextValue( vtkVariant((double)(mPVVolume)) );
	dataSet->InsertNextValue( vtkVariant((double)(mNecRegVolume)) );
	dataSet->InsertNextValue( vtkVariant((double)(mAnalysedDatasetVolume)) );
	dataSet->InsertNextValue( vtkVariant((double)(mNumberLobules)) );
	dataSet->InsertNextValue( vtkVariant(mNumAnalysedNonProlifHepaticNuclei) );
	dataSet->InsertNextValue( vtkVariant(mNumAnalysedProlifHepaticNuclei) );
	dataSet->InsertNextValue( vtkVariant(mNumAnalysedNonProlifNonHepaticNuclei) );
	dataSet->InsertNextValue( vtkVariant(mNumAnalysedProlifNonHepaticNuclei) );
	dataSet->SetName(mDataSetID.c_str());

	readTable->AddColumn(dataSet);

	readTable->Dump(25);

	vtkSmartPointer<vtkDelimitedTextWriter> writer = vtkSmartPointer<vtkDelimitedTextWriter>::New();
	writer->SetFileName((mPath + "../" + "Nuclei_Output_file_1.txt").c_str());
	writer->SetFieldDelimiter(" ");
	writer->SetInput(readTable);
	writer->Update();
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::WriteNucleiInformationFile(std::string filenamePostfix, std::vector<NucleiAnalysisContainer> nuclei)
{
    std::string suffixVolORArea = "area";
    if(ImageDimension == 3) suffixVolORArea = "volume";

    std::fstream file1, file2, tempfile;
    bool hasHeader = false;

    file1.open((mPath + "../Nuclei_Output_file" + filenamePostfix + ".txt").c_str(), fstream::in);

    std::string line;
    getline(file1, line);

    if(line.find("dataset")!=std::string::npos) {
        hasHeader = true;

        tempfile.open((mPath + "../Nuclei_Output_file_tempfile" + filenamePostfix + ".txt").c_str(), fstream::out);
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

        std::remove((mPath + "../Nuclei_Output_file" + filenamePostfix + ".txt").c_str());
        std::rename((mPath + "../Nuclei_Output_file_tempfile" + filenamePostfix + ".txt").c_str(), (mPath + "../Nuclei_Output_file" + filenamePostfix + ".txt").c_str());
    }
    else
        file1.close();

    if(!hasHeader)
        file2.open((mPath + "../Nuclei_Output_file" + filenamePostfix + ".txt").c_str(), fstream::out);
    else
        file2.open((mPath + "../Nuclei_Output_file" + filenamePostfix + ".txt").c_str(), fstream::out | fstream::app);

    file2.flags(fstream::left | fstream::scientific);
    if(!hasHeader) {
        file2.width(40);
        file2 << "dataset";
        file2.width(10);
        file2 << "lobuleLabel";
        file2.width(10);
        file2 << "label";
        file2.width(20);
        file2 << suffixVolORArea;
//        file2.width(20);
//        file2 << "roundness";
//        file2.width(20);
//        file2 << "elongation";
        file2.width(20);
        file2 << "isProliferating";
        file2.width(20);
        file2 << "relativeLobularPos";
//        file2.width(20);
//        file2 << "xCentroid";
//        file2.width(20);
//        file2 << "yCentroid";
//        if(ImageDimension == 3) {
//            file2.width(20);
//            file2 << "zCentroid";
//        }
//        file2.width(20);
//        file2 << "numPxlOnBorder";
        file2.width(20);
        file2 << "distToCV";
        file2.width(20);
        file2 << "distToPV" << std::endl;
    }

    for(unsigned int i=0; i<nuclei.size(); i++) {
        file2.width(40);
        file2 << mDataSetID;
        file2.width(10);
        file2 << nuclei[i].lobuleLabel;
        file2.width(10);
        file2 << nuclei[i].label;
        file2.width(20);
        file2 << nuclei[i].volume;
//        file2.width(20);
//        file2 << nuclei[i].roundness;
//        file2.width(20);
//        file2 << nuclei[i].elongation;
        file2.width(20);
        file2 << nuclei[i].isProliferating;
        file2.width(20);
        file2 << nuclei[i].relativeLobularPos;
//        file2.width(20);
//        file2 << nuclei[i].x;
//        file2.width(20);
//        file2 << nuclei[i].y;
//        if(ImageDimension == 3) {
//            file2.width(20);
//            file2 << nuclei[i].z;
//        }
//        file2.width(20);
//        file2 << nuclei[i].numPxlOnBorder;
        file2.width(20);
        file2 << nuclei[i].distToCV;
        file2.width(20);
        file2 << nuclei[i].distToPV << std::endl;
    }

    file2.close();
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::WriteLobuleInformationFile(std::string filenamePostfix, std::vector<LobulesAnalysisContainer> lobules)
{
    std::string suffixVolORArea = "area";
    if(ImageDimension == 3) suffixVolORArea = "volume";

    std::fstream file1, file2, tempfile;
    bool hasHeader = false;

    file1.open((mPath + "../Lobule_Output_file" + filenamePostfix + ".txt").c_str(), fstream::in);

    std::string line;
    getline(file1, line);

    if(line.find("dataset")!=std::string::npos) {
        hasHeader = true;

        tempfile.open((mPath + "../Lobule_Output_file_tempfile" + filenamePostfix + ".txt").c_str(), fstream::out);
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

        std::remove((mPath + "../Lobule_Output_file" + filenamePostfix + ".txt").c_str());
        std::rename((mPath + "../Lobule_Output_file_tempfile" + filenamePostfix + ".txt").c_str(), (mPath + "../Lobule_Output_file" + filenamePostfix + ".txt").c_str());
    }
    else
        file1.close();

    if(!hasHeader)
        file2.open((mPath + "../Lobule_Output_file" + filenamePostfix + ".txt").c_str(), fstream::out);
    else
        file2.open((mPath + "../Lobule_Output_file" + filenamePostfix + ".txt").c_str(), fstream::out | fstream::app);

    file2.flags(fstream::left | fstream::scientific);
    if(!hasHeader) {
        file2.width(40);
        file2 << "dataset";
        file2.width(10);
        file2 << "label";
        file2.width(20);
        file2 << suffixVolORArea;
        for(unsigned int i=0; i<mNumLobuleSlices; i++) {
            file2.width(20);
            file2 << suffixVolORArea << "Slice_" << i;
        }
        file2.width(30);
        file2 << "numCV";
        file2.width(30);
        file2 << "numPV";
        file2.width(30);
        file2 << suffixVolORArea << "NecReg" << std::endl;
    }

    for(unsigned int i=0; i<lobules.size(); i++) {
        file2.width(40);
        file2 << mDataSetID;
        file2.width(10);
        file2 << lobules[i].label;
        file2.width(20);
        file2 << lobules[i].volume;
        for(unsigned int j=0; j<lobules[i].volumePerRelativePositionSlice.size(); j++) {
            file2.width(20);
            file2 << lobules[i].volumePerRelativePositionSlice[j];
        }
        file2.width(30);
        file2 << lobules[i].numCV;
        file2.width(30);
        file2 << lobules[i].numPV;
        file2.width(30);
        file2 << lobules[i].volumeNecReg << std::endl;
    }

    file2.close();
}


template< unsigned int VImageDimension > void AnalyzeNucleiFilter< VImageDimension >::Update()
{
    ParseParameterContext();

    std::string timeStamp;

    QDateTime time = QDateTime::currentDateTime();
    timeStamp = time.toString("hh:mm:ss dd.MM.yyyy").toStdString();

    std::cout << "Analyze nuclei pipeline: " << std::endl;
    std::cout << " path: " << mPath << std::endl;
    std::cout << " extension: " << mFilenameExtension << std::endl;

    std::cout << " mDataSetFullFilenameNonHepNuc: " << mDataSetFullFilenameNonHepNuc << std::endl;
    std::cout << " mDataSetFullFilenameHepNuc: " << mDataSetFullFilenameHepNuc << std::endl;
    std::cout << " mDataSetFullFilenameCV: " << mDataSetFullFilenameCV << std::endl;
    std::cout << " mDataSetFullFilenamePV: " << mDataSetFullFilenamePV << std::endl;
    std::cout << " mDataSetFullFilenameTissue: " << mDataSetFullFilenameTissue << std::endl;

    CollectBasicImageInformation();

    if(mCVFileExists) {
        this->ReadImageLayer(mDataSetFullFilenameCV, mCVBin, mSpacing, mCentralVeinThreshold);
        this->BuildDistanceMap(mCVBin, mCVDistMap, mSpacing, true);

        mCVVolume = this->ComputeVolumeOfRegions(mDataSetFullFilenameCV, mSpacing, mVoxelVolume, mCentralVeinThreshold);
        if(!mTissueFileExists)  mEffectiveDatasetVolume -= mCVVolume;
    }

    if(mPVFileExists) {
        this->ReadImageLayer(mDataSetFullFilenamePV, mPVBin, mSpacing, mPortalVeinThreshold);
        this->BuildDistanceMap(mPVBin, mPVDistMap, mSpacing, true);

        mPVVolume = this->ComputeVolumeOfRegions(mDataSetFullFilenamePV, mSpacing, mVoxelVolume, mPortalVeinThreshold);
        if(!mTissueFileExists)  mEffectiveDatasetVolume -= mPVVolume;
    }

    if(mNecRegFileExists) {
        mNecRegVolume = this->ComputeVolumeOfRegions(mDataSetFullFilenameNecReg, mSpacing, mVoxelVolume, mNecRegThreshold);
        if(mNecRegVolume<20*mVoxelVolume)   mNecRegVolume = 0;                              //necessary because of stupid text in Kai's image data
        if(!mTissueFileExists)              mEffectiveDatasetVolume -= mNecRegVolume;

        this->ReadImageLayer(mDataSetFullFilenameNecReg, mNecRegBin, mSpacing, mNecRegThreshold);
        this->BuildDistanceMap(mNecRegBin, mNecRegDistMap, mSpacing, true);
    }

    if(mLobulesFileExists) {
        this->ReadImageLayer(mDataSetFullFilenameLobules, mLobuleBin, mSpacing, mLobulesThreshold);

        if(mNecRegFileExists) {
            typename AddCImageFilterType::Pointer addLobuleNecReg = AddCImageFilterType::New();
            addLobuleNecReg->ReleaseDataFlagOn();
            addLobuleNecReg->SetInput1(mNecRegBin);
            addLobuleNecReg->SetInput2(mLobuleBin);
            addLobuleNecReg->Update();

            mLobuleBin = addLobuleNecReg->GetOutput();
            mLobuleBin->DisconnectPipeline();
        }

        if(mCVFileExists) {
            typename AddCImageFilterType::Pointer addLobuleCV = AddCImageFilterType::New();
            addLobuleCV->ReleaseDataFlagOn();
            addLobuleCV->SetInput1(mCVBin);
            addLobuleCV->SetInput2(mLobuleBin);
            addLobuleCV->Update();

            mLobuleBin = addLobuleCV->GetOutput();
            mLobuleBin->DisconnectPipeline();
        }
        mLobuleBin->SetSpacing(mSpacing);

        this->BuildDistanceMap(mLobuleBin, mLobuleDistMap, mSpacing, false);
        this->ReadImageLayer(mDataSetFullFilenameLobules, mLobuleBin, mSpacing, mLobulesThreshold); //Reset lobule bin to original lobule segmentation with necReg / cv hole
    }

    if(mTissueFileExists) {         //otherwise mEffectiveDatasetVolume defined by dataset volume in CollectBasicImageInformation()
        mTissueVolume = this->ComputeVolumeOfRegions(mDataSetFullFilenameTissue, mSpacing, mVoxelVolume, mTissueThreshold);
        mVoidVolume = this->ComputeVolumeOfRegions(mDataSetFullFilenameTissue, mSpacing, mVoxelVolume, mVoidThreshold);
        mEffectiveDatasetVolume = mTissueVolume;
    }
    else if(mLobulesFileExists) {   //in case no tissue specified, but lobule mask
        mTissueVolume = this->ComputeVolumeOfRegions(mDataSetFullFilenameLobules, mSpacing, mVoxelVolume, mLobulesThreshold);
        mVoidVolume = this->ComputeVolumeOfRegions(mDataSetFullFilenameTissue, mSpacing, mVoxelVolume, mVoidThreshold);
        mEffectiveDatasetVolume = mTissueVolume;
    }

    if(mLobulesFileExists) {
        if(mNecRegFileExists && mNecRegVolume>0)
            ComputeRelativeLobularPositions(mNecRegDistMap, mLobuleDistMap);
        else if(mCVFileExists)
            ComputeRelativeLobularPositions(mCVDistMap, mLobuleDistMap);
    }

    CollectBasicNucleiInformation();
    CollectBasicLobuleInformation();

    WriteBasicInformationFile();
    WriteLobuleInformationFile("", mLobules);
    WriteNucleiInformationFile("_hep", mHepaticNuclei);
    WriteNucleiInformationFile("_nonHep", mNonHepaticNuclei);

    if(mSaveAsGraph)
        WriteAsGraph(mPath + "hepNucleiAnalysisGraph.txt", mHepaticNuclei);

    WriteLogFile(timeStamp);

    mCVBin->ReleaseData();
    mCVBin = NULL;

    mCVDistMap->ReleaseData();
    mCVDistMap = NULL;

    mPVBin->ReleaseData();
    mPVBin = NULL;

    mPVDistMap->ReleaseData();
    mPVDistMap = NULL;

    mNecRegBin->ReleaseData();
    mNecRegBin = NULL;

    mNecRegDistMap->ReleaseData();
    mNecRegDistMap = NULL;

    mLobuleBin->ReleaseData();
    mLobuleBin = NULL;

    mLobuleDistMap->ReleaseData();
    mLobuleDistMap = NULL;

    mRelativeLobularPositionMap->ReleaseData();
    mRelativeLobularPositionMap = NULL;
}
