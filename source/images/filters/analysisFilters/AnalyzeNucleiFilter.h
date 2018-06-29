///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeNucleiFilter.h                                                //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-08                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ANALYZENUCLEIFILTER_H_
#define ANALYZENUCLEIFILTER_H_

#include "BasePipeline.h"

#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"

#include <string>
#include <vector>

#include "../../../tools/input/FilenameParser.h"


struct LobulesAnalysisContainer
{
    LobulesAnalysisContainer() {
        label = 0;
        volume = -1;
        volumeNecReg = -1;
        numCV = 0;
        numPV = 0;
    }
    unsigned long int label;
    double volume;
    double volumeNecReg;
    std::vector<double> volumePerRelativePositionSlice;
    unsigned int numCV;
    unsigned int numPV;
};

struct NucleiAnalysisContainer
{
    NucleiAnalysisContainer() {
        lobuleLabel = 0;
        label = 0;
        isProliferating = false;
        volume = -1;
        roundness = -1;
        elongation = -1;
        distToCV = -1;
        distToPV = -1;
        distToNecReg = -1;
        distToLobularBoundary = -1;
        numPxlOnBorder = -1;
        relativeLobularPos = -1.;

        x = 0;
        y = 0;
        z = 0;
    }

    unsigned long int lobuleLabel;
    unsigned long int label;
    bool isProliferating;
    double volume;
    double roundness;
    double elongation;
    double numPxlOnBorder;
    double relativeLobularPos;

    double x;
    double y;
    double z;

    double distToCV;
    double distToPV;
    double distToNecReg;
    double distToLobularBoundary;
};

template< unsigned int VImageDimension > class AnalyzeNucleiFilter : public BasePipeline<VImageDimension>
{
public:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

protected:
    typedef itk::Point<double, VImageDimension>                                     PointType;
    typedef typename BasePipeline<VImageDimension>::CScalarPixelType                CScalarPixelType;
    typedef typename BasePipeline<VImageDimension>::FScalarPixelType                FScalarPixelType;
    typedef typename BasePipeline<VImageDimension>::LScalarPixelType                LScalarPixelType;

    typedef typename BasePipeline<VImageDimension>::CScalarVoImageType              CScalarImageType;
    typedef typename BasePipeline<VImageDimension>::FScalarVoImageType              FScalarImageType;
    typedef typename BasePipeline<VImageDimension>::LScalarVoImageType              LScalarImageType;

    typedef typename BasePipeline<VImageDimension>::CScalarImagePointerType         CScalarImagePointerType;
    typedef typename BasePipeline<VImageDimension>::FScalarImagePointerType         FScalarImagePointerType;

    typedef typename BasePipeline<VImageDimension>::ScalarVoReaderType              CScalarImageReaderType;

    typedef itk::AddImageFilter<CScalarImageType, CScalarImageType, CScalarImageType>       AddCImageFilterType;
    typedef itk::AddImageFilter<FScalarImageType, FScalarImageType, FScalarImageType>       AddFImageFilterType;
    typedef itk::DivideImageFilter<FScalarImageType, FScalarImageType, FScalarImageType>    DivideImageFilterType;
    typedef typename BasePipeline<VImageDimension>::ImageToShapeLabelMapFilterType          ImageToShapeLabelMapFilterType;
    typedef typename BasePipeline<VImageDimension>::LabelMapToLabelImageFilterType          LabelMapToLabelImageFilterType;
    typedef typename BasePipeline<VImageDimension>::ThresholdFilterType1                    ThresholdFilterType;

public:
	AnalyzeNucleiFilter();
	virtual ~AnalyzeNucleiFilter();

	void Update();

protected:
	void ParseParameterContext();
	void WriteLogFile(std::string timeStamp);
	
	void CollectBasicImageInformation();
	void CollectBasicNucleiInformation();
	void CollectBasicLobuleInformation();

	void WriteBasicInformationFile();
	void WriteNucleiInformationFile(std::string filenamePostfix, std::vector<NucleiAnalysisContainer> nuclei);
	void WriteLobuleInformationFile(std::string filenamePostfix, std::vector<LobulesAnalysisContainer> lobules);

	void WriteAsGraph(std::string filename, std::vector<NucleiAnalysisContainer> nuclei);

	void ComputeRelativeLobularPositions(FScalarImagePointerType cvDistMap, FScalarImagePointerType lobuleDistMap);

	std::string mDataSetID;
	std::string mPath;
	std::string mFilenameExtension;

	std::string mDataSetFullFilenameNonHepNuc;
	std::string mDataSetFullFilenameHepNuc;
	std::string mDataSetFullFilenameProlifNuc;
	std::string mDataSetFullFilenameCV;
	std::string mDataSetFullFilenamePV;
	std::string mDataSetFullFilenameTissue;
	std::string mDataSetFullFilenameNecReg;
	std::string mDataSetFullFilenameLobules;

	QFileInfo mInfoFullFilename;

	int mNonProlifNucleiThreshold;
	int mProlifNucleiThreshold;
	int mCentralVeinThreshold;
	int mPortalVeinThreshold;
	int mTissueThreshold;
	int mNecRegThreshold;
	int mVoidThreshold;
	int mLobulesThreshold;

	bool mNonHepNucFileExists;
	bool mHepNucFileExists;
	bool mCVFileExists;
	bool mPVFileExists;
	bool mTissueFileExists;
	bool mNecRegFileExists;
	bool mLobulesFileExists;

	unsigned int mNumLobuleSlices;

	CScalarImagePointerType mCVBin;
	CScalarImagePointerType mPVBin;
	CScalarImagePointerType mNecRegBin;
	CScalarImagePointerType mLobuleBin;

	FScalarImagePointerType mCVDistMap;
	FScalarImagePointerType mPVDistMap;
	FScalarImagePointerType mNecRegDistMap;
	FScalarImagePointerType mLobuleDistMap;

	FScalarImagePointerType mRelativeLobularPositionMap;

	typename CScalarImageType::SpacingType mSpacing;
	typename CScalarImageType::SizeType mSize;

	long double mVoxelVolume;
	long double mDatasetVolume;
	long double mCVVolume;
    long double mPVVolume;
	long double mTissueVolume;
	long double mNecRegVolume;
	long double mVoidVolume;
	long double mEffectiveDatasetVolume;
	long double mAnalysedDatasetVolume;
	long double mNumberLobules;

	int mNumOverallNonProlifHepaticNuclei;
	int mNumOverallProlifHepaticNuclei;
	int mNumOverallNonProlifNonHepaticNuclei;
	int mNumOverallProlifNonHepaticNuclei;
	int mNumAnalysedNonProlifHepaticNuclei;
	int mNumAnalysedProlifHepaticNuclei;
	int mNumAnalysedNonProlifNonHepaticNuclei;
	int mNumAnalysedProlifNonHepaticNuclei;

	std::vector<LobulesAnalysisContainer> mLobules;
	std::vector<NucleiAnalysisContainer> mHepaticNuclei;
	std::vector<NucleiAnalysisContainer> mNonHepaticNuclei;

	bool mSaveAsGraph;
};

#include "AnalyzeNucleiFilter.tpp"

#endif /* ANALYZENUCLEIFILTER_H_ */
