///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeCellsFilter.h                                                 //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-05-03                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ANALYZECELLSFILTER_H_
#define ANALYZECELLSFILTER_H_

#include <string>
#include <vector>

#include <QFileInfo>
#include <QString>

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkNeighborhoodIterator.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include "../../../tools/input/FilenameParser.h"


class CSParameterContext;

struct CellAnalysisContainer
{
    CellAnalysisContainer() {
        mLabel = 0;
        mX = 0;
        mY = 0;
        mZ = 0;

        mVolume = 0;
        mContactAreaWithCells = 0;
        mContactAreaWithCellsRatio = 0;
        mContactAreaWithSinusoids = 0;
        mContactAreaWithSinusoidsRatio = 0;
        mContactAreaWithBile = 0;
        mContactAreaWithBileRatio = 0;
        mContactAreaWithDataSetBorder = 0;
        mContactAreaWithDataSetBorder = 0;
        mContactAreaWithDataSetBorderRatio = 0;
        mContactAreaNil = 0;
        mContactAreaNilRatio = 0;
        mNumberNuclei = 0;
        mDistToCV = 0;
        mDistToPV = 0;
    }

    unsigned long int mLabel;

    double mX;
    double mY;
    double mZ;

    double mVolume;
    double mContactAreaWithCells;
    double mContactAreaWithCellsRatio;
    double mContactAreaWithSinusoids;
    double mContactAreaWithSinusoidsRatio;
    double mContactAreaWithBile;
    double mContactAreaWithBileRatio;
    double mContactAreaWithDataSetBorder;
    double mContactAreaWithDataSetBorderRatio;
    double mContactAreaNil;
    double mContactAreaNilRatio;

    int mNumberNuclei;
    std::vector<double> mNucleiVolume;

    double mDistToCV;
    double mDistToPV;
};



class AnalyzeCellsFilter
{
protected:
    typedef unsigned char                                                                   CScalarPixelType;
    typedef float                                                                           FScalarPixelType;
    typedef long                                                                            LScalarPixelType;
    typedef itk::Point<double, 3>                                                           PointType;

    typedef itk::Image<CScalarPixelType, 3>                                                 CScalarVoImageType;
    typedef itk::Image<FScalarPixelType, 3>                                                 FScalarVoImageType;
    typedef itk::Image<LScalarPixelType, 3>                                                 LScalarVoImageType;

    typedef itk::ImageFileWriter<CScalarVoImageType>                                        CScalarVoWriterType;

    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                                               ImageToShapeLabelMapFilterType;
    typedef itk::ConstantBoundaryCondition<LScalarVoImageType>                                                      BoundaryConditionType;
    typedef itk::ImageFileReader<CScalarVoImageType>                                                                ScalarVoReaderType;
    typedef itk::LabelMapToLabelImageFilter<ImageToShapeLabelMapFilterType::OutputImageType, LScalarVoImageType>    LabelMapToLabelImageFilterType;
    typedef itk::MinimumMaximumImageCalculator<CScalarVoImageType>                                                  CMinMaxCalculatorType;
    typedef itk::NeighborhoodIterator<LScalarVoImageType, BoundaryConditionType>                                    NeighborhoodIteratorType;
    typedef itk::SignedMaurerDistanceMapImageFilter<CScalarVoImageType, FScalarVoImageType>                         SignedMaurerDistanceMapImageFilterType;

public:
    AnalyzeCellsFilter();
    virtual ~AnalyzeCellsFilter();

    void SetParameterContext(CSParameterContext *paramContext)
    {
        mpParamContext = paramContext;
    };

    void Update();
protected:
    void ParseParameterContext();
    void PopulateCellAnalysisContainer();
    void WriteDataFile();

    void WriteAsGraph(std::string filename);

    CSParameterContext *mpParamContext;

    std::string mDataSetID;

	QString mDataSetFullFilenameCell;
    QFileInfo mInfoDataSetFullFilenameCell;
    std::string mDataSetPathCell;
    std::string mDataSetNameCell;
    std::string mDataSetFileExtensionCell;

	QString mDataSetFullFilenameNuclei;
    QFileInfo mInfoDataSetFullFilenameNuclei;
    std::string mDataSetPathNuclei;
    std::string mDataSetNameNuclei;
    std::string mDataSetFileExtensionNuclei;

	QString mDataSetFullFilenameSinusoid;
    QFileInfo mInfoDataSetFullFilenameSinusoid;
    std::string mDataSetPathSinusoid;
    std::string mDataSetNameSinusoid;
    std::string mDataSetFileExtensionSinusoid;

	QString mDataSetFullFilenameBile;
    QFileInfo mInfoDataSetFullFilenameBile;
    std::string mDataSetPathBile;
    std::string mDataSetNameBile;
    std::string mDataSetFileExtensionBile;

	QString mDataSetFullFilenameCV;
    QFileInfo mInfoDataSetFullFilenameCV;
    std::string mDataSetPathCV;
    std::string mDataSetNameCV;
    std::string mDataSetFileExtensionCV;

	QString mDataSetFullFilenamePV;
    QFileInfo mInfoDataSetFullFilenamePV;
    std::string mDataSetPathPV;
    std::string mDataSetNamePV;
    std::string mDataSetFileExtensionPV;

	bool mHasCV;
	bool mHasPV;

    CScalarVoImageType::SpacingType mSpacing;
    double mVoxelVolume;

    std::map<unsigned long int, CellAnalysisContainer> mCells;

    bool mSaveAsGraph;
};

#endif /* ANALYZECELLSFILTER_H_ */
