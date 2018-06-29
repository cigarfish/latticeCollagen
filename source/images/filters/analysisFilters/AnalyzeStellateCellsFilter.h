///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AnalyzeStellateCellsFilter.h                                         //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-01-30                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ANALYZESTELLATECELLSFILTER_H_
#define ANALYZESTELLATECELLSFILTER_H_

#include <string>
#include <vector>

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageFileReader.h"
#include "itkLabelMap.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include "../../../tools/input/FilenameParser.h"

class CSParameterContext;

struct StellateCell
{
    StellateCell() {
        mLabel = 0;
        distToCV = 0;
        distToPV = 0;
    }

    unsigned long int mLabel;
    double distToCV;
    double distToPV;
};


class AnalyzeStellateCellsFilter
{
protected:
    typedef unsigned char                               CScalarPixelType;
    typedef float                                       FScalarPixelType;

    typedef itk::Image<CScalarPixelType, 3>             CScalarVoImageType;
    typedef itk::Image<FScalarPixelType, 3>             FScalarVoImageType;

    typedef itk::ImageRegionConstIterator<FScalarVoImageType>   ImageIteratorType;

    typedef itk::ImageFileReader<CScalarVoImageType>            ScalarVoReaderType;

    typedef itk::BinaryImageToShapeLabelMapFilter<CScalarVoImageType>                           ImageToShapeLabelMapFilterType;
    typedef itk::SignedMaurerDistanceMapImageFilter<CScalarVoImageType, FScalarVoImageType>     SignedMaurerDistanceMapImageFilterType;

public:
    AnalyzeStellateCellsFilter();
    virtual ~AnalyzeStellateCellsFilter();

    void SetParameterContext(CSParameterContext *paramContext)
    {
        m_paramContext = paramContext;
    };

    void Update();

protected:
    void ParseParameterContext();

    void CollectBasicImageInformation();
    void CollectBasicNucleiInformation(bool withCV, bool withPV);
    void CollectBasicStellateCellInformation();
    void CollectBasicHepatocyteInformation();

    void WriteBasicInformationFile();
    void WriteBasicNucleiInformationFile();
    void WriteBasicHepatocyteInformationFile();

    long double ComputeVolumeOfRegions(std::string path, std::string filename, std::string ext);

    CSParameterContext *m_paramContext;

    std::string m_dataSetID;

    std::string m_dataSetFullFilenameNonHepNuc;
    std::string m_dataSetPathNonHepNuc;
    std::string m_dataSetNameNonHepNuc;
    std::string m_dataSetFileExtensionNonHepNuc;

    std::string m_dataSetFullFilenameHepNuc;
    std::string m_dataSetPathHepNuc;
    std::string m_dataSetNameHepNuc;
    std::string m_dataSetFileExtensionHepNuc;

    std::string m_dataSetFullFilenameStCNuc;
    std::string m_dataSetPathStCNuc;
    std::string m_dataSetNameStCNuc;
    std::string m_dataSetFileExtensionStCNuc;

    std::string m_dataSetFullFilenameStCBody;
    std::string m_dataSetPathStCBody;
    std::string m_dataSetNameStCBody;
    std::string m_dataSetFileExtensionStCBody;

    std::string m_dataSetFullFilenameHepCell;
    std::string m_dataSetPathHepCell;
    std::string m_dataSetNameHepCell;
    std::string m_dataSetFileExtensionHepCell;

    std::string m_dataSetFullFilenameCV;
    std::string m_dataSetPathCV;
    std::string m_dataSetNameCV;
    std::string m_dataSetFileExtensionCV;

    std::string m_dataSetFullFilenamePV;
    std::string m_dataSetPathPV;
    std::string m_dataSetNamePV;
    std::string m_dataSetFileExtensionPV;

    long int m_num_dim;
    long int *m_dim;
    CScalarVoImageType::SpacingType m_spacing;

    unsigned int m_numRings;
    unsigned char m_normVeinDistRingStepSize;
    std::map<unsigned int, double> m_normCVDistRings;
    std::map<unsigned int, double> m_normPVDistRings;

    long double m_voxel_volume;
    long double m_dataset_volume;
    long double m_cv_volume;
    long double m_pv_volume;
    long double m_effective_dataset_volume;

    double m_volume_stellateCellBodies;
    double m_volume_stellateCellBodiesPerEffVol;

    int m_num_hepaticNuclei;
    int m_num_nonHepaticNuclei;
    int m_num_stellateCellNuclei;
    int m_num_stellateCellBodies;
    int m_num_hepatocytes;

    double m_num_hepaticNucleiPerEffVol;
    double m_num_nonHepaticNucleiPerEffVol;
    double m_num_stellateCellNucleiPerEffVol;
    double m_num_stellateCellBodiesPerEffVol;
    double m_num_hepatocytesPerEffVol;

    std::map<unsigned long int, StellateCell> m_stellateCells;
    std::vector<double> m_hepatocyteVolume;
};

#endif /* ANALYZESTELLATECELLSFILTER_H_ */
