///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ObjectBasedSegmentation.h                                            //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-07-31                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef OBJECTBASEDSEGMENTATION_H_
#define OBJECTBASEDSEGMENTATION_H_

#include "LabelMapGraphBasePipeline.h"


template< unsigned int VImageDimension > class ObjectBasedSegmentation : public LabelMapGraphBasePipeline<VImageDimension>
{
public:
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

    static const int NumberClasses = 6;
    static const int NumberSegmentationModes = 2;
    static const int NumberClassifiers = 3;

    enum Class {
        BACKGROUND,
        NUCLEUS,
        SINUSOID,
        BILE,
        VEIN,
        NECROTIC_REGION
    };

    enum SegmentationMode {
        LabelObjects,
        LabelObjectPairs
    };

    enum Classifier {
        WithoutSegmentation,
        KNNBasedMerging,
        PythonSVMBasedMerging
    };

    static const std::string ClassName[];
    static const std::string SegmentationModeName[];
    static const std::string ClassifierName[];

private:
    typedef LabelMapGraphBasePipeline<VImageDimension>                  Superclass;

    typedef typename Superclass::SpacingType                            SpacingType;

    typedef typename Superclass::LabelMapType                           LabelMapType;
    typedef typename Superclass::UnaryFeatureType                       UnaryFeatureType;
    typedef typename Superclass::BinaryFeatureType                      BinaryFeatureType;

    typedef typename Superclass::LabelImageToGraphFilterType            LabelImageToGraphFilterType;

public:
    ObjectBasedSegmentation();
    virtual ~ObjectBasedSegmentation();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);
    void WriteDataSetSummary();

    void ActivateFeaturesAtLabelMap();

    //Todo: preliminary name: is supposed to be the method, that handles labelobject / pair classification, incorporation into error function + doing whatever it takes to merge everything
    void ClassifyLabelObjects();
    void ClassifierBasedMerging();

    std::string mTrainingDatasetFilename;
    std::string mTrainingDatabaseFilename;

    bool* mWithFeature;

    bool mKeepOnlyTargetClassLabelObjects;

    SegmentationMode mSegmentationMode;
    Class mSegmentationTargetClass;
    Classifier mSegmentationClassifier;

    bool mSaveEverything;                   //1 - save everything, 0 - save only essentials
    std::string mLogFilenameSave;
    std::string mFilenameSave;
};

#include "ObjectBasedSegmentation.tpp"

#endif /* OBJECTBASEDSEGMENTATION_H_ */
