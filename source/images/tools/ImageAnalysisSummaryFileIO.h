///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ImageAnalysisSummaryFileIO.h                                         //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-10-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef IMAGEANALYSISSUMMARYFILEIO_H_
#define IMAGEANALYSISSUMMARYFILEIO_H_

#include <fstream>
#include <string>


enum FileType
{
    BileSegmentationBin,
    BileSegmentationOverlay,
    BileSkeleton,
    BileGraph,
    SinusoidSegmentationBin,
    SinusoidSegmentationOverlay,
    SinusoidSkeleton,
    SinusoidGraph,
    NonHepNucleiSegmentationBin,
    NonHepNucleiSegmentationOverlay,
    HepNucleiSegmentationBin,
    HepNucleiSegmentationOverlay,
    CentralVeinSegmentationBin,
    CentralVeinSegmentationOverlay,
    PortalVeinSegmentationBin,
    PortalVeinSegmentationOverlay,
    CentralVeinMaskSegmentationBin,
    CentralVeinMaskSegmentationOverlay,
    PortalVeinMaskSegmentationBin,
    PortalVeinMaskSegmentationOverlay,
    NecroticRegionSegmentationBin,
    NecroticRegionSegmentationDMsOverlay,
    NecroticRegionSegmentationDPPIVOverlay,
    CellShapeBin,
    CellShapeOverlay,
    LobuleShapeBin,
    LobuleShapeOverlay
};


class ImageAnalysisSummaryFileIO
{
public:
    ImageAnalysisSummaryFileIO();
    virtual ~ImageAnalysisSummaryFileIO();

    //add entry in file list, params: path - path to file list file, fileName - full filename including path and extension of file that will be added in list
    static void AddEntry(FileType type, std::string path, std::string fileName);
    static std::string GetEntry(FileType type, std::string dataSetFileListPath);
private:
    static std::string GetFileTypeString(FileType type);
    static void TestForAndResolveConflict(std::string path, std::string id);
};


#endif /* IMAGEANALYSISSUMMARYFILEIO_H_ */
