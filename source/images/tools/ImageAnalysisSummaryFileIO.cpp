///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  ImageAnalysisSummaryFileWriter.cpp                                   //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-10-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "ImageAnalysisSummaryFileIO.h"

#include <sstream>
#include <iostream>


ImageAnalysisSummaryFileIO::ImageAnalysisSummaryFileIO()
{
    // TODO Auto-generated constructor stub
}


ImageAnalysisSummaryFileIO::~ImageAnalysisSummaryFileIO()
{
    // TODO Auto-generated destructor stub
}


std::string ImageAnalysisSummaryFileIO::GetFileTypeString(FileType type)
{
    std::string id;

    switch(type)
    {
    case BileSegmentationBin:
        id = "BileSegmentationBinary";
        break;
    case BileSegmentationOverlay:
        id = "BileSegmentationOverlay";
        break;
    case BileSkeleton:
        id = "BileSkeleton";
        break;
    case BileGraph:
        id = "BileGraph";
        break;
    case SinusoidSegmentationBin:
        id = "SinusoidSegmentationBinary";
        break;
    case SinusoidSegmentationOverlay:
        id = "SinusoidSegmentationOverlay";
        break;
    case SinusoidSkeleton:
        id = "SinusoidSkeleton";
        break;
    case SinusoidGraph:
        id = "SinusoidGraph";
        break;
    case HepNucleiSegmentationBin:
        id = "HepNucleiSegmentationBinary";
        break;
    case HepNucleiSegmentationOverlay:
        id = "HepNucleiSegmentationOverlay";
        break;
    case NonHepNucleiSegmentationBin:
        id = "NonHepNucleiSegmentationBinary";
        break;
    case NonHepNucleiSegmentationOverlay:
        id = "NonHepNucleiSegmentationOverlay";
        break;
    case CentralVeinSegmentationBin:
        id = "CentralVeinSegmentationBinary";
        break;
    case CentralVeinSegmentationOverlay:
        id = "CentralVeinSegmentationOverlay";
        break;
    case PortalVeinSegmentationBin:
        id = "PortalVeinSegmentationBinary";
        break;
    case PortalVeinSegmentationOverlay:
        id = "PortalVeinSegmentationOverlay";
        break;
    case CentralVeinMaskSegmentationBin:
        id = "CentralVeinMaskSegmentationBinary";
        break;
    case CentralVeinMaskSegmentationOverlay:
        id = "CentralVeinMaskSegmentationOverlay";
        break;
    case PortalVeinMaskSegmentationBin:
        id = "PortalVeinMaskSegmentationBinary";
        break;
    case PortalVeinMaskSegmentationOverlay:
        id = "PortalVeinMaskSegmentationOverlay";
        break;
    case NecroticRegionSegmentationBin:
        id = "NecroticRegionSegmentationBin";
        break;
    case NecroticRegionSegmentationDMsOverlay:
        id = "NecroticRegionSegmentationDMsOverlay";
        break;
    case NecroticRegionSegmentationDPPIVOverlay:
        id = "NecroticRegionSegmentationDPPIVOverlay";
        break;
    case CellShapeBin:
        id = "CellShapeBinary";
        break;
    case CellShapeOverlay:
        id = "CellShapeDPPIVOverlay";
        break;
    case LobuleShapeBin:
        id = "LobuleShapeBinary";
        break;
    case LobuleShapeOverlay:
        id = "LobuleShapeOverlay";
        break;
    default:
        id = "NoValidType";
        break;
    }

    return id;
}


std::string ImageAnalysisSummaryFileIO::GetEntry(FileType type, std::string dataSetFileListPath)
{
    std::string id, idTargetPath;

    id = GetFileTypeString(type);
    std::cout << "I look for " << id << std::endl;

    std::fstream file;
    std::string line;
    std::stringstream num;
    num << id;

    file.open(dataSetFileListPath.c_str(), std::fstream::in);

    while(!file.eof() && file.is_open()) {
        getline(file, line);
//        std::cout << "Got Line " << line << " from dataSetFileList file." << std::endl;
        std::stringstream ssLine;
        ssLine << line;

        std::string firstWord;
        ssLine >> firstWord;
//        std::cout << "Got First Word " << firstWord << std::endl;

        if(firstWord.compare(num.str())==0) {
            std::cout << "Gotcha, this line starts with " << num.str() << std::endl;
            ssLine >> idTargetPath;
            break;
        }
    }
    file.close();

    return idTargetPath;
}


void ImageAnalysisSummaryFileIO::AddEntry(FileType type, std::string path, std::string fileName)
{
    std::string id;

    id = GetFileTypeString(type);

    TestForAndResolveConflict(path, id);

    std::fstream file;
    file.open((path + "dataSetFileList.ias").c_str(), std::ios::out | std::ios::app);
    file.flags(std::ios::left);

    file.width(50);
    file << id;
    file.width(100);
    file << fileName;

    file.close();
}


void ImageAnalysisSummaryFileIO::TestForAndResolveConflict(std::string path, std::string id)
{
    std::fstream file1, file2;
    std::string line;
    std::stringstream num;
    num << id;

    file1.open((path + "dataSetFileList.ias").c_str(), std::fstream::in);
    file2.open((path + "_temp_dataSetFileList.ias").c_str(), std::fstream::out);
    file2.flags(std::ios::left);

    while(!file1.eof() && file1.is_open()) {
        getline(file1, line);

        if(!line.empty()) {
//            std::cout << "Got Line " << line << " from dataSetFileList file." << std::endl;
            std::stringstream ssLine;
            ssLine << line;

            std::string firstWord;
            ssLine >> firstWord;
//            std::cout << "Got First Word " << firstWord << std::endl;

            if(firstWord.compare(num.str())!=0) {
//                std::cout << "This line didn't start with " << num.str() << std::endl;
                file2 << line << "\n";
            }
        }
    }
    file1.close();
    file2.close();

    std::remove((path + "dataSetFileList.ias").c_str());
    std::rename((path + "_temp_dataSetFileList.ias").c_str(), (path + "dataSetFileList.ias").c_str());
}
