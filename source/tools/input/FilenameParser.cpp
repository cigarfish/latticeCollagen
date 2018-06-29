///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  FilenameParser.cpp                                                   //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  ??                                                                   //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "FilenameParser.h"

#include <stdlib.h>


FilenameParser::FilenameParser()
{
    // TODO Auto-generated constructor stub
}


FilenameParser::~FilenameParser()
{
    // TODO Auto-generated destructor stub
}


//! Method to decompose a string into directory(optional), filename and file extension(mandatory) sub-strings
bool FilenameParser::ParseFilename(std::string fullFilename, std::string &dir, std::string &filename, std::string &ext)
{
    size_t lastSlash, lastDot;

    lastSlash = fullFilename.find_last_of("/\\");
    if(lastSlash == std::string::npos)
        lastSlash = -1;
    dir = fullFilename.substr(0, lastSlash+1);

    lastDot = fullFilename.find_last_of(".");
    if(lastDot == std::string::npos)
        return false;
    ext = fullFilename.substr(lastDot);

    int nameLength = fullFilename.length() - dir.length() - ext.length();
    if(nameLength == 0 || dir.compare("/path/to/")==0)
        return false;
    filename = fullFilename.substr(lastSlash+1, nameLength);

    return true;
}


bool FilenameParser::IsSequence(std::string filename, std::string &prefix, int &startIndex)
{
    std::string suffix;
    char *ep;
    bool foundIndex = false;
    startIndex = 0;
    int length = 0;
    for(unsigned int i=1; i<filename.length(); i++) {
        suffix = filename.substr(filename.length()-i, i);
        int j = strtol(suffix.c_str(),&ep,10);

        if(*ep!=NULL)
            break;
        else {
            foundIndex = true;
            startIndex = j;
            length++;
        }
    }
    prefix = filename.substr(0, filename.length()-length);

    return foundIndex;
}
