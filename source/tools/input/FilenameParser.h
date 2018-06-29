///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  FilenameParser.h                                                     //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  ??                                                                   //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef FILENAMEPARSER_H_
#define FILENAMEPARSER_H_

#include <string>


/*!
  \brief Simple file name parser class that decomposes a file name into fillowing substrings:
  (optional) path (/../), the actual file name and a file extension (.x*)
*/
class FilenameParser
{
public:
    FilenameParser();
    virtual ~FilenameParser();

    static bool ParseFilename(std::string fullFilename, std::string &dir, std::string &filename, std::string &ext);
    static bool IsSequence(std::string filename, std::string &prefix, int &startIndex);
};

#endif /* FILENAMEPARSER_H_ */
