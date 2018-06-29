////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  build_information.h                                           //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2016-07-29 12:57:10                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef BUILD_INFORMATION_H
#define BUILD_INFORMATION_H

#include <iostream>
#include <iomanip>
#include <sstream>

#include "../commit.h"


void buildInformation( char *programName )
{
    // output build information
    std::cout << std::endl;
    std::cout << std::setfill('#') << std::setw(80) << '\n';
    std::cout << std::setfill(' ');

    std::cout << "# Starting " << programName << "," << std::endl;
    //std::cout << "#  built from commit " << CS_BUILD_INFO_COMMIT << std::endl;
    //std::cout << "#          on branch " << CS_BUILD_INFO_BRANCH << std::endl;

#if defined(TAINTED_BY_DIFF)
    std::cout << "#" << std::endl;
    std::cout << "# This commit was modified by local changes:" << std::endl;
    std::cout << std::setfill('#') << std::setw(80) << '\n';
    std::cout << std::setfill(' ')
              << "#" << std::endl;
    std::cout << buildInfo_diff << "#" << std::endl;

#endif

    std::cout << std::setfill('#') << std::setw(80) << '\n';
    std::cout << std::setfill(' ') << std::endl;

};

#endif //BUILD_INFORMATION_H
