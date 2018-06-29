///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  About.cpp                                                            //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-06-10                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "About.h"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <sstream>

#include "Version.h"


// Constructor
About::About(QWidget *parent) : QWidget(parent)
{
    // Setting up the Qt Designer code
    setupUi(this);

    std::stringstream versionNumber;
    std::stringstream developerEnum;

#ifdef CS_TI_QUANT_ONLY
    versionNumber << "TiQuant " << VERSION_TIQUANT_MAJOR << "." << VERSION_TIQUANT_MINOR << "." << VERSION_TIQUANT_PATCH << std::endl;

#elif CS_TI_SIM_ONLY
    versionNumber << "TiSim " << VERSION_TISIM_MAJOR << "." << VERSION_TISIM_MINOR << "." << VERSION_TISIM_PATCH << std::endl;

#else
    versionNumber << "CellSys " << VERSION_CELLSYS_MAJOR << "." << VERSION_CELLSYS_MINOR << "." << VERSION_CELLSYS_PATCH << std::endl;

#endif

    labelRelease->setText(QString::fromStdString(versionNumber.str()));

#ifdef CS_TI_QUANT_ONLY
    developerEnum << "Developers (listed alphabetically): Drasdo, D., Friebel, A., Hoehme, S., Johann, T.";

#elif CS_TI_SIM_ONLY
    developerEnum << "Developers (listed alphabetically): Drasdo, D., Hoehme, S., Johann, T., Neitsch, J.";

#else
    developerEnum << "Developers (listed alphabetically): Drasdo, D., Friebel, A., Hoehme, S., Johann, T., Neitsch, J.";

#endif

    labelDevelopers->setText(QString::fromStdString(developerEnum.str()));

    std::ifstream file("./LICENSE.txt");
    std::ostringstream licenseText;
    copy(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), std::ostreambuf_iterator<char>(licenseText));

    licenseBrowser->setText(QString::fromStdString(licenseText.str()));
}


About::~About()
{
}

