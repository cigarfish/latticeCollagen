////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  ModelElementHollowSphere.h                                    //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-09-09 12:04:28                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef MODEL_ELEMENT_HOLLOWSPHERE_H
#define MODEL_ELEMENT_HOLLOWSPHERE_H

#include "ModelElement.h"

class ModelElementSphere;

#include <vector>
#include <sstream>

namespace H5 { class CompType; };


class ModelElementHollowSphere : public ModelElement
{
public:
    ModelElementHollowSphere( double x, double y, double z );

    BoundingBox * boundingBox();

    void HDF5DataFormat( H5::CompType & );
    H5::CompType ParseHDF5DataFormat( H5::CompType & inputType,
                              std::stringstream & errors,
                              std::stringstream & warnings );


    // Method to adapt the radius according to pressure [ =mLastForceAbsolute/(2*PI*r^2) ]
    void Update();

    std::vector<ModelElement *> mInteractingElements;

    double mRadius;
    double mPressure;
    double mStressResponse;
};

#endif // MODEL_ELEMENT_HOLLOWSPHERE_H
