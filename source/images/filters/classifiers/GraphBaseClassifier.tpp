///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphBaseClassifier.tpp                                              //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-03-23                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "GraphBaseClassifier.h"

#include <iostream>
#include <sstream>
#include <typeinfo>



template< unsigned int VImageDimension > GraphBaseClassifier< VImageDimension >::GraphBaseClassifier()
{
    mTargetClass = ObjectBasedSegmentation<VImageDimension>::NUCLEUS;
}


template< unsigned int VImageDimension > GraphBaseClassifier< VImageDimension >::~GraphBaseClassifier()
{
    // TODO Auto-generated destructor stub
}


template< unsigned int VImageDimension > void GraphBaseClassifier< VImageDimension >::Init()
{
    mNumActiveUnaryFeatures = mLabelMapGraph->GetNumberActiveUnaryFeatures();
    mNumActiveUnaryFeatureComponents = mLabelMapGraph->GetNumberActiveUnaryFeatureComponents();
    mNumActiveBinaryFeatures = mLabelMapGraph->GetNumberActiveBinaryFeatures();
}
