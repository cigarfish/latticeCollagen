///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphEdgeBaseClassifier.tpp                                          //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2014-02-19                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "GraphEdgeBaseClassifier.h"

#include <iostream>
#include <sstream>
#include <typeinfo>



template< unsigned int VImageDimension > GraphEdgeBaseClassifier< VImageDimension >::GraphEdgeBaseClassifier()
{
    this->mTargetClass = ObjectBasedSegmentation<VImageDimension>::NUCLEUS;
}


template< unsigned int VImageDimension > GraphEdgeBaseClassifier< VImageDimension >::~GraphEdgeBaseClassifier()
{
    // TODO Auto-generated destructor stub
}

