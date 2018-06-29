///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  GraphVertexBaseClassifier.tpp                                        //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2015-03-23                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "GraphVertexBaseClassifier.h"

#include <iostream>
#include <sstream>
#include <typeinfo>



template< unsigned int VImageDimension > GraphVertexBaseClassifier< VImageDimension >::GraphVertexBaseClassifier()
{
    this->mTargetClass = ObjectBasedSegmentation<VImageDimension>::NUCLEUS;
}


template< unsigned int VImageDimension > GraphVertexBaseClassifier< VImageDimension >::~GraphVertexBaseClassifier()
{
    // TODO Auto-generated destructor stub
}

