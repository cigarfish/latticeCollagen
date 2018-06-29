///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  MultiChannelBinaryThresholdImageFunction.tpp                         //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-03-16                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "MultiChannelBinaryThresholdImageFunction.h"

namespace itk
{

template <class TInputImage, class TCoordRep> MultiChannelBinaryThresholdImageFunction<TInputImage,TCoordRep>::MultiChannelBinaryThresholdImageFunction()
{
    m_Lower = NumericTraits<PixelType>::NonpositiveMin();
    m_Upper = NumericTraits<PixelType>::max();
}


/**
 * Values greater than or equal to the value are inside
 */
template <class TInputImage, class TCoordRep> void MultiChannelBinaryThresholdImageFunction<TInputImage,TCoordRep>::ThresholdAbove(PixelType thresh)
{
    if (m_Lower != thresh || m_Upper != NumericTraits<PixelType>::max())
    {
        m_Lower = thresh;
        m_Upper = NumericTraits<PixelType>::max();
        this->Modified();
    }
}


/**
 * The values less than or equal to the value are inside
 */
template <class TInputImage, class TCoordRep> void MultiChannelBinaryThresholdImageFunction<TInputImage,TCoordRep>::ThresholdBelow(PixelType thresh)
{
    if (m_Lower != NumericTraits<PixelType>::NonpositiveMin() || m_Upper != thresh)
    {
        m_Lower = NumericTraits<PixelType>::NonpositiveMin();
        m_Upper = thresh;
        this->Modified();
    }
}


/**
 * The values less than or equal to the value are inside
 */
template <class TInputImage, class TCoordRep> void MultiChannelBinaryThresholdImageFunction<TInputImage,TCoordRep>::ThresholdBetween(PixelType lower, PixelType upper)
{
    if (m_Lower != lower || m_Upper != upper)
    {
        m_Lower = lower;
        m_Upper = upper;
        this->Modified();
    }
}


template <class TInputImage, class TCoordRep> void MultiChannelBinaryThresholdImageFunction<TInputImage,TCoordRep>::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf( os, indent );

    os << indent << "Lower: " << m_Lower << std::endl;
    os << indent << "Upper: " << m_Upper << std::endl;
}

} // end namespace itk

