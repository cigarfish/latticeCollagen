///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  MultiChannelBinaryThresholdImageFunction.h                           //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-03-16                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef MULTICHANNELBINARYTHRESHOLDIMAGEFUNCTION_H_
#define MULTICHANNELBINARYTHRESHOLDIMAGEFUNCTION_H_

#include "itkImageFunction.h"

namespace itk
{

template <class TInputImage, class TCoordRep = float> class MultiChannelBinaryThresholdImageFunction : public ImageFunction<TInputImage,bool,TCoordRep>
{
public:
    /** Standard class typedefs. */
    typedef MultiChannelBinaryThresholdImageFunction    Self;
    typedef ImageFunction<TInputImage,bool,TCoordRep>   Superclass;
    typedef SmartPointer<Self>                          Pointer;
    typedef SmartPointer<const Self>                    ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(MultiChannelBinaryThresholdImageFunction, ImageFunction);

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** InputImageType typedef support. */
    typedef typename Superclass::InputImageType InputImageType;

    /** Typedef to describe the type of pixel. */
    typedef typename TInputImage::PixelType PixelType;

    /** Dimension underlying input image. */
    itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

    /** Point typedef support. */
    typedef typename Superclass::PointType PointType;

    /** Index typedef support. */
    typedef typename Superclass::IndexType IndexType;

    /** ContinuousIndex typedef support. */
    typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

    /** BinaryThreshold the image at a point position
     *
     * Returns true if the image intensity at the specified point position
     * satisfies the threshold criteria.  The point is assumed to lie within
     * the image buffer.
     *
     * ImageFunction::IsInsideBuffer() can be used to check bounds before
     * calling the method. */

    virtual bool Evaluate(const PointType& point) const
    {
        IndexType index;
        this->ConvertPointToNearestIndex( point, index );
        return ( this->EvaluateAtIndex( index ) );
    }

    /** BinaryThreshold the image at a continuous index position
     *
     * Returns true if the image intensity at the specified point position
     * satisfies the threshold criteria.  The point is assumed to lie within
     * the image buffer.
     *
     * ImageFunction::IsInsideBuffer() can be used to check bounds before
     * calling the method. */
    virtual bool EvaluateAtContinuousIndex(const ContinuousIndexType & index) const
    {
        IndexType nindex;

        this->ConvertContinuousIndexToNearestIndex(index, nindex);
        return this->EvaluateAtIndex(nindex);
    }

    /** BinaryThreshold the image at an index position.
     *
     * Returns true if the image intensity at the specified point position
     * satisfies the threshold criteria.  The point is assumed to lie within
     * the image buffer.
     *
     * ImageFunction::IsInsideBuffer() can be used to check bounds before
     * calling the method. */
    virtual bool EvaluateAtIndex( const IndexType & index ) const
    {
        PixelType value = this->GetInputImage()->GetPixel(index);

        for(unsigned int i=0; i<PixelType::Dimension; i++)
        {
            if( value[i]<m_Lower[i] || m_Upper[i]<value[i])
                return false;
        }
        return true;
    }

    /** Get the lower threshold value. */
    itkGetConstReferenceMacro(Lower,PixelType);

    /** Get the upper threshold value. */
    itkGetConstReferenceMacro(Upper,PixelType);

    /** Values greater than or equal to the value are inside. */
    void ThresholdAbove(PixelType thresh);

    /** Values less than or equal to the value are inside. */
    void ThresholdBelow(PixelType thresh);

    /** Values that lie between lower and upper inclusive are inside. */
    void ThresholdBetween(PixelType lower, PixelType upper);

protected:
    MultiChannelBinaryThresholdImageFunction();
    virtual ~MultiChannelBinaryThresholdImageFunction(){};
    void PrintSelf(std::ostream& os, Indent indent) const;

private:
    MultiChannelBinaryThresholdImageFunction( const Self& ); //purposely not implemented
    void operator=( const Self& ); //purposely not implemented

    PixelType m_Lower;
    PixelType m_Upper;
};

} // end namespace itk

#include "MultiChannelBinaryThresholdImageFunction.tpp"

#endif /* MULTICHANNELBINARYTHRESHOLDIMAGEFUNCTION_H_ */
