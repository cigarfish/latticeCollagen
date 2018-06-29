///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  MultiChannelBinaryThreshodImageFilter.h                              //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-10-28                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef MULTICHANNELBINARYTHRESHOLDIMAGEFILTER_H_
#define MULTICHANNELBINARYTHRESHOLDIMAGEFILTER_H_

#include "itkUnaryFunctorImageFilter.h"


namespace itk
{

namespace Functor
{

/*!
  \brief A Functor for binary thresholding on multiple channels of an image with pixel type vector
*/
template<class TInput, class TOutput> class MultiChannelBinaryThreshold
{
public:
    MultiChannelBinaryThreshold()
    {
        m_LowerThreshold = NumericTraits<TInput>::NonpositiveMin();
        m_UpperThreshold = NumericTraits<TInput>::max();
        m_OutsideValue   = NumericTraits<TOutput>::Zero;
        m_InsideValue    = NumericTraits<TOutput>::max();
    }
    ~MultiChannelBinaryThreshold() {};

    void SetLowerThreshold(const TInput &thresh)
    {
        m_LowerThreshold = thresh;
    }
    void SetUpperThreshold(const TInput &thresh)
    {
        m_UpperThreshold = thresh;
    }
    void SetInsideValue(const TOutput &value)
    {
        m_InsideValue = value;
    }
    void SetOutsideValue(const TOutput &value)
    {
        m_OutsideValue = value;
    }

    bool operator!=( const MultiChannelBinaryThreshold & other ) const
    {
        if( m_LowerThreshold!=other.m_LowerThreshold || m_UpperThreshold!=other.m_UpperThreshold ||
                m_InsideValue!=other.m_InsideValue || m_OutsideValue!=other.m_OutsideValue )
            return true;
        return false;
    }

    bool operator==( const MultiChannelBinaryThreshold & other ) const
    {
        return !(*this!=other);
    }

    inline TOutput operator()( const TInput &A ) const
    {
        for(unsigned int i=0; i<TInput::Dimension; i++)
        {
            if( A[i]<m_LowerThreshold[i] || m_UpperThreshold[i]<A[i])
                return m_OutsideValue;
        }
        return m_InsideValue;
    }

private:
    TInput      m_LowerThreshold;
    TInput      m_UpperThreshold;
    TOutput     m_InsideValue;
    TOutput     m_OutsideValue;
};
}   //end namespace Functor


/*!
  \brief A Filter for binary thresholding on multiple channels of an image with pixel type vector

  The filter expects as input images templated over pixels of type itk::Vector<Data-Type-Which-Implements-Greater-Lesser-Equal-Operators> and outputs a binary image.
  \n With a vector of same dimensionality as the pixel type vectors lower and upper thresholds for each channel can be set.
  \n For the output image the foreground and background values can be specified (with SetInsideValue / SetOutsideValue).
  \n If the values of all channels of a pixel are within the corresponding lower and upper threshold the pixel is set to foreground otherwise to background.
*/
template <class TInputImage, class TOutputImage> class MultiChannelBinaryThresholdImageFilter :
    public UnaryFunctorImageFilter< TInputImage, TOutputImage, Functor::MultiChannelBinaryThreshold<typename TInputImage::PixelType, typename TOutputImage::PixelType> >
{
public:
    /** Standard class typedefs. */
    typedef MultiChannelBinaryThresholdImageFilter                                  Self;
    typedef UnaryFunctorImageFilter< TInputImage, TOutputImage,
            Functor::MultiChannelBinaryThreshold<typename TInputImage::PixelType,
            typename TOutputImage::PixelType> >                                     Superclass;
    typedef SmartPointer<Self>                                                      Pointer;
    typedef SmartPointer<const Self>                                                ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(MultiChannelBinaryThresholdImageFilter, UnaryFunctorImageFilter);

    /** Pixel types. */
    typedef typename TInputImage::PixelType  InputPixelType;
    typedef typename TOutputImage::PixelType OutputPixelType;

    /** Type of DataObjects to use for scalar inputs */
    typedef SimpleDataObjectDecorator<InputPixelType> InputPixelObjectType;

    /** Set the "outside" pixel value. The default value
     * NumericTraits<OutputPixelType>::Zero. */
    itkSetMacro(OutsideValue,OutputPixelType);

    /** Get the "outside" pixel value. */
    itkGetConstReferenceMacro(OutsideValue,OutputPixelType);

    /** Set the "inside" pixel value. The default value
     * NumericTraits<OutputPixelType>::max() */
    itkSetMacro(InsideValue,OutputPixelType);

    /** Get the "inside" pixel value. */
    itkGetConstReferenceMacro(InsideValue,OutputPixelType);

    /** Set the thresholds. The default lower threshold
     * is NumericTraits<InputPixelType>::NonpositiveMin(). The default upper
     * threshold is NumericTraits<InputPixelType>::max. An execption is thrown
     * if the lower threshold is greater than the upper threshold. */
    virtual void SetUpperThreshold(const InputPixelType threshold);
    virtual void SetUpperThresholdInput( const InputPixelObjectType *);
    virtual void SetLowerThreshold(const InputPixelType threshold);
    virtual void SetLowerThresholdInput( const InputPixelObjectType *);

    /** Get the threshold values. */
    virtual InputPixelType GetUpperThreshold() const;
    virtual InputPixelObjectType *GetUpperThresholdInput();
    virtual const InputPixelObjectType *GetUpperThresholdInput() const;
    virtual InputPixelType GetLowerThreshold() const;
    virtual InputPixelObjectType *GetLowerThresholdInput();
    virtual const InputPixelObjectType *GetLowerThresholdInput() const;

#ifdef ITK_USE_CONCEPT_CHECKING
    itkConceptMacro( InputHasNumericTraitsCheck, (Concept::HasNumericTraits<typename TInputImage::PixelType::ValueType>) );
    itkConceptMacro( OutputHasNumericTraitsCheck,(Concept::HasNumericTraits<typename TOutputImage::PixelType>) );
    itkConceptMacro( InputPixelTypeComparable, (Concept::Comparable<typename TInputImage::PixelType::ValueType>) );
#endif

protected:
    MultiChannelBinaryThresholdImageFilter();
    virtual ~MultiChannelBinaryThresholdImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;

    /** This method is used to set the state of the filter before
     * multi-threading. */
    virtual void BeforeThreadedGenerateData();

private:
    MultiChannelBinaryThresholdImageFilter(const Self&);    //purposely not implemented
    void operator=(const Self&);                            //purposely not implemented

    OutputPixelType     m_InsideValue;
    OutputPixelType     m_OutsideValue;
};

} // end namespace itk

#include "MultiChannelBinaryThresholdImageFilter.tpp"

#endif /* MULTICHANNELBINARYTHRESHOLDIMAGEFILTER_H_ */
