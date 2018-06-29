///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  MultiChannelBinaryThreshodImageFilter.tpp                            //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-10-28                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "MultiChannelBinaryThresholdImageFilter.h"

#include "itkConceptChecking.h"
#include "itkSimpleDataObjectDecorator.h"


namespace itk
{

/*!
  \brief Constructor
*/
//template <class TInputImage, class TOutputImage> MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::MultiChannelBinaryThresholdImageFilter()
template<class TInputImage, class TOutputImage>MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::MultiChannelBinaryThresholdImageFilter()
{
    m_OutsideValue   = NumericTraits<OutputPixelType>::Zero;
    m_InsideValue    = NumericTraits<OutputPixelType>::max();

    // We are going to create the object with a few default inputs to
    // hold the threshold values.

    typename InputPixelObjectType::Pointer lower=InputPixelObjectType::New();
    lower->Set( NumericTraits<InputPixelType>::NonpositiveMin() );
    this->ProcessObject::SetNthInput( 1, lower );

    typename InputPixelObjectType::Pointer upper=InputPixelObjectType::New();
    upper->Set( NumericTraits<InputPixelType>::max() );
    this->ProcessObject::SetNthInput( 2, upper );
}


/*!
  \brief Set the lower threshold
  \param pixel type (vector) of input image, specifies lower threshold for every channel of the vector
*/
template <class TInputImage, class TOutputImage> void MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::SetLowerThreshold(const InputPixelType threshold)
{
    // first check to see if anything changed
    typename InputPixelObjectType::Pointer lower=this->GetLowerThresholdInput();
    if (lower && lower->Get() == threshold)
        return;

    // create a data object to use as the input and to store this
    // threshold. we always create a new data object to use as the input
    // since we do not want to change the value in any current input
    // (the current input could be the output of another filter or the
    // current input could be used as an input to several filters)
    lower = InputPixelObjectType::New();
    this->ProcessObject::SetNthInput(1, lower);

    lower->Set(threshold);
    this->Modified();
}


template <class TInputImage, class TOutputImage> void MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::SetLowerThresholdInput(const InputPixelObjectType *input)
{
    if (input != this->GetLowerThresholdInput())
    {
        this->ProcessObject::SetNthInput(1, const_cast<InputPixelObjectType*>(input));
        this->Modified();
    }
}


/*!
  \brief Get lower threshold vector
*/
template <class TInputImage, class TOutputImage> typename MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::InputPixelType
MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::GetLowerThreshold() const
{
    typename InputPixelObjectType::Pointer lower = const_cast<Self*>(this)->GetLowerThresholdInput();

    return lower->Get();
}


template <class TInputImage, class TOutputImage> typename MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::InputPixelObjectType *
MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::GetLowerThresholdInput()
{
    typename InputPixelObjectType::Pointer lower = static_cast<InputPixelObjectType *>(this->ProcessObject::GetInput(1));
    if (!lower)
    {
        // no input object available, create a new one and set it to the
        // default threshold
        lower = InputPixelObjectType::New();
        lower->Set( NumericTraits<InputPixelType>::NonpositiveMin() );
        this->ProcessObject::SetNthInput( 1, lower );
    }

    return lower;
}


template <class TInputImage, class TOutputImage> const typename MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::InputPixelObjectType *
MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::GetLowerThresholdInput() const
{
    typename InputPixelObjectType::Pointer lower = const_cast<InputPixelObjectType*>( static_cast<const InputPixelObjectType *>(this->ProcessObject::GetInput(1)) );

    if (!lower)
    {
        // no input object available, create a new one and set it to the
        // default threshold
        lower = InputPixelObjectType::New();
        lower->Set( NumericTraits<InputPixelType>::NonpositiveMin() );
        const_cast<Self*>(this)->ProcessObject::SetNthInput( 1, lower );
    }

    return lower;
}


/*!
  \brief Set the upper threshold
  \param pixel type (vector) of input image, specifies upper threshold for every channel of the vector
*/
template <class TInputImage, class TOutputImage> void MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::SetUpperThreshold(const InputPixelType threshold)
{
    // first check to see if anything changed
    typename InputPixelObjectType::Pointer upper=this->GetUpperThresholdInput();
    if (upper && upper->Get() == threshold)
        return;

    // create a data object to use as the input and to store this
    // threshold. we always create a new data object to use as the input
    // since we do not want to change the value in any current input
    // (the current input could be the output of another filter or the
    // current input could be used as an input to several filters)
    upper = InputPixelObjectType::New();
    this->ProcessObject::SetNthInput(2, upper);

    upper->Set(threshold);
    this->Modified();
}


template <class TInputImage, class TOutputImage> void MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::SetUpperThresholdInput(const InputPixelObjectType *input)
{
    if (input != this->GetUpperThresholdInput())
    {
        this->ProcessObject::SetNthInput(2, const_cast<InputPixelObjectType*>(input));
        this->Modified();
    }
}


/*!
  \brief Get upper threshold vector
*/
template <class TInputImage, class TOutputImage> typename MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::InputPixelType
MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::GetUpperThreshold() const
{
    typename InputPixelObjectType::Pointer upper = const_cast<Self*>(this)->GetUpperThresholdInput();

    return upper->Get();
}


template <class TInputImage, class TOutputImage> typename MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::InputPixelObjectType *
MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::GetUpperThresholdInput()
{
    typename InputPixelObjectType::Pointer upper = static_cast<InputPixelObjectType *>(this->ProcessObject::GetInput(2));
    if (!upper)
    {
        // no input object available, create a new one and set it to the
        // default threshold
        upper = InputPixelObjectType::New();
        upper->Set( NumericTraits<InputPixelType>::max() );
        this->ProcessObject::SetNthInput( 2, upper );
    }

    return upper;
}


template <class TInputImage, class TOutputImage> const typename MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::InputPixelObjectType *
MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::GetUpperThresholdInput() const
{
    typename InputPixelObjectType::Pointer upper = const_cast<InputPixelObjectType*>( static_cast<const InputPixelObjectType *>(this->ProcessObject::GetInput(2)) );

    if (!upper)
    {
        // no input object available, create a new one and set it to the
        // default threshold
        upper = InputPixelObjectType::New();
        upper->Set( NumericTraits<InputPixelType>::max() );
        const_cast<Self*>(this)->ProcessObject::SetNthInput( 2, upper );
    }

    return upper;
}


/*!
  \brief Print self method, prints Outside/Inside values and Lower/Upper thresholds
*/
template <class TInputImage, class TOutputImage> void MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "OutsideValue: " << static_cast<typename NumericTraits<OutputPixelType>::PrintType>(m_OutsideValue) << std::endl;
    os << indent << "InsideValue: " << static_cast<typename NumericTraits<OutputPixelType>::PrintType>(m_InsideValue) << std::endl;
    os << indent << "LowerThreshold: " << static_cast<typename NumericTraits<InputPixelType>::PrintType>(this->GetLowerThreshold()) << std::endl;
    os << indent << "UpperThreshold: " << static_cast<typename NumericTraits<InputPixelType>::PrintType>(this->GetUpperThreshold()) << std::endl;
}


template <class TInputImage, class TOutputImage> void MultiChannelBinaryThresholdImageFilter<TInputImage, TOutputImage>::BeforeThreadedGenerateData()
{
    // set up the functor values
    typename InputPixelObjectType::Pointer lowerThreshold = this->GetLowerThresholdInput();
    typename InputPixelObjectType::Pointer upperThreshold = this->GetUpperThresholdInput();

    for(unsigned int i=0; i<TInputImage::PixelType::Dimension; i++)
    {
        if(lowerThreshold->Get()[i] > upperThreshold->Get()[i])
            itkExceptionMacro(<<"A lower threshold cannot be greater than an upper threshold, in no dimension of the vector.");
    }

    // Setup up the functor
    this->GetFunctor().SetLowerThreshold( lowerThreshold->Get() );
    this->GetFunctor().SetUpperThreshold( upperThreshold->Get() );

    this->GetFunctor().SetInsideValue( m_InsideValue );
    this->GetFunctor().SetOutsideValue( m_OutsideValue );

}

} // end namespace itk

