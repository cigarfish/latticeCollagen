///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelMapGraphToLabelImageFilter.tpp                                  //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-09-30                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "LabelMapGraphToLabelImageFilter.h"

#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

template< class TInputImage, class TOutputImage > LabelMapGraphToLabelImageFilter< TInputImage, TOutputImage >::LabelMapGraphToLabelImageFilter()
{
    mProgress = NULL;
}


template< class TInputImage, class TOutputImage > LabelMapGraphToLabelImageFilter< TInputImage, TOutputImage >::~LabelMapGraphToLabelImageFilter()
{
    // be sure that the progress reporter has been destroyed
    if(mProgress != NULL)
    {
        delete mProgress;
    }
}


template< class TInputImage, class TOutputImage > void LabelMapGraphToLabelImageFilter< TInputImage, TOutputImage >::GenerateInputRequestedRegion()
{
    // call the superclass' implementation of this method
    Superclass::GenerateInputRequestedRegion();

    // We need all the input.
    InputImagePointer input = const_cast< InputImageType* >(this->GetInput());

    if(!input)
        return;

    input->SetRequestedRegion(input->GetLargestPossibleRegion());
}


template< class TInputImage, class TOutputImage > void LabelMapGraphToLabelImageFilter< TInputImage, TOutputImage >::EnlargeOutputRequestedRegion(DataObject *)
{
    this->GetOutput()->SetRequestedRegion(this->GetOutput()->GetLargestPossibleRegion());
}


template< class TInputImage, class TOutputImage > void LabelMapGraphToLabelImageFilter< TInputImage, TOutputImage >::AfterThreadedGenerateData()
{
  // destroy progress reporter
  delete mProgress;
  mProgress = NULL;
}


template< class TInputImage, class TOutputImage > void LabelMapGraphToLabelImageFilter< TInputImage, TOutputImage >::ThreadedGenerateData(const OutputImageRegionType &, ThreadIdType itkNotUsed(threadId))
{
    while(true)
    {
        // first lock the mutex
        mLabelObjectContainerLock->Lock();

        if(mLabelObjectIterator.IsAtEnd())
        {
            // no more objects. Release the lock and return
            mLabelObjectContainerLock->Unlock();
            return;
        }

        // get the label object
        LabelObjectType *labelObject = mLabelObjectIterator.GetLabelObject();

        // increment the iterator now, so it will not be invalidated if the object
        // is destroyed
        ++mLabelObjectIterator;

        // pretend one more object is processed, even if it will be done later, to
        // simplify the lock management
        mProgress->CompletedPixel();

        // unlock the mutex, so the other threads can get an object
        mLabelObjectContainerLock->Unlock();

        // and run the user defined method for that object
        this->ThreadedProcessLabelObject(labelObject);
    }
}


template< class TInputImage, class TOutputImage > void LabelMapGraphToLabelImageFilter< TInputImage, TOutputImage >::BeforeThreadedGenerateData()
{
    OutputImageType      *output = this->GetOutput();
    const InputImageType *input = this->GetInput();

    output->FillBuffer(input->GetBackgroundValue());

    // initialize the iterator
    mLabelObjectIterator =  typename InputImageType::Iterator(this->GetLabelMap());
    //  m_LabelObjectIterator = typename InputImageType::Iterator(this->GetLabelMap());

    // and the mutex
    mLabelObjectContainerLock = FastMutexLock::New();

    // be sure that the previous progress reporter has been destroyed
    if(mProgress != NULL)
    {
        delete mProgress;
    }
    // initialize the progress reporter
    mProgress = new ProgressReporter(this, 0, this->GetLabelMap()->GetNumberOfLabelObjects());
}


template< class TInputImage, class TOutputImage > void LabelMapGraphToLabelImageFilter< TInputImage, TOutputImage >::ThreadedProcessLabelObject(LabelObjectType *labelObject)
{
    const typename LabelObjectType::LabelType &label = labelObject->GetLabel();
    typename LabelObjectType::ConstIndexIterator it(labelObject);
    while(!it.IsAtEnd())
    {
        this->GetOutput()->SetPixel(it.GetIndex(), label);
        ++it;
    }
}

} // end namespace itk

