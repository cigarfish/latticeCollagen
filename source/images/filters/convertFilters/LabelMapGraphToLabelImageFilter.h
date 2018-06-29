///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelMapGraphToLabelImageFilter.h                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-09-30                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LABELMAPGRAPHTOLABELIMAGEFILTER_H_
#define LABELMAPGRAPHTOLABELIMAGEFILTER_H_

#include "itkImageToImageFilter.h"

//#include "itkNumericTraits.h"
//#include "itkProgressReporter.h"
//#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

template< class TInputImage, class TOutputImage > class ITK_EXPORT LabelMapGraphToLabelImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
    /** Standard class typedefs. */
    typedef LabelMapGraphToLabelImageFilter                 Self;
    typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef SmartPointer< Self >                            Pointer;
    typedef SmartPointer< const Self >                      ConstPointer;

    /** Some convenient typedefs. */
    typedef TInputImage                              InputImageType;
    typedef TOutputImage                             OutputImageType;
    typedef typename InputImageType::Pointer         InputImagePointer;
    typedef typename InputImageType::ConstPointer    InputImageConstPointer;
    typedef typename InputImageType::RegionType      InputImageRegionType;
    typedef typename InputImageType::PixelType       InputImagePixelType;
    typedef typename InputImageType::LabelObjectType LabelObjectType;

    typedef typename OutputImageType::Pointer      OutputImagePointer;
    typedef typename OutputImageType::ConstPointer OutputImageConstPointer;
    typedef typename OutputImageType::RegionType   OutputImageRegionType;
    typedef typename OutputImageType::PixelType    OutputImagePixelType;
    typedef typename OutputImageType::IndexType    IndexType;

    /** ImageDimension constants */
    itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

    /** Standard New method. */
    itkNewMacro(Self);

    /** Runtime information support. */
    itkTypeMacro(LabelMapGraphToLabelImageFilter, ImageToImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
    itkConceptMacro( SameDimensionCheck, ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
#endif

    void GenerateInputRequestedRegion();
    void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));

protected:
    LabelMapGraphToLabelImageFilter();
    ~LabelMapGraphToLabelImageFilter();

    virtual void BeforeThreadedGenerateData();
    virtual void AfterThreadedGenerateData();
    virtual void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, ThreadIdType threadId);
    virtual void ThreadedProcessLabelObject(LabelObjectType *labelObject);

    virtual InputImageType* GetLabelMap() { return static_cast< InputImageType* >( const_cast< DataObject* >(this->ProcessObject::GetInput(0)) ); };

    typename FastMutexLock::Pointer mLabelObjectContainerLock;

private:
    LabelMapGraphToLabelImageFilter(const Self &);  //purposely not implemented
    void operator=(const Self &);                   //purposely not implemented

    typename InputImageType::Iterator mLabelObjectIterator;
    ProgressReporter *mProgress;
}; // end of class

} // end namespace itk

#include "LabelMapGraphToLabelImageFilter.tpp"

#endif /* LABELMAPGRAPHTOLABELIMAGEFILTER_H_ */
