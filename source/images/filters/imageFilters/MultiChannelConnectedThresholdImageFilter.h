///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  MultiChannelConnectedThresholdImageFilter.h                          //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-03-16                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef MULTICHANNELCONNECTEDTHRESHOLDIMAGEFILTER_H_
#define MULTICHANNELCONNECTEDTHRESHOLDIMAGEFILTER_H_

#include "itkImageToImageFilter.h"


namespace itk
{

/*!
  \brief Label pixels that are connected to a seed and lie within a range of values.

  ConnectedThresholdImageFilter labels pixels with ReplaceValue that are connected to an initial Seed AND lie within a Lower and Upper threshold range.
  Same as itk::ConnectedThresholdImageFilter but works also on images with vector pixel values:
  The filter expects as input images templated over pixels of type itk::Vector<Data-Type-Which-Implements-Greater-Lesser-Equal-Operators> and outputs a binary image.
  \n With a vector of same dimensionality as the pixel type vectors lower and upper thresholds for each channel can be set.
  \n For the output image a foreground value can be specified (with SetReplaceValue).
  \n If the values of all channels of a pixel are within the corresponding lower and upper threshold the pixel is set to replace value.

   RegionGrowingSegmentation
*/
template <class TInputImage, class TOutputImage> class MultiChannelConnectedThresholdImageFilter : public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef MultiChannelConnectedThresholdImageFilter       Self;
    typedef ImageToImageFilter<TInputImage,TOutputImage>    Superclass;
    typedef SmartPointer<Self>                              Pointer;
    typedef SmartPointer<const Self>                        ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(MultiChannelConnectedThresholdImageFilter, ImageToImageFilter);

    typedef TInputImage                           InputImageType;
    typedef typename InputImageType::Pointer      InputImagePointer;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename InputImageType::RegionType   InputImageRegionType;
    typedef typename InputImageType::PixelType    InputImagePixelType;
    typedef typename InputImageType::IndexType    IndexType;
    typedef typename InputImageType::SizeType     SizeType;

    typedef TOutputImage                          OutputImageType;
    typedef typename OutputImageType::Pointer     OutputImagePointer;
    typedef typename OutputImageType::RegionType  OutputImageRegionType;
    typedef typename OutputImageType::PixelType   OutputImagePixelType;

    void PrintSelf ( std::ostream& os, Indent indent ) const;


    /** Set seed point. */
    void SetSeed(const IndexType &seed);
    void AddSeed(const IndexType &seed);

    /** Clear the seed list. */
    void ClearSeeds ();

    /** Set/Get value to replace thresholded pixels. Pixels that lie *
     *  within Lower and Upper (inclusive) will be replaced with this
     *  value. The default is 1. */
    itkSetMacro(ReplaceValue, OutputImagePixelType);
    itkGetConstMacro(ReplaceValue, OutputImagePixelType);

    /** Type of DataObjects to use for scalar inputs */
    typedef SimpleDataObjectDecorator<InputImagePixelType> InputPixelObjectType;

    /** Set Upper and Lower Threshold inputs as values */
    virtual void SetUpper( InputImagePixelType );
    virtual void SetLower( InputImagePixelType );

    /** Set Threshold inputs that are connected to the pipeline */
    virtual void SetUpperInput( const InputPixelObjectType *);
    virtual void SetLowerInput( const InputPixelObjectType *);

    /** Get Upper and Lower Threshold inputs as values */
    virtual InputImagePixelType GetUpper() const;
    virtual InputImagePixelType GetLower() const;

    /** Get Threshold inputs that are connected to the pipeline */
    virtual InputPixelObjectType * GetUpperInput();
    virtual InputPixelObjectType * GetLowerInput();

    /** Image dimension constants */
    itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

#ifdef ITK_USE_CONCEPT_CHECKING
    /** Begin concept checking */
    itkConceptMacro( OutputEqualityComparableCheck, (Concept::EqualityComparable<OutputImagePixelType>) );
    itkConceptMacro( InputEqualityComparableCheck, (Concept::EqualityComparable<InputImagePixelType>) );
//    itkConceptMacro( InputHasNumericTraitsCheck, (Concept::HasNumericTraits<InputImagePixelType::ValueType>) );
//    itkConceptMacro( InputPixelTypeComparable, (Concept::Comparable<InputImagePixelType::ValueType>) );
    itkConceptMacro( OutputHasNumericTraitsCheck, (Concept::HasNumericTraits<OutputImagePixelType>) );
    itkConceptMacro( OutputOStreamWritableCheck, (Concept::OStreamWritable<OutputImagePixelType>) );
    /** End concept checking */
#endif

    /** Face connectivity is 4 connected in 2D, 6  connected in 3D, 2*n   in ND
     *  Full connectivity is 8 connected in 2D, 26 connected in 3D, 3^n-1 in ND
     *  Default is to use FaceConnectivity. */
    typedef enum { FaceConnectivity, FullConnectivity } ConnectivityEnumType;

protected:
    MultiChannelConnectedThresholdImageFilter();
    virtual ~MultiChannelConnectedThresholdImageFilter() {};
    std::vector<IndexType> m_SeedList;
    InputImagePixelType    m_Lower;
    InputImagePixelType    m_Upper;
    OutputImagePixelType   m_ReplaceValue;

    // Override since the filter needs all the data for the algorithm
    void GenerateInputRequestedRegion();

    // Override since the filter produces the entire dataset
    void EnlargeOutputRequestedRegion(DataObject *output);

    void GenerateData();

    // Type of connectivity to use.
    ConnectivityEnumType m_Connectivity;

private:
    MultiChannelConnectedThresholdImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#include "MultiChannelConnectedThresholdImageFilter.tpp"

#endif /* MULTICHANNELCONNECTEDTHRESHOLDIMAGEFILTER_H_ */
