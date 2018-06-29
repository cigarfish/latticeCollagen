///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CavityFillingImageFilter.h                                           //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef CAVITYFILLINGIMAGEFILTER_H_
#define CAVITYFILLINGIMAGEFILTER_H_

#include "itkImageToImageFilter.h"

#include "../../../tools/random/Random.h"


namespace itk
{

/** \class CavityFillingImageFilter
 */
template< class TInputImage, class TOutputImage > class ITK_EXPORT CavityFillingImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  /** Standard class typedefs. */
  typedef CavityFillingImageFilter                                Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType >   Superclass;
  typedef SmartPointer< Self >                                    Pointer;
  typedef SmartPointer< const Self >                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CavityFillingImageFilter, ImageToImageFilter);

  /** Image typedef support. */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType      InputSizeType;
  typedef typename InputImageType::SizeValueType SizeValueType;
  typedef typename InputImageType::SpacingType   InputSpacingType;

  typedef ConstantBoundaryCondition< TInputImage >                              ConstBoundaryConditionType;

  typedef ConstNeighborhoodIterator<TInputImage, ConstBoundaryConditionType>    InputNeighborhoodIteratorType;
  typedef NeighborhoodIterator<TOutputImage>                                    OutputNeighborhoodIteratorType;

  /** Returns the number of pixels that changed when the filter was executed. */
  itkGetConstReferenceMacro(NumberOfPixelsChanged, SizeValueType);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( IntConvertibleToInputCheck, ( Concept::Convertible< int, InputPixelType > ) );
  itkConceptMacro( UnsignedIntConvertibleToInputCheck, ( Concept::Convertible< unsigned int, InputPixelType > ) );
  /** End concept checking */
#endif

  void SetRadius(int radius) { m_radius = radius; };
  void SetNumberIntersectionRays(int numRays) { m_samplePoints = numRays; };
  void SetSpacing(InputSpacingType spacing) { m_spacing = spacing; };
  void SetLowerThreshold(double lower_threshold) { m_lower_threshold = lower_threshold; };
  void SetUpperThreshold(double upper_threshold) { m_upper_threshold = upper_threshold; };

protected:
  CavityFillingImageFilter();
  virtual ~CavityFillingImageFilter();
  void PrintSelf(std::ostream & os, Indent indent) const;

  virtual void BeforeThreadedGenerateData();
  void ThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread, ThreadIdType threadId);
  virtual void AfterThreadedGenerateData();

private:
  CavityFillingImageFilter(const Self &);   //purposely not implemented
  void operator=(const Self &);             //purposely not implemented

  Random *m_random;

  SizeValueType m_NumberOfPixelsChanged;
  int m_radius;
  int m_samplePoints;
  InputSpacingType m_spacing;
  double m_lower_threshold;
  double m_upper_threshold;

  // Auxiliary array for multi-threading
  Array< SizeValueType > m_Count;
};

} // end namespace itk

#include "CavityFillingImageFilter.tpp"

#endif /* CAVITYFILLINGIMAGEFILTER_H_ */
