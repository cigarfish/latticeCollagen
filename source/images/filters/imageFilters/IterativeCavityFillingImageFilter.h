///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  IterativeCavityFillingImageFilter.h                                  //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-13                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ITERATIVECAVITYFILLINGIMAGEFILTER_H_
#define ITERATIVECAVITYFILLINGIMAGEFILTER_H_

#include "itkImageToImageFilter.h"

#include "CavityFillingImageFilter.h"


namespace itk
{

template<class TImage> class ITK_EXPORT IterativeCavityFillingImageFilter : public ImageToImageFilter<TImage, TImage>
{
public:

  /** Convenient typedefs for simplifying declarations. */
  typedef TImage InputImageType;
  typedef TImage OutputImageType;

  /** Standard class typedefs. */
  typedef IterativeCavityFillingImageFilter                     Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(IterativeCavityFillingImageFilter, ImageToImageFilter);

  /** Type of the internal Voting filter that is going to be executed
    iteratively */
  typedef CavityFillingImageFilter<InputImageType, OutputImageType> CavityFilterType;

  /** Image typedef support. */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType     InputSizeType;
  typedef typename InputImageType::SpacingType  InputSpacingType;

  /** Maximum number of iterations. This filter is executed iteratively as
   * long as at least one pixel has changed in a previous iteration, or until
   * the maximum number of iterations has been reached. */
  itkGetConstReferenceMacro(MaximumNumberOfIterations, unsigned int);
  itkSetMacro(MaximumNumberOfIterations, unsigned int);

  void SetRadius(int radius) { m_radius = radius; };
  void SetSpacing(InputSpacingType spacing) { m_spacing = spacing; };
  void SetLowerThreshold(double threshold) { m_lower_threshold = threshold; };
  void SetUpperThreshold(double threshold) { m_upper_threshold = threshold; };

  /** Number of iterations executed at any given time. This is useful at the
   * end of the execution in order to verify how many iterations were
   * performed.  */
  itkGetConstReferenceMacro(CurrentNumberOfIterations, unsigned int);
  itkSetMacro(CurrentNumberOfIterations, unsigned int);

  /** Returns the number of pixels that changed when the filter was executed. */
  itkGetConstReferenceMacro(NumberOfPixelsChanged, unsigned int);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputEqualityComparableCheck, ( Concept::EqualityComparable< InputPixelType > ) );
  itkConceptMacro( InputOStreamWritableeCheck, ( Concept::OStreamWritable< InputPixelType > ) );
  /** End concept checking */
#endif
protected:
  IterativeCavityFillingImageFilter();
  virtual ~IterativeCavityFillingImageFilter();
  void PrintSelf(std::ostream & os, Indent indent) const;

  void GenerateData();

private:
  IterativeCavityFillingImageFilter(const Self &);  //purposely not implemented
  void operator=(const Self &);                     //purposely not


  SizeValueType m_NumberOfPixelsChanged;
  int m_radius;
  InputSpacingType m_spacing;
  double m_lower_threshold;
  double m_upper_threshold;

  unsigned int m_MaximumNumberOfIterations;
  unsigned int m_CurrentNumberOfIterations;
};

} // end namespace itk

#include "IterativeCavityFillingImageFilter.tpp"

#endif /* ITERATIVECAVITYFILLINGIMAGEFILTER_H_ */
