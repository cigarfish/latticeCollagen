///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelImageToLabelMapGraphFilter.h                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-09-04                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LabelImageToLabelMapGraphFilter_H_
#define LabelImageToLabelMapGraphFilter_H_

#include <itkLabelImageToLabelMapFilter.h>

#include "../../tools/FeatureLabelObject.h"
#include "../../tools/LabelMapGraph.h"

namespace itk
{

template< class TOriginalImage, class TInputImage, class TOutputImage = LabelMapGraph< TOriginalImage, TInputImage, FeatureLabelObject< typename TInputImage::PixelType, ::itk::GetImageDimension< TInputImage >::ImageDimension > > >
class ITK_EXPORT LabelImageToLabelMapGraphFilter : public LabelImageToLabelMapFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef LabelImageToLabelMapGraphFilter                           Self;
  typedef LabelImageToLabelMapFilter< TInputImage, TOutputImage >   Superclass;
  typedef SmartPointer< Self >                                      Pointer;
  typedef SmartPointer< const Self >                                ConstPointer;

  /** Some convenient typedefs. */
  typedef TOriginalImage                        OriginalImageType;
  typedef TInputImage                           InputImageType;
  typedef TOutputImage                          OutputImageType;

  typedef typename OriginalImageType::Pointer   OriginalImagePointer;

  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputImagePixelType;
  typedef typename InputImageType::IndexType    IndexType;

  typedef typename OutputImageType::Pointer         OutputImagePointer;
  typedef typename OutputImageType::ConstPointer    OutputImageConstPointer;
  typedef typename OutputImageType::RegionType      OutputImageRegionType;
  typedef typename OutputImageType::PixelType       OutputImagePixelType;
  typedef typename OutputImageType::LabelObjectType LabelObjectType;
  typedef typename LabelObjectType::LengthType      LengthType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(LabelImageToLabelMapGraphFilter, LabelImageToLabelMapFilter);

  /**
   * Set/Get the value used as "background" in the output image.
   * Defaults to NumericTraits<PixelType>::NonpositiveMin().
   */
  void SetBackgroundValue(OutputImagePixelType b) { mBackgroundValue = b; };
  OutputImagePixelType GetBackgroundValue() { return mBackgroundValue; };

  void SetOriginalImage(OriginalImagePointer image) { mOriginalImage = image; };
  void SetLabelImage(InputImagePointer image) { mLabelImage = image; };
  void UsePrecalculatedGraph(vtkSmartPointer<vtkUndirectedGraph>  graph);

#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro( SameDimensionCheck, ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
#endif

protected:
  LabelImageToLabelMapGraphFilter();
  ~LabelImageToLabelMapGraphFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** LabelImageToLabelMapFilter needs the entire input be
   * available. Thus, it needs to provide an implementation of
   * GenerateInputRequestedRegion(). */
  void GenerateInputRequestedRegion();
  /** LabelImageToLabelMapFilter will produce the entire output. */
  void EnlargeOutputRequestedRegion( DataObject *itkNotUsed(output) );

  virtual void GenerateData();

private:
  LabelImageToLabelMapGraphFilter(const Self &);   //purposely not implemented
  void operator=(const Self &);                         //purposely not implemented

  OutputImagePixelType mBackgroundValue;

  OriginalImagePointer mOriginalImage;
  InputImagePointer mLabelImage;
  vtkSmartPointer<vtkMutableUndirectedGraph>  mObjectGraph;

  bool usePrecalculatedGraph;
}; // end of class

} // end namespace itk

#include "LabelImageToLabelMapGraphFilter.tpp"

#endif /* LabelImageToLabelMapGraphFilter_H_ */
