///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelShapeKeepNObjectsNextToPosImageFilter.h                         //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-06-28                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LABELSHAPEKEEPNOBJECTSNEXTTOPOSIMAGEFILTER_H_
#define LABELSHAPEKEEPNOBJECTSNEXTTOPOSIMAGEFILTER_H_

#include "itkInPlaceLabelMapFilter.h"


namespace itk
{

template<class TImage> class ITK_EXPORT LabelShapeKeepNObjectsNextToPosImageFilter : public InPlaceLabelMapFilter<TImage>
{
public:
    /** Standard class typedefs. */
    typedef LabelShapeKeepNObjectsNextToPosImageFilter  Self;
    typedef InPlaceLabelMapFilter<TImage>               Superclass;
    typedef SmartPointer<Self>                          Pointer;
    typedef SmartPointer<const Self>                    ConstPointer;

    /** Some convenient typedefs. */
    typedef TImage                              ImageType;
    typedef typename ImageType::Pointer         ImagePointer;
    typedef typename ImageType::ConstPointer    ImageConstPointer;
    typedef typename ImageType::PixelType       PixelType;
    typedef typename ImageType::IndexType       IndexType;
    typedef typename ImageType::LabelObjectType LabelObjectType;

    typedef typename LabelObjectType::AttributeType AttributeType;

    /** ImageDimension constants */
    itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);

    /** Standard New method. */
    itkNewMacro(Self);

    /** Runtime information support. */
    itkTypeMacro(LabelShapeKeepNObjectsNextToPosImageFilter, InPlaceLabelMapFilter);

    //Some more convenient typedefs
    typedef Point<double, ImageDimension>       PointType;

    /**
     * Set/Get the value used as "background" in the output image.
     * Defaults to NumericTraits<PixelType>::NonpositiveMin().
     */
    itkSetMacro(BackgroundValue, PixelType);
    itkGetConstMacro(BackgroundValue, PixelType);

    /**
     * Set/Get the number of objects to keep
     */
    itkGetConstMacro(NumberOfObjects, SizeValueType);
    itkSetMacro(NumberOfObjects, SizeValueType);

    /**
     * Set/Get the ordering of the objects. By default, the ones with the
     * highest value are kept. Turming ReverseOrdering to true make this filter
     * keep the objects with the smallest values
     */
    itkGetConstMacro(ReverseOrdering, bool);
    itkSetMacro(ReverseOrdering, bool);
    itkBooleanMacro(ReverseOrdering);

    /**
     * Set/Get the position
     */
    itkGetConstMacro(Position, PointType);
    itkSetMacro(Position, PointType);

protected:
  LabelShapeKeepNObjectsNextToPosImageFilter();
  ~LabelShapeKeepNObjectsNextToPosImageFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  void GenerateData();

private:
  LabelShapeKeepNObjectsNextToPosImageFilter(const Self &);  //purposely not implemented
  void operator=(const Self &);                              //purposely not implemented

  PixelType             m_BackgroundValue;
  SizeValueType         m_NumberOfObjects;
  bool                  m_ReverseOrdering;
  PointType             m_Position;
}; // end of class

} // end namespace itk

#include "LabelShapeKeepNObjectsNextToPosImageFilter.tpp"

#endif /* LABELSHAPEKEEPNOBJECTSNEXTTOPOSIMAGEFILTER_H_ */
