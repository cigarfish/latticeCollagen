///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LineConstIteratorDerivative.h                                        //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-13                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LINECONSTITERATORDERIVATIVE_H_
#define LINECONSTITERATORDERIVATIVE_H_


namespace itk
{
/**
 * \class LineConstIterator
 * \brief An iterator that walks a Bresenham line through an ND image
 *        with read-only access to pixels.
 *
 * Same as in itk, but without warning after it leaves image region
 *
 */
template< class TImage >
class ITK_EXPORT LineConstIteratorDerivative
{
public:
    /** Standard class typedefs. */
    typedef LineConstIteratorDerivative Self;

    /** Dimension of the image that the iterator walks.  This constant is needed so
     * that functions that are templated over image iterator type (as opposed to
     * being templated over pixel type and dimension) can have compile time
     * access to the dimension of the image that the iterator walks. */
    itkStaticConstMacro(ImageIteratorDimension, unsigned int, TImage::ImageDimension);

    /** Index typedef support. */
    typedef typename TImage::IndexType      IndexType;

    /** Offset typedef support. */
    typedef typename TImage::OffsetType      OffsetType;

    /** Size typedef support. */
    typedef typename TImage::SizeType      SizeType;

    /** Region typedef support */
    typedef typename TImage::RegionType RegionType;

    /** Spacing typedef support */
    typedef typename TImage::SpacingType SpacingType;

    /** Origin typedef support */
    typedef typename TImage::PointType PointType;

    /** Image typedef support. */
    typedef TImage ImageType;

    /** PixelContainer typedef support. Used to refer to the container for
     * the pixel data. While this was already typdef'ed in the superclass,
     * it needs to be redone here for this subclass to compile properly with gcc. */
    typedef typename TImage::PixelContainer  PixelContainer;
    typedef typename PixelContainer::Pointer PixelContainerPointer;

    /** Internal Pixel Type */
    typedef typename TImage::InternalPixelType InternalPixelType;

    /** External Pixel Type */
    typedef typename TImage::PixelType PixelType;

    /**  Accessor type that convert data between internal and external
     *  representations. */
    typedef typename TImage::AccessorType AccessorType;

    /** Run-time type information (and related methods). */
    itkTypeMacro(LineConstIteratorDerivative, None);

    /** Get the dimension (size) of the index. */
    static unsigned int GetImageIteratorDimension()
    {
        return TImage::ImageDimension;
    }

    /** Get the index. This provides a read only reference to the index. */
    const IndexType GetIndex()
    {
        return m_CurrentImageIndex;
    }

    /** Get the pixel value */
    const PixelType & Get(void) const
    {
        return m_Image->GetPixel(m_CurrentImageIndex);
    }

    /** Is the iterator at the end of the line? */
    bool IsAtEnd()
    {
        return m_IsAtEnd;
    }

    /** Move an iterator to the beginning of the line. */
    void GoToBegin();

    /** Walk forward along the line to the next index in the image. */
    void operator++();

    /** operator= is provided to make sure the handle to the image is properly
     * reference counted. */
    Self & operator=(const Self & it);

    /** Constructor establishes an iterator to walk along a line */
    LineConstIteratorDerivative(const ImageType *imagePtr, const IndexType & firstIndex, const IndexType & lastIndex);

    /** Default Destructor. */
    virtual ~LineConstIteratorDerivative() {}
protected: //made protected so other iterators can access
    /** Smart pointer to the source image. */
    typename ImageType::ConstWeakPointer m_Image;

    /** Region type to iterate over. */
    RegionType m_Region;

    /** Is the iterator at the end of its walk? */
    bool m_IsAtEnd;

    /** Start, end and current ND index position in the image of the line */
    IndexType m_CurrentImageIndex;
    IndexType m_StartIndex;
    IndexType m_LastIndex;
    IndexType m_EndIndex;  // one past the end of the line in the m_MainDirection

    /** Variables that drive the Bresenham-Algorithm */
    // The dimension with the largest difference between start and end
    unsigned int m_MainDirection;

    // Accumulated error for the other dimensions
    IndexType m_AccumulateError;

    // Increment for the error for each step. Two times the difference between
    // start and end
    IndexType m_IncrementError;

    // If enough is accumulated for a dimension, the index has to be
    // incremented. Will be the number of pixels in the line
    IndexType m_MaximalError;

    // Direction of increment. -1 or 1
    IndexType m_OverflowIncrement;

    // After an overflow, the accumulated error is reduced again. Will be
    // two times the number of pixels in the line
    IndexType m_ReduceErrorAfterIncrement;
};

} // end namespace itk

#include "LineConstIteratorDerivative.tpp"

#endif /* LINECONSTITERATORDERIVATIVE_H_ */
