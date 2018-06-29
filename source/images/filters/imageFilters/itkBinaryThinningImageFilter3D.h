///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  itkBinaryThinningImageFilter3D.h                                     //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2011-10-28                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//     code by Hanno Homann, Oxford University, Wolfson Medical Vision Lab, UK,      //
//     published in the Insight Journal, and can be downloaded here:                 //
//     http://hdl.handle.net/1926/1292                                               //
//                                                                                   //
//     License information:                                                          //
//     The work made available by the Insight Journal is distributed under the       //
//     Creative Commons Attribution License Version 3.0. Legal details are           //
//     available at http://creativecommons.org/licenses/by/3.0/legalcode             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef __itkBinaryThinningImageFilter3D_h
#define __itkBinaryThinningImageFilter3D_h

#include <itkImageToImageFilter.h>


namespace itk
{

/*!
 \brief This filter computes a one-pixel-wide skeleton of a 3D input image.

 \n This class is parametrized over the type of the input image and the type of the output image.

 \n The input is assumed to be a binary image. All non-zero valued voxels are set to 1 internally to simplify the computation. The filter will produce a skeleton of the object.
 The output background values are 0, and the foreground values are 1.

 \n A 26-neighbourhood configuration is used for the foreground and a 6-neighbourhood configuration for the background. Thinning is performed symmetrically in order to
 guarantee that the skeleton lies medial within the object.


 \n
 \author Hanno Homann, Oxford University, Wolfson Medical Vision Lab, UK.
 \n Code was published in the Insight Journal, and can be downloaded here: http://hdl.handle.net/1926/1292.
 \n License information:
 \n The work made available by the Insight Journal is distributed under the Creative Commons Attribution License Version 3.0. Legal details are available at http://creativecommons.org/licenses/by/3.0/legalcode.

 \n This filter is a parallel thinning algorithm and is an implementation of the algorithm described in:
 \n T.C. Lee, R.L. Kashyap, and C.N. Chu.
 \n Building skeleton models via 3-D medial surface/axis thinning algorithms.
 \n Computer Vision, Graphics, and Image Processing, 56(6):462--478, 1994.
*/

template <class TInputImage,class TOutputImage> class BinaryThinningImageFilter3D : public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef BinaryThinningImageFilter3D    Self;
    typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro( BinaryThinningImageFilter3D, ImageToImageFilter );

    /** Type for input image. */
    typedef   TInputImage       InputImageType;

    /** Type for output image: Skelenton of the object.  */
    typedef   TOutputImage      OutputImageType;

    /** Type for the region of the input image. */
    typedef typename InputImageType::RegionType RegionType;

    /** Type for the index of the input image. */
    typedef typename RegionType::IndexType  IndexType;

    /** Type for the pixel type of the input image. */
    typedef typename InputImageType::PixelType InputImagePixelType ;

    /** Type for the pixel type of the input image. */
    typedef typename OutputImageType::PixelType OutputImagePixelType ;

    /** Type for the size of the input image. */
    typedef typename RegionType::SizeType SizeType;

    /** Pointer Type for input image. */
    typedef typename InputImageType::ConstPointer InputImagePointer;

    /** Pointer Type for the output image. */
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Boundary condition type for the neighborhood iterator */
    typedef ConstantBoundaryCondition< TInputImage > ConstBoundaryConditionType;

    /** Neighborhood iterator type */
    typedef NeighborhoodIterator<TInputImage, ConstBoundaryConditionType> NeighborhoodIteratorType;

    /** Neighborhood type */
    typedef typename NeighborhoodIteratorType::NeighborhoodType NeighborhoodType;

    /** Get Skelenton by thinning image. */
    OutputImageType * GetThinning(void);

    /** ImageDimension enumeration   */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension );
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TOutputImage::ImageDimension );

#ifdef ITK_USE_CONCEPT_CHECKING
    /** Begin concept checking */
    itkConceptMacro(SameDimensionCheck,
                    (Concept::SameDimension<InputImageDimension, 3>));
    itkConceptMacro(SameTypeCheck,
                    (Concept::SameType<InputImagePixelType, OutputImagePixelType>));
    itkConceptMacro(InputAdditiveOperatorsCheck,
                    (Concept::AdditiveOperators<InputImagePixelType>));
    itkConceptMacro(InputConvertibleToIntCheck,
                    (Concept::Convertible<InputImagePixelType, int>));
    itkConceptMacro(IntConvertibleToInputCheck,
                    (Concept::Convertible<int, InputImagePixelType>));
    itkConceptMacro(InputIntComparableCheck,
                    (Concept::Comparable<InputImagePixelType, int>));
    /** End concept checking */
#endif

protected:
    BinaryThinningImageFilter3D();
    virtual ~BinaryThinningImageFilter3D() {};
    void PrintSelf(std::ostream& os, Indent indent) const;

    /** Compute thinning Image. */
    void GenerateData();

    /** Prepare data. */
    void PrepareData();

    /**  Compute thinning Image. */
    void ComputeThinImage();

    /**  isEulerInvariant [Lee94] */
    bool isEulerInvariant(NeighborhoodType neighbors, int *LUT);
    void fillEulerLUT(int *LUT);
    /**  isSimplePoint [Lee94] */
    bool isSimplePoint(NeighborhoodType neighbors);
    /**  Octree_labeling [Lee94] */
    void Octree_labeling(int octant, int label, int *cube);


private:
    BinaryThinningImageFilter3D(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

}; // end of BinaryThinningImageFilter3D class

} //end namespace itk

#include "itkBinaryThinningImageFilter3D.tpp"

#endif
