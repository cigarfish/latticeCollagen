///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SLICImageFilter.h                                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-07-03                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SLICIMAGEFILTER_H_
#define SLICIMAGEFILTER_H_

#include "itkImageToImageFilter.h"


namespace itk
{

template<class TInputImage, class TLabelImage> class ITK_EXPORT SLICImageFilter : public ImageToImageFilter<TInputImage, TLabelImage>
{
public:
    /**
     * Standard class typedefs
     */
    typedef SLICImageFilter                                 Self;
    typedef ImageToImageFilter< TInputImage, TLabelImage >  Superclass;
    typedef SmartPointer< Self >                            Pointer;
    typedef SmartPointer< const Self >                      ConstPointer;

    typedef TInputImage                                         InputImageType;
    typedef typename InputImageType::Pointer                    InputImagePointer;
    typedef typename InputImageType::ConstPointer               InputConstImagePointer;
    typedef typename InputImageType::PixelType                  InputPixelType;
    typedef typename InputImageType::RegionType                 InputImageRegionType;
    typedef typename InputImageType::IndexType                  InputIndexType;
    typedef typename InputImageType::SizeType                   InputSizeType;
    typedef typename InputImageType::SpacingType                InputSpacingType;
    typedef ImageRegionConstIterator<InputImageType>            InputIteratorType;

    typedef TLabelImage                                         LabelImageType;
    typedef typename LabelImageType::Pointer                    LabelImagePointer;
    typedef typename LabelImageType::PixelType                  LabelPixelType;
    typedef typename LabelImageType::RegionType                 LabelImageRegionType;
    typedef typename LabelImageType::IndexType                  LabelIndexType;
    typedef typename LabelImageType::SizeType                   LabelSizeType;

    itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int, TLabelImage::ImageDimension);

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SLICImageFilter, ImageToImageFilter);

    void SetSuperpixelSpacing(InputSpacingType voxelSpacing, InputSpacingType superpixelSpacing) { mSuperpixelBySize = true; mVoxelSpacing = voxelSpacing; mSuperpixelSpacing = superpixelSpacing; }
    void SetNumberSuperPixel(unsigned long numSuperpixel) { mSuperpixelBySize = false; mNumberSuperPixels = numSuperpixel; }
    void SetPerturbSeeds(bool perturbSeeds) { mPerturbSeeds = perturbSeeds; }

    unsigned long GetNumberSuperPixel() { return mNumberSuperPixels; }

protected:
    typedef double                                                          PixelDoubleType;
    typedef long int                                                        HelperPixelLongType;
    typedef itk::RGBPixel<PixelDoubleType>                                  LABPixelType;

    typedef itk::Image<PixelDoubleType, InputImageDimension>                DoubleScalarImageType;
    typedef itk::Image<LABPixelType, InputImageDimension>                   LABScalarImageType;
    typedef itk::Image<HelperPixelLongType, InputImageDimension>            HelperScalarVoImageType;

    typedef ImageRegionIterator<LABScalarImageType>                         LABIteratorType;

    typedef typename DoubleScalarImageType::IndexType                       DoubleScalarIndexType;
    typedef typename LABScalarImageType::IndexType                          LABIndexType;
    typedef typename HelperScalarVoImageType::IndexType                     HelperIndexType;

    SLICImageFilter()
    {
        mPerturbSeeds = true;
        mNumberSuperPixels = 500;

        mImageVolume = 0;
    }

    virtual ~SLICImageFilter(){}
    void PrintSelf(std::ostream & os, Indent indent) const;

    void RGB2XYZ(InputPixelType& sRGB, LABPixelType& sXYZ);
    void RGB2LAB(InputPixelType& sRGB, LABPixelType& sLAB);
    void DoRGBtoLABConversion();

    void DetectLABEdges();
    void PerturbSeeds(std::vector<LABPixelType>& kseedsPxl, std::vector<LABIndexType>& kseedsIdx);
    void EnforceLabelConnectivity(int& numlabels);

    void PerformSuperpixelSegmentation_VariableSandM(std::vector<LABPixelType>& kseedsPxl, std::vector<LABIndexType>& kseedsIdx, const int& STEP, const int& NUMITR);

    void GetLABXYZSeeds_ForGivenK(std::vector<LABPixelType>& kseedsPxl, std::vector<LABIndexType>& kseedsIdx);

    void PerformSLICO_ForGivenK(const double& m);

    void GenerateData();

private:
    SLICImageFilter(const Self &);       //purposely not implemented
    void operator=(const Self &);        //purposely not implemented

    typename LABScalarImageType::Pointer mLABInputImage;
    typename DoubleScalarImageType::Pointer mEdgeImage;
    typename LabelImageType::Pointer mLabelImage;
    typename HelperScalarVoImageType::Pointer mEnforceConnectivityImage;

    bool mPerturbSeeds;
    bool mSuperpixelBySize;
    InputSpacingType mVoxelSpacing;
    InputSpacingType mSuperpixelSpacing;
    unsigned long mNumberSuperPixels;

    InputSizeType mInputSize;
    unsigned int mImageVolume;
};

} // end namespace itk

#include "SLICImageFilter.tpp"

#endif /* SLICIMAGEFILTER_H_ */
