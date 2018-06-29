///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  FastCLAHEImageFilter.h                                               //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-06-18                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef FASTCLAHE_H_
#define FASTCLAHE_H_

#include "itkBoxImageFilter.h"


namespace itk
{

template<class TImageType> class ITK_EXPORT FastCLAHEImageFilter : public BoxImageFilter< TImageType, TImageType >
{
public:
    /**
     * Standard class typedefs
     */
    typedef FastCLAHEImageFilter                         Self;
    typedef ImageToImageFilter< TImageType, TImageType > Superclass;
    typedef SmartPointer< Self >                         Pointer;
    typedef SmartPointer< const Self >                   ConstPointer;

    itkStaticConstMacro(ImageDimension, unsigned int, TImageType::ImageDimension);

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(FastCLAHEImageFilter, ImageToImageFilter);

    /** Image type typedef support. */
    typedef TImageType                      ImageType;

    typedef typename ImageType::RegionType  ImageRegionType;

    typedef typename ImageType::SizeType    ImageSizeType;
    typedef typename ImageType::PixelType   PixelType;

    void SetClipLimit(double cLi)
    {
        m_clipLimit = cLi;
    }

    double GetClipLimit()
    {
        return m_clipLimit;
    }

    void SetStepSizeRadius(unsigned long ssr)
    {
        m_stepSizeRadius = ssr;
    }

    unsigned long GetStepSizeRadius()
    {
        return m_stepSizeRadius;
    }

protected:
    FastCLAHEImageFilter()
    {
        min = itk::NumericTraits<PixelType>::min();
        max = itk::NumericTraits<PixelType>::max();

        m_stepSizeRadius = 1;
        m_clipLimit = 0.1;
        m_nrGreyLevels = 256;
        this->SetRadius(5);
    }

    virtual ~FastCLAHEImageFilter(){}
    void PrintSelf(std::ostream & os, Indent indent) const;

    /**
     * Standard pipeline method
     */
//    void GenerateData();

    void ThreadedGenerateData(const ImageRegionType & outputRegionForThread, ThreadIdType threadId);

private:
    FastCLAHEImageFilter(const Self &);       //purposely not
    // implemented
    void operator=(const Self &);           //purposely not

    void ClipHistogram(std::map<PixelType,unsigned  long> &histo);

    PixelType min;
    PixelType max;

    unsigned long m_stepSizeRadius;
    unsigned long m_nrGreyLevels;
    unsigned long m_histoVolume;
    double m_clipLimit;
    unsigned int m_scaledClipLimit;
};

} // end namespace itk

#include "FastCLAHEImageFilter.tpp"

#endif /* FASTAHE_H_ */
