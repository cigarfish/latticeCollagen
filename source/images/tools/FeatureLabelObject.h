///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  FeatureLabelObject.h                                                 //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-09-11                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef FEATURELABELOBJECT_H_
#define FEATURELABELOBJECT_H_

#include <itkLabelObject.h>

#include <itkHistogram.h>
#include <itkVectorContainer.h>


namespace itk {

template< class TLabel, unsigned int VImageDimension > class ITK_EXPORT FeatureLabelObject : public LabelObject< TLabel, VImageDimension >
{
public:
    /** Standard class typedefs */
    typedef FeatureLabelObject                     Self;
    typedef LabelObject< TLabel, VImageDimension > Superclass;
    typedef typename Superclass::LabelObjectType   LabelObjectType;
    typedef SmartPointer< Self >                   Pointer;
    typedef SmartPointer< const Self >             ConstPointer;
    typedef WeakPointer< const Self >              ConstWeakPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(FeatureLabelObject, LabelObject);

    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

    typedef TLabel LabelType;
    typedef typename Superclass::IndexType                                          IndexType;
    typedef typename Superclass::LineType                                           LineType;
    typedef typename Superclass::LengthType                                         LengthType;
    typedef typename Superclass::AttributeType                                      AttributeType;
    typedef Point< double, VImageDimension >                                        CentroidType;
    typedef itk::FixedArray<int, itkGetStaticConstMacro(ImageDimension)*2>          BoundingBoxType;

    typedef double                                                                                  MeasurementType;
    typedef itk::Statistics::Histogram<MeasurementType, itk::Statistics::DenseFrequencyContainer2>  HistogramType;

    typedef itk::VectorContainer<unsigned char, double>                                             FeatureValueVectorType;

    static AttributeType GetAttributeFromName(const std::string &s);
    static std::string GetNameFromAttribute(const AttributeType &a);

    //Shape feature Set-/Getter
    double GetPhysicalSize() const { return mPhysicalSize; }
    void SetPhysicalSize(const double v) { mPhysicalSize = v; }

    SizeValueType GetNumberOfPixels() const { return mNumberOfPixels; }
    void SetNumberOfPixels(const SizeValueType &v) { mNumberOfPixels = v; };

    CentroidType GetCoordAcc() const { return mCoordAcc; };
    void SetCoordAcc(const CentroidType coordAcc) { mCoordAcc = coordAcc; };
    CentroidType GetCentroid() const { return mCentroid; };
    void SetCentroid(const CentroidType centroid) { mCentroid = centroid; };
    BoundingBoxType GetBoundingBox() const { return mBoundingBox; };
    void SetBoundingBox(const BoundingBoxType boundingBox) { mBoundingBox = boundingBox; };

    //Color feature Set-/Getter
    double GetRedAcc() const { return mRedAcc; };
    void SetRedAcc(const double v) { mRedAcc = v; };
    double GetRedMean() const { return mRedMean; };
    void SetRedMean(const double v) { mRedMean = v; };

    double GetGreenAcc() const { return mGreenAcc; };
    void SetGreenAcc(const double v) { mGreenAcc = v; };
    double GetGreenMean() const { return mGreenMean; };
    void SetGreenMean(const double v) { mGreenMean = v; };

    double GetBlueAcc() const { return mBlueAcc; };
    void SetBlueAcc(const double v) { mBlueAcc = v; };
    double GetBlueMean() const { return mBlueMean; };
    void SetBlueMean(const double v) { mBlueMean = v; };

    HistogramType::Pointer GetRedHistogram() const { return mRedHistogram; };
    void SetRedHistogram(const HistogramType::Pointer h) { mRedHistogram = h; };

    HistogramType::Pointer GetGreenHistogram() const { return mGreenHistogram; };
    void SetGreenHistogram(const HistogramType::Pointer h) { mGreenHistogram = h; };

    HistogramType::Pointer GetBlueHistogram() const { return mBlueHistogram; };
    void SetBlueHistogram(const HistogramType::Pointer h) { mBlueHistogram = h; };

    HistogramType::Pointer GetRGBHistogram() const { return mRGBHistogram; };
    void SetRGBHistogram(const HistogramType::Pointer h) { mRGBHistogram = h; };

    //TODO: prelim
    FeatureValueVectorType::Pointer GetTextureFeatures() const { return mTextureFeatures; };
    void SetTextureFeatures(FeatureValueVectorType::Pointer h) { mTextureFeatures = h; };

    virtual void CopyAttributesFrom(const LabelObjectType *lo);

protected:
  FeatureLabelObject();
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  FeatureLabelObject(const Self &); //purposely not implemented
  void operator=(const Self &);     //purposely not implemented

  //Shape Features
  SizeValueType     mNumberOfPixels;
  double            mPhysicalSize;
  CentroidType      mCoordAcc;          //mCentroid = mCoordAcc / mNumberOfPixels
  CentroidType      mCentroid;
  BoundingBoxType   mBoundingBox;

  //Color Features
  double        mRedAcc;            //mRedMean = mRedAcc / mNumberOfPixels
  double        mRedMean;
  double        mGreenAcc;          //mGreenMean = mGreenAcc / mNumberOfPixels
  double        mGreenMean;
  double        mBlueAcc;           //mBlueMean = mBlueAcc / mNumberOfPixels
  double        mBlueMean;

  HistogramType::Pointer mRedHistogram;
  HistogramType::Pointer mGreenHistogram;
  HistogramType::Pointer mBlueHistogram;
  HistogramType::Pointer mRGBHistogram;
  //TODO: prelim
  FeatureValueVectorType::Pointer mTextureFeatures;
};

} // end namespace itk

#include "FeatureLabelObject.tpp"

#endif
