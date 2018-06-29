///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  FeatureLabelObject.tpp                                               //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-09-11                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "FeatureLabelObject.h"

#include "itkLabelObjectLine.h"


namespace itk {

template< class TLabel, unsigned int VImageDimension > FeatureLabelObject< TLabel, VImageDimension >::FeatureLabelObject()
{
    mNumberOfPixels = 0;
    mPhysicalSize = 0;
    mCoordAcc.Fill(0);
    mCentroid.Fill(0);
    mBoundingBox.Fill(0);

    mRedAcc = 0;
    mRedMean = 0;
    mGreenAcc = 0;
    mGreenMean = 0;
    mBlueAcc = 0;
    mBlueMean = 0;

    mRedHistogram = HistogramType::New();
    mGreenHistogram = HistogramType::New();
    mBlueHistogram = HistogramType::New();
    mRGBHistogram = HistogramType::New();
    mTextureFeatures = FeatureValueVectorType::New();
}


template< class TLabel, unsigned int VImageDimension > void FeatureLabelObject< TLabel, VImageDimension >::CopyAttributesFrom(const LabelObjectType *lo)
{
    Superclass::CopyAttributesFrom(lo);

    // copy the data of the current type if possible
    const Self *src = dynamic_cast< const Self * >( lo );
    if ( src == NULL )
    {
        return;
    }
    mNumberOfPixels = src->mNumberOfPixels;
    mPhysicalSize = src->mPhysicalSize;
    mCoordAcc = src->mCoordAcc;
    mCentroid = src->mCentroid;
    mBoundingBox = src->mBoundingBox;

    mRedAcc = src->mRedAcc;
    mRedMean = src->mRedMean;
    mGreenAcc = src->mGreenAcc;
    mGreenMean = src->mGreenMean;
    mBlueAcc = src->mBlueAcc;
    mBlueMean = src->mBlueMean;

    mRedHistogram = src->mRedHistogram;
    mGreenHistogram = src->mGreenHistogram;
    mBlueHistogram = src->mBlueHistogram;
    mRGBHistogram = src->mRGBHistogram;
    mTextureFeatures = src->mTextureFeatures;
}


template< class TLabel, unsigned int VImageDimension > void FeatureLabelObject< TLabel, VImageDimension >::PrintSelf(std::ostream & os, Indent indent) const
{
   Superclass::PrintSelf(os, indent);

   os << indent << "Shape feaures: " << std::endl;
   os << indent << "NumberOfPixels: " << mNumberOfPixels << std::endl;
   os << indent << "PhysicalSize: " << mPhysicalSize << std::endl;
   os << indent << "CoordAcc: " << mCoordAcc << std::endl;
   os << indent << "Centroid: " << mCentroid << std::endl;
   os << indent << "BoundingBox: " << mBoundingBox << std::endl;
   os << indent << std::endl;

   os << indent << "Color feaures: " << std::endl;
   os << indent << "RedAcc: " << mRedAcc << std::endl;
   os << indent << "RedMean: " << mRedMean << std::endl;
   os << indent << "GreenAcc: " << mGreenAcc << std::endl;
   os << indent << "GreenMean: " << mGreenMean << std::endl;
   os << indent << "BlueAcc: " << mBlueAcc << std::endl;
   os << indent << "BlueMean: " << mBlueMean << std::endl;

   os << indent << "Red Histogram: " << mRedHistogram << std::endl;
   os << indent << "Green Histogram: " << mGreenHistogram << std::endl;
   os << indent << "Blue Histogram: " << mBlueHistogram << std::endl;
   os << indent << "RGB Histogram: " << mRGBHistogram << std::endl;

   os << indent << "Texture Features: " << mTextureFeatures << std::endl;
}

} // end namespace itk

