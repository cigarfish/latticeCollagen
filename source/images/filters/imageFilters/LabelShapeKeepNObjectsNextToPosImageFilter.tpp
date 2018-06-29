///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelShapeKeepNObjectsNextToPosImageFilter.tpp                       //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-06-28                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "LabelShapeKeepNObjectsNextToPosImageFilter.h"

#include "itkLabelMapUtilities.h"
#include "itkShapeLabelObject.h"
#include "itkShapeLabelObjectAccessors.h"
#include "itkShapeKeepNObjectsLabelMapFilter.h"


namespace itk
{

template<class TImage> LabelShapeKeepNObjectsNextToPosImageFilter<TImage>::LabelShapeKeepNObjectsNextToPosImageFilter()
{
    m_BackgroundValue = NumericTraits<PixelType>::NonpositiveMin();
    m_NumberOfObjects = 1;
    m_ReverseOrdering = false;

    for(unsigned int i=0; i<ImageDimension; i++)
        m_Position[i] = 95.0;                                                                                               //TODO: largestpossibleregion based...
    m_Position[2] = 52;

    // create the output image for the removed objects
    this->SetNumberOfRequiredOutputs(2);
    this->SetNthOutput(1, static_cast<TImage*>(this->MakeOutput(1).GetPointer()) );
}


template<class TImage> void LabelShapeKeepNObjectsNextToPosImageFilter<TImage>::GenerateData()
{
    // Allocate the output
    this->AllocateOutputs();

    ImageType *output = this->GetOutput();
    ImageType *output2 = this->GetOutput(1);

    output2->SetBackgroundValue(output->GetBackgroundValue());

    typedef std::multimap<double, LabelObjectType*> DistMapType;
    DistMapType distMap;

    for(unsigned int i=0; i<output->GetNumberOfLabelObjects(); i++) {
        //compute dist to centroid
        double dist = m_Position.EuclideanDistanceTo(output->GetNthLabelObject(i)->GetCentroid());

        //push to map
        distMap.insert(std::pair<double, LabelObjectType*>(dist, output->GetNthLabelObject(i)));
    }

    typename DistMapType::iterator it;

    int counter = m_NumberOfObjects;

    if(m_ReverseOrdering) {
        it = distMap.end();

        while(counter>0 && it!=distMap.begin()) {
            output2->AddLabelObject(it->second);
            output->RemoveLabelObject(it->second);
            it--;
            counter--;
        }
    }
    else {
        it = distMap.begin();

        while(counter>0 && it!=distMap.end()) {
            output2->AddLabelObject(it->second);
            output->RemoveLabelObject(it->second);
            it++;
            counter--;
        }
    }
}


template<class TImage> void LabelShapeKeepNObjectsNextToPosImageFilter<TImage>::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);

    os << indent << "BackgroundValue: " << static_cast< typename NumericTraits<PixelType>::PrintType >(m_BackgroundValue) << std::endl;
    os << indent << "NumberOfObjects: "  << m_NumberOfObjects << std::endl;
    os << indent << "ReverseOrdering: "  << m_ReverseOrdering << std::endl;
    os << indent << "Pos: (";
    for(unsigned int i=0; i<ImageDimension-1; i++)
        os << indent << m_Position[i] << ",";
    os << indent << m_Position[ImageDimension-1] << ")" << std::endl;
}

} // end namespace itk

