///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CavityFillingImageFilter.tpp                                         //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-12                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "CavityFillingImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkProgressReporter.h"

#include "../../../model/BasicDatatypes/Vector.h"
#include "../../tools/LineConstIteratorDerivative.h"

#include <algorithm>
#include <math.h>
#include <vector>


namespace itk
{

template< class TInputImage, class TOutputImage > CavityFillingImageFilter< TInputImage, TOutputImage >::CavityFillingImageFilter()
{
    m_random = new Random();
    m_random->Init();

    m_NumberOfPixelsChanged = 0;
    m_radius = 1;
    m_samplePoints = 20;
    m_lower_threshold = NumericTraits< InputPixelType >::Zero;
    m_upper_threshold = NumericTraits< InputPixelType >::max();
    m_spacing[0] = 1.;
    m_spacing[1] = 1.;
    m_spacing[2] = 1.;
}


template< class TInputImage, class TOutputImage > CavityFillingImageFilter< TInputImage, TOutputImage >::~CavityFillingImageFilter()
{
    delete m_random;
}


template <class TInputImage, class TOutputImage> void CavityFillingImageFilter<TInputImage, TOutputImage>::BeforeThreadedGenerateData()
{
    typename InputImageType::ConstPointer input    = this->GetInput();
    typename OutputImageType::Pointer output       = this->GetOutput();

    output->SetRegions(input->GetLargestPossibleRegion());
    output->Allocate();
    output->FillBuffer(NumericTraits< InputPixelType >::Zero);

    unsigned int numberOfThreads = this->GetNumberOfThreads();
    this->m_Count.SetSize(numberOfThreads);
    for(unsigned int i = 0; i<numberOfThreads; i++)
        this->m_Count[i] = 0;

    std::cout << "Cavity filter here: number of threads = " << numberOfThreads << std::endl;
}


template< class TInputImage, class TOutputImage > void CavityFillingImageFilter<TInputImage, TOutputImage>::ThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread, ThreadIdType threadId)
{
    typename InputImageType::ConstPointer input = this->GetInput();
    typename OutputImageType::Pointer output = this->GetOutput();

    ConstBoundaryConditionType boundaryCondition;
    boundaryCondition.SetConstant(NumericTraits< InputPixelType >::max());

    typename InputNeighborhoodIteratorType::RadiusType radius;
    radius.Fill(m_radius);
    InputNeighborhoodIteratorType it(radius, input, outputRegionForThread);
    OutputNeighborhoodIteratorType ot(radius, output, outputRegionForThread);
    it.SetBoundaryCondition(boundaryCondition);

    unsigned int numberOfPixelsChanged = 0;

    // Define points
    typedef typename InputNeighborhoodIteratorType::OffsetType OffsetType;
    OffsetType *directions = new OffsetType[m_samplePoints];

    double scaledRadius = m_radius*m_spacing[0];

    for(int i=0; i<m_samplePoints; i++) {
        double x,y,z;
        m_random->GetRandomUnitVector(&x, &y, &z);

        directions[i][0] = std::ceil(x*scaledRadius/m_spacing[0]);
        directions[i][1] = std::ceil(y*scaledRadius/m_spacing[1]);
        directions[i][2] = std::ceil(z*scaledRadius/m_spacing[2]);
//        std::cout << "point " << i << " on sphere = " << directions[i] << std::endl;
    }
	
    // Loop through the image.
    for(it.GoToBegin(), ot.GoToBegin(); !it.IsAtEnd(); ++it, ++ot)
    {
        typename InputNeighborhoodIteratorType::IndexType centerIndex = it.GetIndex();

//        if(centerIndex[0]==0 && centerIndex[1] == 0 && centerIndex[2]==0)
//            std::cout << "For neighborhood at pixel " << centerIndex << std::endl;

        if(it.GetCenterPixel()>0) {
            ot.SetCenterPixel(NumericTraits<OutputPixelType>::max());
        }
        else {
            int hit = 0;
            int numOutside = 0;

            //Loop through the different directions
            for(int i=0; i<m_samplePoints; i++)
            {
                typename InputNeighborhoodIteratorType::IndexType directionIndex = it.GetIndex(directions[i]);

//                if(centerIndex[0]==0 && centerIndex[1] == 0 && centerIndex[2]==0)
//                    std::cout << "Build line iterator to " << directionIndex << std::endl;

                itk::LineConstIteratorDerivative<InputImageType> lt(input, centerIndex, directionIndex);
                lt.GoToBegin();

                //Loop through the line
                while(!lt.IsAtEnd())
                {
                    if(lt.Get()>0) {
//                        if(centerIndex[0]==0 && centerIndex[1] == 0 && centerIndex[2]==0)
//                            std::cout << "For index " << centerIndex << " we have hit at " << lt.GetIndex() << std::endl;
                        hit++;
                        break;
                    }

                    ++lt;
                }
                if(!input->GetLargestPossibleRegion().IsInside(lt.GetIndex()))
                    numOutside++;
            }
            double hitQuot = (double)hit/(double)(m_samplePoints-numOutside) * (1.-(double)numOutside/(double)m_samplePoints);
//            double hitQuot = (double)hit/(double)(m_samplePoints-numOutside);

//            if(centerIndex[0]==0 && centerIndex[1] == 0 && centerIndex[2]==0)
//                std::cout << "For index " << centerIndex << " we have hit = " << hit << " & numOutside=" << numOutside << " & hitQuot = " << hitQuot << std::endl;

            if(m_lower_threshold <= hitQuot && hitQuot <= m_upper_threshold) {
                ot.SetCenterPixel(NumericTraits<OutputPixelType>::max());
                numberOfPixelsChanged++;
            }
        }
    } // end image iteration loop

    delete directions;

    this->m_Count[threadId] = numberOfPixelsChanged;
}


template<class TInputImage, class TOutputImage> void CavityFillingImageFilter<TInputImage, TOutputImage>::AfterThreadedGenerateData()
{
    this->m_NumberOfPixelsChanged = 0;

    for ( unsigned int t = 0; t < m_Count.GetSize(); t++ )
        this->m_NumberOfPixelsChanged += this->m_Count[t];
}


template< class TInputImage, class TOutput > void CavityFillingImageFilter< TInputImage, TOutput >::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    os << indent << "Radius: " << m_radius << std::endl;
    os << indent << "Lower Threshold: " << m_lower_threshold << std::endl;
    os << indent << "Upper Threshold: " << m_upper_threshold << std::endl;
}
} // end namespace itk

