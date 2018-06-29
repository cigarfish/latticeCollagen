///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  IterativeCavityFillingImageFilter.tpp                                //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-11-13                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "IterativeCavityFillingImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"

#include <vector>
#include <algorithm>


namespace itk
{

template<class TInputImage> IterativeCavityFillingImageFilter< TInputImage >::IterativeCavityFillingImageFilter()
{
    m_radius = 1;
    m_MaximumNumberOfIterations = 8;
    m_CurrentNumberOfIterations = 0;
    m_NumberOfPixelsChanged = 0;
    m_spacing = new double[3];
    m_spacing[0] = 1;
    m_spacing[1] = 1;
    m_spacing[2] = 1;

    m_lower_threshold = NumericTraits< InputPixelType >::Zero;
    m_upper_threshold = NumericTraits< InputPixelType >::max();
}


template<class TInputImage> IterativeCavityFillingImageFilter< TInputImage >::~IterativeCavityFillingImageFilter()
{
    //nothing to do here
}


template<class TInputImage> void IterativeCavityFillingImageFilter<TInputImage>::GenerateData()
{
    typename InputImageType::ConstPointer input  = this->GetInput();

    m_NumberOfPixelsChanged = 0;

    typename CavityFilterType::Pointer filter = CavityFilterType::New();
    filter->SetRadius( m_radius );
    filter->SetLowerThreshold( m_lower_threshold );
    filter->SetUpperThreshold( m_upper_threshold );
    filter->SetSpacing(m_spacing);

    m_CurrentNumberOfIterations = 0;

    typename OutputImageType::Pointer output;

//    ProgressReporter progress(this, 0, m_MaximumNumberOfIterations);

    while(m_CurrentNumberOfIterations < m_MaximumNumberOfIterations)
    {
        filter->SetInput(input);
        filter->Update();

        m_CurrentNumberOfIterations++;
//        progress.CompletedPixel();   // not really a pixel but an iteration
//        this->InvokeEvent( IterationEvent() );

        const unsigned int numberOfPixelsChangedInThisIteration = filter->GetNumberOfPixelsChanged();
        m_NumberOfPixelsChanged += numberOfPixelsChangedInThisIteration;

        output = filter->GetOutput();
        output->DisconnectPipeline();
        input = output;
        if(numberOfPixelsChangedInThisIteration == 0)
        {
            break;
        }
    }
    std::cout << "finished cavity filter iterations" << std::endl;
    this->GraftOutput(output);
    std::cout << "grafted cavity output" << std::endl;
}

/**
 * Standard "PrintSelf" method
 */
template<class TInputImage> void IterativeCavityFillingImageFilter< TInputImage >::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    os << indent << "Radius: " << m_radius << std::endl;
    os << indent << "Lower Threshold: " << m_lower_threshold << std::endl;
    os << indent << "Upper Threshold: " << m_upper_threshold << std::endl;

    os << indent << "Maximum Number of Iterations : " << m_MaximumNumberOfIterations << std::endl;
    os << indent << "Current Number of Iterations : " << m_CurrentNumberOfIterations << std::endl;
    os << indent << "Number of Pixels Changed     : " << m_NumberOfPixelsChanged << std::endl;
}

} // end namespace itk

