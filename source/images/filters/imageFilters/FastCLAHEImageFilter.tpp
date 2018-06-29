///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  FastCLAHEImageFilter.tpp                                             //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2012-06-18                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "FastCLAHEImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkProgressReporter.h"

#include "vnl/vnl_math.h"

#include <map>
#include <set>


namespace itk
{

template<class TImageType> void FastCLAHEImageFilter<TImageType>::ClipHistogram(std::map<PixelType, unsigned long> &histo)
{
    typedef std::map<PixelType, unsigned long> HistoMapType;
    typename HistoMapType::iterator itMap;

    unsigned long nrExcess, oldNrExcess, binIncr, upper, step_size;
    long binExcess;

    nrExcess = 0;

    itMap = histo.begin();
    while(itMap!=histo.end())
    {
        binExcess = itMap->second - m_scaledClipLimit;
        if(binExcess > 0)
            nrExcess += binExcess;
        itMap++;
    }
//    std::cout << "clip limit = " << m_scaledClipLimit << std::endl;
//    std::cout << "1. excess = " << nrExcess << std::endl;

    binIncr = nrExcess / m_nrGreyLevels;
    upper = m_scaledClipLimit - binIncr;

    for(long i=min; i<(max+1); i++)
    {
        itMap = histo.find(i);

        if(itMap != histo.end())
        {
            if(histo[i] > m_scaledClipLimit)
                histo[i] = m_scaledClipLimit;
            else {
                if(histo[i] > upper) {
                    nrExcess -= (m_scaledClipLimit - histo[i]);
                    histo[i] = m_scaledClipLimit;
                }
                else {
                    nrExcess -= binIncr;
                    histo[i] += binIncr;
                }
            }
        }
        else {
            nrExcess -= binIncr;
            histo.insert(typename HistoMapType::value_type(i, binIncr));
        }
    }
//    std::cout << "2. excess = " << nrExcess << std::endl;

    oldNrExcess = nrExcess+1;
    while(nrExcess>0 && nrExcess<oldNrExcess)
    {
        oldNrExcess = nrExcess;

        step_size = m_nrGreyLevels / nrExcess;
        if(step_size<1)
            step_size = 1;

        for(long i=min; (nrExcess>0)&&(i<max+1); i+=step_size)
        {
            itMap = histo.find(i);

            if(itMap != histo.end()) {
                if(histo[i] < m_scaledClipLimit) {
                    histo[i]++;
                    nrExcess--;
                }
            }
            else {
                histo.insert(typename HistoMapType::value_type(i, 1));
                nrExcess--;
            }
        }
    }
//    std::cout << "3. excess = " << nrExcess << std::endl;
}


template<class TInputImage> void FastCLAHEImageFilter<TInputImage>::ThreadedGenerateData(const ImageRegionType & outputRegionForThread, ThreadIdType threadId)
{
    typename ImageType::ConstPointer input = this->GetInput();
    typename ImageType::Pointer output = this->GetOutput();

    // Allocate the output
    this->AllocateOutputs();

    int win_size = 2*m_stepSizeRadius+1;

    itk::Size<ImageDimension> imageSize;
    itk::Size<ImageDimension> winRadius;
    for(unsigned int i=0; i<ImageDimension; i++) {
        imageSize[i] = output->GetLargestPossibleRegion().GetSize()[i]; //outputRegionForThread.GetSize(i);
        winRadius[i] = m_stepSizeRadius;
    }

    //number of pixels/voxel in window
    m_histoVolume = 1;
    for(unsigned int i=0; i<ImageDimension; i++)
        m_histoVolume = m_histoVolume * (2 * this->GetRadius()[i] + 1);

    m_scaledClipLimit = m_clipLimit * m_histoVolume;

    bool calcNewHisto = true;

    bool isOnGrid[ImageDimension], isInLastRow[ImageDimension], isOneMissing[ImageDimension];

    for(unsigned int i=0; i<ImageDimension; i++) {
        int rest = (imageSize[i] - (m_stepSizeRadius + 1)) % win_size;

        if( (rest>0) && (rest<(m_stepSizeRadius+1)) )
            isOneMissing[i] = true;
        else
            isOneMissing[i] = false;
    }

    ImageRegionConstIterator<ImageType> itInput(input, input->GetRequestedRegion());

    // Setup for processing the image
    ZeroFluxNeumannBoundaryCondition<ImageType> nbc;
    ConstNeighborhoodIterator<ImageType>        hit;
    ConstNeighborhoodIterator<ImageType>        bit;

    // Find the data-set boundary "faces"
    NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<ImageType> bC;

    typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<ImageType>::FaceListType faceList;
    faceList = bC(input, outputRegionForThread, this->GetRadius());

    typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<ImageType>::FaceListType::iterator fit;

    // Map stores (number of pixel)/(window size) for each gray value.
    typedef std::map<PixelType, unsigned long> HistoMapType, CDFMapType;
    HistoMapType histo, cdf;
    double cdfMin, quot, s;
    typename HistoMapType::iterator itMap;


    // Process each faces.  These are N-d regions which border
    // the edge of the buffer
//    int fC = 0;
    for(fit=faceList.begin(); fit!=faceList.end(); ++fit)
    {
//        std::cout << "face = " << fC << std::endl;
//        fC++;
//        std::cout << "start index: " << fit->GetIndex() << " size: " << fit->GetSize() << std::endl;

        // Create a neighborhood iterator for the window
        bit = ConstNeighborhoodIterator<ImageType>(winRadius, input, *fit);
        bit.OverrideBoundaryCondition(&nbc);
        bit.GoToBegin();

        // Create a neighborhood iterator for the histogram
        hit = ConstNeighborhoodIterator<ImageType>(this->GetRadius(), input, *fit);
        hit.OverrideBoundaryCondition(&nbc);
        hit.GoToBegin();


        // iterate over the region for this face
        while(!bit.IsAtEnd())
        {
            //determine if it is necessary to compute a new histogram
            calcNewHisto = true;
            for(unsigned int i=0; i<ImageDimension; i++) {
                if(bit.GetIndex()[i] % win_size != 0)
                    isOnGrid[i] = false;
                else
                    isOnGrid[i] = true;

                if(bit.GetIndex()[i] == (int)(imageSize[i]-1))
                    isInLastRow[i] = true;
                else
                    isInLastRow[i] = false;


                if(!((isOneMissing[i] && isInLastRow[i]) || isOnGrid[i])) {
                    calcNewHisto = false;
                    break;
                }
            }

            //if we need a new histogram
            if(calcNewHisto)
            {
//                std::cout << "Calculate histogram for (" << hit.GetIndex()[0] << ", " << hit.GetIndex()[1] << ", " << hit.GetIndex()[2] << ") " << std::endl;

                //CLAHE algorithm

                PixelType f;

                //Compute the histogram
                histo.clear();
                for(unsigned int i=0; i<hit.Size(); ++i)
                {
                    f = hit.GetPixel(i);
                    itMap = histo.find(f);
                    if(itMap != histo.end())
                        histo[f] = histo[f] + 1;
                    else
                        histo.insert(typename HistoMapType::value_type(f, 1));
                }

                //clip the histogram
                ClipHistogram(histo);

//                itMap = histo.begin();
//                while(itMap!=histo.end()) {
//                    std::cout << "histo[" << (int)(itMap->first) << "] = " << itMap->second << std::endl;
//                    itMap++;
//                }

                //compute the cumulative distribution function
                cdf.clear();
                itMap = histo.begin();
                cdf[itMap->first] = itMap->second;
                if(itMap->first == min)
                    cdfMin = itMap->second;
                else
                    cdfMin = 0;
                f = itMap->first;
                itMap++;

                while(itMap!=histo.end())
                {
                    cdf[itMap->first] = cdf[f] + itMap->second;
                    f = itMap->first;
                    itMap++;
                }

//                itMap = cdf.begin();
//                while(itMap!=cdf.end()) {
//                    std::cout << "cdf[" << (int)(itMap->first) << "] = " << itMap->second << std::endl;
//                    itMap++;
//                }
//                itMap = cdf.end();
//                itMap--;
//                std::cout << "cdf[" << (int)(itMap->first) << "] = " << itMap->second << std::endl;

                quot = (double)(max)/(double)(m_histoVolume-cdfMin);
                s = (double)(cdfMin)*quot;

                //compute the equalized intensity values for all pixels within the window, based on this cdf
                //IMPORTANT NOTE: Use iterators here for output image, using direct pixel access via output.Set method results in crash during filter destructor call
                typename ImageType::RegionType outRegion = output->GetLargestPossibleRegion();
                typename ImageType::RegionType bitRegion = bit.GetBoundingBoxAsImageRegion();
                typename ImageType::RegionType iscRegion = outRegion;
                iscRegion.Crop(bitRegion);

                ImageRegionIterator<ImageType> itOut(output, iscRegion);
                for(unsigned int i=0; i<bit.Size(); ++i) {
                    bool isInBounds;
                    double px = bit.GetPixel(i, isInBounds);

                    if(isInBounds) {
                        itOut.Set( (PixelType)(cdf[px]*quot - s) );
                        ++itOut;
                    }
                }
            }

            // move the neighborhood
            ++bit;
            ++hit;
        }
    }
}


template<class TImageType> void FastCLAHEImageFilter<TImageType>::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

} // end namespace itk

