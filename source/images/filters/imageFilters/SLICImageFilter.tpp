///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  SLICImageFilter.tpp                                                  //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-07-03                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "SLICImageFilter.h"

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include <algorithm>
#include <assert.h>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


namespace itk
{

template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::RGB2XYZ(InputPixelType& sRGB, LABPixelType& sXYZ)
{
    double R = sRGB[0]/255.0;
    double G = sRGB[1]/255.0;
    double B = sRGB[2]/255.0;

    double r, g, b;

    if(R <= 0.04045)    r = R/12.92;
    else                r = pow((R+0.055)/1.055,2.4);
    if(G <= 0.04045)    g = G/12.92;
    else                g = pow((G+0.055)/1.055,2.4);
    if(B <= 0.04045)    b = B/12.92;
    else                b = pow((B+0.055)/1.055,2.4);

    sXYZ[0] = r*0.4124564 + g*0.3575761 + b*0.1804375;
    sXYZ[1] = r*0.2126729 + g*0.7151522 + b*0.0721750;
    sXYZ[2] = r*0.0193339 + g*0.1191920 + b*0.9503041;
}


template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::RGB2LAB(InputPixelType& sRGB, LABPixelType& sLAB)
{
    //------------------------
    // sRGB to XYZ conversion
    //------------------------
    LABPixelType sXYZ;
    RGB2XYZ(sRGB, sXYZ);

    //------------------------
    // XYZ to LAB conversion
    //------------------------
    double epsilon = 0.008856;  //actual CIE standard
    double kappa   = 903.3;     //actual CIE standard

    double Xr = 0.950456;   //reference white
    double Yr = 1.0;        //reference white
    double Zr = 1.088754;   //reference white

    double xr = sXYZ[0]/Xr;
    double yr = sXYZ[1]/Yr;
    double zr = sXYZ[2]/Zr;

    double fx, fy, fz;
    if(xr > epsilon)    fx = pow(xr, 1.0/3.0);
    else                fx = (kappa*xr + 16.0)/116.0;
    if(yr > epsilon)    fy = pow(yr, 1.0/3.0);
    else                fy = (kappa*yr + 16.0)/116.0;
    if(zr > epsilon)    fz = pow(zr, 1.0/3.0);
    else                fz = (kappa*zr + 16.0)/116.0;

    sLAB[0] = 116.0*fy-16.0;
    sLAB[1] = 500.0*(fx-fy);
    sLAB[2] = 200.0*(fy-fz);
}


template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::DoRGBtoLABConversion()
{
    InputConstImagePointer input  = this->GetInput();

    InputIteratorType itIn(input, input->GetLargestPossibleRegion());
    LABIteratorType itOut(mLABInputImage, mLABInputImage->GetLargestPossibleRegion());
    itIn.GoToBegin();
    itOut.GoToBegin();

    while(!itIn.IsAtEnd() || !itOut.IsAtEnd())
    {
        InputPixelType iPx = itIn.Get();
        LABPixelType oPx;
        RGB2LAB(iPx, oPx);

        itOut.Set(oPx);

        ++itIn;
        ++itOut;
    }
}


template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::DetectLABEdges()
{
    mEdgeImage = DoubleScalarImageType::New();
    mEdgeImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    mEdgeImage->Allocate();

    itk::ImageRegionIterator<DoubleScalarImageType> itEdge(mEdgeImage, mEdgeImage->GetLargestPossibleRegion());
    itEdge.GoToBegin();

    typename LABScalarImageType::SizeType radius;
    radius.Fill(1);

    itk::ConstNeighborhoodIterator<LABScalarImageType> itLAB(radius, mLABInputImage, mLABInputImage->GetLargestPossibleRegion());
    itLAB.GoToBegin();

    typedef typename ConstNeighborhoodIterator<LABScalarImageType>::OffsetType OffsetType;
    OffsetType *off = new OffsetType[6];

    for(unsigned int i=0; i<InputImageDimension; i++) {
        off[2*i].Fill(0);
        off[2*i+1].Fill(0);

        off[2*i][i] = -1;
        off[2*i+1][i] = 1;
    }

    while(!itLAB.IsAtEnd() || !itEdge.IsAtEnd())
    {
        double dx(0), dy(0), dz(0);

        dx = pow(itLAB.GetPixel(off[0])[0]-itLAB.GetPixel(off[1])[0], 2) + pow(itLAB.GetPixel(off[0])[1]-itLAB.GetPixel(off[1])[1], 2) + pow(itLAB.GetPixel(off[0])[2]-itLAB.GetPixel(off[1])[2], 2);
        dy = pow(itLAB.GetPixel(off[2])[0]-itLAB.GetPixel(off[3])[0], 2) + pow(itLAB.GetPixel(off[2])[1]-itLAB.GetPixel(off[3])[1], 2) + pow(itLAB.GetPixel(off[2])[2]-itLAB.GetPixel(off[3])[2], 2);
        if(InputImageDimension == 3)
            dz = pow(itLAB.GetPixel(off[4])[0]-itLAB.GetPixel(off[5])[0], 2) + pow(itLAB.GetPixel(off[4])[1]-itLAB.GetPixel(off[5])[1], 2) + pow(itLAB.GetPixel(off[4])[2]-itLAB.GetPixel(off[5])[2], 2);

        itEdge.Set(dx + dy + dz);

        ++itLAB;
        ++itEdge;
    }
    if(off) delete [] off;
}


template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::PerturbSeeds(std::vector<LABPixelType>& kseedsPxl, std::vector<LABIndexType>& kseedsIdx)
{
    typename DoubleScalarImageType::SizeType radius;
    radius.Fill(1);

    itk::ConstNeighborhoodIterator<DoubleScalarImageType> it(radius, mEdgeImage, mEdgeImage->GetLargestPossibleRegion());
    int neighborhoodSize = 1;
    for(unsigned int i=0; i<InputImageDimension; i++)
        neighborhoodSize *= it.GetNeighborhood().GetSize()[i];

    int numseeds = kseedsIdx.size();
    for(int n=0; n<numseeds; n++)
    {
        DoubleScalarIndexType storedIdx = kseedsIdx[n];
        it.SetLocation(kseedsIdx[n]);

        for(int i=0; i<neighborhoodSize; i++)
        {
            bool inBounds;
            PixelDoubleType pxl = it.GetPixel(i, inBounds);

            if(inBounds && pxl < mEdgeImage->GetPixel(storedIdx))
                storedIdx = it.GetIndex(i);
        }

        if(storedIdx != kseedsIdx[n])
        {
            kseedsIdx[n] = storedIdx;
            kseedsPxl[n] = mLABInputImage->GetPixel(kseedsIdx[n]);
        }
    }
}


template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::EnforceLabelConnectivity(int& numlabels)
{
    const int SUPSZ = mImageVolume/mNumberSuperPixels;

    mEnforceConnectivityImage->FillBuffer(-1);

    itk::Offset<InputImageDimension> *off = new itk::Offset<InputImageDimension>[6];
    for(unsigned int i=0; i<InputImageDimension; i++) {
        off[2*i].Fill(0);
        off[2*i+1].Fill(0);

        off[2*i][i] = -1;
        off[2*i+1][i] = 1;
    }
    int nOff = InputImageDimension*2;

    int label = 0;
    int adjlabel = 0;

    LabelIndexType* vec = new LabelIndexType[mImageVolume];
    LabelIndexType idxO;

    itk::ImageRegionIteratorWithIndex<HelperScalarVoImageType> it(mEnforceConnectivityImage, mEnforceConnectivityImage->GetLargestPossibleRegion());

    while(!it.IsAtEnd())
    {
        idxO = it.GetIndex();

        if(0>it.Value())
        {
            mEnforceConnectivityImage->SetPixel(idxO, label);

            vec[0] = it.GetIndex();

            for(int n=0; n<nOff; n++)
            {
                LabelIndexType e = vec[0] + off[n];
                if(mEnforceConnectivityImage->GetLargestPossibleRegion().IsInside(e))
                {
                    idxO = e;

                    if(mEnforceConnectivityImage->GetPixel(idxO) >= 0)
                        adjlabel = mEnforceConnectivityImage->GetPixel(idxO);
                }
            }

            int count = 1;
            for(int c=0; c<count; c++)
            {
                for(int n=0; n<nOff; n++)
                {
                    LabelIndexType e = vec[c] + off[n];

                    if(mEnforceConnectivityImage->GetLargestPossibleRegion().IsInside(e))
                    {
                        if( 0>mEnforceConnectivityImage->GetPixel(e) && mLabelImage->GetPixel(idxO)==mLabelImage->GetPixel(e) )
                        {
                            vec[count] = e;
                            mEnforceConnectivityImage->SetPixel(e, label);
                            count++;
                        }
                    }

                }
            }

            if(count <= SUPSZ >> 2)
            {
                for(int c=0; c<count; c++)
                    mEnforceConnectivityImage->SetPixel(vec[c], adjlabel);
                label--;
            }
            label++;
        }
        ++it;
    }

    numlabels = label;

    if(vec) delete [] vec;
    if(off) delete [] off;

    itk::ImageRegionConstIteratorWithIndex<HelperScalarVoImageType> itCopy(mEnforceConnectivityImage, mEnforceConnectivityImage->GetLargestPossibleRegion());
    while(!itCopy.IsAtEnd())
    {
        mLabelImage->SetPixel(itCopy.GetIndex(), itCopy.Get()+1);
        ++itCopy;
    }
}


template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::PerformSuperpixelSegmentation_VariableSandM(std::vector<LABPixelType>& kseedsPxl,
        std::vector<LABIndexType>& kseedsIdx, const int& STEP, const int& NUMITR)
{
    const int numk = kseedsIdx.size();
    int numitr = 0;

    //----------------
    int offset = STEP;
    if(STEP < 10) offset = STEP*1.5;
    //----------------

    int maxxyzfac=0;
    if(InputImageDimension==2)      maxxyzfac = STEP*STEP;
    else if(InputImageDimension==3) maxxyzfac = STEP*STEP*STEP;

    std::vector<LABPixelType> sigmalab(numk, 0);
    std::vector<LABPixelType> sigmaxyz(numk, 0);    //use pixel (instead of index) type to have double array
    std::vector<int> clustersize(numk, 0);
    std::vector<double> inv(numk, 0);               //to store 1/clustersize[k] values
    std::vector<double> maxlab(numk, 10*10);        //THIS IS THE VARIABLE VALUE OF M, just start with 10
    std::vector<double> maxxyz(numk, maxxyzfac);    //THIS IS THE VARIABLE VALUE OF M, just start with 10

    typename DoubleScalarImageType::Pointer distxyz = DoubleScalarImageType::New();
    distxyz->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    distxyz->Allocate();
    distxyz->FillBuffer(DBL_MAX);
    typename DoubleScalarImageType::Pointer distlab = DoubleScalarImageType::New();
    distlab->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    distlab->Allocate();
    distlab->FillBuffer(DBL_MAX);
    typename DoubleScalarImageType::Pointer distvec = DoubleScalarImageType::New();
    distvec->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    distvec->Allocate();
    distvec->FillBuffer(DBL_MAX);

    double invxyzwt = 1.0/(STEP*STEP);           //NOTE: this is different from how usual SLIC/LKM works

    typename LABScalarImageType::SizeType radius;
    radius.Fill(offset);

    itk::ConstNeighborhoodIterator<LABScalarImageType> itLABNeigh(radius, mLABInputImage, mLABInputImage->GetLargestPossibleRegion());
    itk::ImageRegionIteratorWithIndex<LABScalarImageType> itLAB(mLABInputImage, mLABInputImage->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<LabelImageType> itLI(mLabelImage, mLabelImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator<DoubleScalarImageType> itDLAB(distlab, distlab->GetLargestPossibleRegion());
    itk::ImageRegionIterator<DoubleScalarImageType> itDXYZ(distxyz, distxyz->GetLargestPossibleRegion());

    int neighborhoodSize = 1;
    for(unsigned int i=0; i<InputImageDimension; i++)
        neighborhoodSize *= itLABNeigh.GetNeighborhood().GetSize()[i];

    while(numitr < NUMITR)
    {
        numitr++;
        distvec->FillBuffer(DBL_MAX);

        for(int n=0; n<numk; n++)
        {
            itLABNeigh.SetLocation(kseedsIdx[n]);

            for(int i=0; i<neighborhoodSize; i++)
            {
                bool inBounds;
                LABPixelType pxl = itLABNeigh.GetPixel(i, inBounds);

                if(inBounds) {
                    LABIndexType idx = itLABNeigh.GetIndex(i);

                    double d = 0;
                    for(unsigned int j=0; j<pxl.GetNumberOfComponents(); j++)
                        d += pow(pxl[j] - kseedsPxl[n][j], 2);
                    distlab->SetPixel(idx, d);

                    d = 0;
                    for(unsigned int j=0; j<InputImageDimension; j++)
                        d += pow((double)(idx[j] - kseedsIdx[n][j]), 2);
                    distxyz->SetPixel(idx, d);

                    //------------------------------------------------------------------------
                    double dist = distlab->GetPixel(idx)/maxlab[n] + distxyz->GetPixel(idx)*invxyzwt;     //only varying m, prettier superpixels
                    //------------------------------------------------------------------------

                    if(dist < distvec->GetPixel(idx))
                    {
                        distvec->SetPixel(idx, dist);
                        mLabelImage->SetPixel(idx, n);
                    }
                }
            }
        }

        //-----------------------------------------------------------------
        // Assign the max color distance for a cluster
        //-----------------------------------------------------------------
        if(0 == numitr)
        {
            maxlab.assign(numk,1);
            maxxyz.assign(numk,1);
        }

        itLI.GoToBegin();
        itDLAB.GoToBegin();
        itDXYZ.GoToBegin();

        while(!itLI.IsAtEnd())
        {
            if(maxlab[ itLI.Value() ] < itDLAB.Value())
                maxlab[ itLI.Value() ] = itDLAB.Value();
            if(maxxyz[ itLI.Value() ] < itDXYZ.Value())
                maxxyz[ itLI.Value() ] = itDXYZ.Value();

            ++itLI;
            ++itDLAB;
            ++itDXYZ;
        }

        //-----------------------------------------------------------------
        // Recalculate the centroid and store in the seed values
        //-----------------------------------------------------------------
        sigmalab.assign(numk, 0);
        sigmaxyz.assign(numk, 0);
        clustersize.assign(numk, 0);

        itLI.GoToBegin();
        itLAB.GoToBegin();

        while(!itLI.IsAtEnd())
        {
            LabelPixelType pxl = itLI.Value();
            assert(itLI.Value() >= 0);

            sigmalab[ itLI.Value() ] += itLAB.Value();
            for(unsigned int i=0; i<InputImageDimension; i++)
                sigmaxyz[ itLI.Value() ][i] += itLAB.GetIndex()[i];

            clustersize[pxl]++;

            ++itLI;
            ++itLAB;
        }

        for(int k=0; k<numk; k++)
        {
//            assert(clustersize[k] > 0);
            if(clustersize[k] <= 0)
                clustersize[k] = 1;

            inv[k] = 1.0/double(clustersize[k]);    //computing inverse now to multiply, than divide later
        }

        for(int k=0; k<numk; k++)
        {
            kseedsPxl[k] = sigmalab[k]*inv[k];         //pxl = pxl*scalar
            for(unsigned int i=0; i<InputImageDimension; i++)
                kseedsIdx[k][i] = sigmaxyz[k][i]*inv[k];
        }
    }
}


template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::GetLABXYZSeeds_ForGivenK(std::vector<LABPixelType>& kseedsPxl, std::vector<LABIndexType>& kseedsIdx)
{
    std::fstream file;

//    file.open("/scratch/friebel/TestData/SLIC/debugFile.txt", fstream::out);
//    file.flags(fstream::left | fstream::scientific);

    if(InputImageDimension==2) {
        double step[2];
        if(mSuperpixelBySize) {
            step[0] = mSuperpixelSpacing[0] / mVoxelSpacing[0];
            step[1] = mSuperpixelSpacing[1] / mVoxelSpacing[1];
        }
        else {
            step[0] = sqrt(double(mImageVolume)/double(mNumberSuperPixels));
            step[1] = step[0];
        }

        int off[2];
        off[0] = step[0]/2;
        off[1] = step[1]/2;

        LABIndexType idx;

//        file << "mVoxelSpacing = " << mVoxelSpacing << ", mSuperpixelSpacing = " << mSuperpixelSpacing << std::endl;
//        file << "step = " << step[0] << ", " << step[1] << " off = " << off[0] << ", " << off[1] << std::endl;

        for(unsigned int y=0; y<mInputSize[1]; y++)
        {
            unsigned int Y = y*step[1] + off[1];
            if(Y > mInputSize[1]-1)
                break;

            for(unsigned int x=0; x<mInputSize[0]; x++)
            {
                unsigned int X = x*step[0] + off[0];
                if(X > mInputSize[0]-1)
                    break;

                idx[0] = X; idx[1] = Y;
//                file << "idx = " << idx << std::endl;

                kseedsPxl.push_back(mLABInputImage->GetPixel(idx));
                kseedsIdx.push_back(idx);
            }
        }
    }
    else {
        double step[3];
        if(mSuperpixelBySize) {
            step[0] = mSuperpixelSpacing[0] / mVoxelSpacing[0];
            step[1] = mSuperpixelSpacing[1] / mVoxelSpacing[1];
            step[2] = mSuperpixelSpacing[2] / mVoxelSpacing[2];
        }
        else {
            step[0] = pow(double(mImageVolume)/double(mNumberSuperPixels), 1./3.);
            step[1] = step[0];
            step[2] = step[0];
        }

        int off[3];
        off[0] = step[0]/2;
        off[1] = step[1]/2;
        off[2] = step[2]/2;

//        file << "mVoxelSpacing = " << mVoxelSpacing << ", mSuperpixelSpacing = " << mSuperpixelSpacing << std::endl;
//        file << "step = " << step[0] << ", " << step[1] << ", " << step[2] << " off = " << off[0] << ", " << off[1] << ", " << off[2] << std::endl;

        LABIndexType idx;

        for(unsigned int z=0; z<mInputSize[2]; z++)
        {
            unsigned int Z = z*step[2] + off[2];
            if(Z > mInputSize[2]-1)
                break;

            for(unsigned int y=0; y<mInputSize[1]; y++)
            {
                unsigned int Y = y*step[1] + off[1];
                if(Y > mInputSize[1]-1)
                    break;

                for(unsigned int x=0; x<mInputSize[0]; x++)
                {
                    unsigned int X = x*step[0] + off[0];
                    if(X > mInputSize[0]-1)
                        break;

                    idx[0] = X; idx[1] = Y; idx[2] = Z;
//                    file << "idx = " << idx << std::endl;

                    kseedsPxl.push_back(mLABInputImage->GetPixel(idx));
                    kseedsIdx.push_back(idx);
                }
            }
        }
    }
//    file.close();

    if(mSuperpixelBySize)
        mNumberSuperPixels = kseedsIdx.size();

    if(mPerturbSeeds)
        PerturbSeeds(kseedsPxl, kseedsIdx);
}


template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::PerformSLICO_ForGivenK(const double& m)
{
    std::vector<LABPixelType> kseedsPxl(0);
    std::vector<LABIndexType> kseedsIdx(0);

    mLabelImage = LabelImageType::New();
    mLabelImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    mLabelImage->Allocate();

    mEnforceConnectivityImage = HelperScalarVoImageType::New();
    mEnforceConnectivityImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    mEnforceConnectivityImage->Allocate();

    mLABInputImage = LABScalarImageType::New();
    mLABInputImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    mLABInputImage->Allocate();

    DoRGBtoLABConversion();
    std::cout << "DoRGBtoLABConversion finished" << std::endl;

    if(mPerturbSeeds)
        DetectLABEdges();
    std::cout << "DetectLABEdges finished" << std::endl;

    GetLABXYZSeeds_ForGivenK(kseedsPxl, kseedsIdx);
    std::cout << "GetLABXYSeeds_ForGivenK finished" << std::endl;

    int STEP = 0;
    if(InputImageDimension == 2)
        STEP = sqrt(double(mImageVolume)/double(mNumberSuperPixels)) + 2.0;   //adding a small value in the even the STEP size is too small.
    else if(InputImageDimension == 3)
        STEP = pow(double(mImageVolume)/double(mNumberSuperPixels), 1./3.) + 2.0;

    PerformSuperpixelSegmentation_VariableSandM(kseedsPxl, kseedsIdx, STEP, 10);
    std::cout << "PerformSuperpixelSegmentation_VariableSandM finished" << std::endl;

    int numlabels = kseedsIdx.size();
    EnforceLabelConnectivity(numlabels);
    std::cout << "EnforceLabelConnectivity finished" << std::endl;

    std::cout << "numlabels = " << numlabels << std::endl;
}


template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::GenerateData()
{
    InputConstImagePointer inputImage = this->GetInput();
    mInputSize = inputImage->GetLargestPossibleRegion().GetSize();

    this->AllocateOutputs();
    this->GetOutput()->FillBuffer(0);

    LabelImagePointer outputImage = this->GetOutput();

    mImageVolume = 1;
    for(unsigned int i=0; i<InputImageDimension; i++)
        mImageVolume *= mInputSize[i];

    PerformSLICO_ForGivenK(10.0);

    outputImage->Graft(mLabelImage);
}


template<class TInputImage, class TLabelImage> void SLICImageFilter<TInputImage, TLabelImage>::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

} // end namespace itk

