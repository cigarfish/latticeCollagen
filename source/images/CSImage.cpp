///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSImage.cpp                                                          //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2011-10-20 14:44:46                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "CSImage.h"

#include "vtkImageData.h"
#include "itkImageFileReader.h"
#include "vtkImageFlip.h"

#include <iostream>


CSImage::CSImage()
{
  mpITKImage  =   itkImageType::New();
  mpVTKImage  =   vtkImageData::New();
  mpConnector = itkToVTKFilter::New();
}


CSImage::~CSImage()
{
  mpConnector->Delete();
  mpVTKImage->Delete();
  mpITKImage->Delete();
}


void
CSImage::loadFromFile(const char * fileName)
{
  // read in the image to mpITKImage
  itk::ImageFileReader<itkImageType>::Pointer reader = itk::ImageFileReader<itkImageType>::New();
  reader->SetFileName(fileName);
  reader->Update();

  // itkImageType::Pointer saveITKImage = NULL;

  // if (mpITKImage)
  //   {
  //     saveITKImage = mpITKImage;
  //   }

  mpITKImage = reader->GetOutput();

  //reader->Delete();

  // if (saveITKImage)
  //   saveITKImage->Delete();

}


void
CSImage::loadFromFile( const std::string fileName )
{loadFromFile(fileName.c_str());}


void
CSImage::setImageData( vtkImageData * newImage )
{
  // TODO:
  // transform into itkImageType... in later revisions
  // for now:

  if (mpITKImage)
    {
      mpITKImage->Delete();
      mpITKImage = NULL;
    }
  
  if (!mpVTKImage)
    mpVTKImage = vtkImageData::New();

  mpVTKImage->DeepCopy(newImage);
}


void
CSImage::setImageData( itkImageType * /* newImage */)
{
}


void
CSImage::addITKFilter( CSITKImageFilter * itkFilter )
{
  mITKFilters.push_back(itkFilter);
}


void
CSImage::addVTKFilter( CSVTKImageFilter * vtkFilter )
{
  mVTKFilters.push_back(vtkFilter);
}


void
CSImage::update()
{
  if (mpITKImage)
    {
      // apply itk filters (first, only one)
      if ( mITKFilters.size() )
	{
	  CSITKImageFilter * iter = mITKFilters[0];
	  iter->SetInput(mpITKImage);

	  // convert to vtkImageData
	  mpConnector->SetInput(iter->GetOutput());
	}
      else
	{
	  mpConnector->SetInput(mpITKImage);
	}

      // flip the y axis:
      vtkImageFlip * flipper = vtkImageFlip::New();
      flipper->SetInput( mpConnector->GetOutput() );
      flipper->SetFilteredAxis(1);

      mpVTKImage = flipper->GetOutput();

      // itmImageToVTKImageFilter does, unfortunately, not take care
      // of the image size, so we now have to set it manually:
      itkImageType::RegionType wholeImageRegion = mpITKImage->GetLargestPossibleRegion();
      itkImageType::SizeType imageSize = wholeImageRegion.GetSize();

      int width  = imageSize[0];
      int height = imageSize[1];

      std::cout << "itk image size is:  " << width << "x" << height << std::endl;

      mpVTKImage->SetExtent(0,width-1,0,height-1,0,0);

      std::cout << "vtk image size is:  " << mpVTKImage->GetDimensions()[0] << "x" << mpVTKImage->GetDimensions()[1] << std::endl;
    }

  if ( mVTKFilters.size() )
    {
      // apply vtk filters (first, only one)
      CSVTKImageFilter * vter = mVTKFilters[0];

      vter->SetInput(mpVTKImage);

      // the result:
      mpVTKImage = vter->GetOutput();
    }

}
