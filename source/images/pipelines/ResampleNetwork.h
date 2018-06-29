/*
 * ResampleNetwork.h
 *
 *  Created on: Nov 20, 2012
 *      Author: friebel
 */

#ifndef RESAMPLENETWORK_H_
#define RESAMPLENETWORK_H_

#include <string>
#include <vector>

#include "itkAddImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkUndirectedGraph.h>

//TODO: get rid of defines, only for fast prototyping phase
#define dim 3
#define rad 5.0

#define resampleInputPath "/home/ls1/friebel/Workspace/Data/Confocal/TestData/ResampleNetwork/3D/041/original_graph/"
#define resampleOutputPath "/home/ls1/friebel/Workspace/Data/Confocal/TestData/ResampleNetwork/3D/041/spacing_5micron/"

#define useRandomCells 1

#define loadIntensityMap 0

#define xRes 1024
#define yRes 1024
#define zRes 125
#define xScale 0.620996
#define yScale 0.620996
#define zScale 0.54

#define voxPerCell  53.0 //53 for 5micron, 19 for 7micron rescaling


class ResampleNetwork
{
	typedef itk::RGBPixel<unsigned char>        RGBPixelType;

	typedef itk::Image<unsigned char, dim>      ImageType;
	typedef itk::Image<int, dim>				IntImageType;
	typedef itk::Image<unsigned short, dim>     ShoImageType;
	typedef itk::Image<float, dim>              FloImageType;
	typedef itk::Image<RGBPixelType, dim>       RGBImageType;

    typedef itk::ImageFileReader<ImageType>     ReaderType;
    typedef itk::ImageFileReader<ShoImageType>  ShoReaderType;
    typedef itk::ImageFileWriter<ImageType>     WriterType;
    typedef itk::ImageFileWriter<ShoImageType>  ShoWriterType;
	typedef itk::ImageFileWriter<RGBImageType>  RGBWriterType;

	typedef itk::AddImageFilter<FloImageType, FloImageType, FloImageType>                                       AddImageFilterType;
	typedef itk::AddImageFilter<ImageType, IntImageType, ShoImageType>                                          AddImageFilterType2;
	typedef itk::DanielssonDistanceMapImageFilter<ImageType, FloImageType>										DanielssonDistanceMapImageFilterType;
	typedef itk::DivideImageFilter<FloImageType, FloImageType, FloImageType>                                    DivideImageFilterType;
	typedef itk::InvertIntensityImageFilter<ImageType>															InvertIntensityImageFilterType;
	typedef itk::InvertIntensityImageFilter<IntImageType>														IInvertIntensityImageFilterType;
	typedef itk::MaskImageFilter<IntImageType, ImageType, IntImageType>                                         MaskImageFilterType;
	typedef itk::MorphologicalWatershedImageFilter<FloImageType, IntImageType>									MorphoWatershedImageFilterType;
	typedef itk::LabelImageToShapeLabelMapFilter<MorphoWatershedImageFilterType::OutputImageType>               LabelImageToShapeLabelMapFilterType;
	typedef itk::LabelMapToLabelImageFilter<LabelImageToShapeLabelMapFilterType::OutputImageType, IntImageType>	LabelMapToLabelImageFilterType;
	typedef itk::LabelOverlayImageFilter<ImageType, IntImageType, RGBImageType>                    				LabelOverlayImageFilterType;
	typedef itk::RescaleIntensityImageFilter<IntImageType, ShoImageType>                                        RescaleImageFilterType;
	typedef itk::RescaleIntensityImageFilter<ShoImageType, IntImageType>                                        RescaleImageFilterType2;
	typedef itk::RescaleIntensityImageFilter<FloImageType, ShoImageType>                                        RescaleImageFilterType3;

public:
    ResampleNetwork();
    virtual ~ResampleNetwork();

    void Update();
private:
    void LoadGraphs();
    void LoadSegmentationImages();
    void RescaleGraphsToNewGrid();
    void AddRandomCellsToNewGrid();
	void RescaleCellsToNewGrid();
    void SaveRescaledNetworkAsImage();
    void SaveRescaledNetworkAsText();
    void SaveRescaledCellsAsImage();
	void SaveRescaledCellsAsText();

    std::vector< vtkSmartPointer<vtkUndirectedGraph> > m_graphs;
    std::vector< vtkSmartPointer<vtkUndirectedGraph> > m_rescaledGraphs;
    ImageType::Pointer m_rescaledNetworkImage;
	IntImageType::Pointer m_rescaledCellsImage;
	ImageType::Pointer m_dppivChannel;
	ImageType::Pointer m_nucleiImage;
	ImageType::Pointer m_sinusoidImage;
	ImageType::Pointer m_cvImage;
	ImageType::Pointer m_pvImage;
	IntImageType::Pointer m_cellImage;

    std::string m_graphPath;
    std::string m_graphFilename;
    std::string m_graphFileExtension;

    std::string m_cellImagePath;
    std::string m_cellImageFilename;
    std::string m_cellImageFileExtension;

	std::string m_dppivImagePath;
    std::string m_dppivImageFilename;
    std::string m_dppivImageFileExtension;

	std::string m_rescaledImagePath;
	std::string m_rescaledImageFileExtension;

    double m_rad;

    long int m_num_dim;
    long int *m_oldDim;
    long int *m_newDim;
    double *m_oldVoxSpacing;
    double *m_oldVoxSpacingMiddlePointOffset;
    double *m_newVoxSpacing;
    double *m_newVoxSpacingMiddlePointOffset;
};

#endif /* RESAMPLENETWORK_H_ */
