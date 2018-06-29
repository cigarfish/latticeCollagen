/*
 * ResampleNetwork.cpp
 *
 *  Created on: Nov 20, 2012
 *      Author: friebel
 */

#include "ResampleNetwork.h"

#include <iostream>
#include <fstream>

#if (ITK_VERSION_MAJOR >= 4)
#include <itkTIFFImageIO.h>
#endif

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMap.h"
#include "itkShapeLabelObject.h"
//#include "itkImageRandomNonRepeatingIteratorWithIndex.h"

#include <vtkGraphReader.h>
#include <vtkGraphWriter.h>
#include "vtkPoints.h"



ResampleNetwork::ResampleNetwork()
{
    // TODO Auto-generated constructor stub
    m_graphPath = resampleInputPath;
    m_graphFilename = "pruned_graph";
    m_graphFileExtension = ".txt";

    m_cellImagePath = resampleInputPath;
    m_cellImageFilename = "nuclei_step6_bin";
    m_cellImageFileExtension = ".tif";

	m_dppivImagePath = resampleInputPath;
    m_dppivImageFilename = "C2";
    m_dppivImageFileExtension = ".tif";

	m_rescaledImagePath = resampleOutputPath;
	m_rescaledImageFileExtension = ".tif";

    m_rad = rad;

    m_num_dim = dim;
    m_oldDim = new long int[dim];
    m_newDim = new long int[dim];

    m_oldDim[0] = xRes;
    m_oldDim[1] = yRes;
#if(dim==3)
    m_oldDim[2] = zRes;
#endif


    m_oldVoxSpacing = new double[dim];
    m_newVoxSpacing = new double[dim];

    m_oldVoxSpacing[0] = xScale;
    m_oldVoxSpacing[1] = yScale;
#if(dim==3)
    m_oldVoxSpacing[2] = zScale;
#endif

    m_oldVoxSpacingMiddlePointOffset = new double[dim];
    m_oldVoxSpacingMiddlePointOffset[0] = m_oldVoxSpacing[0]/2.0;
    m_oldVoxSpacingMiddlePointOffset[1] = m_oldVoxSpacing[1]/2.0;
#if(dim==3)
    m_oldVoxSpacingMiddlePointOffset[2] = m_oldVoxSpacing[2]/2.0;
#endif

    m_newVoxSpacing[0] = m_rad;
    m_newVoxSpacing[1] = m_rad;
#if(dim==3)
    m_newVoxSpacing[2] = m_rad;
#endif

    m_newVoxSpacingMiddlePointOffset = new double[dim];
    m_newVoxSpacingMiddlePointOffset[0] = m_newVoxSpacing[0]/2.0;
    m_newVoxSpacingMiddlePointOffset[1] = m_newVoxSpacing[1]/2.0;
#if(dim==3)
    m_newVoxSpacingMiddlePointOffset[2] = m_newVoxSpacing[2]/2.0;
#endif
}


ResampleNetwork::~ResampleNetwork()
{
    // TODO Auto-generated destructor stub
}


void ResampleNetwork::LoadGraphs()
{
    vtkSmartPointer<vtkGraphReader> reader = vtkSmartPointer<vtkGraphReader>::New();

    int i=0;
    while(true) {
        std::stringstream name;
        name << m_graphPath << m_graphFilename << i << m_graphFileExtension;
        std::cout << "load graph file " << name.str() << std::endl;

        reader->SetFileName(name.str().c_str());
        int s = reader->OpenVTKFile();
        if(s != 0) {
            reader->Update();
            vtkSmartPointer<vtkUndirectedGraph> undirectedGraph = vtkSmartPointer<vtkUndirectedGraph>::New();
            reader->GetOutput()->ToUndirectedGraph(undirectedGraph);

            m_graphs.push_back(undirectedGraph);

            reader->CloseVTKFile();
        }
        else
            break;

        i++;
    }
}


void ResampleNetwork::LoadSegmentationImages()
{
	ImageType::SpacingType s;
	for(int i=0; i<dim; i++)
		s = m_oldVoxSpacing;

#if(!useRandomCells)
	ReaderType::Pointer readerCells = ReaderType::New();
	readerCells->SetFileName(m_cellImagePath + m_cellImageFilename + m_cellImageFileExtension);
	readerCells->ReleaseDataBeforeUpdateFlagOn();
	readerCells->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerCells->SetImageIO( itk::TIFFImageIO::New() );
#endif
	readerCells->Update();

	m_nucleiImage = readerCells->GetOutput();
	m_nucleiImage->DisconnectPipeline();
	m_nucleiImage->SetSpacing(s);


	ReaderType::Pointer readerDPPIV = ReaderType::New();
	readerDPPIV->SetFileName(m_dppivImagePath + m_dppivImageFilename + m_dppivImageFileExtension);
	readerDPPIV->ReleaseDataBeforeUpdateFlagOn();
	readerDPPIV->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerDPPIV->SetImageIO( itk::TIFFImageIO::New() );
#endif
	readerDPPIV->Update();

	m_dppivChannel = readerDPPIV->GetOutput();
	m_dppivChannel->DisconnectPipeline();
	m_dppivChannel->SetSpacing(s);
#endif

    ReaderType::Pointer readerSinusoids = ReaderType::New();
    readerSinusoids->SetFileName(m_cellImagePath + "sinus_step4_bin" + m_cellImageFileExtension);
#if (ITK_VERSION_MAJOR >= 4)
    readerSinusoids->SetImageIO( itk::TIFFImageIO::New() );
#endif
    readerSinusoids->Update();

    m_sinusoidImage = readerSinusoids->GetOutput();
    m_sinusoidImage->DisconnectPipeline();
    m_sinusoidImage->SetSpacing(s);


    ReaderType::Pointer readerCV = ReaderType::New();
    readerCV->SetFileName(m_rescaledImagePath + "vein_central_bin" + m_rescaledImageFileExtension);
#if (ITK_VERSION_MAJOR >= 4)
    readerCV->SetImageIO( itk::TIFFImageIO::New() );
#endif
    readerCV->Update();

    m_cvImage = readerCV->GetOutput();
    m_cvImage->DisconnectPipeline();
    m_cvImage->SetSpacing(s);


    ReaderType::Pointer readerPV = ReaderType::New();
    readerPV->SetFileName(m_rescaledImagePath + "vein_portal_bin" + m_rescaledImageFileExtension);
#if (ITK_VERSION_MAJOR >= 4)
    readerPV->SetImageIO( itk::TIFFImageIO::New() );
#endif
    readerPV->Update();

    m_pvImage = readerPV->GetOutput();
    m_pvImage->DisconnectPipeline();
    m_pvImage->SetSpacing(s);
}


void ResampleNetwork::RescaleGraphsToNewGrid()
{
    //***Initialize the new image********************************
    ImageType::IndexType start;
    start.Fill(0);

    ImageType::SizeType size;
    size[0] = std::ceil((m_oldVoxSpacing[0] * m_oldDim[0]) / m_newVoxSpacing[0]);
    size[1] = std::ceil((m_oldVoxSpacing[1] * m_oldDim[1]) / m_newVoxSpacing[1]);
#if(dim==3)
    size[2] = std::ceil((m_oldVoxSpacing[2] * m_oldDim[2]) / m_newVoxSpacing[2]);
#endif

    ImageType::RegionType region(start, size);

    m_rescaledNetworkImage = ImageType::New();
    m_rescaledNetworkImage->SetRegions(region);
    m_rescaledNetworkImage->Allocate();
    m_rescaledNetworkImage->FillBuffer(0);
    //***********************************************************

    //***Resample graph to new grid******************************
    std::map<vtkIdType, ImageType::IndexType> *vertexIndexToNewPos;
    vertexIndexToNewPos = new std::map<vtkIdType, ImageType::IndexType>[m_graphs.size()];

    for(unsigned int i=0; i<m_graphs.size(); i++) {
        std::cout << "processing graph " << i << std::endl;
        vtkSmartPointer<vtkPoints> vertexPos = m_graphs[i]->GetPoints();

        for(int j=0; j<m_graphs[i]->GetNumberOfVertices(); j++) {
            double *pos;
            pos = vertexPos->GetPoint(j);

            ImageType::IndexType newPos;
            for(unsigned int k=0; k<m_num_dim; k++) {
                pos[k] *= m_oldVoxSpacing[k];
                newPos[k] = std::floor((pos[k] + m_oldVoxSpacingMiddlePointOffset[k])/m_newVoxSpacing[k]);
            }

            vertexIndexToNewPos[i][j] = newPos;
            m_rescaledNetworkImage->SetPixel(newPos, 255);
        }
    }
    //***********************************************************

    //***Build 6-neighborhood-connections between connected vertices,if missing********************
#if(dim==3)
	//3D connection building
    for(unsigned int i=0; i<m_graphs.size(); i++) {
        for(unsigned int j=0; j<m_graphs[i]->GetNumberOfVertices(); j++) {
            ImageType::IndexType thisIdx = vertexIndexToNewPos[i][j];

            for(unsigned int k=0; k<m_graphs[i]->GetDegree(j); k++) {
                ImageType::IndexType neighborIdx = vertexIndexToNewPos[i][m_graphs[i]->GetOutEdge(j, k).Target];

                ImageType::OffsetType offsets[3];

                offsets[0][0] = thisIdx[0]-neighborIdx[0];  offsets[0][1] = 0;                          offsets[0][2] = 0;
                offsets[1][0] = 0;                          offsets[1][1] = thisIdx[1]-neighborIdx[1];  offsets[1][2] = 0;
                offsets[2][0] = 0;                          offsets[2][1] = 0;                          offsets[2][2] = thisIdx[2]-neighborIdx[2];


                ImageType::IndexType sixNeighborsOfV1[3];         //thisIdx = V1
                sixNeighborsOfV1[0] = thisIdx-offsets[0];
                sixNeighborsOfV1[1] = thisIdx-offsets[1];
                sixNeighborsOfV1[2] = thisIdx-offsets[2];

                ImageType::IndexType sixNeighborsOfV2[3];         //neighborIdx = V2
                sixNeighborsOfV2[0] = neighborIdx+offsets[0];
                sixNeighborsOfV2[1] = neighborIdx+offsets[1];
                sixNeighborsOfV2[2] = neighborIdx+offsets[2];

                std::vector<bool> isThereSixNeighborsOfV1;
                std::vector<bool> isThereSixNeighborsOfV2;

                for(unsigned int l=0; l<3; l++) {
                    (m_rescaledNetworkImage->GetPixel( sixNeighborsOfV1[l] )==255) ?
                            isThereSixNeighborsOfV1.push_back(true) : isThereSixNeighborsOfV1.push_back(false);
                }
                for(unsigned int l=0; l<3; l++) {
                    (m_rescaledNetworkImage->GetPixel( sixNeighborsOfV2[l] )==255) ?
                            isThereSixNeighborsOfV2.push_back(true) : isThereSixNeighborsOfV2.push_back(false);
                }

                bool isThereSixNeighborhoodConnectionFromV1ToV2 = false;
                for(unsigned int l=0; l<isThereSixNeighborsOfV1.size(); l++) {
                    if(isThereSixNeighborsOfV1[l]) {
                        for(unsigned int m=0; m<isThereSixNeighborsOfV2.size(); m++) {
                            if(isThereSixNeighborsOfV2[m]) {
                                ImageType::OffsetType o = sixNeighborsOfV1[l] - sixNeighborsOfV2[m];
                                int accO = std::abs(o[0]) + std::abs(o[1]) + std::abs(o[2]);
                                if(accO==1 || accO==0) {
                                    isThereSixNeighborhoodConnectionFromV1ToV2 = true;
                                    break;
                                }
                            }
                        }
                    }
                    if(isThereSixNeighborhoodConnectionFromV1ToV2)
                        break;
                }
                if(!isThereSixNeighborhoodConnectionFromV1ToV2) {
                    bool pathClosed = false;
                    for(unsigned int l=0; l<isThereSixNeighborsOfV1.size(); l++) {
                        if(isThereSixNeighborsOfV1[l]) {
                            ImageType::OffsetType o = {{sixNeighborsOfV1[l][0]-neighborIdx[0], 0, 0}};
                            int accO = std::abs(o[0]) + std::abs(o[1]) + std::abs(o[2]);
                            if(accO==0) {
								o[0] = 0;
								o[1] = sixNeighborsOfV1[l][1]-neighborIdx[1];
								o[2] = 0;
							}
                            m_rescaledNetworkImage->SetPixel(neighborIdx + o, 255);

                            pathClosed = true;
                            break;
                        }
                    }
                    if(!pathClosed) {
                        for(unsigned int l=0; l<isThereSixNeighborsOfV2.size(); l++) {
                            if(isThereSixNeighborsOfV2[l]) {
                                ImageType::OffsetType o = {{sixNeighborsOfV2[l][0]-thisIdx[0], 0, 0}};
                                int accO = std::abs(o[0]) + std::abs(o[1]) + std::abs(o[2]);
                                if(accO==0) {
									o[0] = 0;
									o[1] = sixNeighborsOfV2[l][1]-thisIdx[1];
									o[2] = 0;
								}
                                m_rescaledNetworkImage->SetPixel(thisIdx + o, 255);

                                pathClosed = true;
                                break;
                            }
                        }
                    }
                    if(!pathClosed) {
                        m_rescaledNetworkImage->SetPixel(sixNeighborsOfV1[0], 255);
                        m_rescaledNetworkImage->SetPixel(sixNeighborsOfV2[1], 255);
                    }
//                    std::cout << "between vertex at " << thisIdx[0] << ", " << thisIdx[1] << " and the vertex at " << neighborIdx[0] << ", " << neighborIdx[1] << " is no connection" << std::endl;
                }
//                else
//                    std::cout << "between vertex at " << thisIdx[0] << ", " << thisIdx[1] << " and the vertex at " << neighborIdx[0] << ", " << neighborIdx[1] << " is a connection" << std::endl;
            }
        }
    }
#endif

#if(dim==2)
    //2D-connection building
    for(unsigned int i=0; i<m_graphs.size(); i++) {
        for(unsigned int j=0; j<m_graphs[i]->GetNumberOfVertices(); j++) {
            ImageType::IndexType thisIdx = vertexIndexToNewPos[i][j];

            for(unsigned int k=0; k<m_graphs[i]->GetDegree(j); k++) {
                ImageType::IndexType neighborIdx = vertexIndexToNewPos[i][m_graphs[i]->GetOutEdge(j, k).Target];

				ImageType::OffsetType offsets[2];

                offsets[0][0] = thisIdx[0]-neighborIdx[0]; offsets[0][1] = 0;
                offsets[1][0] = 0; offsets[1][1] = thisIdx[1]-neighborIdx[1];

                bool sixNeighborhoodConnection = false;
                for(unsigned int l=0; l<2; l++) {
                    if(m_rescaledNetworkImage->GetPixel( neighborIdx+offsets[l] )==255) {
                        sixNeighborhoodConnection = true;
                        break;
                    }
                }
                if(sixNeighborhoodConnection)
                    std::cout << "between vertex at " << thisIdx[0] << ", " << thisIdx[1] << " and the vertex at " << neighborIdx[0] << ", " << neighborIdx[1] << " is a connection" << std::endl;
                else {
                    m_rescaledNetworkImage->SetPixel(neighborIdx+offsets[0], 255);
                    std::cout << "between vertex at " << thisIdx[0] << ", " << thisIdx[1] << " and the vertex at " << neighborIdx[0] << ", " << neighborIdx[1] << " is no connection" << std::endl;
                }
            }
        }
    }
#endif
	//*********************************************************************************************
}


void ResampleNetwork::RescaleCellsToNewGrid()
{
	//**GENERATE*CELL*SEGMENTATION*****************************************************************
#if(!loadIntensityMap)
	//----------DISTANCE-MAP---------------------------------------------------------------------------------------------------
	DanielssonDistanceMapImageFilterType::Pointer distanceMapNuc = DanielssonDistanceMapImageFilterType::New();
	distanceMapNuc->SetInput(m_nucleiImage);

    DanielssonDistanceMapImageFilterType::Pointer distanceMapSin = DanielssonDistanceMapImageFilterType::New();
    distanceMapSin->SetInput(m_sinusoidImage);

    AddImageFilterType::Pointer add = AddImageFilterType::New();
    add->SetInput1(distanceMapNuc->GetOutput());
    add->SetInput2(distanceMapSin->GetOutput());

    DivideImageFilterType::Pointer divide = DivideImageFilterType::New();
    divide->SetInput1(distanceMapNuc->GetOutput());
    divide->SetInput2(add->GetOutput());


    RescaleImageFilterType3::Pointer distMapRescaler = RescaleImageFilterType3::New();
    distMapRescaler->ReleaseDataFlagOn();
    distMapRescaler->SetInput(divide->GetOutput());

    ShoWriterType::Pointer writer5 = ShoWriterType::New();
    writer5->ReleaseDataFlagOn();
    writer5->SetFileName(m_rescaledImagePath + "combDistMap" + m_rescaledImageFileExtension);
    writer5->SetInput(distMapRescaler->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer5->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer5->Update();
    //-------------------------------------------------------------------------------------------------------------------------
#else
    ShoReaderType::Pointer readerIntensityMap = ShoReaderType::New();
    readerIntensityMap->SetFileName(m_rescaledImagePath + "invertIntensity" + m_rescaledImageFileExtension);
    readerIntensityMap->ReleaseDataBeforeUpdateFlagOn();
    readerIntensityMap->UseStreamingOn();
#if (ITK_VERSION_MAJOR >= 4)
    readerIntensityMap->SetImageIO( itk::TIFFImageIO::New() );
#endif

    RescaleImageFilterType2::Pointer intensityMapRescaler = RescaleImageFilterType2::New();
    intensityMapRescaler->ReleaseDataFlagOn();
    intensityMapRescaler->SetInput(readerIntensityMap->GetOutput());
    intensityMapRescaler->Update();

    IntImageType::Pointer rescaledImage = intensityMapRescaler->GetOutput();
    rescaledImage->DisconnectPipeline();

    ImageType::SpacingType s;
    for(int i=0; i<dim; i++)
        s = m_oldVoxSpacing;
    rescaledImage->SetSpacing(s);
#endif

	//----------WATERSHED-----------------------------------------------------------------------------------------------------
	MorphoWatershedImageFilterType::Pointer morphWatershed = MorphoWatershedImageFilterType::New();
	morphWatershed->SetLevel(0.20);
	morphWatershed->FullyConnectedOn();
	morphWatershed->ReleaseDataFlagOn();
#if(!loadIntensityMap)
	morphWatershed->SetInput(divide->GetOutput());
#else
	morphWatershed->SetInput(rescaledImage);
#endif
	//-------------------------------------------------------------------------------------------------------------------------


	InvertIntensityImageFilterType::Pointer invertIntensity1 = InvertIntensityImageFilterType::New();
	invertIntensity1->SetInput(m_sinusoidImage);

	MaskImageFilterType::Pointer maskImage1 = MaskImageFilterType::New();
    maskImage1->ReleaseDataFlagOn();
    maskImage1->SetInput1(morphWatershed->GetOutput());
    maskImage1->SetInput2(invertIntensity1->GetOutput());


    InvertIntensityImageFilterType::Pointer invertIntensity2 = InvertIntensityImageFilterType::New();
    invertIntensity2->SetInput(m_cvImage);

    MaskImageFilterType::Pointer maskImage2 = MaskImageFilterType::New();
    maskImage2->ReleaseDataFlagOn();
    maskImage2->SetInput1(maskImage1->GetOutput());
    maskImage2->SetInput2(invertIntensity2->GetOutput());


    InvertIntensityImageFilterType::Pointer invertIntensity3 = InvertIntensityImageFilterType::New();
    invertIntensity3->SetInput(m_pvImage);

    MaskImageFilterType::Pointer maskImage3 = MaskImageFilterType::New();
    maskImage3->ReleaseDataFlagOn();
    maskImage3->SetInput1(maskImage2->GetOutput());
    maskImage3->SetInput2(invertIntensity3->GetOutput());

	//----------TO-LABEL-MAP-AND-BACK-----------------------------------------------------------------------------------------
	LabelImageToShapeLabelMapFilterType::Pointer watershedImageToLabelMap = LabelImageToShapeLabelMapFilterType::New();
	watershedImageToLabelMap->ReleaseDataFlagOn();
	watershedImageToLabelMap->SetInput(maskImage3->GetOutput());

    LabelMapToLabelImageFilterType::Pointer watershedLabelMapToImage = LabelMapToLabelImageFilterType::New();
    watershedLabelMapToImage->ReleaseDataFlagOn();
    watershedLabelMapToImage->SetInput(watershedImageToLabelMap->GetOutput());
    watershedLabelMapToImage->Update();

	m_cellImage = watershedLabelMapToImage->GetOutput();
	m_cellImage->DisconnectPipeline();

	RescaleImageFilterType::Pointer rescaler2 = RescaleImageFilterType::New();
	rescaler2->ReleaseDataFlagOn();
	rescaler2->SetInput(m_cellImage);

	ShoWriterType::Pointer writer4 = ShoWriterType::New();
	writer4->ReleaseDataFlagOn();
	writer4->SetFileName(m_rescaledImagePath + "cellLabelMap" + m_rescaledImageFileExtension);
	writer4->SetInput(rescaler2->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
	writer4->SetImageIO( itk::TIFFImageIO::New() );
#endif
	writer4->Update();
	//-------------------------------------------------------------------------------------------------------------------------

    //----------VISULALIZATION-ONLY-FILTER---NOT-NEEDED-BY-OTHER-PIPELINE-FILTERS----------------------------------------------
    LabelOverlayImageFilterType::Pointer watershedOverlayImageFilter = LabelOverlayImageFilterType::New();
    watershedOverlayImageFilter->SetOpacity(0.6);
    watershedOverlayImageFilter->ReleaseDataFlagOn();
    watershedOverlayImageFilter->SetInput(m_dppivChannel);
    watershedOverlayImageFilter->SetLabelImage(m_cellImage);

	RGBWriterType::Pointer writer = RGBWriterType::New();
	writer->ReleaseDataFlagOn();
	writer->SetFileName(m_rescaledImagePath + "cell_segmentation" + m_rescaledImageFileExtension);
	writer->SetInput(watershedOverlayImageFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
	writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
	writer->Update();
	//*********************************************************************************************

	//**ALLOCATE*CELL*IMAGE************************************************************************
    ImageType::IndexType startCellImage;
    startCellImage.Fill(0);

    ImageType::SizeType sizeCellImage;
    sizeCellImage[0] = std::ceil((m_oldVoxSpacing[0] * m_oldDim[0]) / m_newVoxSpacing[0]);
    sizeCellImage[1] = std::ceil((m_oldVoxSpacing[1] * m_oldDim[1]) / m_newVoxSpacing[1]);
#if(dim==3)
    sizeCellImage[2] = std::ceil((m_oldVoxSpacing[2] * m_oldDim[2]) / m_newVoxSpacing[2]);
#endif

    ImageType::RegionType region(startCellImage, sizeCellImage);

    m_rescaledCellsImage = IntImageType::New();
    m_rescaledCellsImage->SetRegions(region);
    m_rescaledCellsImage->Allocate();
    m_rescaledCellsImage->FillBuffer(0);

    m_newDim[0] = sizeCellImage[0];
    m_newDim[1] = sizeCellImage[1];
    m_newDim[2] = sizeCellImage[2];
	//*********************************************************************************************

	//**RESCALE*CELLS*TO*NEW*LATTICE***************************************************************
    std::cout << "Now the rescaling begins.." << std::endl;

	ImageType::SizeType sizeBackup;
	sizeBackup[0] = std::ceil(m_newVoxSpacing[0] / m_oldVoxSpacing[0]);
	sizeBackup[1] = std::ceil(m_newVoxSpacing[1] / m_oldVoxSpacing[1]);
	sizeBackup[2] = std::ceil(m_newVoxSpacing[2] / m_oldVoxSpacing[2]);

	for(unsigned int i=0; i<m_newDim[0]; i++) {
		for(unsigned int j=0; j<m_newDim[1]; j++) {
			for(unsigned int k=0; k<m_newDim[2]; k++) {
			    ImageType::SizeType size;
			    size = sizeBackup;

				ImageType::IndexType idxToSet;
                idxToSet[0] = i;
                idxToSet[1] = j;
                idxToSet[2] = k;

				if(m_rescaledNetworkImage->GetPixel(idxToSet)!=255) {
					ImageType::IndexType start;
					start[0] = std::floor((m_newVoxSpacing[0] * i) / m_oldVoxSpacing[0]);
					start[1] = std::floor((m_newVoxSpacing[1] * j) / m_oldVoxSpacing[1]);
					start[2] = std::floor((m_newVoxSpacing[2] * k) / m_oldVoxSpacing[2]);

					for(int m=0; m<dim; m++)
					    if(start[m]+size[m]>=m_oldDim[m]) size[m] = m_oldDim[m]-1-start[m];

					ImageType::RegionType region(start, size);

					std::cout << "Lattice voxel " << i << ", " << j << ", " << k << std::endl;
					std::cout << "covers in original image region at " << start[0] << ", " << start[1] << ", " << start[2] << " of size " << size[0] << ", " << size[1] << ", " << size[2] << std::endl;

					typedef itk::ImageRegionConstIterator<IntImageType> IteratorType;
					IteratorType it(m_cellImage, region);

					std::map<int, int> labelIDToFrequency;
					for(it.GoToBegin(); !it.IsAtEnd(); ++it) {
						if(labelIDToFrequency.count(it.Get()) > 0)
							labelIDToFrequency[it.Get()]++;
						else
							labelIDToFrequency[it.Get()] = 0;

					}
					int max = 0;
					int highscoreID;
					for(std::map<int,int>::iterator iter = labelIDToFrequency.begin(); iter != labelIDToFrequency.end(); iter++) {
						std::cout << "cellID " << (*iter).first << " with freq " << (*iter).second << " was found" << std::endl;
						if((*iter).second>max)
							highscoreID = (*iter).first;
					}
					std::cout << "and thus this lattice voxel is stained " << highscoreID << std::endl;
					m_rescaledCellsImage->SetPixel(idxToSet, highscoreID);
				}
			}
		}
	}
	//*********************************************************************************************
}


void ResampleNetwork::AddRandomCellsToNewGrid()
{
    //allocate cell image
    ImageType::IndexType startCellImage;
    startCellImage.Fill(0);

    ImageType::SizeType sizeCellImage;
    sizeCellImage[0] = std::ceil((m_oldVoxSpacing[0] * m_oldDim[0]) / m_newVoxSpacing[0]);
    sizeCellImage[1] = std::ceil((m_oldVoxSpacing[1] * m_oldDim[1]) / m_newVoxSpacing[1]);
#if(dim==3)
    sizeCellImage[2] = std::ceil((m_oldVoxSpacing[2] * m_oldDim[2]) / m_newVoxSpacing[2]);
#endif

    ImageType::RegionType region(startCellImage, sizeCellImage);

    m_rescaledCellsImage = IntImageType::New();
    m_rescaledCellsImage->SetRegions(region);
    m_rescaledCellsImage->Allocate();
    m_rescaledCellsImage->FillBuffer(0);

    ImageType::Pointer randomNucleiImage = ImageType::New();
    randomNucleiImage->SetRegions(region);
    randomNucleiImage->Allocate();
    randomNucleiImage->FillBuffer(0);

    m_newDim[0] = sizeCellImage[0];
    m_newDim[1] = sizeCellImage[1];
    m_newDim[2] = sizeCellImage[2];
    //----------------------------

    //compute volume of data stack
    double dataSetVol = 1;

    dataSetVol *= m_newDim[0];
    dataSetVol *= m_newDim[1];
    dataSetVol *= m_newDim[2];

    std::cout << "data set volume " << dataSetVol << std::endl;
    //----------------------------

    //compute volume of sinusoidal network
    typedef itk::BinaryImageToShapeLabelMapFilter<ImageType> ImageToShapeLabelMapFilterType;

    ImageToShapeLabelMapFilterType::Pointer imageToShaLabMapFilter = ImageToShapeLabelMapFilterType::New();
    imageToShaLabMapFilter->SetInput(m_rescaledNetworkImage);
    imageToShaLabMapFilter->Update();

    int sinusoidVol = 0;

    int numLabObjects = imageToShaLabMapFilter->GetOutput()->GetNumberOfLabelObjects();
    for(int i=0; i<numLabObjects; ++i) {
        imageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->Optimize();
        sinusoidVol += imageToShaLabMapFilter->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixels();
    }
    std::cout << "sinusoid volume " << sinusoidVol << std::endl;
    //------------------------------

    int numCells = std::ceil((double)(dataSetVol - sinusoidVol)/voxPerCell);      //for 7 micron scaling, one cell ~19 lattice voxel
    std::cout << "num cells = " << numCells << std::endl;


    for(int i=0; i<numCells; i++) {
        ImageType::IndexType randIdx;
        randIdx[0] = rand() % m_newDim[0];
        randIdx[1] = rand() % m_newDim[1];
        randIdx[2] = rand() % m_newDim[2];
        if(m_rescaledNetworkImage->GetPixel(randIdx)!=255)
            randomNucleiImage->SetPixel(randIdx, 255);
        else
            i--;
    }

//    typedef itk::ImageRandomNonRepeatingIteratorWithIndex<ImageType> RandomIteratorType;
//
//    RandomIteratorType randIter(randomNucleiImage, randomNucleiImage->GetLargestPossibleRegion());
//    randIter.SetNumberOfSamples(numCells);
//    randIter.ReinitializeSeed();
//    randIter.GoToBegin();
//
//    while(!randIter.IsAtEnd()) {
//        randIter.Set(255);
//        ++randIter;
//    }

    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(m_rescaledImagePath + "randNuclei" + m_rescaledImageFileExtension);
    writer->SetInput(randomNucleiImage);
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();

    DanielssonDistanceMapImageFilterType::Pointer distanceMapNuc = DanielssonDistanceMapImageFilterType::New();
    distanceMapNuc->SetInput(randomNucleiImage);

    DanielssonDistanceMapImageFilterType::Pointer distanceMapSin = DanielssonDistanceMapImageFilterType::New();
    distanceMapSin->SetInput(m_rescaledNetworkImage);

    AddImageFilterType::Pointer add = AddImageFilterType::New();
    add->SetInput1(distanceMapNuc->GetOutput());
    add->SetInput2(distanceMapSin->GetOutput());

    DivideImageFilterType::Pointer divide = DivideImageFilterType::New();
    divide->SetInput1(distanceMapNuc->GetOutput());
    divide->SetInput2(add->GetOutput());

    RescaleImageFilterType3::Pointer distMapRescaler = RescaleImageFilterType3::New();
    distMapRescaler->ReleaseDataFlagOn();
    distMapRescaler->SetInput(divide->GetOutput());

    ShoWriterType::Pointer writer5 = ShoWriterType::New();
    writer5->ReleaseDataFlagOn();
    writer5->SetFileName(m_rescaledImagePath + "combDistMap" + m_rescaledImageFileExtension);
    writer5->SetInput(distMapRescaler->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer5->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer5->Update();

    MorphoWatershedImageFilterType::Pointer morphWatershed = MorphoWatershedImageFilterType::New();
    morphWatershed->SetLevel(0.20);
    morphWatershed->FullyConnectedOn();
    morphWatershed->ReleaseDataFlagOn();
    morphWatershed->SetInput(divide->GetOutput());
    morphWatershed->MarkWatershedLineOff();


    InvertIntensityImageFilterType::Pointer invertIntensity1 = InvertIntensityImageFilterType::New();
    invertIntensity1->SetInput(m_rescaledNetworkImage);

    MaskImageFilterType::Pointer maskImage1 = MaskImageFilterType::New();
    maskImage1->ReleaseDataFlagOn();
    maskImage1->SetInput1(morphWatershed->GetOutput());
    maskImage1->SetInput2(invertIntensity1->GetOutput());

    ReaderType::Pointer readerCV = ReaderType::New();
    readerCV->SetFileName(m_rescaledImagePath + "vein_central_bin" + m_rescaledImageFileExtension);
#if (ITK_VERSION_MAJOR >= 4)
    readerCV->SetImageIO( itk::TIFFImageIO::New() );
#endif

    InvertIntensityImageFilterType::Pointer invertIntensity2 = InvertIntensityImageFilterType::New();
    invertIntensity2->SetInput(readerCV->GetOutput());

    MaskImageFilterType::Pointer maskImage2 = MaskImageFilterType::New();
    maskImage2->ReleaseDataFlagOn();
    maskImage2->SetInput1(maskImage1->GetOutput());
    maskImage2->SetInput2(invertIntensity2->GetOutput());

    ReaderType::Pointer readerPV = ReaderType::New();
    readerPV->SetFileName(m_rescaledImagePath + "vein_portal_bin" + m_rescaledImageFileExtension);
#if (ITK_VERSION_MAJOR >= 4)
    readerPV->SetImageIO( itk::TIFFImageIO::New() );
#endif

    InvertIntensityImageFilterType::Pointer invertIntensity3 = InvertIntensityImageFilterType::New();
    invertIntensity3->SetInput(readerPV->GetOutput());

    MaskImageFilterType::Pointer maskImage3 = MaskImageFilterType::New();
    maskImage3->ReleaseDataFlagOn();
    maskImage3->SetInput1(maskImage2->GetOutput());
    maskImage3->SetInput2(invertIntensity3->GetOutput());
    maskImage3->Update();

    m_rescaledCellsImage = maskImage3->GetOutput();
    m_rescaledCellsImage->DisconnectPipeline();

    AddImageFilterType2::Pointer add2 = AddImageFilterType2::New();
    add2->SetInput1(m_rescaledNetworkImage);
    add2->SetInput2(m_rescaledCellsImage);

    ShoWriterType::Pointer writer1 = ShoWriterType::New();
    writer1->SetFileName(m_rescaledImagePath + "rescaledNetworkAndCells" + m_rescaledImageFileExtension);
    writer1->SetInput(add2->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer1->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer1->Update();
}


void ResampleNetwork::SaveRescaledNetworkAsImage()
{
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(m_rescaledImagePath + "rescaledNetwork" + m_rescaledImageFileExtension);
    writer->SetInput(m_rescaledNetworkImage);
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();
}


void ResampleNetwork::SaveRescaledNetworkAsText()
{
    itk::ImageRegionIterator<ImageType> iter(m_rescaledNetworkImage, m_rescaledNetworkImage->GetLargestPossibleRegion());

    std::fstream file;
    file.open((m_rescaledImagePath + "resampledGraph" + m_graphFileExtension).c_str(), std::fstream::out);

    while(!iter.IsAtEnd()) {
        if(iter.Get()==255)
            file << iter.GetIndex()[0] << ", " << iter.GetIndex()[1] << ", " << iter.GetIndex()[2] << std::endl;

        ++iter;
    }
    file.close();
}


void ResampleNetwork::SaveRescaledCellsAsImage()
{
    RescaleImageFilterType::Pointer rescaleFilter = RescaleImageFilterType::New();
    rescaleFilter->ReleaseDataFlagOn();
    rescaleFilter->SetInput(m_rescaledCellsImage);

    ShoWriterType::Pointer writer = ShoWriterType::New();
    writer->SetFileName(m_rescaledImagePath + "rescaledCells" + m_rescaledImageFileExtension);
    writer->SetInput(rescaleFilter->GetOutput());
#if (ITK_VERSION_MAJOR >= 4)
    writer->SetImageIO( itk::TIFFImageIO::New() );
#endif
    writer->Update();
}


void ResampleNetwork::SaveRescaledCellsAsText()
{
    itk::ImageRegionIterator<IntImageType> iter(m_rescaledCellsImage, m_rescaledCellsImage->GetLargestPossibleRegion());

    std::fstream file;
    file.open((m_rescaledImagePath + "resampledCells" + m_graphFileExtension).c_str(), std::fstream::out);

    while(!iter.IsAtEnd()) {
        if(iter.Get()!=0)
			file << iter.GetIndex()[0] << ", " << iter.GetIndex()[1] << ", " << iter.GetIndex()[2] << "       " << iter.Get() << std::endl;

        ++iter;
    }
    file.close();
}


void ResampleNetwork::Update()
{
    LoadGraphs();
    LoadSegmentationImages();
    RescaleGraphsToNewGrid();
#if(useRandomCells)
    AddRandomCellsToNewGrid();
#else
	RescaleCellsToNewGrid();
#endif
    SaveRescaledNetworkAsImage();
#if(dim==3)
        SaveRescaledNetworkAsText();
		SaveRescaledCellsAsText();
#endif
	SaveRescaledCellsAsImage();
}

