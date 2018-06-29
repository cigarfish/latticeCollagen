///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  FileFormatConverter.h                                                //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-08-01                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef FILEFORMATCONVERTER_H_
#define FILEFORMATCONVERTER_H_

#include <sstream>

#include "itkAddImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkCastImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkGroupSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMaskedSpatialObjectToImageFilter.h"
#include "itkMaskNegatedImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkTubeSpatialObject.h"
#include "itkVesselTubeSpatialObject.h"
#include "itkVTKImageImport.h"

#include "vtkImageData.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include "BasePipeline.h"


class FileFormatConverter : public BasePipeline<3>
{
public:
    enum InputData {
        Graph_Mevis,
        Graph_Dresden,
        Triangualation_Dresden,
        Bit8_RGB_2D_Image,
        Bit8_RGB_3D_Image
    };

    enum OutputData {
        Graph_TiQuant,
        Binary_Image,
        Bit8_Greyscale_2D_Image,
        Bit8_Greyscale_3D_Image
    };

    static const unsigned int NumberInputFileFormats = 5;
    static const unsigned int NumberOutputFileFormats = 4;

    static const std::string InputDataName[];
    static const std::string OutputDataName[];

    static const std::vector<FileFormatConverter::OutputData> AvailableConversions(FileFormatConverter::InputData input);

private:
    typedef unsigned char                           CPixelType;
    typedef float                                   FPixelType;
    typedef unsigned int                            IPixelType;
    typedef unsigned long                           LPixelType;
    typedef itk::RGBPixel<CPixelType>               RGBCPixelType;

    typedef itk::Image<CPixelType, 2>               C2DImageType;
    typedef itk::Image<CPixelType, 3>               C3DImageType;
    typedef itk::Image<FPixelType, 3>               FImageType;
    typedef itk::Image<IPixelType, 3>               IImageType;
    typedef itk::Image<LPixelType, 3>               LImageType;
    typedef itk::Image<RGBCPixelType, 2>            RGBC2DImageType;
    typedef itk::Image<RGBCPixelType, 3>            RGBC3DImageType;
    typedef itk::Point<double, 3>                   PointType;

    typedef itk::TubeSpatialObject<3>               TubeType;
    typedef itk::TubeSpatialObjectPoint<3>          TubePointType;
    typedef TubePointType::CovariantVectorType      VectorType;
    typedef itk::VesselTubeSpatialObject<3>         VesselTubeType;
    typedef itk::VesselTubeSpatialObjectPoint<3>    VesselTubePointType;
    typedef itk::GroupSpatialObject<3>              GroupType;
    typedef GroupType::TransformType                TransformType;

    typedef itk::ImageFileReader<RGBC2DImageType>   RGBC2DReaderType;
    typedef itk::ImageFileReader<C3DImageType>      C3DReaderType;
    typedef itk::ImageFileWriter<C2DImageType>      C2DWriterType;
    typedef itk::ImageFileWriter<C3DImageType>      C3DWriterType;
    typedef itk::ImageFileWriter<IImageType>        IWriterType;

    typedef itk::AddImageFilter<C3DImageType, C3DImageType>                                                 AddImageFilterType;
    typedef itk::BinaryImageToShapeLabelMapFilter<C3DImageType>                                             ImageToShapeLabelMapFilterType;
    typedef itk::CastImageFilter<LImageType, IImageType>                                                    CastLIImageFilterType;
    typedef itk::ConstantBoundaryCondition<C3DImageType>                                                    BoundaryConditionType;
    typedef itk::LabelMapToLabelImageFilter<ImageToShapeLabelMapFilterType::OutputImageType, LImageType>    LabelMapToLabelImageFilterType;
    typedef itk::MaskedSpatialObjectToImageFilter<GroupType, C3DImageType>                                  SpatialObjectToImageFilterType;
    typedef itk::MaskNegatedImageFilter<C3DImageType, C3DImageType>                                         MaskImageFilterType;
    typedef itk::NeighborhoodIterator<C3DImageType, BoundaryConditionType>                                  NeighborhoodIteratorType;
    typedef itk::SignedMaurerDistanceMapImageFilter<C3DImageType, FImageType>                               SignedMaurerDistanceMapImageFilterType;
    typedef itk::VTKImageImport<C3DImageType>                                                               VTKImageImportType;

    typedef typename ImageToShapeLabelMapFilterType::OutputImageType                                        LabelMapType;
    typedef typename LabelMapType::LabelObjectType                                                          LabelObjectType;

public:
    FileFormatConverter();
    virtual ~FileFormatConverter();

    void Update();

private:
    void ParseParameterContext();
    void WriteLogFile(std::string timeStamp);


    static vtkSmartPointer<vtkPolyData> ReadSTLFile(std::string path);
    static vtkSmartPointer<vtkMutableDirectedGraph> ReadAndConvertMevisGraphToVTKGraph(std::string path);
    static std::vector< vtkSmartPointer<vtkMutableUndirectedGraph> > ReadAndConvertDresdenGraphToVTKGraph(std::string path);

    static vtkSmartPointer<vtkImageData> ConvertTriangulationToImageData(vtkSmartPointer<vtkPolyData> polyData);
    static void ConvertVtkToItkImageData(vtkSmartPointer<vtkImageData> vtkImageData, C3DImageType::Pointer &itkImageData);
    static void ConvertRGBToGreyscaleImageData(RGBC2DImageType::Pointer &rgbImageData, C2DImageType::Pointer &scalarImageData, int mode);
    static void ConvertRGBToGreyscaleImageData(RGBC3DImageType::Pointer &rgbImageData, C3DImageType::Pointer &scalarImageData, int mode);

    static void ReadGraph(std::string fullFilename, std::vector< vtkSmartPointer<vtkMutableUndirectedGraph> > &graphs);
    static void SaveGraph(std::string path, vtkSmartPointer<vtkMutableDirectedGraph> graph);
    static void SaveGraph(std::string path, vtkSmartPointer<vtkMutableUndirectedGraph> graph, double spacingX = 1, double spacingY = 1, double spacingZ = 1);
    static void ReadImage(std::string path, RGBC2DImageType::Pointer &image);
    static void SaveImage(std::string path, C2DImageType::Pointer image);
    static void SaveImage(std::string path, C3DImageType::Pointer image);

    static void ComputeBounds(vtkSmartPointer<vtkGraph> graph, double*& boundsBoth, double*& offset);
    static void ComputeBounds(vtkSmartPointer<vtkGraph> graph1, vtkSmartPointer<vtkGraph> graph2, double*& boundsBoth, double*& offset);

public:
    //Deprecated: not supported by new input -> output scheme
    static void ReadMevisDataAndConvertToImageData();
    //Deprecated: not supported by new input -> output scheme
    static void ConvertTreeToImage(vtkSmartPointer<vtkMutableDirectedGraph> graph, std::string pathEndSections, std::string pathRestSections, double *bounds, double *offset);
    //Deprecated: not supported by new input -> output scheme
    static void ComposeSegmentationFilesForInterfaceToStefan();
    //Deprecated: not supported by new input -> output scheme
    static void ComposeSegmentationFilesForInterfaceToStefan2();
    //Deprecated: not supported by new input -> output scheme
    static void ComposeGraphFromDresdenFormat();

private:
    std::string mInputFilename;
    std::string mPath;
    std::string mFilename;
    std::string mFilenameExtension;

    QFileInfo mInfoFullFilename;

    InputData mInputFileFormat;
    OutputData mOutputFileFormat;
};

#endif /* MEVISFILEFORMATCONVERTER_H_ */
