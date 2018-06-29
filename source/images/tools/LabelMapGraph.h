///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelMapGraph.h                                                      //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-09-04                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LABELMAPGRAPH_H_
#define LABELMAPGRAPH_H_

#include <itkImageBase.h>
#include <itkLabelGeometryImageFilter.h>
#include <itkMaskedImageToHistogramFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkRGBToLuminanceImageFilter.h>
#include <itkScalarImageToTextureFeaturesFilter.h>

#include <vtkMutableUndirectedGraph.h>
#include <vtkSmartPointer.h>


namespace itk
{

template< class TOriginalImage, class TLabelImage, class TLabelObject > class ITK_EXPORT LabelMapGraph : public ImageBase< ::itk::GetImageDimension<TLabelObject>::ImageDimension >
{
public:
    static const int NumberUnaryFeatures = 9;
    static const int NumberBinaryFeatures = 3;

    enum UnaryFeatures {
        MEAN_RED_INTENSITY,
        MEAN_GREEN_INTENSITY,
        MEAN_BLUE_INTENSITY,
        RED_HISTOGRAM,
        GREEN_HISTOGRAM,
        BLUE_HISTOGRAM,
        RGB_HISTOGRAM,
        TEXTURE_FEATURES,
        NUMBER_OF_PIXELS
    };

    enum BinaryFeatures {
        MEAN_RED_SIMILARITY,
        MEAN_GREEN_SIMILARITY,
        MEAN_BLUE_SIMILARITY
    };

    enum FeatureState {
        IsActive,
        IsInactive,
        IsNotAttached
    };

    static const std::string UnaryFeatureName[];
    static const std::string BinaryFeatureName[];

    static int UnaryFeatureComponents[];

    void AttachUnaryFeature(UnaryFeatures f);
    void ActivateUnaryFeature(UnaryFeatures f);
    void DeactivateUnaryFeature(UnaryFeatures f);
    void AttachBinaryFeature(BinaryFeatures f);
    void ActivateBinaryFeature(BinaryFeatures f);
    void DeactivateBinaryFeature(BinaryFeatures f);
    FeatureState GetUnaryFeatureState(UnaryFeatures f)      { return mUnaryFeatureState[f]; };
    FeatureState GetBinaryFeatureState(BinaryFeatures f)    { return mBinaryFeatureState[f]; };

    unsigned int GetNumberActiveUnaryFeatures()             { return mActiveUnaryFeatures; };
    unsigned int GetNumberActiveUnaryFeatureComponents()    { return mActiveUnaryFeatureComponents; };
    unsigned int GetNumberActiveBinaryFeatures()            { return mActiveBinaryFeatures; };

    static const int NUMBER_HISTOGRAM_BINS = 4;

private:
    FeatureState mUnaryFeatureState[NumberUnaryFeatures];
    FeatureState mBinaryFeatureState[NumberBinaryFeatures];

    unsigned int mActiveUnaryFeatures;
    unsigned int mActiveUnaryFeatureComponents;
    unsigned int mActiveBinaryFeatures;

public:
    /** Standard class typedefs */
    typedef LabelMapGraph                                                           Self;
    typedef ImageBase< ::itk::GetImageDimension< TLabelObject >::ImageDimension >   Superclass;
    typedef SmartPointer< Self >                                                    Pointer;
    typedef SmartPointer< const Self >                                              ConstPointer;
    typedef WeakPointer< const Self >                                               ConstWeakPointer;

    itkNewMacro(Self);
    itkTypeMacro(LabelMapGraph, ImageBase);

    typedef TLabelObject                        LabelObjectType;
    typedef typename LabelObjectType::Pointer   LabelObjectPointerType;
    typedef typename Superclass::SizeValueType  SizeValueType;
    typedef SizeValueType                       LengthType;
    typedef typename LabelObjectType::LabelType LabelType;
    typedef LabelType                           PixelType;

    itkStaticConstMacro(ImageDimension, unsigned int, LabelObjectType::ImageDimension);

    typedef std::vector< LabelType >                LabelVectorType;
    typedef std::vector< LabelObjectPointerType >   LabelObjectVectorType;
    typedef typename Superclass::IndexType          IndexType;
    typedef typename Superclass::OffsetType         OffsetType;
    typedef typename Superclass::SizeType           SizeType;
    typedef typename Superclass::DirectionType      DirectionType;
    typedef typename Superclass::RegionType         RegionType;
    typedef typename Superclass::SpacingType        SpacingType;
    typedef typename Superclass::PointType          PointType;
    typedef typename Superclass::OffsetValueType    OffsetValueType;

    typedef TOriginalImage                              OriginalImageType;
    typedef typename OriginalImageType::Pointer         OriginalImagePointer;
    typedef typename OriginalImageType::ConstPointer    OriginalImageConstPointer;
    typedef typename OriginalImageType::PixelType       OriginalImagePixelType;
    typedef typename OriginalImagePixelType::ComponentType  OriginalImagePixelComponentType;
    typedef typename OriginalImageType::IndexType       OriginalImageIndexType;
    typedef typename OriginalImageType::SizeType        OriginalImageSizeType;
    typedef typename OriginalImageType::RegionType      OriginalImageRegionType;
    typedef typename OriginalImageType::OffsetType      OriginalImageOffsetType;

    typedef TLabelImage                                 LabelImageType;
    typedef typename LabelImageType::Pointer            LabelImagePointer;
    typedef typename LabelImageType::ConstPointer       LabelImageConstPointer;
    typedef typename LabelImageType::PixelType          LabelImagePixelType;
    typedef typename LabelImageType::IndexType          LabelImageIndexType;
    typedef typename LabelImageType::SizeType           LabelImageSizeType;
    typedef typename LabelImageType::RegionType         LabelImageRegionType;
    typedef typename LabelImageType::OffsetType         LabelImageOffsetType;

    typedef itk::Image<LabelImagePixelType, ImageDimension> GrayscaleImageType;
    typedef typename GrayscaleImageType::Pointer            GrayscaleImagePointer;


    void SetOriginalImage(OriginalImagePointer image) { mOriginalImage = image; };
    void SetLabelImage(LabelImagePointer image) { mLabelImage = image; };
    void SetGraph(vtkSmartPointer<vtkMutableUndirectedGraph>  graph);
    static bool GraphHasFeature(vtkSmartPointer<vtkUndirectedGraph> graph, UnaryFeatures f);
    static bool GraphHasFeature(vtkSmartPointer<vtkUndirectedGraph> graph, BinaryFeatures f);

    vtkSmartPointer<vtkMutableUndirectedGraph> GetLabelMapGraph() { return mObjectGraph; };
    LabelImagePointer GetLabelImage() { return mLabelImage; };

    void SetRegions(const RegionType & region)
    {
        this->SetLargestPossibleRegion(region);
        this->SetBufferedRegion(region);
        this->SetRequestedRegion(region);
    }

    void SetRegions(const SizeType & size)
    {
        RegionType region; region.SetSize(size);

        this->SetLargestPossibleRegion(region);
        this->SetBufferedRegion(region);
        this->SetRequestedRegion(region);
    }

    virtual void Initialize();
    virtual void Allocate();
    virtual void Graft(const DataObject *data);

    LabelObjectType * GetLabelObject(const LabelType & label);
    const LabelObjectType * GetLabelObject(const LabelType & label) const;
    LabelObjectType * GetLabelObject(const IndexType & idx) const;
    LabelObjectVectorType GetLabelObjects() const;

    LabelObjectType * GetNthLabelObject(const SizeValueType & pos);
    const LabelObjectType * GetNthLabelObject(const SizeValueType & pos) const;

    LabelVectorType GetLabels() const;

    bool HasLabel(const LabelType label) const;
    typename Self::SizeValueType GetNumberOfLabelObjects() const { return mLabelObjectContainer.size(); };

    void AddLabelObject(LabelObjectType *labelObject);
    void AddLabelVertex(const IndexType &idx, const LabelType &label);
    void AddLabelEdge(const IndexType &idx, const LabelType &label);
    vtkIdType LazyMergeLabelObjects(vtkIdType o1, vtkIdType o2);
    vtkIdType MergeLabelObjects(vtkIdType o1, vtkIdType o2);
    vtkIdType LazyMergeLabelObjects(LabelObjectPointerType o1, LabelObjectPointerType o2);
    int LazyCollapseLabelObjects(std::set<vtkIdType> &labelsToCollapse);
    vtkIdType MergeLabelObjects(LabelObjectPointerType o1, LabelObjectPointerType o2);
    void PushLabelObject(LabelObjectType *labelObject);
    void RemoveLabelObject(LabelObjectType *labelObject);
    void RemoveLabelVertex(const LabelType &label);
    void RemoveLabel(const LabelType &label);
    void RemoveLabels(std::set<vtkIdType> labels);
    void KeepLabelsRemoveRest(std::set<vtkIdType> labelsToKeep);

    //Doesn't reset LabelImage, only LabelObjects and LabelGraph are discarded
    void ClearLabels();

    void InitializeAllLabelObjects(bool fromPrecalculatedGraph);

    void PrintLabelObjects(std::ostream & os) const;
    void PrintLabelObjects() const { this->PrintLabelObjects(std::cerr); };

    const LabelType& GetPixel(const IndexType &idx) const;
    void SetPixel(const IndexType & idx, const LabelType &label);
    void AddPixel(const IndexType & idx, const LabelType &label);
    void RemovePixel(const IndexType & idx, const LabelType &label);

    void LazySetLine(const IndexType &idx, const LengthType &length, const LabelType &label);
    void SetLine(const IndexType &idx, const LengthType &length, const LabelType &label);

    void SetBackgroundValue(LabelType bv) { mBackgroundValue = bv; };
    const LabelType GetBackgroundValue() const;

    void Optimize();

    /** ConstIterator
     */
    class ConstIterator
    {
    public:
        ConstIterator() {}

        ConstIterator(const Self *lm)
        {
            m_Begin = lm->mLabelObjectContainer.begin();
            m_End = lm->mLabelObjectContainer.end();
            m_Iterator = m_Begin;
        }

        ConstIterator(const ConstIterator &iter)
        {
            m_Iterator = iter.m_Iterator;
            m_Begin = iter.m_Begin;
            m_End = iter.m_End;
        }

        ConstIterator& operator=(const ConstIterator &iter)
        {
            m_Iterator = iter.m_Iterator;
            m_Begin = iter.m_Begin;
            m_End = iter.m_End;
            return *this;
        }

        const LabelObjectType* GetLabelObject() const
        {
            return m_Iterator->second;
        }

        const LabelType& GetLabel() const
        {
            return m_Iterator->first;
        }

        ConstIterator operator++(int)
        {
            ConstIterator tmp = *this;
            ++(*this);
            return tmp;
        }

        ConstIterator& operator++()
        {
            ++m_Iterator;
            return *this;
        }

        bool operator==(const ConstIterator & iter) const
        {
            return m_Iterator == iter.m_Iterator && m_Begin == iter.m_Begin && m_End == iter.m_End;
        }

        bool operator!=(const ConstIterator & iter) const
        {
            return !( *this == iter );
        }

        void GoToBegin()
        {
            m_Iterator = m_Begin;
        }

        bool IsAtEnd() const
        {
            return m_Iterator == m_End;
        }

    private:
        typedef typename std::map<LabelType, LabelObjectPointerType>::const_iterator InternalIteratorType;
        InternalIteratorType m_Iterator;
        InternalIteratorType m_Begin;
        InternalIteratorType m_End;
    };

    /** Iterator
     */
    class Iterator
    {
    public:

        Iterator() {}

        Iterator(Self *lm)
        {
            m_Begin = lm->mLabelObjectContainer.begin();
            m_End = lm->mLabelObjectContainer.end();
            m_Iterator = m_Begin;
        }

        Iterator(const Iterator &iter)
        {
            m_Iterator = iter.m_Iterator;
            m_Begin = iter.m_Begin;
            m_End = iter.m_End;
        }

        Iterator& operator=(const Iterator &iter)
        {
            m_Iterator = iter.m_Iterator;
            m_Begin = iter.m_Begin;
            m_End = iter.m_End;
            return *this;
        }

        LabelObjectType* GetLabelObject()
        {
            return m_Iterator->second;
        }

        const LabelType& GetLabel() const
        {
            return m_Iterator->first;
        }

        Iterator operator++(int)
        {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        Iterator& operator++()
        {
            ++m_Iterator;
            return *this;
        }

        bool operator==(const Iterator &iter) const
        {
            return m_Iterator == iter.m_Iterator && m_Begin == iter.m_Begin && m_End == iter.m_End;
        }

        bool operator!=(const Iterator &iter) const
        {
            return !( *this == iter );
        }

        void GoToBegin()
        {
            m_Iterator = m_Begin;
        }

        bool IsAtEnd() const
        {
            return m_Iterator == m_End;
        }

    private:
        typedef typename std::map<LabelType, LabelObjectPointerType>::iterator InternalIteratorType;
        InternalIteratorType m_Iterator;
        InternalIteratorType m_Begin;
        InternalIteratorType m_End;

        friend class LabelMapGraph;
    };

protected:
  LabelMapGraph();
  virtual ~LabelMapGraph() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
    LabelMapGraph(const Self &);        //purposely not implemented
    void operator=(const Self &);       //purposely not implemented

    double ComputeSimilarity(double a, double b);
    void ComputeRedSimilarity(vtkIdType edge, vtkIdType source, vtkIdType target);
    void ComputeGreenSimilarity(vtkIdType edge, vtkIdType source, vtkIdType target);
    void ComputeBlueSimilarity(vtkIdType edge, vtkIdType source, vtkIdType target);

    void BuildLabelMapGraphAnnotations();
    void PreparationsForInitialization();
    void InitializeUnaryLabelObjectFeatures(LabelObjectPointerType &lo);
    void InitializePrecalculatedUnaryLabelObjectFeatures(LabelObjectPointerType &lo);
    void InitializeBinaryLabelObjectFeatures(LabelObjectPointerType &lo);
    void UpdateUnaryLabelObjectFeatures(LabelObjectPointerType &lo1, LabelObjectPointerType &lo2);      //update lo1 features preceding the merging of lo2 into lo1
    void UpdateBinaryLabelObjectFeatures(LabelObjectPointerType &lo);

    /** the LabelObject container type */
    typedef std::map<LabelType, LabelObjectPointerType>         LabelObjectContainerType;
    typedef typename LabelObjectContainerType::iterator         LabelObjectContainerIterator;
    typedef typename LabelObjectContainerType::const_iterator   LabelObjectContainerConstIterator;

    typedef itk::RGBToLuminanceImageFilter<OriginalImageType, GrayscaleImageType>   RGBToGrayscaleFilterType;
    typedef itk::MinimumMaximumImageCalculator<GrayscaleImageType>                  ImageCalculatorFilterType;

    typedef itk::LabelGeometryImageFilter<LabelImageType, GrayscaleImageType>       LabelGeometryFilterType;

    typedef itk::Statistics::MaskedImageToHistogramFilter<OriginalImageType, LabelImageType>    HistogramFilterType;
    typedef typename HistogramFilterType::HistogramMeasurementVectorType                        HistogramMeasurementVectorType;
    typedef typename HistogramFilterType::HistogramSizeType                                     HistogramSizeType;
    typedef typename HistogramFilterType::HistogramType                                         HistogramType;
    typedef typename HistogramType::FrequencyContainerType                                      FrequencyContainerType;
    typedef typename FrequencyContainerType::AbsoluteFrequencyType                              FrequencyType;

    typedef itk::Statistics::ScalarImageToTextureFeaturesFilter<GrayscaleImageType> TextureFeaturesFilterType;
    typedef typename TextureFeaturesFilterType::FeatureValueVector                  FeatureVectorType;


    void AddPixel(const LabelObjectContainerIterator& it, const IndexType& idx, const LabelType& iLabel);
    void RemovePixel(const LabelObjectContainerIterator& it, const IndexType& idx, bool iEmitModifiedEvent);

    OriginalImagePointer                        mOriginalImage;
    GrayscaleImagePointer                       mOriginalGrayscaleImage;
    LabelImagePointer                           mLabelImage;

    vtkSmartPointer<vtkMutableUndirectedGraph>  mObjectGraph;
    LabelObjectContainerType                    mLabelObjectContainer;

    LabelType                                   mBackgroundValue;
    typename GrayscaleImageType::PixelType      mGrayscaleMinPixelValue;
    typename GrayscaleImageType::PixelType      mGrayscaleMaxPixelValue;

    double  mVoxelVolume;
};

}   //end namespace

#include "LabelMapGraph.tpp"

#endif /* LABELMAPGRAPH_H_ */
