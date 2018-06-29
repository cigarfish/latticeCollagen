///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelMapGraph.tpp                                                    //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-09-04                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "LabelMapGraph.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageBase.h"
#include "itkWeakPointer.h"

#include <vtkFloatArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkOutEdgeIterator.h>
#include <vtkPoints.h>
#include <vtkUnsignedLongArray.h>

#include "GraphAnnotationHelper.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#define DEBUG_MODE true


namespace itk
{

template< class TOriginalImage, class TLabelImage, class TLabelObject > const std::string LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::UnaryFeatureName[] = {
        "Mean red intensity",
        "Mean green intensity",
        "Mean blue intensity",
        "Red Histogram",
        "Green Histogram",
        "Blue Histogram",
        "RGB Histogram",
        "Texture features",
        "Number of pixels",
        "" // be sure this empty string is the last entry
};


template< class TOriginalImage, class TLabelImage, class TLabelObject > const std::string LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::BinaryFeatureName[] = {
        "Mean red similarity",
        "Mean green similarity",
        "Mean blue similarity",
        "" // be sure this empty string is the last entry
};


template< class TOriginalImage, class TLabelImage, class TLabelObject > int LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::UnaryFeatureComponents[] = {
        1,
        1,
        1,
        NUMBER_HISTOGRAM_BINS,
        NUMBER_HISTOGRAM_BINS,
        NUMBER_HISTOGRAM_BINS,
        NUMBER_HISTOGRAM_BINS*NUMBER_HISTOGRAM_BINS*NUMBER_HISTOGRAM_BINS,
        6,
        1
};


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::AttachUnaryFeature(UnaryFeatures f)
{
    if(mUnaryFeatureState[f] != IsActive) {
        mActiveUnaryFeatures++;
        mActiveUnaryFeatureComponents += UnaryFeatureComponents[f];
        mUnaryFeatureState[f] = IsActive;
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::ActivateUnaryFeature(UnaryFeatures f)
{
    if(mUnaryFeatureState[f] == IsNotAttached) {
        std::cout << "Error: failed to activate unary feature " << UnaryFeatureName[f] << std::endl;
        return;                                                             //can't activate features that are not attached to the graph
    }

    if(mUnaryFeatureState[f] == IsInactive) {
        mActiveUnaryFeatures++;
        mActiveUnaryFeatureComponents += UnaryFeatureComponents[f];
        mUnaryFeatureState[f] = IsActive;
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::DeactivateUnaryFeature(UnaryFeatures f)
{
    if(mUnaryFeatureState[f] == IsNotAttached)
        return;                                                             //can't deactivate features that are not attached to the graph

    if(mUnaryFeatureState[f] == IsActive) {
        mActiveUnaryFeatures--;
        mActiveUnaryFeatureComponents -= UnaryFeatureComponents[f];
        mUnaryFeatureState[f] = IsInactive;
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::AttachBinaryFeature(BinaryFeatures f)
{
    if(mBinaryFeatureState[f] != IsActive) {
        mActiveBinaryFeatures++;
        mBinaryFeatureState[f] = IsActive;
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::ActivateBinaryFeature(BinaryFeatures f)
{
    if(mBinaryFeatureState[f] == IsNotAttached) {
        std::cout << "Error: failed to activate binary feature " << BinaryFeatureName[f] << std::endl;
        return;                                                             //can't activate features that are not attached to the graph
    }

    if(mBinaryFeatureState[f] == IsInactive) {
        mActiveBinaryFeatures++;
        mBinaryFeatureState[f] = IsActive;
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::DeactivateBinaryFeature(BinaryFeatures f)
{
    if(mBinaryFeatureState[f] == IsNotAttached)
        return;                                                             //can't deactivate features that are not attached to the graph

    if(mBinaryFeatureState[f] == IsActive) {
        mActiveBinaryFeatures--;
        mBinaryFeatureState[f] = IsInactive;
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LabelMapGraph()
{
    for(int i=0; i<NumberUnaryFeatures; i++)
        mUnaryFeatureState[i] = IsNotAttached;

    mActiveUnaryFeatures = 0;
    mActiveUnaryFeatureComponents = 0;

    for(int i=0; i<NumberBinaryFeatures; i++)
        mBinaryFeatureState[i] = IsNotAttached;

    mActiveBinaryFeatures = 0;

    mVoxelVolume = 0.;
    mBackgroundValue = NumericTraits< LabelType >::Zero;

    mObjectGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    mObjectGraph->SetPoints(points);
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIds = vtkSmartPointer<vtkUnsignedLongArray>::New();
    pedigreeIds->Initialize();
    pedigreeIds->SetName("Pedigree IDs");
    mObjectGraph->GetVertexData()->SetPedigreeIds(pedigreeIds);

    this->Initialize();         //calls clear labels, nothing more
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);

    os << indent << "BackgroundValue: " << static_cast< typename NumericTraits< LabelType >::PrintType >( mBackgroundValue ) << std::endl;
    os << indent << "LabelObjectContainer: " << &mLabelObjectContainer << std::endl;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > bool LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GraphHasFeature(vtkSmartPointer<vtkUndirectedGraph> graph, UnaryFeatures f)
{
    std::stringstream featureName;
    featureName << UnaryFeatureName[f];
    if(UnaryFeatureComponents[f] > 1)
        featureName << 0;

    vtkDataSetAttributes* vd = graph->GetVertexData();
    for(int i = 0; i < vd->GetNumberOfArrays(); ++i) {
        if(vd->GetArrayName(i) && vd->GetArrayName(i) == '\0')
            continue;
        else {
            std::string arrayName = vd->GetArrayName(i);

            if(arrayName.compare(featureName.str()) == 0) {
                return true;
            }
        }
    }
    return false;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > bool LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GraphHasFeature(vtkSmartPointer<vtkUndirectedGraph> graph, BinaryFeatures f)
{
    vtkDataSetAttributes* ed = graph->GetEdgeData();
    for(int i=0; i<ed->GetNumberOfArrays(); ++i) {
        if(ed->GetArrayName(i) && ed->GetArrayName(i) == '\0')
            continue;
        else {
            std::string arrayName = ed->GetArrayName(i);

            if(arrayName.compare(BinaryFeatureName[f]) == 0) {
                return true;
            }
        }
    }
    return false;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::SetGraph(vtkSmartPointer<vtkMutableUndirectedGraph>  graph)
{
    mObjectGraph = graph;

    for(int i=0; i<NumberUnaryFeatures; i++)
        if(Self::GraphHasFeature(mObjectGraph, static_cast<UnaryFeatures>(i)))
            mUnaryFeatureState[i] = IsInactive;

    for(int j=0; j<NumberBinaryFeatures; j++)
        if(Self::GraphHasFeature(mObjectGraph, static_cast<BinaryFeatures>(j)))
            mBinaryFeatureState[j] = IsInactive;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::Initialize()
{
    //TODO: Algorithm to construct LabelObjects & Graph based on LabelImage

    this->ClearLabels();
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::Allocate()
{
    this->Initialize();
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::Graft(const DataObject *data)
{
    // call the superclass' implementation
    Superclass::Graft(data);

    if(data) {
        // Attempt to cast data to an Image
        const Self *imgData;

        try
        {
            imgData = dynamic_cast< const Self * >( data );
        }
        catch ( ... )
        {
            return;
        }

        if(imgData)
        {
            // Now copy anything remaining that is needed
            mLabelObjectContainer   = imgData->mLabelObjectContainer;
            mLabelImage             = imgData->mLabelImage;
            mOriginalImage          = imgData->mOriginalImage;
            mObjectGraph            = imgData->mObjectGraph;
            mBackgroundValue        = imgData->mBackgroundValue;
        }
        else
        {
            // pointer could not be cast back down
            itkExceptionMacro( << "itk::Image::Graft() cannot cast " << typeid( data ).name() << " to " << typeid( const Self * ).name() );
        }
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > typename LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LabelObjectType* LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GetLabelObject(const LabelType & label)
{
    if(mBackgroundValue == label)
    {
        itkExceptionMacro(<< "Label " << static_cast< typename NumericTraits< LabelType >::PrintType >(label) << " is the background label.");
    }
    LabelObjectContainerIterator it = mLabelObjectContainer.find(label);
    if(it == mLabelObjectContainer.end())
    {
        itkExceptionMacro(<< "No label object with label " << static_cast< typename NumericTraits< LabelType >::PrintType >(label) << ".");
    }

    return it->second;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > const typename LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LabelObjectType* LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GetLabelObject(const LabelType & label) const
{
    if(mBackgroundValue == label)
    {
        itkExceptionMacro(<< "Label " << static_cast< typename NumericTraits< LabelType >::PrintType >(label) << " is the background label.");
    }
    LabelObjectContainerConstIterator it = mLabelObjectContainer.find(label);
    if(it == mLabelObjectContainer.end())
    {
        itkExceptionMacro(<< "No label object with label " << static_cast< typename NumericTraits< LabelType >::PrintType >(label) << ".");
    }

    return it->second;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > bool LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::HasLabel(const LabelType label) const
{
    return mLabelObjectContainer.find(label) != mLabelObjectContainer.end();
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > const typename LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LabelType& LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GetPixel(const IndexType & idx) const
{
    return mLabelImage->GetPixel(idx);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > const typename LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LabelType LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GetBackgroundValue() const
{
    return mBackgroundValue;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > typename LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LabelObjectType* LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GetNthLabelObject(const SizeValueType & pos)
{
    SizeValueType i = 0;

    for(LabelObjectContainerIterator it = mLabelObjectContainer.begin(); it != mLabelObjectContainer.end(); it++)
    {
        if(i == pos)
        {
            return it->second;
        }
        i++;
    }
    itkExceptionMacro(<< "Can't access to label object at position " << pos << ". The label map has only " << this->GetNumberOfLabelObjects() << " label objects registered.");
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > const typename LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LabelObjectType* LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GetNthLabelObject(const SizeValueType & pos) const
{
    SizeValueType i = 0;

    for(LabelObjectContainerConstIterator it = mLabelObjectContainer.begin(); it != mLabelObjectContainer.end(); it++)
    {
        if(i == pos)
        {
            return it->second;
        }
        i++;
    }
    itkExceptionMacro(<< "Can't access to label object at position " << pos << ". The label map has only " << this->GetNumberOfLabelObjects() << " label objects registered.");
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::AddLabelVertex(const IndexType &idx, const LabelType &label)
{
    if(label != mBackgroundValue && mObjectGraph->FindVertex(label) == -1) {
        vtkSmartPointer<vtkPoints> points = mObjectGraph->GetPoints();

        vtkIdType v = mObjectGraph->AddVertex(label);
        mObjectGraph->GetVertexData()->GetPedigreeIds()->DataChanged();     //for vtk5.10 compatibility

        points->InsertPoint(v, idx[0], idx[1], idx[2]);
        mObjectGraph->SetPoints(points);
        mObjectGraph->Modified();
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::AddLabelEdge(const IndexType &idx, const LabelType &label)
{
    // Define offsets
    int numOffInDimMode = 4;
    if(ImageDimension == 3) numOffInDimMode = 6;
    OffsetType *off;
    if(ImageDimension == 2) {
        off = new OffsetType[4];        // face connected
        off[0][0] = -1; off[0][1] =  0; // 0 - NM
        off[1][0] =  1, off[1][1] =  0; // 1 - SM
        off[2][0] =  0; off[2][1] =  1; // 2 - MS
        off[3][0] =  0, off[3][1] = -1; // 3 - MN
    }
    else if(ImageDimension == 3) {
        off = new OffsetType[6];            // face connected
        off[0][0] = -1; off[0][1] =  0; off[0][2] =  0;    // 0 - NMM
        off[1][0] =  1; off[1][1] =  0; off[1][2] =  0;    // 1 - SMM
        off[2][0] =  0; off[2][1] =  1; off[2][2] =  0;    // 2 - MSM
        off[3][0] =  0; off[3][1] = -1; off[3][2] =  0;    // 3 - MNM
        off[4][0] =  0; off[4][1] =  0; off[4][2] = -1;    // 4 - MMN
        off[5][0] =  0; off[5][1] =  0; off[5][2] =  1;    // 5 - MMS
    }

    SizeType radius;
    radius.Fill(1);

    for(unsigned int i=0; i<numOffInDimMode; i++) {
        IndexType neiIdx = idx+off[i];

        if(mLabelImage->GetLargestPossibleRegion().IsInside(neiIdx) && mLabelImage->GetPixel(neiIdx) != label && mLabelImage->GetPixel(neiIdx) != mBackgroundValue) {
            vtkIdType v1 = mObjectGraph->FindVertex(label);
            vtkIdType v2 = mObjectGraph->FindVertex(mLabelImage->GetPixel(neiIdx));

            if(v2 != -1 && mObjectGraph->GetEdgeId(v1, v2) == -1) {
                mObjectGraph->AddEdge(v1, v2);
                mObjectGraph->Modified();
            }
        }
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::RemoveLabelVertex(const LabelType &label)
{
    mObjectGraph->RemoveVertex(mObjectGraph->FindVertex(label));
    mObjectGraph->Squeeze();
    mObjectGraph->Modified();
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::SetPixel(const IndexType &idx, const LabelType &iLabel)
{
    bool newLabel = true; // or can be initialized by ( iLabel == m_BackgroundValue )

    LabelObjectContainerIterator it = mLabelObjectContainer.begin();

    while(it != mLabelObjectContainer.end())
    {
        // increment the iterator before removing the pixel because
        // RemovePixel() can remove the object and thus invalidate the
        // iterator
        if(it->first != iLabel)
        {
            LabelObjectContainerIterator tempIt = it;
            ++it;
            bool emitModifiedEvent = (iLabel == mBackgroundValue);
            this->RemovePixel(tempIt, idx, emitModifiedEvent);
        }
        else
        {
            newLabel = false;
            this->AddPixel(it, idx, iLabel);
            ++it;
        }
    }
    if(newLabel)
    {
        this->AddPixel(mLabelObjectContainer.end(), idx, iLabel);
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::AddPixel(const IndexType &idx, const LabelType &label)
{
    if(label == mBackgroundValue)
    {
        // just do nothing
        return;
    }

    LabelObjectContainerIterator it = mLabelObjectContainer.find(label);
    this->AddPixel(it, idx, label);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::AddPixel(const LabelObjectContainerIterator &it, const IndexType &idx, const LabelType &label)
{
    if(label == mBackgroundValue)
    {
        // just do nothing
        return;
    }

    if(it != mLabelObjectContainer.end())
    {
        // the label already exist - add the pixel to it
        ( *it ).second->AddIndex(idx);
        this->Modified();

        mLabelImage->SetPixel(idx, label);
        this->AddLabelEdge(idx, label);
    }
    else
    {
        // the label does not exist yet - create a new one
        LabelObjectPointerType labelObject = LabelObjectType::New();
        labelObject->SetLabel(label);
        labelObject->AddIndex(idx);
        // Modified() is called in AddLabelObject()
        this->AddLabelObject(labelObject);
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::RemovePixel(const LabelObjectContainerIterator &it, const IndexType &idx, bool iEmitModifiedEvent)
{
    if(it != mLabelObjectContainer.end())
    {
        // the label already exist - add the pixel to it
        if(it->second->RemoveIndex(idx))
        {
            if(it->second->Empty())
            {
                this->RemoveLabelObject(it->second);
            }
            if(iEmitModifiedEvent)
            {
                this->Modified();
            }
            mLabelImage->SetPixel(idx, mBackgroundValue);
        }
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::RemovePixel(const IndexType &idx, const LabelType &label)
{
    if(label == mBackgroundValue)
    {
        // just do nothing
        return;
    }

    LabelObjectContainerIterator it = mLabelObjectContainer.find(label);

    bool emitModifiedEvent = true;
    RemovePixel(it, idx, emitModifiedEvent);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::SetLine(const IndexType &idx, const LengthType &length, const LabelType &label)
{
    if(label == mBackgroundValue)
    {
        // just do nothing
        return;
    }

    LabelObjectContainerIterator it = mLabelObjectContainer.find(label);

    if(it != mLabelObjectContainer.end())
    {
        // the label already exist - add the pixel to it
        (*it).second->AddLine(idx, length);
        this->Modified();

        for(unsigned int pixelId = 0; pixelId < length; pixelId++) {
            OffsetType offLine;
            offLine[0] = pixelId;
            offLine[1] = 0;
            if(ImageDimension==3)
                offLine[2] = 0;
            IndexType lidx = idx + offLine;
            mLabelImage->SetPixel(lidx, (*it).second->GetLabel());

            this->AddLabelEdge(lidx, (*it).second->GetLabel());
        }
    }
    else
    {
        // the label does not exist yet - create a new one
        LabelObjectPointerType labelObject = LabelObjectType::New();
        labelObject->SetLabel(label);
        labelObject->AddLine(idx, length);
        // Modified() is called in AddLabelObject()
        this->AddLabelObject(labelObject);
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LazySetLine(const IndexType &idx, const LengthType &length, const LabelType &label)
{
    if(label == mBackgroundValue)
    {
        // just do nothing
        return;
    }

    LabelObjectContainerIterator it = mLabelObjectContainer.find(label);

    if(it != mLabelObjectContainer.end())
    {
        // the label already exist - add the pixel to it
        (*it).second->AddLine(idx, length);
        this->Modified();

        for(unsigned int pixelId = 0; pixelId < length; pixelId++) {
            OffsetType off;
            off[0] = pixelId;
            off[1] = 0;
            if(ImageDimension==3)
                off[2] = 0;
            IndexType lidx = idx + off;
            mLabelImage->SetPixel(lidx, (*it).second->GetLabel());
        }
    }
    else
    {
        // the label does not exist yet - create a new one
        LabelObjectPointerType labelObject = LabelObjectType::New();
        labelObject->SetLabel(label);
        labelObject->AddLine(idx, length);
        // Modified() is called in AddLabelObject()
        this->AddLabelObject(labelObject);
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > typename LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LabelObjectType* LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GetLabelObject(const IndexType &idx) const
{
    for(LabelObjectContainerConstIterator it = mLabelObjectContainer.begin(); it != mLabelObjectContainer.end(); it++)
    {
        if(it->second->HasIndex(idx))
        {
            return it->second.GetPointer();
        }
    }
    itkExceptionMacro(<< "No label object at index " << idx << ".");
    //   return NULL;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::AddLabelObject(LabelObjectType *labelObject)
{
    itkAssertOrThrowMacro((labelObject != NULL), "Input LabelObject can't be Null");

    mLabelObjectContainer[labelObject->GetLabel()] = labelObject;
    this->Modified();

    for(unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++)
        mLabelImage->SetPixel(labelObject->GetIndex(pixelId), labelObject->GetLabel());

    this->AddLabelVertex(labelObject->GetIndex(0), labelObject->GetLabel());
    for(unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++)
        this->AddLabelEdge(labelObject->GetIndex(pixelId), labelObject->GetLabel());
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::PushLabelObject(LabelObjectType *labelObject)
{
    itkAssertOrThrowMacro((labelObject != NULL), "Input LabelObject can't be Null");

    if(mLabelObjectContainer.empty())
    {
        if(mBackgroundValue == 0)
        {
            labelObject->SetLabel(1);
        }
        else
        {
            labelObject->SetLabel(0);
        }
    }
    else
    {
        LabelType lastLabel = mLabelObjectContainer.rbegin()->first;
        LabelType firstLabel = mLabelObjectContainer.begin()->first;
        if(lastLabel != NumericTraits< LabelType >::max() && lastLabel + 1 != mBackgroundValue)
        {
            labelObject->SetLabel(lastLabel + 1);
        }
        else if(lastLabel != NumericTraits< LabelType >::max() && lastLabel + 1 != NumericTraits< LabelType >::max() && lastLabel + 2 != mBackgroundValue)
        {
            labelObject->SetLabel(lastLabel + 2);
        }
        else if(firstLabel != NumericTraits< LabelType >::NonpositiveMin() && firstLabel - 1 != mBackgroundValue)
        {
            labelObject->SetLabel(firstLabel - 1);
        }
        else
        {
            // search for an unused label
            LabelType label = firstLabel;
            LabelObjectContainerConstIterator it;
            for(it = mLabelObjectContainer.begin(); it != mLabelObjectContainer.end(); it++, label++)
            {
                assert( (it->second.IsNotNull()) );
                if(label == mBackgroundValue)
                {
                    label++;
                }
                if(label != it->first)
                {
                    labelObject->SetLabel(label);
                    break;
                }
            }
            if(label == lastLabel)
            {
                itkExceptionMacro(<< "Can't push the label object: the label map is full.");
            }
        }
    }
    // modified is called in AddLabelObject()
    this->AddLabelObject(labelObject);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::RemoveLabelObject(LabelObjectType *labelObject)
{
    itkAssertOrThrowMacro( (labelObject != NULL), "Input LabelObject can't be Null" );
    // modified is called in RemoveLabel()
    this->RemoveLabel(labelObject->GetLabel());
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::RemoveLabel(const LabelType &label)
{
    if(mBackgroundValue == label)
    {
        itkExceptionMacro(<< "Label " << static_cast< typename NumericTraits< LabelType >::PrintType >(label) << " is the background label.");
    }

    LabelObjectContainerIterator it = mLabelObjectContainer.find(label);

    for(unsigned int pixelId = 0; pixelId < (*it).second->Size(); pixelId++)
        mLabelImage->SetPixel((*it).second->GetIndex(pixelId), mBackgroundValue);

    mLabelObjectContainer.erase(label);
    this->Modified();

    this->RemoveLabelVertex(label);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::RemoveLabels(std::set<vtkIdType> labels)
{
    for(std::set<vtkIdType>::iterator it=labels.begin(); it!=labels.end(); ++it)
        RemoveLabel(*it);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::KeepLabelsRemoveRest(std::set<vtkIdType> labelsToKeep)
{
    std::set<vtkIdType> labelsToRemove;

    for(LabelObjectContainerIterator it=mLabelObjectContainer.begin(); it!=mLabelObjectContainer.end(); ++it)
        if(labelsToKeep.count(it->first)==0 && it->first != mBackgroundValue)
            labelsToRemove.insert(it->first);

    RemoveLabels(labelsToRemove);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::ClearLabels()
{
    if(!mLabelObjectContainer.empty())
    {
        mLabelObjectContainer.clear();
        this->Modified();

        mObjectGraph->Initialize();
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > typename LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LabelVectorType LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GetLabels() const
{
    LabelVectorType res;

    res.reserve(this->GetNumberOfLabelObjects());
    for(LabelObjectContainerConstIterator it = mLabelObjectContainer.begin(); it != mLabelObjectContainer.end(); it++)
    {
        res.push_back(it->first);
    }
    return res;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > typename LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LabelObjectVectorType LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::GetLabelObjects() const
{
    LabelObjectVectorType res;

    res.reserve(this->GetNumberOfLabelObjects());
    for(LabelObjectContainerConstIterator it = mLabelObjectContainer.begin(); it != mLabelObjectContainer.end(); it++)
    {
        res.push_back(it->second);
    }
    return res;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::PrintLabelObjects(std::ostream &os) const
{
    for(LabelObjectContainerConstIterator it = mLabelObjectContainer.begin(); it != mLabelObjectContainer.end(); it++)
    {
        assert( ( it->second.IsNotNull() ) );
        it->second->Print(os);
        os << std::endl;
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::Optimize()
{
    for(LabelObjectContainerConstIterator it = mLabelObjectContainer.begin(); it != mLabelObjectContainer.end(); it++ )
    {
        assert( ( it->second.IsNotNull() ) );
        it->second->Optimize();
    }
    this->Modified();
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > vtkIdType LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::MergeLabelObjects(vtkIdType o1, vtkIdType o2)
{
    return MergeLabelObjects(mLabelObjectContainer[o1], mLabelObjectContainer[o2]);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > vtkIdType LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LazyMergeLabelObjects(vtkIdType o1, vtkIdType o2)
{
    return LazyMergeLabelObjects(mLabelObjectContainer[o1], mLabelObjectContainer[o2]);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > vtkIdType LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::MergeLabelObjects(LabelObjectPointerType o1, LabelObjectPointerType o2)
{
    if(o1->Size() < o2->Size()) {
        LabelObjectPointerType temp = o1;
        o1 = o2;
        o2 = temp;
    }
    vtkIdType removedLabel = o2->GetLabel();

    for(unsigned int i=0; i<o2->GetNumberOfLines(); i++)
        this->SetLine(o2->GetLine(i).GetIndex(), o2->GetLine(i).GetLength(), o1->GetLabel());
    o1->Optimize();

    this->UpdateUnaryLabelObjectFeatures(o1, o2);
    this->UpdateBinaryLabelObjectFeatures(o1);

    o2->Clear();
    this->RemoveLabelObject(o2);

    return removedLabel;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > vtkIdType LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LazyMergeLabelObjects(LabelObjectPointerType o1, LabelObjectPointerType o2)
{
    if(o1->Size() < o2->Size()) {
        LabelObjectPointerType temp = o1;
        o1 = o2;
        o2 = temp;
    }
    vtkIdType removedLabel = o2->GetLabel();

    for(unsigned int i=0; i<o2->GetNumberOfLines(); i++)
        this->LazySetLine(o2->GetLine(i).GetIndex(), o2->GetLine(i).GetLength(), o1->GetLabel());
    o1->Optimize();

    vtkIdType sourceO1 = mObjectGraph->FindVertex(o1->GetLabel());
    vtkIdType sourceO2 = mObjectGraph->FindVertex(o2->GetLabel());
    vtkSmartPointer<vtkOutEdgeIterator> it = vtkSmartPointer<vtkOutEdgeIterator>::New();
    mObjectGraph->GetOutEdges(sourceO2, it);

    while(it->HasNext()) {
        vtkOutEdgeType e = it->Next();
        vtkIdType targetO2 = e.Target;

        if(targetO2 != sourceO1 && mObjectGraph->GetEdgeId(sourceO1, targetO2) == -1) {
            mObjectGraph->AddEdge(sourceO1, targetO2);
            mObjectGraph->Modified();
        }
    }

    this->UpdateUnaryLabelObjectFeatures(o1, o2);
    this->UpdateBinaryLabelObjectFeatures(o1);

    o2->Clear();
    this->RemoveLabelObject(o2);

    return removedLabel;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > int LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::LazyCollapseLabelObjects(std::set<vtkIdType> &labelsToCollapse)
{
    int mergeOpsOverall = 0;
    int mergeOpsLastIteration = 1;

    std::set<vtkIdType> isolatedLabels = labelsToCollapse;

    while(mergeOpsLastIteration > 0) {
        mergeOpsLastIteration = 0;

        std::map<vtkIdType, vtkIdType> labelPairsToCollapse;
        std::set<vtkIdType> labelsUsedThisIteration;
        std::set<vtkIdType> labelsDeletedThisIteration;

        for(std::set<vtkIdType>::iterator itOuter=labelsToCollapse.begin(); itOuter!=labelsToCollapse.end(); ++itOuter) {
            if(labelsUsedThisIteration.count(*itOuter)!=0)
                continue;

            std::set<vtkIdType>::iterator itInner=itOuter;
            itInner++;
            for(itInner; itInner!=labelsToCollapse.end(); ++itInner) {
                if(labelsUsedThisIteration.count(*itOuter)==0 && labelsUsedThisIteration.count(*itInner)==0 && mObjectGraph->GetEdgeId(mObjectGraph->FindVertex(*itOuter), mObjectGraph->FindVertex(*itInner))!=-1) {
                    labelsUsedThisIteration.insert(*itOuter);
                    labelsUsedThisIteration.insert(*itInner);
                    isolatedLabels.erase(*itOuter);
                    isolatedLabels.erase(*itInner);
                    labelPairsToCollapse.insert(std::pair<vtkIdType, vtkIdType>(*itOuter, *itInner));
                    break;
                }
            }
        }
        std::cout << "collected all pairs: " << labelPairsToCollapse.size() << std::endl;

        for(std::map<vtkIdType, vtkIdType>::iterator it=labelPairsToCollapse.begin(); it!=labelPairsToCollapse.end(); ++it) {
            labelsDeletedThisIteration.insert( LazyMergeLabelObjects(it->first, it->second) );
            mergeOpsLastIteration++;
        }
        mergeOpsOverall += mergeOpsLastIteration;

        for(std::set<vtkIdType>::iterator it=labelsDeletedThisIteration.begin(); it!=labelsDeletedThisIteration.end(); ++it) {
            labelsToCollapse.erase(*it);
            isolatedLabels.erase(*it);
        }

        for(std::set<vtkIdType>::iterator it=isolatedLabels.begin(); it!=isolatedLabels.end(); ++it)
            labelsToCollapse.erase(*it);

        std::cout << "finished all merges: " << mergeOpsLastIteration << std::endl;
    }
    labelsToCollapse = isolatedLabels;

    return mergeOpsOverall;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > double LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::ComputeSimilarity(double a, double b)
{
    double similarity = 1. - (abs(a - b) / (double)NumericTraits<OriginalImagePixelComponentType>::max());
    return similarity;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::ComputeRedSimilarity(vtkIdType edge, vtkIdType source, vtkIdType target)
{
    double similarity = ComputeSimilarity(mLabelObjectContainer[source]->GetRedMean(), mLabelObjectContainer[target]->GetRedMean());
    vtkFloatArray::SafeDownCast(mObjectGraph->GetEdgeData()->GetArray("Mean red similarity"))->InsertValue(edge, similarity);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::ComputeGreenSimilarity(vtkIdType edge, vtkIdType source, vtkIdType target)
{
    double similarity = ComputeSimilarity(mLabelObjectContainer[source]->GetGreenMean(), mLabelObjectContainer[target]->GetGreenMean());
    vtkFloatArray::SafeDownCast(mObjectGraph->GetEdgeData()->GetArray("Mean green similarity"))->InsertValue(edge, similarity);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::ComputeBlueSimilarity(vtkIdType edge, vtkIdType source, vtkIdType target)
{
    double similarity = ComputeSimilarity(mLabelObjectContainer[source]->GetBlueMean(), mLabelObjectContainer[target]->GetBlueMean());
    vtkFloatArray::SafeDownCast(mObjectGraph->GetEdgeData()->GetArray("Mean blue similarity"))->InsertValue(edge, similarity);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::UpdateUnaryLabelObjectFeatures(LabelObjectPointerType &lo1, LabelObjectPointerType &lo2)
{
//    std::cout << "UpdateUnaryLabelObjectFeatures: lo1 size = " << lo1->GetNumberOfPixels() << ", lo2 size = " << lo2->GetNumberOfPixels() << ", lo1 blueAcc = " << lo1->GetBlueAcc() << ", lo2 blueAcc = " << lo2->GetBlueAcc() << ", lo1 blueMean = " << lo1->GetBlueMean() << ", lo2 blueMean = " << lo2->GetBlueMean() << std::endl;

    SizeValueType nbOfPixels = lo1->GetNumberOfPixels() + lo2->GetNumberOfPixels();
    double physicalSize = nbOfPixels * mVoxelVolume;

    lo1->SetNumberOfPixels(nbOfPixels);
    lo1->SetPhysicalSize(physicalSize);

    itk::Point<double, ImageDimension> coordAcc1 = lo1->GetCoordAcc();
    itk::Point<double, ImageDimension> coordAcc2 = lo2->GetCoordAcc();
    itk::FixedArray<int, itkGetStaticConstMacro(ImageDimension)*2> bbox1 = lo1->GetBoundingBox();
    itk::FixedArray<int, itkGetStaticConstMacro(ImageDimension)*2> bbox2 = lo2->GetBoundingBox();

    itk::Point<double, ImageDimension> coordAcc;
    ContinuousIndex<double, ImageDimension> centroid;
    itk::FixedArray<int, itkGetStaticConstMacro(ImageDimension)*2> bbox;

    int loc = 0;
    for(unsigned int i=0; i<ImageDimension; i++) {
        coordAcc[i] = coordAcc1[i] + coordAcc2[i];
        centroid[i] = coordAcc[i] / nbOfPixels;

        bbox[loc]   = std::min(bbox1[loc],   bbox2[loc]);
        bbox[loc+1] = std::max(bbox1[loc+1], bbox2[loc+1]);
        loc+=2;
    }

    lo1->SetCoordAcc(coordAcc);
    lo1->SetCentroid(centroid);
    lo1->SetBoundingBox(bbox);

    vtkIdType vId = mObjectGraph->FindVertex(lo1->GetLabel());
    if(ImageDimension==2)
        mObjectGraph->GetPoints()->SetPoint(vId, lo1->GetCentroid()[0], lo1->GetCentroid()[1], 0);
    else if(ImageDimension==3)
        mObjectGraph->GetPoints()->SetPoint(vId, lo1->GetCentroid()[0], lo1->GetCentroid()[1], lo1->GetCentroid()[2]);

    if(mUnaryFeatureState[MEAN_RED_INTENSITY] == IsActive) {
        double redAcc = lo1->GetRedAcc() + lo2->GetRedAcc();
        double redMean = redAcc / (double)nbOfPixels;

        lo1->SetRedAcc(redAcc);
        lo1->SetRedMean(redMean);
        vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[MEAN_RED_INTENSITY].c_str()))->InsertValue(vId, lo1->GetRedMean());
    }
    if(mUnaryFeatureState[MEAN_GREEN_INTENSITY] == IsActive) {
        double greenAcc = lo1->GetGreenAcc() + lo2->GetGreenAcc();
        double greenMean = greenAcc / (double)nbOfPixels;

        lo1->SetGreenAcc(greenAcc);
        lo1->SetGreenMean(greenMean);
        vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[MEAN_GREEN_INTENSITY].c_str()))->InsertValue(vId, lo1->GetGreenMean());
    }
    if(mUnaryFeatureState[MEAN_BLUE_INTENSITY] == IsActive) {
        double blueAcc = lo1->GetBlueAcc() + lo2->GetBlueAcc();
        double blueMean = blueAcc / (double)nbOfPixels;

        lo1->SetBlueAcc(blueAcc);
        lo1->SetBlueMean(blueMean);
        vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[MEAN_BLUE_INTENSITY].c_str()))->InsertValue(vId, lo1->GetBlueMean());
    }
    if(mUnaryFeatureState[RED_HISTOGRAM] == IsActive) {
        typename HistogramType::Pointer h1 = lo1->GetRedHistogram();
        typename HistogramType::Pointer h2 = lo2->GetRedHistogram();

        for(int i=0; i<UnaryFeatureComponents[RED_HISTOGRAM]; i++) {
            h1->IncreaseFrequency( i, h2->GetFrequency(i) );

            std::stringstream arrayName;
            arrayName << UnaryFeatureName[RED_HISTOGRAM] << i;
            vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->InsertValue(vId, (double)h1->GetFrequency(i) / (double)nbOfPixels);
        }
        lo1->SetRedHistogram(h1);
    }
    if(mUnaryFeatureState[GREEN_HISTOGRAM] == IsActive) {
        typename HistogramType::Pointer h1 = lo1->GetGreenHistogram();
        typename HistogramType::Pointer h2 = lo2->GetGreenHistogram();

        for(int i=0; i<UnaryFeatureComponents[GREEN_HISTOGRAM]; i++) {
            h1->IncreaseFrequency( i, h2->GetFrequency(i) );

            std::stringstream arrayName;
            arrayName << UnaryFeatureName[GREEN_HISTOGRAM] << i;
            vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->InsertValue(vId, (double)h1->GetFrequency(i) / (double)nbOfPixels);
        }
        lo1->SetGreenHistogram(h1);
    }
    if(mUnaryFeatureState[BLUE_HISTOGRAM] == IsActive) {
        typename HistogramType::Pointer h1 = lo1->GetBlueHistogram();
        typename HistogramType::Pointer h2 = lo2->GetBlueHistogram();

        for(int i=0; i<UnaryFeatureComponents[BLUE_HISTOGRAM]; i++) {
            h1->IncreaseFrequency( i, h2->GetFrequency(i) );

            std::stringstream arrayName;
            arrayName << UnaryFeatureName[BLUE_HISTOGRAM] << i;
            vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->InsertValue(vId, (double)h1->GetFrequency(i) / (double)nbOfPixels);
        }
        lo1->SetBlueHistogram(h1);
    }
    if(mUnaryFeatureState[RGB_HISTOGRAM] == IsActive) {
        typename HistogramType::Pointer h1 = lo1->GetRGBHistogram();
        typename HistogramType::Pointer h2 = lo2->GetRGBHistogram();

        for(int i=0; i<UnaryFeatureComponents[RGB_HISTOGRAM]; i++) {
            h1->IncreaseFrequency( i, h2->GetFrequency(i) );

            std::stringstream arrayName;
            arrayName << UnaryFeatureName[RGB_HISTOGRAM] << i;
            vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->InsertValue(vId, (double)h1->GetFrequency(i) / (double)nbOfPixels);
        }
        lo1->SetRGBHistogram(h1);
    }
    if(mUnaryFeatureState[TEXTURE_FEATURES] == IsActive) {
        typename GrayscaleImageType::RegionType region;
        typename GrayscaleImageType::SizeType size;
        typename GrayscaleImageType::IndexType index;

        int loc = 0;
        for (unsigned int dim=0; dim<ImageDimension; dim++) {
            index[dim] = (long int)lo1->GetBoundingBox()[loc];                               //bbox min
            size[dim] = (long unsigned int)lo1->GetBoundingBox()[loc+1] - index[dim] + 1;    //bbox max - min + 1
            loc += 2;
        }
        region.SetSize(size);
        region.SetIndex(index);

        mLabelImage->SetRequestedRegion(region);
        mOriginalGrayscaleImage->SetRequestedRegion(region);

        typename TextureFeaturesFilterType::Pointer textureFilter = TextureFeaturesFilterType::New();
        textureFilter->SetInput(mOriginalGrayscaleImage);
        textureFilter->SetMaskImage(mLabelImage);
        textureFilter->SetInsidePixelValue((unsigned int)lo1->GetLabel());
        textureFilter->SetPixelValueMinMax(mGrayscaleMinPixelValue, mGrayscaleMaxPixelValue);
        textureFilter->FastCalculationsOn();
        textureFilter->Update();

        typename FeatureVectorType::Pointer featureVector = textureFilter->GetFeatureMeans();
        lo1->SetTextureFeatures(featureVector);

        for(int i=0; i<UnaryFeatureComponents[TEXTURE_FEATURES]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[TEXTURE_FEATURES] << i;
            vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->InsertValue(vId, lo1->GetTextureFeatures()->ElementAt(i));
        }
    }
    if(mUnaryFeatureState[NUMBER_OF_PIXELS] == IsActive)
        vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[NUMBER_OF_PIXELS].c_str()))->InsertValue(vId, lo1->GetNumberOfPixels());

//    std::cout << "UpdateUnaryLabelObjectFeatures: updated lo1 size = " << lo1->GetNumberOfPixels() << ", updated lo1 blueAcc = " << lo1->GetBlueAcc() << ", updated lo1 blueMean = " << lo1->GetBlueMean() << std::endl;
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::UpdateBinaryLabelObjectFeatures(LabelObjectPointerType &lo)
        {
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIds = vtkUnsignedLongArray::SafeDownCast(mObjectGraph->GetVertexData()->GetPedigreeIds());

    vtkSmartPointer<vtkOutEdgeIterator> it = vtkSmartPointer<vtkOutEdgeIterator>::New();
    mObjectGraph->GetOutEdges(mObjectGraph->FindVertex(lo->GetLabel()), it);

    while(it->HasNext()) {
        vtkOutEdgeType e = it->Next();
        vtkIdType source = lo->GetLabel();
        vtkIdType target = pedigreeIds->GetValue(e.Target);

        if(mBinaryFeatureState[MEAN_RED_SIMILARITY] == IsActive)      ComputeRedSimilarity(e.Id, source, target);
        if(mBinaryFeatureState[MEAN_GREEN_SIMILARITY] == IsActive)    ComputeGreenSimilarity(e.Id, source, target);
        if(mBinaryFeatureState[MEAN_BLUE_SIMILARITY] == IsActive)     ComputeBlueSimilarity(e.Id, source, target);
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::InitializeUnaryLabelObjectFeatures(LabelObjectPointerType &lo)
{
    SizeValueType                               nbOfPixels = 0;
    double                                      physicalSize;
    ContinuousIndex< double, ImageDimension >   centroid;
    centroid.Fill(0);
    double                                      accRed = 0.;
    double                                      accGreen = 0.;
    double                                      accBlue = 0.;

    nbOfPixels = lo->Size();
    physicalSize = nbOfPixels * mVoxelVolume;

    lo->SetNumberOfPixels(nbOfPixels);
    lo->SetPhysicalSize(physicalSize);

    for(unsigned int pixelId=0; pixelId<nbOfPixels; pixelId++)
    {
        for(unsigned int j=0; j<ImageDimension; j++)
            centroid[j] += lo->GetIndex(pixelId)[j];

        if(mUnaryFeatureState[MEAN_RED_INTENSITY] == IsActive)    accRed += mOriginalImage->GetPixel(lo->GetIndex(pixelId)).GetRed();
        if(mUnaryFeatureState[MEAN_GREEN_INTENSITY] == IsActive)  accGreen += mOriginalImage->GetPixel(lo->GetIndex(pixelId)).GetGreen();
        if(mUnaryFeatureState[MEAN_BLUE_INTENSITY] == IsActive)   accBlue += mOriginalImage->GetPixel(lo->GetIndex(pixelId)).GetBlue();
    }
    lo->SetCoordAcc(centroid);

    const unsigned int MeasurementVectorSize = 3;               // RGB
    HistogramSizeType histogramSize( MeasurementVectorSize );
    histogramSize.Fill(0);

    HistogramMeasurementVectorType lowerBound(MeasurementVectorSize);
    lowerBound.Fill(itk::NumericTraits<OriginalImagePixelComponentType>::min() - 0.5);

    HistogramMeasurementVectorType upperBound(MeasurementVectorSize);
    upperBound.Fill(itk::NumericTraits<OriginalImagePixelComponentType>::max() + 0.5) ;

    typename HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
    histogramFilter->SetInput(mOriginalImage);
    histogramFilter->SetMaskImage(mLabelImage);
    histogramFilter->SetMaskValue(lo->GetLabel());
    histogramFilter->SetHistogramBinMinimum(lowerBound);
    histogramFilter->SetHistogramBinMaximum(upperBound);
    histogramFilter->SetAutoMinimumMaximum(false);
    histogramFilter->SetMarginalScale(10);                      // Required (could this be set in the filter?)

    if(mUnaryFeatureState[RED_HISTOGRAM] == IsActive) {
        histogramSize[0] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Red   channel
        histogramSize[1] = 1;                                   // number of bins for the Green channel
        histogramSize[2] = 1;                                   // number of bins for the Blue  channel

        histogramFilter->SetHistogramSize(histogramSize);
        histogramFilter->Update();

        const typename HistogramType::Pointer histogram = histogramFilter->GetOutput();
        histogram->DisconnectPipeline();

        lo->SetRedHistogram(histogram);
    }
    if(mUnaryFeatureState[GREEN_HISTOGRAM] == IsActive) {
        histogramSize[0] = 1;                                   // number of bins for the Red   channel
        histogramSize[1] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Green channel
        histogramSize[2] = 1;                                   // number of bins for the Blue  channel

        histogramFilter->SetHistogramSize(histogramSize);
        histogramFilter->Update();

        const typename HistogramType::Pointer histogram = histogramFilter->GetOutput();
        histogram->DisconnectPipeline();

        lo->SetGreenHistogram(histogram);
    }
    if(mUnaryFeatureState[BLUE_HISTOGRAM] == IsActive) {
        histogramSize[0] = 1;                                   // number of bins for the Red   channel
        histogramSize[1] = 1;                                   // number of bins for the Green channel
        histogramSize[2] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Blue  channel

        histogramFilter->SetHistogramSize(histogramSize);
        histogramFilter->Update();

        const typename HistogramType::Pointer histogram = histogramFilter->GetOutput();
        histogram->DisconnectPipeline();

        lo->SetBlueHistogram(histogram);
    }
    if(mUnaryFeatureState[RGB_HISTOGRAM] == IsActive) {
        histogramSize[0] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Red   channel
        histogramSize[1] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Green channel
        histogramSize[2] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Blue  channel

        histogramFilter->SetHistogramSize(histogramSize);
        histogramFilter->Update();

        const typename HistogramType::Pointer histogram = histogramFilter->GetOutput();
        histogram->DisconnectPipeline();

        lo->SetRGBHistogram(histogram);
    }
    if(mUnaryFeatureState[TEXTURE_FEATURES] == IsActive) {
        typename GrayscaleImageType::RegionType region;
        typename GrayscaleImageType::SizeType size;
        typename GrayscaleImageType::IndexType index;

        int loc = 0;
        for (unsigned int dim=0; dim<ImageDimension; dim++) {
            index[dim] = (long int)lo->GetBoundingBox()[loc];                               //bbox min
            size[dim] = (long unsigned int)lo->GetBoundingBox()[loc+1] - index[dim] + 1;    //bbox max - min + 1
            loc += 2;
        }
        region.SetSize(size);
        region.SetIndex(index);

        mLabelImage->SetRequestedRegion(region);
        mOriginalGrayscaleImage->SetRequestedRegion(region);

        typename TextureFeaturesFilterType::Pointer textureFilter = TextureFeaturesFilterType::New();
        textureFilter->SetInput(mOriginalGrayscaleImage);
        textureFilter->SetMaskImage(mLabelImage);
        textureFilter->SetInsidePixelValue((unsigned int)lo->GetLabel());
        textureFilter->SetPixelValueMinMax(mGrayscaleMinPixelValue, mGrayscaleMaxPixelValue);
        textureFilter->FastCalculationsOn();
        textureFilter->Update();

        typename FeatureVectorType::Pointer featureVector = textureFilter->GetFeatureMeans();
        lo->SetTextureFeatures(featureVector);
    }

    for(unsigned int i=0; i<ImageDimension; i++)
        centroid[i] /= (double)nbOfPixels;
    lo->SetCentroid(centroid);


    vtkIdType vId = mObjectGraph->FindVertex(lo->GetLabel());
    if(ImageDimension==2)
        mObjectGraph->GetPoints()->SetPoint(vId, lo->GetCentroid()[0], lo->GetCentroid()[1], 0);
    else if(ImageDimension==3)
        mObjectGraph->GetPoints()->SetPoint(vId, lo->GetCentroid()[0], lo->GetCentroid()[1], lo->GetCentroid()[2]);

    if(mUnaryFeatureState[MEAN_RED_INTENSITY] == IsActive) {
        lo->SetRedAcc(accRed);
        lo->SetRedMean(accRed / (double)nbOfPixels);
        vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[MEAN_RED_INTENSITY].c_str()))->InsertValue(vId, lo->GetRedMean());
    }
    if(mUnaryFeatureState[MEAN_GREEN_INTENSITY] == IsActive) {
        lo->SetGreenAcc(accGreen);
        lo->SetGreenMean(accGreen / (double)nbOfPixels);
        vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[MEAN_GREEN_INTENSITY].c_str()))->InsertValue(vId, lo->GetGreenMean());
    }
    if(mUnaryFeatureState[MEAN_BLUE_INTENSITY] == IsActive) {
        lo->SetBlueAcc(accBlue);
        lo->SetBlueMean(accBlue / (double)nbOfPixels);
        vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[MEAN_BLUE_INTENSITY].c_str()))->InsertValue(vId, lo->GetBlueMean());
    }
    if(mUnaryFeatureState[RED_HISTOGRAM] == IsActive) {
        for(int i=0; i<UnaryFeatureComponents[RED_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[RED_HISTOGRAM] << i;
            vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->InsertValue(vId, (double)lo->GetRedHistogram()->GetFrequency(i) / (double)nbOfPixels);
        }
    }
    if(mUnaryFeatureState[GREEN_HISTOGRAM] == IsActive) {
        for(int i=0; i<UnaryFeatureComponents[GREEN_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[GREEN_HISTOGRAM] << i;
            vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->InsertValue(vId, (double)lo->GetGreenHistogram()->GetFrequency(i) / (double)nbOfPixels);
        }
    }
    if(mUnaryFeatureState[BLUE_HISTOGRAM] == IsActive) {
        for(int i=0; i<UnaryFeatureComponents[BLUE_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[BLUE_HISTOGRAM] << i;
            vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->InsertValue(vId, (double)lo->GetBlueHistogram()->GetFrequency(i) / (double)nbOfPixels);
        }
    }
    if(mUnaryFeatureState[RGB_HISTOGRAM] == IsActive) {
        for(int i=0; i<UnaryFeatureComponents[RGB_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[RGB_HISTOGRAM] << i;
            vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->InsertValue(vId, (double)lo->GetRGBHistogram()->GetFrequency(i) / (double)nbOfPixels);
        }
    }
    if(mUnaryFeatureState[TEXTURE_FEATURES] == IsActive) {
        for(int i=0; i<UnaryFeatureComponents[TEXTURE_FEATURES]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[TEXTURE_FEATURES] << i;
            vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->InsertValue(vId, lo->GetTextureFeatures()->ElementAt(i));
        }
    }
    if(mUnaryFeatureState[NUMBER_OF_PIXELS] == IsActive)
        vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[NUMBER_OF_PIXELS].c_str()))->InsertValue(vId, lo->GetNumberOfPixels());
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::InitializePrecalculatedUnaryLabelObjectFeatures(LabelObjectPointerType &lo)
{
    SizeValueType                               nbOfPixels = 0;
    double                                      physicalSize;
    ContinuousIndex< double, ImageDimension >   centroid, centroidAcc;
    centroid.Fill(0);
    centroidAcc.Fill(0);

    vtkIdType vId = mObjectGraph->FindVertex(lo->GetLabel());

    nbOfPixels = lo->Size();

    lo->SetNumberOfPixels(nbOfPixels);
    lo->SetPhysicalSize(nbOfPixels * mVoxelVolume);

    for(unsigned int i=0; i<ImageDimension; i++)
        centroid[i] = mObjectGraph->GetPoints()->GetPoint(vId)[i];
    lo->SetCentroid(centroid);

    for(unsigned int i=0; i<ImageDimension; i++)
        centroidAcc[i] = centroid[i] * nbOfPixels;
    lo->SetCoordAcc(centroidAcc);

    const unsigned int MeasurementVectorSize = 3;               // RGB
    HistogramSizeType histogramSize( MeasurementVectorSize );
    histogramSize.Fill(0);

    HistogramMeasurementVectorType lowerBound(MeasurementVectorSize);
    lowerBound.Fill(itk::NumericTraits<OriginalImagePixelComponentType>::min() - 0.5);

    HistogramMeasurementVectorType upperBound(MeasurementVectorSize);
    upperBound.Fill(itk::NumericTraits<OriginalImagePixelComponentType>::max() + 0.5);

    if(mUnaryFeatureState[MEAN_RED_INTENSITY] == IsActive) {
        double redMean = vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[MEAN_RED_INTENSITY].c_str()))->GetValue(vId);
        lo->SetRedMean(redMean);
        lo->SetRedAcc(redMean * nbOfPixels);
    }
    if(mUnaryFeatureState[MEAN_GREEN_INTENSITY] == IsActive) {
        double greenMean = vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[MEAN_GREEN_INTENSITY].c_str()))->GetValue(vId);
        lo->SetGreenMean(greenMean);
        lo->SetGreenAcc(greenMean * nbOfPixels);
    }
    if(mUnaryFeatureState[MEAN_BLUE_INTENSITY] == IsActive) {
        double blueMean = vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(UnaryFeatureName[MEAN_BLUE_INTENSITY].c_str()))->GetValue(vId);
        lo->SetBlueMean(blueMean);
        lo->SetBlueAcc(blueMean * nbOfPixels);
    }
    if(mUnaryFeatureState[RED_HISTOGRAM] == IsActive) {
        histogramSize[0] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Red   channel
        histogramSize[1] = 1;                                   // number of bins for the Green channel
        histogramSize[2] = 1;                                   // number of bins for the Blue channel

        typename HistogramType::Pointer histogram = HistogramType::New();
        histogram->SetMeasurementVectorSize(MeasurementVectorSize);
        histogram->Initialize(histogramSize, lowerBound, upperBound );

        for(int i=0; i<UnaryFeatureComponents[RED_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[RED_HISTOGRAM] << i;
            histogram->SetFrequency(i, static_cast<FrequencyType>(nbOfPixels * vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->GetValue(vId)) );
        }
        lo->SetRedHistogram(histogram);

        std::cout << "Red Histogram of LO " << lo->GetLabel() << " init discrepancy: total frequ = " << histogram->GetTotalFrequency() << " vs. nbOfPixels = " << nbOfPixels << std::endl;
    }
    if(mUnaryFeatureState[GREEN_HISTOGRAM] == IsActive) {
        histogramSize[0] = 1;                                   // number of bins for the Red   channel
        histogramSize[1] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Green channel
        histogramSize[2] = 1;                                   // number of bins for the Blue channel

        typename HistogramType::Pointer histogram = HistogramType::New();
        histogram->SetMeasurementVectorSize(MeasurementVectorSize);
        histogram->Initialize(histogramSize, lowerBound, upperBound );

        for(int i=0; i<UnaryFeatureComponents[GREEN_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[GREEN_HISTOGRAM] << i;
            histogram->SetFrequency(i, static_cast<FrequencyType>(nbOfPixels * vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->GetValue(vId)) );
        }
        lo->SetGreenHistogram(histogram);
    }
    if(mUnaryFeatureState[BLUE_HISTOGRAM] == IsActive) {
        histogramSize[0] = 1;                                   // number of bins for the Red   channel
        histogramSize[1] = 1;                                   // number of bins for the Green channel
        histogramSize[2] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Blue channel

        typename HistogramType::Pointer histogram = HistogramType::New();
        histogram->SetMeasurementVectorSize(MeasurementVectorSize);
        histogram->Initialize(histogramSize, lowerBound, upperBound );

        for(int i=0; i<UnaryFeatureComponents[BLUE_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[BLUE_HISTOGRAM] << i;
            histogram->SetFrequency(i, static_cast<FrequencyType>(nbOfPixels * vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->GetValue(vId)) );
        }
        lo->SetBlueHistogram(histogram);
    }
    if(mUnaryFeatureState[RGB_HISTOGRAM] == IsActive) {
        histogramSize[0] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Red   channel
        histogramSize[1] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Green channel
        histogramSize[2] = NUMBER_HISTOGRAM_BINS;               // number of bins for the Blue channel

        typename HistogramType::Pointer histogram = HistogramType::New();
        histogram->SetMeasurementVectorSize(MeasurementVectorSize);
        histogram->Initialize(histogramSize, lowerBound, upperBound );

        for(int i=0; i<UnaryFeatureComponents[RGB_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[RGB_HISTOGRAM] << i;
            histogram->SetFrequency(i, static_cast<FrequencyType>(nbOfPixels * vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->GetValue(vId)) );
        }
        lo->SetRGBHistogram(histogram);
    }
    if(mUnaryFeatureState[TEXTURE_FEATURES] == IsActive) {
        typename FeatureVectorType::Pointer featureVector = FeatureVectorType::New();
        featureVector->Reserve(UnaryFeatureComponents[TEXTURE_FEATURES]);

        for(int i=0; i<UnaryFeatureComponents[TEXTURE_FEATURES]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[TEXTURE_FEATURES] << i;
            featureVector->SetElement(i, vtkFloatArray::SafeDownCast(mObjectGraph->GetVertexData()->GetArray(arrayName.str().c_str()))->GetValue(vId) );
        }
        lo->SetTextureFeatures(featureVector);
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::InitializeBinaryLabelObjectFeatures(LabelObjectPointerType &lo)
{
    vtkSmartPointer<vtkUnsignedLongArray> pedigreeIds = vtkUnsignedLongArray::SafeDownCast(mObjectGraph->GetVertexData()->GetPedigreeIds());

    vtkSmartPointer<vtkOutEdgeIterator> it = vtkSmartPointer<vtkOutEdgeIterator>::New();
    mObjectGraph->GetOutEdges(mObjectGraph->FindVertex(lo->GetLabel()), it);

    while(it->HasNext()) {
        vtkOutEdgeType e = it->Next();
        LabelType source = lo->GetLabel();
        LabelType target = pedigreeIds->GetValue(e.Target);

        if(mBinaryFeatureState[MEAN_RED_SIMILARITY] == IsActive)      ComputeRedSimilarity(e.Id, source, target);
        if(mBinaryFeatureState[MEAN_GREEN_SIMILARITY] == IsActive)    ComputeGreenSimilarity(e.Id, source, target);
        if(mBinaryFeatureState[MEAN_BLUE_SIMILARITY] == IsActive)     ComputeBlueSimilarity(e.Id, source, target);

//        if(lo->GetLabel()==1)
//            std::cout << "Similarity Calculation: lo1=" << source << ", lo2=" << target << ", lo1_blueMean=" << mean1 << ", lo2_blueMean=" << mean2 << ", similarity=" << similarity << std::endl;
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::BuildLabelMapGraphAnnotations()
{
    GraphAnnotationHelper anno;
    std::map<vtkIdType, float> emptyFloat;

    if(mUnaryFeatureState[MEAN_RED_INTENSITY] == IsActive)
        anno.AddCustomVertexAnnotation(mObjectGraph, UnaryFeatureName[MEAN_RED_INTENSITY], emptyFloat);
    if(mUnaryFeatureState[MEAN_GREEN_INTENSITY] == IsActive)
        anno.AddCustomVertexAnnotation(mObjectGraph, UnaryFeatureName[MEAN_GREEN_INTENSITY], emptyFloat);
    if(mUnaryFeatureState[MEAN_BLUE_INTENSITY] == IsActive)
        anno.AddCustomVertexAnnotation(mObjectGraph, UnaryFeatureName[MEAN_BLUE_INTENSITY], emptyFloat);
    if(mUnaryFeatureState[RED_HISTOGRAM] == IsActive)
        for(int i=0; i<UnaryFeatureComponents[RED_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[RED_HISTOGRAM] << i;
            anno.AddCustomVertexAnnotation(mObjectGraph, arrayName.str(), emptyFloat);
        }
    if(mUnaryFeatureState[GREEN_HISTOGRAM] == IsActive)
        for(int i=0; i<UnaryFeatureComponents[GREEN_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[GREEN_HISTOGRAM] << i;
            anno.AddCustomVertexAnnotation(mObjectGraph, arrayName.str(), emptyFloat);
        }
    if(mUnaryFeatureState[BLUE_HISTOGRAM] == IsActive)
        for(int i=0; i<UnaryFeatureComponents[BLUE_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[BLUE_HISTOGRAM] << i;
            anno.AddCustomVertexAnnotation(mObjectGraph, arrayName.str(), emptyFloat);
        }
    if(mUnaryFeatureState[RGB_HISTOGRAM] == IsActive)
        for(int i=0; i<UnaryFeatureComponents[RGB_HISTOGRAM]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[RGB_HISTOGRAM] << i;
            anno.AddCustomVertexAnnotation(mObjectGraph, arrayName.str(), emptyFloat);
        }
    if(mUnaryFeatureState[TEXTURE_FEATURES] == IsActive)
        for(int i=0; i<UnaryFeatureComponents[TEXTURE_FEATURES]; i++) {
            std::stringstream arrayName;
            arrayName << UnaryFeatureName[TEXTURE_FEATURES] << i;
            anno.AddCustomVertexAnnotation(mObjectGraph, arrayName.str(), emptyFloat);
        }
    if(mUnaryFeatureState[NUMBER_OF_PIXELS] == IsActive)
        anno.AddCustomVertexAnnotation(mObjectGraph, UnaryFeatureName[NUMBER_OF_PIXELS], emptyFloat);

    if(mBinaryFeatureState[MEAN_RED_SIMILARITY] == IsActive)
        anno.AddCustomEdgeAnnotation(mObjectGraph, BinaryFeatureName[MEAN_RED_SIMILARITY], emptyFloat);
    if(mBinaryFeatureState[MEAN_GREEN_SIMILARITY] == IsActive)
        anno.AddCustomEdgeAnnotation(mObjectGraph, BinaryFeatureName[MEAN_GREEN_SIMILARITY], emptyFloat);
    if(mBinaryFeatureState[MEAN_BLUE_SIMILARITY] == IsActive)
        anno.AddCustomEdgeAnnotation(mObjectGraph, BinaryFeatureName[MEAN_BLUE_SIMILARITY], emptyFloat);
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::PreparationsForInitialization()
{
    std::cout << "PreparationsForInitialization" << std::endl;

    if(mUnaryFeatureState[TEXTURE_FEATURES] == IsActive) {
        typename RGBToGrayscaleFilterType::Pointer rgbToGrayscaleFilter = RGBToGrayscaleFilterType::New();
        rgbToGrayscaleFilter->SetInput(mOriginalImage);
        rgbToGrayscaleFilter->Update();

        mOriginalGrayscaleImage = rgbToGrayscaleFilter->GetOutput();
        mOriginalGrayscaleImage->DisconnectPipeline();

        typename ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
        imageCalculatorFilter->SetImage(mOriginalGrayscaleImage);
        imageCalculatorFilter->Compute();

        mGrayscaleMinPixelValue = imageCalculatorFilter->GetMinimum();
        mGrayscaleMaxPixelValue = imageCalculatorFilter->GetMaximum();


        std::cout << "start compute geometry stuff" << std::endl;
        //TODO: atm only bounding box used; other features later (makes manual centroid computation obsolete)

        typename LabelGeometryFilterType::Pointer labelGeometryImageFilter = LabelGeometryFilterType::New();
        labelGeometryImageFilter->SetInput(mLabelImage);
        labelGeometryImageFilter->SetIntensityInput(mOriginalGrayscaleImage);
        labelGeometryImageFilter->Update();

        typename LabelGeometryFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
        typename LabelGeometryFilterType::LabelsType::iterator allLabelsIt;

        for(allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++)
        {
            typename LabelGeometryFilterType::LabelPixelType labelValue = *allLabelsIt;

            LabelObjectContainerIterator it = mLabelObjectContainer.find((int)labelValue);
            it->second->SetBoundingBox(labelGeometryImageFilter->GetBoundingBox(labelValue));
        }
        std::cout << "end compute geometry stuff" << std::endl;
    }
}


template< class TOriginalImage, class TLabelImage, class TLabelObject > void LabelMapGraph<TOriginalImage, TLabelImage, TLabelObject>::InitializeAllLabelObjects(bool fromPrecalculatedGraph)
{
    if(mObjectGraph->GetNumberOfVertices() != mLabelObjectContainer.size()) {
        std::cout << "LabelMapGraph::InitializeAllLabelObjects(): number vertices (" << mObjectGraph->GetNumberOfVertices() << ") doesn't equal number label objects (" <<
                mLabelObjectContainer.size() << ") -- and that's a problem!" << std::endl;
        exit(EXIT_FAILURE);
    }

    mVoxelVolume = 1;
    for(unsigned int i = 0; i < ImageDimension; i++)
        mVoxelVolume *= mOriginalImage->GetSpacing()[i];

    if(!fromPrecalculatedGraph)
        BuildLabelMapGraphAnnotations();
    PreparationsForInitialization();    //TODO: will be able to compute couple of other shape features,
                                        //other name or save results and do allocation of shape features in the InitializeUnaryLabelObjectFeatures method

    for(LabelObjectContainerIterator it = mLabelObjectContainer.begin(); it != mLabelObjectContainer.end(); it++) {
        if(fromPrecalculatedGraph)
            this->InitializePrecalculatedUnaryLabelObjectFeatures(it->second);
        else
            this->InitializeUnaryLabelObjectFeatures(it->second);
    }

    for(LabelObjectContainerIterator it = mLabelObjectContainer.begin(); it != mLabelObjectContainer.end(); it++) {
        if(!fromPrecalculatedGraph)
            this->InitializeBinaryLabelObjectFeatures(it->second);
    }
}

} // end namespace
