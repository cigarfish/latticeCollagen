///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  LabelImageToLabelMapGraphFilter.tpp                                  //
//                                                                                   //
//     Author:  Adrian Friebel <adrian.friebel@uni-leipzig.de>                       //
//    Created:  2013-09-04                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "LabelImageToLabelMapGraphFilter.h"

#include <itkLabelObject.h>

#include <vtkMutableUndirectedGraph.h>
#include <vtkSmartPointer.h>


namespace itk
{

template< class TOriginalImage, class TInputImage, class TOutputImage > LabelImageToLabelMapGraphFilter< TOriginalImage, TInputImage, TOutputImage >::LabelImageToLabelMapGraphFilter()
{
    mBackgroundValue = NumericTraits< OutputImagePixelType >::NonpositiveMin();
    usePrecalculatedGraph = false;
}


template< class TOriginalImage, class TInputImage, class TOutputImage > void LabelImageToLabelMapGraphFilter< TOriginalImage, TInputImage, TOutputImage >::UsePrecalculatedGraph(vtkSmartPointer<vtkUndirectedGraph>  graph)
{
    usePrecalculatedGraph = true;
    mObjectGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    mObjectGraph->CheckedDeepCopy(graph);
}


template< class TOriginalImage, class TInputImage, class TOutputImage > void LabelImageToLabelMapGraphFilter< TOriginalImage, TInputImage, TOutputImage >::GenerateInputRequestedRegion()
{
    // call the superclass' implementation of this method
    Superclass::GenerateInputRequestedRegion();

    // We need all the input.
    InputImagePointer input = const_cast< InputImageType* >( this->GetInput() );
    if(!input)
    {
        return;
    }
    input->SetRequestedRegion( input->GetLargestPossibleRegion() );
}


template< class TOriginalImage, class TInputImage, class TOutputImage > void LabelImageToLabelMapGraphFilter< TOriginalImage, TInputImage, TOutputImage >::EnlargeOutputRequestedRegion(DataObject *)
{
    this->GetOutput()->SetRequestedRegion( this->GetOutput()->GetLargestPossibleRegion() );
}


template< class TOriginalImage, class TInputImage, class TOutputImage > void LabelImageToLabelMapGraphFilter< TOriginalImage, TInputImage, TOutputImage >::GenerateData()
{
    InputImageConstPointer input  = this->GetInput();

    OutputImagePointer output = this->GetOutput();
    output->SetLabelImage(mLabelImage);
    output->SetOriginalImage(mOriginalImage);
    output->SetBackgroundValue(mBackgroundValue);

    typedef ImageLinearConstIteratorWithIndex< InputImageType > InputLineIteratorType;
    InputLineIteratorType it(input, input->GetLargestPossibleRegion());
    it.SetDirection(0);

    for(it.GoToBegin(); !it.IsAtEnd(); it.NextLine())
    {
        it.GoToBeginOfLine();

        while(!it.IsAtEndOfLine())
        {
            const InputImagePixelType &v = it.Get();

            if(v != static_cast< InputImagePixelType >(mBackgroundValue))
            {
                // We've hit the start of a run
                IndexType idx = it.GetIndex();
                LengthType      length = 1;
                ++it;
                while(!it.IsAtEndOfLine() && it.Get() == v)
                {
                    ++length;
                    ++it;
                }
                // create the run length object to go in the vector
                if(usePrecalculatedGraph)   output->LazySetLine(idx, length, v);
                else                        output->SetLine(idx, length, v);
            }
            else
            {
                // go the the next pixel
                ++it;
            }
        }
    }
    if(usePrecalculatedGraph)   output->SetGraph(mObjectGraph);
    output->Optimize();
}


template< class TOriginalImage, class TInputImage, class TOutputImage > void LabelImageToLabelMapGraphFilter< TOriginalImage, TInputImage, TOutputImage >::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);

    os << indent << "BackgroundValue: " << static_cast< typename NumericTraits< OutputImagePixelType >::PrintType >(mBackgroundValue) << std::endl;
}


} // end namespace itk

