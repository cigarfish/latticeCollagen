////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  ListContainer.h                                               //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2013-10-22 12:22:20                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef CS_LIST_CONTAINER_H
#define CS_LIST_CONTAINER_H

#include <cstdlib>
#include <cassert>
#include <cmath>

template <class T>
class CSListContainer
{
public:
    CSListContainer( const size_t preallocate=1024 )
        : mAllocationSize( preallocate ),
          mDataSize(0),
          mpStart(NULL),
          mpEnd(NULL)
        {
            mpStart = (T *)malloc( mAllocationSize*sizeof(T) );
            mpEnd   = mpStart;
        };

    ~CSListContainer()
        { free( mpStart ); };

    void push_back( T element );

    T * find( T token );

    void clear();

    void swap( CSListContainer<T> & other );

    const T & operator[] (size_t index) const;
    T & operator[](size_t index);

    T * begin() const { return mpStart; };
    T * end()   const { return mpEnd; };

    size_t size() const { return mDataSize; };

    void reallocate( size_t allocationChunk=0 );


private:
    size_t mAllocationSize;
    size_t mDataSize;

    T * mpStart; // pointer to the (1st element of the) data array;
    T * mpEnd;   // pointer to the element after the last of the data array;
    
};


template <class T>
inline
void
CSListContainer<T>::push_back( T element )
{
    *mpEnd = element;

    ++mpEnd; ++mDataSize;

    if ( mDataSize >= mAllocationSize )
        reallocate();
};


template <class T>
inline
T *
CSListContainer<T>::find( T token )
{
    T * iter;
    for ( iter = mpStart; iter != mpEnd; ++iter )
        if ( *iter == token ) break;
    return iter;
}


template <class T>
inline
void
CSListContainer<T>::clear()
{
    // avoiding reallocation and zeroing out
    mpEnd = mpStart;
    mDataSize = 0;
}


template <class T>
inline
void
CSListContainer<T>::swap( CSListContainer<T> & otherContainer )
{
    size_t allocationSize = mAllocationSize;
    size_t dataSize       = mDataSize;
    T * start             = mpStart;
    T * end               = mpEnd;

    mAllocationSize = otherContainer.mAllocationSize;
    mDataSize       = otherContainer.mDataSize;
    mpStart         = otherContainer.mpStart;
    mpEnd           = otherContainer.mpEnd;

    otherContainer.mAllocationSize = allocationSize;
    otherContainer.mDataSize       = dataSize;
    otherContainer.mpStart         = start;
    otherContainer.mpEnd           = end;
}


template <class T>
inline
void
CSListContainer<T>::reallocate( size_t allocationChunk )
{
    if (!allocationChunk)
        allocationChunk = (mAllocationSize<4096) ?
            1024 :
            ((mAllocationSize >> 12) + (mAllocationSize&4095)?1:0) *1024;

    mAllocationSize += allocationChunk;

    mpStart = (T *)realloc( mpStart, mAllocationSize*sizeof(T) );
    mpEnd   = mpStart + mDataSize;
};


template <class T>
inline
const T &
CSListContainer<T>::operator[]( size_t index ) const
{
    assert ( index >=0 && index <= mDataSize );
    return mpStart[index];
};


template <class T>
inline
T &
CSListContainer<T>::operator[]( size_t index )
{
    assert ( index >=0 && index <= mDataSize );
    return mpStart[index];
};


#endif // CS_LIST_CONTAINER_H
