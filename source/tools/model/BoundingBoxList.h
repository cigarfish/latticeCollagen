#ifndef BOUNDINGBOX_LIST_H
#define BOUNDINGBOX_LIST_H

#include "../../model/BasicDatatypes/CSListContainer.h"

#include <iostream>
#include <vector>
#include <stack>
#include <map>

#include "../../model/Elements/ModelElement.h"


class BoundingBoxList
{
public:
    BoundingBoxList( int Dimensions, unsigned long size=1048576 );
    ~BoundingBoxList();

    BoundingBoxList( const BoundingBoxList & other );

    void add( ModelElement * newElement );
    void remove( ModelElement * obsoleteElement );

    //! The size of the element container up to the last non-NULL element.
    //  This method can be used when iterating over all elements.
    unsigned long size() const {return EndIndex();};

    //! The actual number of non-NULL elements.
    //  This number may be smaller than size() - if elements have been removed -
    //  and thus should not be used for iteration over all elements with
    //  element(size_t).
    unsigned long numElements() const {return mNumElements;};

    void reIndex();

    virtual void update();
    void sweepNew();

    void moveY( unsigned long from, unsigned long toAfter );
    void moveToFrontY( unsigned long from );
    void sortAndBinY();

    CSListContainer< unsigned long > &
    intersectingListByIndex( ModelElement * currentObj );

    CSListContainer< ModelElement * >
    intersectingListByElement( ModelElement * currentObj );

    //! The element index of the mEnd token (divided by two)
    unsigned int EndIndex() const {return mEnd>>1;};

    // The general routine for finding the next ModelElement to 'object'.
    // Contains the code in which the 'kind' is Unspecified delegates work
    // to NextNeighborByType() if 'kind' is given otherwise.
    ModelElement * NextNeighbor( ModelElement * object, ModelElement::Type kind =ModelElement::TypeUnspecified );

    ModelElement * NextNeighbor( unsigned long index, ModelElement::Type kind =ModelElement::TypeUnspecified );

    ModelElement * NextNeighborByType( ModelElement * object, ModelElement::Type type );

    // debug routine:
    void output();

    // An iterator-like class used for going through the sorted list of
    // ModelElements.  Mainly for the loop in which interactions between
    // ModelElements are calculated.
    class sorted_iterator
    {
        friend class BoundingBoxList;
    public:
        sorted_iterator( BoundingBoxList * bblist, unsigned long index ) : mIndex(index), mpBBList(bblist) {};
        ModelElement * operator *() { return mpBBList->mppElements[mIndex]; };
        sorted_iterator & operator++()
            {
                while( (mIndex=mpBBList->mpForward[mIndex]) & 1 )
                    if(mIndex>=mpBBList->mEnd)
                        break;
                return *this;
            };
        sorted_iterator & operator--()
            {
                while( (mIndex=mpBBList->mpReverse[mIndex]) & 1 ) {};
                return *this;
            }
        bool operator == (sorted_iterator & rhs) { return (this->mIndex == rhs.mIndex); };
        bool operator != (sorted_iterator & rhs) { return (this->mIndex != rhs.mIndex); };

        unsigned long index() const { return mIndex; };

    private:
        unsigned long mIndex;
        BoundingBoxList *mpBBList;
    };

    // Accessor methods iterator-like interface:
    sorted_iterator begin_sorted() { return sorted_iterator( mIteratorStart ); };
    sorted_iterator & end_sorted() { return mIteratorEnd; };

    class iterator
        {
        public:
            iterator( BoundingBoxList * bblist, unsigned long index=0 )
                : mIndex(index), mpBBList(bblist)
                {
                    if (mIndex >= mpBBList->mEnd)
                    {
                        mIndex = mpBBList->mEnd;
                        return;
                    }
                    // make sure, our mIndex is even
                    if ( mIndex & 1) --mIndex;

                    while (mpBBList->mppElements[mIndex] == NULL)
                    {
                        if (mIndex>=mpBBList->mEnd)
                        {
                            mIndex = mpBBList->mEnd;
                            break;
                        }
                        mIndex += 2;
                    }
                }
            ModelElement * operator *() { return mpBBList->mppElements[mIndex]; };

            iterator & operator++ ()
                {
                    while(mpBBList->mppElements[mIndex+=2]==NULL)
                        if (mIndex>=mpBBList->mEnd) {mIndex=mpBBList->mEnd; break;}
                    return *this;
                };

            iterator & operator-- ()
                {
                    while(mpBBList->mppElements[mIndex-=2]==NULL)
                        if (mIndex<=0) {mIndex=0; break;}
                    return *this;
                };
            bool operator== (iterator rhs) { return  (this->mIndex == rhs.mIndex); };
            bool operator!= (iterator rhs) { return !(this->mIndex == rhs.mIndex); };

            unsigned long index() const {return mIndex;};

        private:
            unsigned long mIndex;
            BoundingBoxList *mpBBList;
        };

        iterator begin() {return iterator(this, 0); };
        iterator end() {return iterator(this, mEnd); };


    ModelElement * element(size_t i) const {return mppElements[i<<1];};
    ModelElement * at( unsigned long i ) const { return mppElements[i]; };

    int getDimensions() const {return mDimensions;};
    void setDimensions(int dimensions) { mDimensions = dimensions; };

    CSListContainer<ModelElement*> GetContainedElements( BoundingBox &, ModelElement::Type=ModelElement::TypeUnspecified );

    const BoundingBox & absoluteLimits() {return mGlobalLimits;};


protected:
    void removeEntry( unsigned long internalIndex );
    void move( unsigned long from, unsigned long toAfter );
    void moveToFront( unsigned long from );
    void insertionSort();
    void removeFromTableau( unsigned long index );


//private:
public:
    BoundingBox mGlobalLimits;

    int mDimensions;
    unsigned long mSize;
    unsigned long mNumElements; // contains the number of ModelElements
    unsigned long mNumEntries;  // contains 2* mNumElements;

    unsigned long mStart;
    unsigned long mEnd;
    unsigned long *mpForward;
    unsigned long *mpReverse;

    unsigned long *mpOpenEndsFwd;
    unsigned long *mpOpenEndsRev;

    // additional y axis aligned limits for subdivision into regions (for 3D)
    unsigned long mStartY;
    unsigned long *mpForwardY;
    unsigned long *mpReverseY;
    double **mppLimitsY;

    // for 3D bounding box, space will be divided into 8 or 64 regions along the
    // y axis with an (ideally) equal number of bounding boxes.  A ModelElement
    // can span several regions, each region a ModelElement belongs to is marked
    // by a bit in the bitfield (unsigned long) corresponding to its
    // mGlobalIndex:
    unsigned long *mpYRegions;


    sorted_iterator mIteratorStart;
    sorted_iterator mIteratorEnd;

    std::stack<unsigned long> nextEmpty;

    CSListContainer<unsigned long> mIntersectingList;

    // the payload
    double ** mppLimits;
    ModelElement **mppElements;
};

    
#endif // BOUNDINGBOX_LIST_H
