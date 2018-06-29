////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  BoundingBoxList.cpp                                           //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2012-11-04 23:24:36                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "BoundingBoxList.h"

#include "../../model/Elements/ModelElement.h"

#include <limits>
#include <iostream>
#include <cstring>


BoundingBoxList::BoundingBoxList( int dimensions, unsigned long size )
    : mDimensions(dimensions),
      mSize(2*size),
      mIteratorStart( this, 0 ),
      mIteratorEnd( this, 0 )
{
    mpForward = new unsigned long[mSize];
    mpReverse = new unsigned long[mSize];

    mppLimits   = new double*[mSize];
    mppElements = new ModelElement*[mSize];

    mStart = 0;
    mEnd = 0;

    mNumEntries  = 0;
    mNumElements = 0;

    mpOpenEndsFwd = new unsigned long[mSize];
    mpOpenEndsRev = new unsigned long[mSize];

    mppLimitsY = new double*[mSize];
    mpForwardY = new unsigned long[mSize];
    mpReverseY = new unsigned long[mSize];
    mpYRegions = new unsigned long[size];

    mStartY = 0;

}


BoundingBoxList::BoundingBoxList( const BoundingBoxList & other )
    : mDimensions( other.mDimensions ),
      mSize( other.mSize ),
      mNumElements( other.mNumElements ),
      mNumEntries( other.mNumEntries ),
      mStart( other.mStart ),
      mEnd( other.mEnd ),
      mStartY( other.mStartY ),
      mIteratorStart( this, other.mIteratorStart.mIndex ),
      mIteratorEnd( this, other.mIteratorStart.mIndex )
{
    mppElements = new ModelElement*[mSize];
    std::copy( other.mppElements,
               other.mppElements + mSize,
               mppElements);

    mppLimits  = new double*[mSize];
    mppLimitsY = new double*[mSize];
    std::copy( other.mppLimits, other.mppLimits + mSize, mppLimits );
    std::copy( other.mppLimitsY, other.mppLimitsY + mSize, mppLimitsY );

    mpForward = new unsigned long[mSize];
    mpReverse = new unsigned long[mSize];
    std::copy( other.mpForward, other.mpForward + mSize, mpForward );
    std::copy( other.mpReverse, other.mpReverse + mSize, mpReverse );

    mpForwardY = new unsigned long[mSize];
    mpReverseY = new unsigned long[mSize];
    mpYRegions = new unsigned long[mSize/2];
    std::copy( other.mpForwardY, other.mpForwardY + mSize, mpForwardY );
    std::copy( other.mpReverseY, other.mpReverseY + mSize, mpReverseY );
    std::copy( other.mpYRegions, other.mpYRegions + mSize/2, mpYRegions );

    // The following is copied for completeness
    mpOpenEndsFwd = new unsigned long[mSize];
    mpOpenEndsRev = new unsigned long[mSize];
    std::copy( other.mpOpenEndsFwd, other.mpOpenEndsFwd + mSize, mpOpenEndsFwd );
    std::copy( other.mpOpenEndsRev, other.mpOpenEndsRev + mSize, mpOpenEndsRev );
}


BoundingBoxList::~BoundingBoxList()
{
    delete[] mpForward;
    delete[] mpReverse;
    delete[] mppLimits;
    delete[] mppElements;

    delete[] mpOpenEndsFwd;
    delete[] mpOpenEndsRev;

    delete[] mppLimitsY;
    delete[] mpForwardY;
    delete[] mpReverseY;
    delete[] mpYRegions;

}


void
BoundingBoxList::add( ModelElement * newElement )
{
    mNumEntries+=2;
    ++mNumElements;

    // TODO:  reallocation
    if ( mNumEntries > mSize )
        std::cout << "overflow\n";

    if ( nextEmpty.empty() )
    {
        newElement->mGlobalIndex = mEnd >>1;

        mppLimits[mEnd] =  & newElement->getBoundingBox()->xmin;
        mppElements[mEnd] = newElement;

        mpReverse[mEnd +1] = mEnd;
        mpForward[mEnd] = mEnd+1;
        ++mEnd;

        mppLimits[mEnd] = & newElement->getBoundingBox()->xmax;
        mppElements[mEnd] = newElement;

        mpReverse[mEnd +1] = mEnd;
        mpForward[mEnd] = mEnd+1;
        ++mEnd;

        mIteratorEnd = sorted_iterator( this, mEnd );

        // mEnd is the same index for x and y since it is the element after
        // the last element in mppElements.
        // Therefore we have to go backwards...
        mpReverseY[mEnd] = mEnd-1;
        mpForwardY[mEnd-1] = mEnd;
        mpReverseY[mEnd-1] = mEnd-2;
        mpForwardY[mEnd-2] = mEnd-1;
        mppLimitsY[mEnd-2] = & newElement->getBoundingBox()->ymin;
        mppLimitsY[mEnd-1] = & newElement->getBoundingBox()->ymax;

    }
    else
    {
        unsigned long insertAt = nextEmpty.top();
        nextEmpty.pop();

        newElement->mGlobalIndex = insertAt >>1;

        mppLimits[insertAt] = & newElement->getBoundingBox()->xmin;
        mppElements[insertAt] = newElement;

        mpForward[insertAt] = insertAt +1;
        mpForward[ mpReverse[mEnd] ] = insertAt;
        mpReverse[insertAt] = mpReverse[mEnd];

        ++insertAt;

        mppLimits[insertAt] = & newElement->getBoundingBox()->xmax;
        mppElements[insertAt] = newElement;

        mpForward[insertAt] = mEnd;
        mpReverse[insertAt] = insertAt -1;
        mpReverse[mEnd] = insertAt;

        if ( mNumEntries == 2 )
            mEnd = 2;

        mpForwardY[ mpReverseY[mEnd] ] = insertAt - 1;
        mpReverseY[ insertAt - 1 ] = mpReverseY[mEnd];

        mpReverseY[mEnd] = insertAt;
        mpForwardY[insertAt] = mEnd;

        mpReverseY[insertAt] = insertAt - 1;
        mpForwardY[insertAt-1] = insertAt;

        mppLimitsY[insertAt-1] = & newElement->getBoundingBox()->ymin;
        mppLimitsY[insertAt]   = & newElement->getBoundingBox()->ymax;

    }
}


void
BoundingBoxList::remove( ModelElement * obsoleteElement )
{
    unsigned long index = obsoleteElement->mGlobalIndex<<1;

    if ( obsoleteElement != mppElements[index] )
        return;

    nextEmpty.push(index);

    removeEntry( index );
    removeEntry( index + 1 );

    --mNumElements;

    if ( !mNumEntries )
    {
        mStart = 0;
        mEnd   = 0;
        mStartY = 0;
        while (!nextEmpty.empty())
            nextEmpty.pop();
    }
}


void
BoundingBoxList::removeEntry( unsigned long internalIndex )
{
    --mNumEntries;

    mppElements[internalIndex] = NULL;

    if ( internalIndex == mStartY )
    {
        mStartY = mpForwardY[mStartY];
        mpReverseY[mStartY] = mStartY;
    }
    else
    {
        mpForwardY[mpReverseY[internalIndex]] = mpForwardY[internalIndex];
        mpReverseY[mpForwardY[internalIndex]] = mpReverseY[internalIndex];
    }

    if ( internalIndex == mStart )
    {
        mStart = *(mpForward + mStart);
        *(mpReverse + mStart) = mStart;
        mIteratorStart = sorted_iterator( this, mStart );
        return;
    }

    if ( internalIndex == mEnd )
    {
        mEnd = *(mpReverse + mEnd);
        *(mpForward + mEnd) = mEnd;
        mIteratorEnd = sorted_iterator(this, mEnd);
        return;
    }

    mpForward[mpReverse[internalIndex]] = mpForward[internalIndex];
    mpReverse[mpForward[internalIndex]] = mpReverse[internalIndex];
}


void
BoundingBoxList::update()
{
    mGlobalLimits = {0.,0.,0.,0.,0.,0.};
    for ( unsigned long i=0; i < mEnd; i+=2 )
        if (mppElements[i])
        {
            mppElements[i]->mIntersectingList.clear();

            if ( mppElements[i]->boundingBox()->xmin < mGlobalLimits.xmin )
                mGlobalLimits.xmin = mppElements[i]->boundingBox()->xmin;
            if ( mppElements[i]->boundingBox()->xmax > mGlobalLimits.xmax )
                mGlobalLimits.xmax = mppElements[i]->boundingBox()->xmax;
            if ( mppElements[i]->boundingBox()->ymin < mGlobalLimits.ymin )
                mGlobalLimits.ymin = mppElements[i]->boundingBox()->ymin;
            if ( mppElements[i]->boundingBox()->ymax > mGlobalLimits.ymax )
                mGlobalLimits.ymax = mppElements[i]->boundingBox()->ymax;
            if ( mppElements[i]->boundingBox()->zmin < mGlobalLimits.zmin )
                mGlobalLimits.zmin = mppElements[i]->boundingBox()->zmin;
            if ( mppElements[i]->boundingBox()->zmax > mGlobalLimits.zmax )
                mGlobalLimits.zmax = mppElements[i]->boundingBox()->zmax;
        }

    //sort
    insertionSort();

    sortAndBinY();
    sweepNew();

}

// Implementation of the algorithm shortly described in paragraph 2 of:
// DJ Tracy, SR Buss, BM Woods: Efficient Large-scale Sweep and Prune Methods
// with AABB Insertion and Removal; IEEE-VR 2009.
void
BoundingBoxList::sweepNew()
{
    unsigned long currentIndex = mStart;

    // linked list for the possible overlapping element's indices:
    unsigned long openEndsStart = mEnd;

    while ( currentIndex != mEnd )
    {
        if ( currentIndex & 1 ) // a max limit
        {
            // remove currentIndex-1 from openEnds linked list:
            if ( openEndsStart == currentIndex-1 )
            {
                openEndsStart = mpOpenEndsFwd[openEndsStart];
                mpOpenEndsRev[ openEndsStart ] = openEndsStart;
            }
            else
            {
                mpOpenEndsRev[ mpOpenEndsFwd[currentIndex-1] ] = mpOpenEndsRev[ currentIndex-1 ];
                mpOpenEndsFwd[ mpOpenEndsRev[currentIndex-1] ] = mpOpenEndsFwd[ currentIndex-1 ];
            }

            // test with all bounding boxes in openEnds
            for ( unsigned long possibleContactIndex = openEndsStart; possibleContactIndex != mEnd; possibleContactIndex = mpOpenEndsFwd[ possibleContactIndex ] )
            {
                if ( (mpYRegions[possibleContactIndex>>1] & mpYRegions[currentIndex>>1]) )
                {
                    if (mDimensions==3)
                    {
                        if ( mppElements[currentIndex]->mBoundingBox.zmin > mppElements[possibleContactIndex]->mBoundingBox.zmax ||
                             mppElements[currentIndex]->mBoundingBox.zmax < mppElements[possibleContactIndex]->mBoundingBox.zmin )
                            continue;
                    }

                    if ( mppElements[currentIndex]->mBoundingBox.ymin > mppElements[possibleContactIndex]->mBoundingBox.ymax ||
                         mppElements[currentIndex]->mBoundingBox.ymax < mppElements[possibleContactIndex]->mBoundingBox.ymin )
                        continue;

                    // we've come this far, the bounding boxes overlap, add to the iterator's mIntersectingList:
                    mppElements[possibleContactIndex]->mIntersectingList.push_back(currentIndex>>1);
                    mppElements[currentIndex]->mIntersectingList.push_back(possibleContactIndex>>1);
                }
            }

        }
        else
        {
            // insert currentIndex into openEnds linked list:
            if ( openEndsStart == mEnd )
            {
                openEndsStart = currentIndex;
            }
            else
            {
                mpOpenEndsFwd[ mpOpenEndsRev[mEnd] ] = currentIndex;
                mpOpenEndsRev[ currentIndex ] = mpOpenEndsRev[mEnd];
            }

            mpOpenEndsRev[ mEnd ] = currentIndex;
            mpOpenEndsFwd[ currentIndex ] = mEnd;
        }

        currentIndex = mpForward[ currentIndex ];
    }
}

void
BoundingBoxList::sortAndBinY()
{
    if ( mNumElements < 2 )
        return;

    memset(mpYRegions, 0, mNumElements*sizeof(unsigned long));

    // initialise:
    // at the start
    unsigned long current = mpForwardY[mStartY];
    unsigned long seeker = mStartY;

    bool insertAtStart = false;

    //mppElements[mStartY]->boundingBox();
    while ( current != mEnd )
    {
        // updating the bounding box is not necessary if this method is called
        // after insertionSort()!
        // mppElements[current]->boundingBox();

        while ( *mppLimitsY[seeker] > *mppLimitsY[current] )
        {
            if ( seeker == mStartY )
            {
                insertAtStart = true;
                break;
            }

            seeker = mpReverseY[seeker];
        }

        unsigned long next =  mpForwardY[current];

        if ( insertAtStart )
        {
            insertAtStart = false;
            moveToFrontY( current );
        }
        else
        {
            if( current != mpForwardY[seeker] )
                moveY( current, seeker );
        }

        current = next;
        seeker = mpReverseY[current];
    }

    // Binning:
    // if mNumElements is high enough, bin into one of 64 compartments else use
    // only 8.  High enough means at least 128 per compartment, i.e. 8192
    // elements.
    // unsigned long elementsPerRegion = ( mNumElements > 8192 ) ? mNumElements>>6 : mNumElements>>3;

    // bin into up to 64 compartments:
    unsigned long elementsPerRegion = (mNumElements>>6) + ( (mNumElements&63) ? 1 : 0 );

    unsigned long elementsInRegion = 0;
    unsigned long currentRegion = 1;

    for ( current = mStartY; current != mEnd; current = mpForwardY[current] )
    {
        if ( elementsInRegion > elementsPerRegion )
        {
            currentRegion <<= 1; // left shift bit
            elementsInRegion = 0;
        }

        // if the element's minimum has been registered in a different region it
        // might be that it even spans more than two region for which we need to
        // fill the bits in its mRegion bit field.
        if ( current & 1 ) // max boundary
        {
            if ( mpYRegions[current>>1] & (currentRegion-1) )
            {
                unsigned long probeRegion = (currentRegion>>1);
                while ( !(probeRegion & mpYRegions[current>>1]) )
                {
                    mpYRegions[current>>1] |= probeRegion;
                    probeRegion >>=1;
                }
            }
        }
        else
            ++elementsInRegion;

        mpYRegions[current>>1] |= currentRegion;
    }
}

void BoundingBoxList::reIndex(){

  for( unsigned long i = 0 ; i < this->mEnd/2 ; i++ ){
    if( this->mppElements[2*i] != NULL )
      this->mppElements[2*i]->mGlobalIndex = i;
  }

}

CSListContainer< unsigned long > &
BoundingBoxList::intersectingListByIndex( ModelElement * currentObj )
{
    // go through mpForward to see what overlaps
    // skipping max entries
    unsigned long currentIndex = currentObj->mGlobalIndex<<1;
    unsigned long intersectingIndex = currentIndex;

    mIntersectingList.clear();

    const BoundingBox * intersectingBB;
    const BoundingBox * currentBB = currentObj->getBoundingBox();

    while ( 1 )
    {
        intersectingIndex = mpForward[intersectingIndex];

        // stop, if we found the xmax of the currentObj
        if ( currentIndex+1 == intersectingIndex )
            break;

        if ( intersectingIndex & 1 )
            continue;

        intersectingBB = mppElements[intersectingIndex]->getBoundingBox();

        if( currentBB->ymin > intersectingBB->ymax || currentBB->ymax < intersectingBB->ymin )
          continue;

        if ( mDimensions == 3 )
        {
            if( currentBB->zmin > intersectingBB->zmax || currentBB->zmax < intersectingBB->zmin )
        	          continue;
        }

        // okay!! overlap!! put into return list
        mIntersectingList.push_back( (unsigned long) intersectingIndex >>1 );
    }

    return mIntersectingList;
}

CSListContainer< ModelElement * >
BoundingBoxList::intersectingListByElement( ModelElement * currentObj )
{
    // go through mpForward to see what overlaps
    // skipping max entries
    unsigned long currentIndex = currentObj->mGlobalIndex<<1;
    unsigned long intersectingIndex = currentIndex;

    const BoundingBox * intersectingBB;
    const BoundingBox * currentBB = currentObj->getBoundingBox();

    CSListContainer< ModelElement * > overlapList(32);

    while ( 1 )
    {
        intersectingIndex = mpForward[intersectingIndex];

        // stop, if we found the xmax of the currentObj
        if ( currentIndex+1 == intersectingIndex )
            break;

        if ( intersectingIndex & 1 )
            continue;

        intersectingBB = mppElements[intersectingIndex]->getBoundingBox();

        if( currentBB->ymin > intersectingBB->ymax || currentBB->ymax < intersectingBB->ymin )
          continue;

        if ( mDimensions == 3 )
        {
            if( currentBB->zmin > intersectingBB->zmax || currentBB->zmax > intersectingBB->zmin )
                continue;
        }

        // okay!! overlap!! put into return list
        overlapList.push_back( mppElements[intersectingIndex] );
    }

    return overlapList;
}


void
BoundingBoxList::move( unsigned long from, unsigned long toAfter )
{
    unsigned long *fromFwdPtr = mpForward + from;
    unsigned long *fromRvsPtr = mpReverse + from;
    unsigned long *toFwdPtr   = mpForward + toAfter;

    if ( from == mStart )
    {
        mStart = *fromFwdPtr;
        mpReverse[*fromFwdPtr] = mStart;
        mIteratorStart = sorted_iterator(this, mStart);
    }
    else
    {
        // short-cut the neighbors
        mpForward[*fromRvsPtr] = *fromFwdPtr;
        mpReverse[*fromFwdPtr] = *fromRvsPtr;
    }

    // relink the pointers from the moved element
    *fromFwdPtr = *toFwdPtr;
    *fromRvsPtr = mpReverse[*toFwdPtr];

    // relink the pointers to the moved element
    mpReverse[*toFwdPtr] = from;
    *toFwdPtr = from;
}

void
BoundingBoxList::moveToFront( unsigned long from )
{
    // unlink
    mpForward[ mpReverse[from] ] = mpForward[from];
    mpReverse[ mpForward[from] ] = mpReverse[from];

    // add to front
    mpForward[from] = mStart;
    mpReverse[mStart] = from;
    mpReverse[from] = from;

    // reset mStart
    mStart = from;
    mIteratorStart = sorted_iterator(this, mStart);
}

void
BoundingBoxList::moveToFrontY( unsigned long from )
{
    // unlink
    mpForwardY[ mpReverseY[from] ] = mpForwardY[from];
    mpReverseY[ mpForwardY[from] ] = mpReverseY[from];

    // add to front
    mpForwardY[from] = mStartY;
    mpReverseY[mStartY] = from;
    mpReverseY[from] = from;

    // reset mStart
    mStartY = from;
}

void
BoundingBoxList::moveY( unsigned long from, unsigned long toAfter )
{
    unsigned long *fromFwdPtr = mpForwardY + from;
    unsigned long *fromRvsPtr = mpReverseY + from;
    unsigned long *toFwdPtr   = mpForwardY + toAfter;

    if ( from == mStartY )
    {
        mStartY = *fromFwdPtr;
        mpReverseY[*fromFwdPtr] = mStartY;
    }
    else
    {
        // short-cut the neighbors
        mpForwardY[*fromRvsPtr] = *fromFwdPtr;
        mpReverseY[*fromFwdPtr] = *fromRvsPtr;
    }

    // relink the pointers from the moved element
    *fromFwdPtr = *toFwdPtr;
    *fromRvsPtr = mpReverseY[*toFwdPtr];

    // relink the pointers to the moved element
    mpReverseY[*toFwdPtr] = from;
    *toFwdPtr = from;
}

void
BoundingBoxList::insertionSort()
{
    if ( !mNumElements )
        return;

    // initialise:
    // at the start
    unsigned long current = mpForward[mStart];
    unsigned long seeker = mStart;

    bool insertAtStart = false;


    mppElements[mStart]->boundingBox();
    while ( current != mEnd )
    {
        mppElements[current]->boundingBox();

        while ( *mppLimits[seeker] > *mppLimits[current] )
        {
            if ( seeker == mStart )
            {
                insertAtStart = true;
                break;
            }

            seeker = mpReverse[seeker];
        }

        unsigned long next =  mpForward[current];

        if ( insertAtStart )
        {
            insertAtStart = false;
            moveToFront( current );
        }
        else
        {
            if( current != mpForward[seeker] )
                move( current, seeker );
        }

        current = next;
        seeker = mpReverse[current];
    }
}


void
BoundingBoxList::output()
{
    unsigned long i = mStart;
    unsigned long* next = mpForward + mStart;
    ModelElement *element = mppElements[mStart];

    while ( i != mEnd )
    {
        std::cout << **(mppLimits+i) << "\t" << element << "\t" << ((i&1)?"max":"min") << std::endl;
        i = *next;
        element = *(mppElements + i);
        next = mpForward + i;
    }

    // std::cout << "reverse\n";

    // i = *(mpReverse + mEnd);
    // while (1)
    // {
    //     element = *(mppElements + i);
    //     std::cout << element << std::endl;
    //     if ( i == mStart )
    //         break;
    //     i = *(mpReverse + i);
    // }
}


// First implementation of a next neighbor search, assuming only spherical
// ModelElements, thus the surface to surface distance can be calculated easily.
ModelElement *
BoundingBoxList::NextNeighbor( ModelElement * object, ModelElement::Type elementType )
{
    if ( elementType )
        return NextNeighborByType( object, elementType );

    const Vector3f thisPosition = object->position;

    const double thisXCoordinate = thisPosition.x;

    double maxX = std::numeric_limits<double>::max();
    double minX = -maxX;

    double distanceRight = std::numeric_limits<double>::max();
    double distanceLeft  = std::numeric_limits<double>::max();

    // only for initialization
    unsigned int toLeft = object->mGlobalIndex<<1;
    unsigned int toRight = toLeft;

    while ( toRight != mEnd )
    {
        unsigned long right = mpForward[toRight];

        if (*mppLimits[right] > maxX )
            break;

        toRight = right;

        if ( toRight == (toLeft+1) )
            continue;

        Vector3f otherPosition = mppElements[toRight]->position;
        // for circular objects, the radius is equal to the difference between centre and l
        double otherRadius = std::fabs(*mppLimits[toLeft] - otherPosition.x);

        otherPosition.Multiply(-1.);
        otherPosition.Add( thisPosition );
        distanceRight = otherPosition.Norm();

        double distanceToSurface = thisXCoordinate + distanceRight - otherRadius;

        maxX = (distanceToSurface < maxX) ? distanceToSurface : maxX;
    }


    while (toLeft != mStart)
    {
        unsigned long left = mpReverse[toLeft];
        if (*mppLimits[left] < minX )
            break;
        toLeft = left;

        Vector3f otherPosition = mppElements[toLeft]->position;
        // for circular objects, the radius is equal to the difference between centre and l
        double otherRadius = std::fabs(*mppLimits[toLeft] - otherPosition.x);

        otherPosition.Multiply(-1.);
        otherPosition.Add( thisPosition );
        distanceLeft = otherPosition.Norm();

        double distanceToSurface = thisXCoordinate - distanceLeft + otherRadius;

        minX = (distanceToSurface < minX) ? minX : distanceToSurface;
    }

    unsigned long found = ( distanceLeft < distanceRight ) ? toLeft : toRight;

    return mppElements[found];
}


ModelElement *
BoundingBoxList::NextNeighbor( unsigned long int index, ModelElement::Type kind )
{
    return NextNeighbor( mppElements[index], kind );
}


ModelElement *
BoundingBoxList::NextNeighborByType( ModelElement * object, ModelElement::Type elementType )
{
    const Vector3f thisPosition = object->position;

    const double thisXCoordinate = thisPosition.x;

    double maxX = std::numeric_limits<double>::max();
    double minX = -maxX;

    double distanceRight = std::numeric_limits<double>::max();
    double distanceLeft  = std::numeric_limits<double>::max();

    // only for initialization
    // search to the right will start from initial left limit
    unsigned int toRight = object->mGlobalIndex<<1;
    // search to the left will start with the initial right limit
    unsigned int toLeft = toRight+1;
    unsigned int closestToRight = toRight;
    unsigned int closestToLeft  = toLeft;

    unsigned long nextRight;
    while ( (nextRight = mpForward[toRight]) != mEnd )
    {
        if (*mppLimits[nextRight] > maxX )
            break;

        toRight = nextRight;

        if ( mppElements[toRight]->mType != elementType )
            continue;

        // ignore this object's other boundary
        if ( toRight == toLeft )
            continue;

        Vector3f otherPosition = mppElements[toRight]->position;
        // for circular objects, the radius is equal to the difference between centre and its (right) BBox limit
        double otherRadius = std::fabs( *mppLimits[toRight|1] - otherPosition.x );

        // re-using otherPosition to calculate the distance:
        otherPosition.Multiply(-1.);
        otherPosition.Add( thisPosition );
        double distanceToSurface = otherPosition.Norm() - otherRadius;

        if ( distanceToSurface < distanceRight )
        {
            closestToRight = toRight;
            distanceRight = distanceToSurface;
            // this is the x coordinate beyond which we don't need to search:
            maxX = thisXCoordinate + distanceToSurface;
        }
    }

    toRight = object->mGlobalIndex<<1;

    unsigned long nextLeft;
    while ( (nextLeft = mpReverse[toLeft]) != mStart)
    {
        if (*mppLimits[nextLeft] < minX )
            break;

        toLeft = nextLeft;

        if ( mppElements[toLeft]->mType != elementType )
            continue;

        // ignore this object's other boundary
        if ( toLeft == toRight )
            continue;

        Vector3f otherPosition = mppElements[toLeft]->position;
        // for circular objects, the radius is equal to the difference between centre and its (right) BBox limit
        double otherRadius = std::fabs(*mppLimits[toLeft|1] - otherPosition.x);

        // re-using otherPosition to calculate the distance:
        otherPosition.Multiply(-1.);
        otherPosition.Add( thisPosition );
        double distanceToSurface = otherPosition.Norm() - otherRadius;

        if ( distanceToSurface < distanceLeft )
        {
            closestToLeft = toLeft;
            distanceLeft = distanceToSurface;
            // this is the x coordinate beyond which we don't need to search:
            minX = thisXCoordinate - distanceToSurface;
        }
    }

    // if there is no object of that type, we have to return NULL:
    if ( (closestToRight == toRight) && (closestToLeft == closestToRight+1) )
        return NULL;

    // the minimum of distanceLeft and distanceRight is the distance of the next neighbor
    // take the corresponding index and return its element.
    unsigned long found = ( distanceLeft < distanceRight ) ? closestToLeft : closestToRight;

    return mppElements[found];
}


CSListContainer<ModelElement*>
BoundingBoxList::GetContainedElements( BoundingBox &box, ModelElement::Type elementType )
{
    CSListContainer<ModelElement*> returnList;

    // search forward to the xmin of box
    unsigned long index = mStart;
    while ( *mppLimits[index] < box.xmin )
    {
        if ( (index = mpForward[index]) == mEnd )
            break;
    }

    if (index == mEnd )
        return returnList;

    // find all limits up to box.xmax
    do
    {
        if ( *mppLimits[index] > box.xmax )
            break;

        // test if found element's other bb limits are in range (one each dim suffices).
        BoundingBox *testBB = mppElements[index]->getBoundingBox();

        if ( testBB->ymin > box.ymax || box.ymin > testBB->ymax )
            continue;

        if (mDimensions==3)
            if ( testBB->zmin > box.zmax || box.zmin > testBB->zmax )
                continue;

        if ( (!elementType) || mppElements[index]->mType==elementType)
            returnList.push_back(mppElements[index]);
    } while ( (index = mpForward[index]) != mEnd );

    return returnList;
}
