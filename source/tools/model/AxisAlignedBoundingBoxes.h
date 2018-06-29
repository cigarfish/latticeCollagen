///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  AxisAlignedBoundingBoxes.h                                           //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-07-23 22:28:14                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef CS_AABB_H
#define CS_AABB_H

#include <algorithm>
#include <map>
#include <vector>
#include <list>

// t1m-debug
#include <iostream>
class ARGBColor;
char colorToChar(ARGBColor &);
// !t1m-debug

template <class Element>
class AxisAlignedBoundingBoxes;


template <class Element>
class PairContainer
{
  friend class AxisAlignedBoundingBoxes<Element>;

 public:
  PairContainer( std::map< Element*, std::vector<Element*> > * sparseMatrix );

  std::pair< Element *, Element * > first();
  std::pair< Element *, Element * > next();
  bool end();

 private:
  typename std::map< Element*, std::vector<Element*> > * mpSparseMatrix;
  typename std::map< Element*, std::vector<Element*> >::iterator rowIterator;
  typename std::vector<Element*>::iterator columnIterator;
};


template <class Element>
class AxisAlignedBoundingBoxes
{
 public:
  AxisAlignedBoundingBoxes(int dimension=3);

  void update();
  void add( Element *newElement );
  void erase( Element *newElement );

  const PairContainer<Element> & returnPairs();

 private:
  void sort( int whichAxis );
  void sortXYZ();

  void buildXPairs();

  void prunePairs();

  PairContainer<Element> * mpPairs;
  typename std::map< Element*, std::vector<Element*> > mSparseMatrix;

  typename std::list<Element *> mObjectsOrdered[3];
  typename std::list<double *> mProjections[3];

  int mDimensions;
};


template <class Element>
inline
std::pair< Element *, Element * >
PairContainer<Element>::first()
{
  // t1m-debug
  // std::cout << "mpSparseMatrix in PairContainer::first:\n";
  // std::map<ModelElement *, std::vector<ModelElement *> >::const_iterator mapIterator
  //   = mpSparseMatrix->begin();
  // while ( mapIterator != mpSparseMatrix->end() )
  //   {
  //     std::cout << colorToChar( mapIterator->first->color ) << " overlaps with\n";
  //     std::vector<ModelElement *>::const_iterator vecIterator = mapIterator->second.begin();
  //     while ( vecIterator != mapIterator->second.end() )
  //       {
  //         std::cout << "    " << colorToChar((*vecIterator)->color) << std::endl;
  //         ++vecIterator;
  //       }
  //     ++mapIterator;
  //   }
  // std::cout << std::endl;
  // !t1m-debug

  rowIterator = mpSparseMatrix->begin();
  columnIterator = rowIterator->second.begin();
  return std::pair< Element *, Element * >(rowIterator->first,
                                           *columnIterator );
}


template <class Element>
inline
std::pair< Element *, Element * >
PairContainer<Element>::next()
{
  // t1m-debug
  // std::cout << "mpSparseMatrix in PairContainer::first:\n";
  // std::map<ModelElement *, std::vector<ModelElement *> >::const_iterator mapIterator
  //   = mpSparseMatrix->begin();
  // while ( mapIterator != mpSparseMatrix->end() )
  //   {
  //     std::cout << colorToChar( mapIterator->first->color ) << " overlaps with\n";
  //     std::vector<ModelElement *>::const_iterator vecIterator = mapIterator->second.begin();
  //     while ( vecIterator != mapIterator->second.end() )
  //       {
  //         std::cout << "    " << colorToChar((*vecIterator)->color) << std::endl;
  //         ++vecIterator;
  //       }
  //     ++mapIterator;
  //   }
  // std::cout << std::endl;
  // !t1m-debug

  if ( ++columnIterator == rowIterator->second.end() )
    {
      ++rowIterator;
      columnIterator = rowIterator->second.begin();
    }

  return std::pair< Element *, Element * >(rowIterator->first,
                                           *columnIterator );
}


template <class Element>
inline
bool
PairContainer<Element>::end()
{ return ( rowIterator == mpSparseMatrix->end() ); }


template <class Element>
AxisAlignedBoundingBoxes<Element>::AxisAlignedBoundingBoxes( int dimensions )
  : mDimensions(dimensions)
{
    mpPairs = new PairContainer<Element>( &mSparseMatrix );
}


template <class Element>
void
AxisAlignedBoundingBoxes<Element>::update()
{
  sortXYZ();
  buildXPairs();
  prunePairs();
}


template <class Element>
void
AxisAlignedBoundingBoxes<Element>::add( Element *newElement )
{
  BoundingBox * newBBox = newElement->boundingBox();

  mProjections[0].push_back(&newBBox->xmin);
  mProjections[0].push_back(&newBBox->xmax);
  mObjectsOrdered[0].push_back(newElement);
  mObjectsOrdered[0].push_back(newElement);

  mProjections[1].push_back(&newBBox->ymin);
  mProjections[1].push_back(&newBBox->ymax);
  mObjectsOrdered[1].push_back(newElement);
  mObjectsOrdered[1].push_back(newElement);

  if ( mDimensions > 2 )
    {
      mProjections[2].push_back(&newBBox->zmin);
      mProjections[2].push_back(&newBBox->zmax);
      mObjectsOrdered[2].push_back(newElement);
      mObjectsOrdered[2].push_back(newElement);
    }
}

template <class Element>
void
AxisAlignedBoundingBoxes<Element>::erase( Element *obsoleteElement )
{
  typedef typename std::list<Element *>::iterator ElementIterator;

  BoundingBox * obsoleteBBox = obsoleteElement->boundingBox();

  int i = 0;
  while ( i < mDimensions )
    {
//      std::list<Element *>::iterator remover =
      ElementIterator remover =
          std::find( mObjectsOrdered[i].begin(), mObjectsOrdered[i].end(), obsoleteElement );
      if ( remover == mObjectsOrdered[i].end() )
        return;
//      std::list<Element *>::iterator oldRemover = remover++;
      ElementIterator oldRemover = remover++;
      mObjectsOrdered[i].erase(oldRemover);

      remover = std::find( remover, mObjectsOrdered[i].end(), obsoleteElement );
      mObjectsOrdered[i].erase(remover);

      ++i;
    }

  std::list<double *>::iterator eraser =
    std::find( mProjections[0].begin(), mProjections[0].end(), &obsoleteBBox->xmin );

  std::list<double *>::iterator oldEraser = eraser++;
  mProjections[0].erase( oldEraser );

  eraser = std::find( eraser, mProjections[0].end(), &obsoleteBBox->xmax );
  mProjections[0].erase( eraser );


  // the y-lists
  eraser = std::find( mProjections[1].begin(), mProjections[1].end(), &obsoleteBBox->ymin );

  oldEraser = eraser++;
  mProjections[1].erase( oldEraser );

  eraser = std::find( eraser, mProjections[1].end(), &obsoleteBBox->ymax );
  mProjections[1].erase( eraser );


  if ( mDimensions > 2 )
    {
      // and the z-lists
        eraser = std::find( mProjections[2].begin(), mProjections[2].end(), &obsoleteBBox->zmin );

      oldEraser = eraser++;
      mProjections[2].erase( oldEraser );

      eraser = std::find( eraser, mProjections[2].end(), &obsoleteBBox->zmax );
      mProjections[2].erase( eraser );
    }
}


// sort keys by value, change the order of 'objects' accordingly
// using insertion sort
template <class Element>
void
AxisAlignedBoundingBoxes<Element>::sort( int xyz )
{
  typedef typename std::list<Element *>::iterator ElementIterator;

  std::list<double *>::iterator progressKeys = mProjections[xyz].begin();
  std::list<double *>::iterator walkKeys = progressKeys++;

//  std::list<Element *>::iterator progressObjects
  ElementIterator progressObjects =
      mObjectsOrdered[xyz].begin();
//  std::list<Element *>::iterator walkObjects = progressObjects++;
  ElementIterator walkObjects = progressObjects++;

  // updating the bounding boxes might be expensive:  do it only for xyz==0
  if (!xyz)
    {
      (*walkObjects)->boundingBox();
      (*progressObjects)->boundingBox();
    }

  std::list<double *>::iterator oldKey;
//  std::list<Element *>::iterator oldObj;
  ElementIterator oldObj;

  bool insertAtBegin = false;

  while ( progressKeys != mProjections[xyz].end() )
    {
      while ( **walkKeys > **progressKeys )
        {
          --walkKeys;
          --walkObjects;

          if ( walkKeys == mProjections[xyz].begin() )
            {
              // if this happens, progressKeys' value is the
              //  yet smallest value in the list:
              if ( **walkKeys > **progressKeys )
                  insertAtBegin = true;

              break;
            }
        }

      // oldKey = progressKey; ++progressKey;
      oldKey = progressKeys++;
      oldObj = progressObjects++;

      // if ( !xyz )
      //   (*progressObjects)->boundingBox();

      if ( insertAtBegin )
        {
          insertAtBegin = false;
          mProjections[xyz].push_front( *oldKey );
          mObjectsOrdered[xyz].push_front( *oldObj );
          mProjections[xyz].erase( oldKey );
          mObjectsOrdered[xyz].erase( oldObj );
        }
      else if ( ++walkKeys != progressKeys )
        {
          ++walkObjects;
          mProjections[xyz].insert( walkKeys, *oldKey );
          mObjectsOrdered[xyz].insert( walkObjects, *oldObj );
          mProjections[xyz].erase( oldKey );
          mObjectsOrdered[xyz].erase( oldObj );
        }

      walkKeys = progressKeys;
      --walkKeys;
      walkObjects = progressObjects;
      --walkObjects;
    }
}


template <class Element>
void
AxisAlignedBoundingBoxes<Element>::sortXYZ()
{
  sort(0);

  sort(1);
  if (mDimensions >2)
    sort(2);
}


template <class Element>
void
AxisAlignedBoundingBoxes<Element>::buildXPairs()
{
  typedef typename std::list<Element *>::iterator ElementIterator;
  typedef typename std::list<Element *>::const_iterator ConstElementIterator;

  std::list<Element *> & inputList = mObjectsOrdered[0];
  std::list<double *> & keyPtrList = mProjections[0];

  std::map< Element *, std::vector<Element *> > & outputMap = mSparseMatrix;
  outputMap.clear();

//  std::list<Element *>::const_iterator objProgressor = inputList.begin();
  ConstElementIterator objProgressor = inputList.begin();

  std::list<double *>::const_iterator keyProgressor = keyPtrList.begin();

//  std::list<Element *>::const_iterator objSeeker;
  ConstElementIterator objSeeker;
  std::list<double *>::const_iterator keySeeker;

  //std::vector<Element *> openIntervals;
  //std::vector<Element *>::iterator openIterator;

  //openIntervals.clear();

  while ( objProgressor != inputList.end() )
    {
      // if this is the second occurance of *objProgressor we skip
      if ( *keyProgressor == &(*objProgressor)->boundingBox()->xmax )
        {
          ++keyProgressor; ++objProgressor;
          continue;
        }

      objSeeker = objProgressor;
      keySeeker = keyProgressor;

      while ( true )
        {
          if ( *(++objSeeker) == *objProgressor )
            break;

          ++keySeeker;

          // if this is the second occurance of *objSeeker, skip
          if ( *keySeeker == &(*objSeeker)->boundingBox()->xmax )
            continue;

          // put *progressor into *seekers intersection list
          std::vector<Element *> & inserter = outputMap[*objProgressor];
          if ( std::find( inserter.begin(), inserter.end(), *objSeeker )
               == inserter.end() )
            inserter.push_back(*objSeeker);
        }

      ++objProgressor;
      ++keyProgressor;
    }

  // t1m-debug
  // std::cout << "mSparseMatrix after buildXPairs:\n";
  // std::map<Element *, std::vector<Element *> >::const_iterator mapIterator
  //   = mSparseMatrix.begin();
  // while ( mapIterator != mSparseMatrix.end() )
  //   {
  //     std::cout << colorToChar( mapIterator->first->color ) << " overlaps with\n";
  //     std::vector<Element *>::const_iterator vecIterator = mapIterator->second.begin();
  //     while ( vecIterator != mapIterator->second.end() )
  //       {
  //         std::cout << "    " << colorToChar((*vecIterator)->color) << std::endl;
  //         ++vecIterator;
  //       }
  //     ++mapIterator;
  //   }
  // std::cout << std::endl;
  // !t1m-debug
}


template <class Element>
void
AxisAlignedBoundingBoxes<Element>::prunePairs()
{
  typedef typename std::list<Element *>::const_iterator ConstElementIterator;

  std::map< Element*, std::vector<Element*> > & compareMap =
    mSparseMatrix;

  std::map< Element*, std::vector<Element*> > outputMap;

  // std::vector<Element *> openIntervals;
  // std::vector<Element *>::iterator openIterator;


  int pass = 0;

  while ( ++pass < mDimensions )
    {
      std::list<Element *> & inputList = mObjectsOrdered[pass];
      std::list<double *> &  keyPtrList = mProjections[pass];

      //std::list<Element *>::const_iterator objProgressor = inputList.begin();
      //std::list<Element *>::const_iterator objSeeker;

      ConstElementIterator objProgressor = inputList.begin();
      ConstElementIterator objSeeker;

      std::list<double *>::const_iterator keyProgressor = keyPtrList.begin();
      std::list<double *>::const_iterator keySeeker;

      double * compareCoord;

      outputMap.clear();

      while ( objProgressor != inputList.end() )
        {
          // if this is the second occurance of *objProgressor we skip
          if ( pass == 1 )
            compareCoord =  &(*objProgressor)->boundingBox()->ymax;
          else
            compareCoord =  &(*objProgressor)->boundingBox()->zmax;

          if ( *keyProgressor == compareCoord )
            {
              ++keyProgressor; ++objProgressor;
              continue;
            }

          objSeeker = objProgressor;
          keySeeker = keyProgressor;

          while ( *(++objSeeker) != *objProgressor )
            {
              ++keySeeker;

              // if this is the second occurance of *objSeeker, skip
              if ( pass == 1 )
                compareCoord =  &(*objSeeker)->boundingBox()->ymax;
              else
                compareCoord =  &(*objSeeker)->boundingBox()->zmax;

              if ( *keySeeker == compareCoord )
                continue;

              // put *progressor into *seekers intersection list
              // only if this pair exists in the compareList;
              std::vector<Element *> & compSeeker = compareMap[*objSeeker];
              std::vector<Element *> & compProgressor
                = compareMap[*objProgressor];
              if ( std::find( compSeeker.begin(), compSeeker.end(), *objProgressor )
                   != compSeeker.end()
                   ||
                   std::find( compProgressor.begin(), compProgressor.end(), *objSeeker )
                   != compProgressor.end()
                   )
                {
                  std::vector<Element *> & inserter = outputMap[*objProgressor];

                  if ( std::find( inserter.begin(), inserter.end(), *objSeeker )
                       == inserter.end() )
                    inserter.push_back(*objSeeker);
                }
            }

          ++objProgressor;
          ++keyProgressor;
        }

      // do
      //   {
      //     openIterator = find( openIntervals.begin(), openIntervals.end(), *progressor );
      //     // if progressor was found in the list of openIntervals, we reached the second
      //     //  entry oy *progressor in inputList, and the interval is closed.
      //     if ( openIterator != openIntervals.end() )
      //       {
      //         openIntervals.erase( openIterator );
      //         continue;
      //       }
      //     else
      //       openIntervals.push_back( *progressor );

      //     // if *progressor is in no pair on the x-axis go to the next element
      //     // if ( find( allOverlapping.begin(), allOverlapping.end(), *progressor ) == allOverlapping.end() )
      //     // continue

      //     seeker = progressor;

      //     while ( *(++seeker) != *progressor )
      //       {
      //         openIterator = find( openIntervals.begin(), openIntervals.end(), *seeker );
      //         if ( openIterator != openIntervals.end() )
      //           continue;

      //         // put *progressor into *seekers intersection list
      //         // only if this pair exists in the x list
      //           std::vector<Element *> compSeeker = compareMap[*seeker];
      //           std::vector<Element *> compProgressor = compareMap[*progressor];
      //           if ( find(compSeeker.begin(),compSeeker.end(),*progressor)
      //                    != compSeeker.end()
      //                ||
      //                find(compProgressor.begin(),compProgressor.end(),*seeker)
      //                    != compProgressor.end() )
      //             {
      //               std::vector<Element *> & inserter = outputMap[*progressor];

      //               if ( find( inserter.begin(), inserter.end(), *seeker )
      //                    == inserter.end() )
      //                 inserter.push_back( *seeker );
      //             }
      //       }
      //   }
      // while ( ++progressor != inputList.end() );

      mSparseMatrix.swap(outputMap);
    }

  // t1m-debug
  // std::cout << "mSparseMatrix after prunePairs:\n";
  // std::map<Element *, std::vector<Element *> >::const_iterator mapIterator
  //   = mSparseMatrix.begin();
  // while ( mapIterator != mSparseMatrix.end() )
  //   {
  //     std::cout << colorToChar( mapIterator->first->color ) << " overlaps with\n";
  //     std::vector<Element *>::const_iterator vecIterator = mapIterator->second.begin();
  //     while ( vecIterator != mapIterator->second.end() )
  //       {
  //         std::cout << "    " << colorToChar((*vecIterator)->color) << std::endl;
  //         ++vecIterator;
  //       }
  //     ++mapIterator;
  //   }
  // std::cout << std::endl;
  // !t1m-debug
}


template <class Element>
const PairContainer<Element> &
AxisAlignedBoundingBoxes<Element>::returnPairs()
{
  update();
  return *mpPairs;
}



template <class Element>
PairContainer<Element>::PairContainer( std::map< Element*, std::vector<Element*> > * sparseMatrix )
  : mpSparseMatrix( sparseMatrix )//,
    // mIterator( sparseMatrix )
{
}

#endif // CS_AABB_H
