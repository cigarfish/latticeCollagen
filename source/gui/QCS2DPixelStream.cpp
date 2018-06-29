#ifndef __BUILD_WINDOWS__

#include "QCS2DPixelStream.h"
#include "2DTools/Pixel.h"

#include <QGraphicsView>
#include <QPair>

#include "QCS2DDisplay.h"


/*!
  \brief Constructor
  \param scene The GraphicsScene to draw on
*/
QCS2DPixelStream::QCS2DPixelStream( QGraphicsScene *scene )
  : mScene(scene),
    mCurrentColor("Black")
{
  mColorQueue.push_back(mCurrentColor);
}


/*!
  \brief Constructor
  \param view The GraphicsView to use for drawing
*/
QCS2DPixelStream::QCS2DPixelStream( QGraphicsView *view )
  : mScene(view->scene()),
    mCurrentColor("Black")
{
  mColorQueue.push_back(mCurrentColor);
}


/*!
  \brief Copy constructor
  \param stream Stream to copy
*/
QCS2DPixelStream::QCS2DPixelStream( QCS2DPixelStream &stream )
  : mScene(stream.mScene),
    mCurrentColor(stream.mCurrentColor)
{
  mPixelBuffer = std::deque<QPoint>(stream.mPixelBuffer);
  mColorQueue  = std::deque<QColor>(stream.mColorQueue);
  mColorChangeQueue = std::deque< std::deque<QPoint>::iterator >(stream.mColorChangeQueue);
}


/*!
  \brief Overloaded left shift operator
  \param point The position of the pixel to add
  \returns The resulting 2DPixelStream
*/
QCS2DPixelStream
QCS2DPixelStream::operator << ( QPoint point )
{
  mPixelBuffer.push_back(point);
  return *this;
}


/*!
  \brief Overloaded left shift operator
  \param point The position of the pixel to add
  \returns The resulting 2DPixelStream
*/
QCS2DPixelStream
QCS2DPixelStream::operator << ( QPair<int,int> point )
{
  mPixelBuffer.push_back(QPoint(point.first,point.second));
  return *this;
}


/*!
  \brief Overloaded left shift operator
  \param point The position of the pixel to add
  \returns The resulting 2DPixelStream
*/
QCS2DPixelStream
QCS2DPixelStream::operator << ( std::pair<int,int> point )
{
  mPixelBuffer.push_back(QPoint(point.first,point.second));
  return *this;
}


/*!
  \brief Overloaded left shift operator
  \param stream A 2DPixelStream to append
  \returns The resulting 2DPixelStream
*/
QCS2DPixelStream
QCS2DPixelStream::operator << ( QCS2DPixelStream stream )
{
  std::deque<QPoint>::iterator it = stream.mPixelBuffer.begin();

  for ( ; it != stream.mPixelBuffer.end(); ++it)
    mPixelBuffer.push_back(*it);

  return *this;
}


/*!
  \brief Overloaded left shift operator to change color
  \param color The color to draw pixels added after using this flavor of operator
  \returns The resulting 2DPixelStream
*/
QCS2DPixelStream
QCS2DPixelStream::operator << (QColor color)
{
  mColorChangeQueue.push_back(mPixelBuffer.end());
  mCurrentColor = color;
  mColorQueue.push_back(mCurrentColor);

  return *this;
}


/*!
  \brief Paint the pixels and forget the pixels added so far

  Pops all pixels after drawing to the GraphicsScene.  Will still leave the color as set before.
*/
void
QCS2DPixelStream::flush()
{
  // finalizing the queue
  mColorChangeQueue.push_back(mPixelBuffer.end());

  std::deque<QPoint>::iterator pointIter = mPixelBuffer.begin();

  while ( mColorChangeQueue.size() )
    {
      QColor color = mColorQueue.front();
      mColorQueue.pop_front();

      std::deque<QPoint>::iterator colChange = *(mColorChangeQueue.begin());

      for ( ; pointIter != colChange; ++pointIter )
	{
	  QPoint point = *pointIter;
	  mScene->addItem(new Pixel(point,color));
	  mPixelBuffer.pop_front();
	}

      mColorChangeQueue.pop_front();
    }
}

#endif // __BUILD_WINDOWS__
