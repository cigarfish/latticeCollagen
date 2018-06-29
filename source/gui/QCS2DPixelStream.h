#ifndef __BUILD_WINDOWS__

#ifndef QCS_2D_PIXELSTREAM_H
#define QCS_2D_PIXELSTREAM_H

#include <deque>
#include <utility>

#include <QPair>
#include <QColor>

class QCS2DDisplay;
class QGraphicsView;
class QGraphicsScene;
class QPoint;
class QColor;


/*!
  \brief Stream imlementation for drawing pixels on a GraphicsView, GraphicsScene or 2DDisplay
*/
class QCS2DPixelStream
{
 public:
  QCS2DPixelStream( QGraphicsScene *scene );
  QCS2DPixelStream( QGraphicsView *view );
  QCS2DPixelStream( QCS2DPixelStream &stream );

  QCS2DPixelStream operator << ( QPoint point );
  QCS2DPixelStream operator << ( QPair<int,int> point );
  QCS2DPixelStream operator << ( std::pair<int,int> point );
  QCS2DPixelStream operator << ( QCS2DPixelStream );
  QCS2DPixelStream operator << ( QColor );

  void flush();

 protected:
  std::deque<QPoint> mPixelBuffer;
  std::deque<QColor> mColorQueue;
  std::deque< std::deque<QPoint>::iterator > mColorChangeQueue;

  QGraphicsScene *mScene;

  QColor mCurrentColor;
};

#endif // QCS_2D_PIXELSTREAM_H

#endif // __BUILD_WINDOWS__
