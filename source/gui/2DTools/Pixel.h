
#ifndef QCS_2DTOOLS_PIXEL_H
#define QCS_2DTOOLS_PIXEL_H

#include "../QCS2DObject.h"

#include <QPainter>


/*!
  \brief A pixel for the 2DDisplay
 */
class Pixel : public QCS2DObject
{
 public:
  Pixel( int x, int y, QColor color=QColor("Black"), QGraphicsItem *parent=0 );
  Pixel( QPointF position, QColor color=QColor("Black"), QGraphicsItem *parent=0 );

  ~Pixel();

  void draw(QPainter *painter, const QStyleOptionGraphicsItem *options, QWidget *parent=0);
  QRectF boundingRect() const;
};

#endif // QCS_2DTOOLS_PIXEL_H
