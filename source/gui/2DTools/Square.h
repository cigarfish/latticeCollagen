
#ifndef QCS_2DTOOLS_SQUARE_H
#define QCS_2DTOOLS_SQUARE_H

#include "../QCS2DObject.h"

class QPointF;
class QColor;
class QPainter;
class QGraphicsItem;
class QStyleOptionGraphicsItem;


/*!
  \brief A 2DObject with the shape of a square
*/
class Square : public QCS2DObject
{
 public:
  Square(QPointF position, qreal sideLength, QGraphicsItem * parent=0);
  ~Square();

  void draw(QPainter *painter, const QStyleOptionGraphicsItem *options, QWidget *parent=0);
  QRectF boundingRect() const;

 private:
  qreal mLengthSide;
};


#endif // QCS_2DTOOLS_SQUARE_H
