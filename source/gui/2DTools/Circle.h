
#ifndef QCS_2DTOOLS_CIRCLE_H
#define QCS_2DTOOLS_CIRCLE_H

#include "../QCS2DObject.h"

class QPointF;
class QColor;
class QPainter;
class QGraphicsItem;
class QStyleOptionGraphicsItem;


/*!
  \brief A 2DObject with the shape of a circle
*/
class Circle : public QCS2DObject
{
 public:
  Circle(QPointF position, qreal radius, QGraphicsItem * parent=0);
  ~Circle();

  void draw(QPainter *painter, const QStyleOptionGraphicsItem *options, QWidget *parent=0);
  QRectF boundingRect() const;

 private:
  qreal mRadius;
};

#endif // QCS_2DTOOLS_CIRCLE_H
