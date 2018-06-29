
#include "../2DTools/Circle.h"

#include <QPainter>


/*!
  \brief Constructor
  \param position The inital position of the circle's center
  \param radius The circle's radius
  \param parent The parent GraphicsItem
*/
Circle::Circle(QPointF position, qreal radius, QGraphicsItem * parent)
  : QCS2DObject(position,parent),
    mRadius(radius)
{}


Circle::~Circle()
{}


/*!
  \brief The draw routine

  Reimlementing the base class's \p draw(..) routine to draw the actual circle.
*/
void
Circle::draw(QPainter *painter, const QStyleOptionGraphicsItem * /*options*/, QWidget * /*parent*/)
{
  painter->setRenderHint(QPainter::Antialiasing);
  painter->setPen(mColor);
  painter->drawEllipse(boundingRect());
}


/*!
  \brief The bounding rectangle.
  \returns The rectangle defined by bounding the circle

  This is used by the 2DDisplay/GraphicsScene to define the position and extension of the drawn shape.
*/
QRectF
Circle::boundingRect() const
{
  return QRectF(mPosition.x()-mRadius, mPosition.y()-mRadius, 2*mRadius, 2*mRadius);
}

