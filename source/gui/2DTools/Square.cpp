
#include "../2DTools/Square.h"

#include <QPainter>


/*!
  \brief Constructor
  \param position The position of the center of the square on the 2DDisplay
  \param sideLength Length of each side of the square
  \param parent The parent GraphicsItem
*/
Square::Square(QPointF position, qreal sideLength, QGraphicsItem * parent)
  : QCS2DObject(position,parent),
    mLengthSide(sideLength)
{}


Square::~Square()
{}


/*!
  \brief The draw routine

  Reimplementing the base class's \p draw(..) routine to draw the actual square.
*/
void
Square::draw(QPainter *painter, const QStyleOptionGraphicsItem * /*options*/, QWidget * /*parent*/)
{
  painter->setRenderHint(QPainter::Antialiasing);
  painter->setPen(mColor);
  painter->drawRect(boundingRect());
}


/*!
  \brief The bounding rectangle
  \return The rect defined by position and extension of the square on the 2DDisplay's GraphicsScene

  This is used by the 2DDisplay/GraphicsScene to define the position and size of the drawn shape.
*/
QRectF
Square::boundingRect() const
{
  return QRectF(mPosition.x()-mLengthSide/2,mPosition.y()-mLengthSide/2,mLengthSide,mLengthSide);
}
