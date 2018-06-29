
#include "../2DTools/Pixel.h"


/*!
  \brief Constructor
  \param x X coordinate of the pixel's position
  \param y Y coordinate of the pixel's position
  \param color Color of the pixel
  \param parent The parent GraphicsItem
*/
Pixel::Pixel( int x, int y, QColor color, QGraphicsItem *parent )
  :QCS2DObject(QPointF(x,y),parent)
{
  setColor(color);
}


/*!
  \brief Constructor
  \param position The position of the pixel
  \param color Color of the pixel
  \param parent The parent GraphicsItem
*/
Pixel::Pixel( QPointF position, QColor color, QGraphicsItem *parent )
  : QCS2DObject(position,parent)
{
  setColor(color);
}


Pixel::~Pixel()
{}


/*!
  \brief The draw routine

  Reimplementing the base class's \p draw() routine to draw the actual pixel.
*/
void
Pixel::draw(QPainter *painter, const QStyleOptionGraphicsItem * /*options*/, QWidget * /*parent*/)
{
  painter->setPen(mColor);
  painter->setBrush(mColor);
  painter->drawRect(boundingRect());
}


/*!
  \brief The bounding rectangle
  \returns A dimensionless rectangle at the pixel's position

  This is used by the 2DDisplay/GraphicsScene to define the position and extension of the drawn shape.
*/
QRectF
Pixel::boundingRect() const
{
  return QRectF(mPosition,mPosition+QPointF(1,1));
}
