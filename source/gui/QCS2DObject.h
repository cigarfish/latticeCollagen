
#ifndef QCS_2D_OBJECT_H
#define QCS_2D_OBJECT_H

#include <QGraphicsItem>

class QPointF;
class QColor;


/*!
  \brief Base class for 2D objects

  The base class for 2D objects to be drawn on the 2DDisplay.
  Declares the necessary member functions for 2D drawing on Qt graphicScenes:\n
  <tt>void paint(..)</tt>\n
  <tt>QRect boundingRect() const</tt>\n
  Overloading of \p paint(..) is handled by calling the helper routine \p draw(..)
 */
class QCS2DObject : public QGraphicsItem
{
 public:
  QCS2DObject(QPointF position, QGraphicsItem *parent=0);

  ~QCS2DObject();

  /*!
    \brief Routine to draw the shape of the 2D object
    \param painter Painter to use for actual drawing (see Qt documentation!)
    \param options Style options
    \param parent The parent widget

    This routine is called by the different flavours of \p paint(), which is used internally to
    draw the GraphicsItems on a GraphicsScene.  The \p painter is used to draw the actual shape.
  */
  virtual void draw(QPainter *painter, const QStyleOptionGraphicsItem *options, QWidget *parent=0) =0;
  /*!
    \brief Return the rectangles defined by the extension of the shape
    \returns The rectangle bounding the object's shape

    Qt uses this function internally to determine position and extension of a GraphicsItem.
    Please reimplement this to return the actual bounding rectangle.
  */
  virtual QRectF boundingRect() const =0;

  void paint(QPainter *painter, QStyleOptionGraphicsItem *options, QWidget *parent=0) {draw(painter,options,parent);};
  void paint(QPainter &painter, QStyleOptionGraphicsItem *options, QWidget *parent=0) {draw(&painter,options,parent);};
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *options, QWidget *parent=0){draw(painter,options,parent);};

  void setColor(QColor color) {mColor = color;};

 protected:
  QPointF mPosition;
  QColor mColor;
};

#endif // QCS_2D_OBJECT_H
