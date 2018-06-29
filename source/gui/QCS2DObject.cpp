
#include "QCS2DObject.h"


/*!
  \brief Constructor
  \param position The initial position of the objects center
  \param parent The object's parent graphics item
*/
QCS2DObject::QCS2DObject(QPointF position, QGraphicsItem *parent)
  :QGraphicsItem(parent),
   mPosition(position),
   mColor()
{}


QCS2DObject::~QCS2DObject()
{}
