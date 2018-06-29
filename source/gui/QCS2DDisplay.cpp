
#include "QCS2DDisplay.h"

#include <QGraphicsView>
#include <QMouseEvent>
#include <QGLContext>

#include "QVTKWidget2.h"
#include "vtkRenderer.h"
#include "vtkGenericOpenGLRenderWindow.h"

#include "QCS2DImageItem.h"

#include "2DTools/Circle.h"
#include "2DTools/Square.h"
#include "2DTools/Pixel.h"

#include "QCS2DPixelStream.h"

#include <math.h>
#include <iostream>


#define max(x,y) (x>y)?(x):(y)


/*!
  \brief Constructor
  \param parent The parent widget
*/
QCS2DDisplay::QCS2DDisplay(QWidget * parent)
  : QGraphicsView(parent),
    lastpos(QPoint(0,0)),
    mCenter(QPointF(0,0)),
    mScale(1.)
{
  mpGLContext = new QGLContext(QGLFormat());
  mpQVTKViewport = new QVTKWidget2(mpGLContext);
  setViewport(mpQVTKViewport);
  setViewportUpdateMode(QGraphicsView::FullViewportUpdate);

  mpArena = new QGraphicsScene();
  setScene(mpArena);

  mpRenderer = vtkRenderer::New();
  mpRenderer->ResetCamera();
  mpQVTKViewport->GetRenderWindow()->AddRenderer(mpRenderer);
  mpQVTKViewport->GetRenderWindow()->SetSwapBuffers(1);
  mpQVTKViewport->setAutoBufferSwap(true);
  mpQVTKViewport->GetRenderWindow()->Render();

  resize(800,600);

  // Circle *circle = new Circle(QPointF(0.,20.),40.);
  // Circle *circle2 = new Circle(QPointF(-30.,-20.),30.);
  // Square *square = new Square(QPointF(0.,20.),80.);
  // Pixel *pixel = new Pixel(QPointF(0.,20.));
  // mpArena->addItem(circle);
  // mpArena->addItem(circle2);
  // mpArena->addItem(square);
  // mpArena->addItem(pixel);

  // circle2->setColor(QColor("Blue"));

  // QCS2DPixelStream stream(mpArena);
  // stream << QColor("Green");
  // for ( int i=0; i<50; ++i )
  //   stream << QPoint(0.,i);
  // stream << QColor("Blue");
  // for ( int i=50; i<100; ++i )
  //   stream << QPoint(0.,i);
  // stream << QColor("Red");
  // for ( int i=100; i<150; ++i )
  //   stream << QPoint(0.,i);

  // stream.flush();

  // setAlignment(0);

  // mpQVTKViewport->update();

  // int newSceneWidth = max( sceneRect().width(), width() );
  // int newSceneHeight = max( sceneRect().height(), height() );

  // QRect newSceneRect = QRect(-newSceneWidth/2,-newSceneHeight/2,newSceneWidth,newSceneHeight);

  // setSceneRect(newSceneRect);
  // mpQVTKViewport->resize(newSceneWidth,newSceneHeight);
}


QCS2DDisplay::~QCS2DDisplay()
{}


/*!
  \brief Add a 2DObject to the scene
*/
void
QCS2DDisplay::addObject(QCS2DObject *object)
{
  mpArena->addItem(object);
}


/*!
  \brief Open and display an image in the current layer
  \param fileName Name of the image file to open
*/
void
QCS2DDisplay::openImage(QString fileName)
{
  QCS2DImageItem * image = new QCS2DImageItem(mpGLContext);
 
  mpArena->addItem(image);

  image->loadFromFile(fileName.toStdString());

  image->setZValue(-10.);

  image->display();
  
  setSceneRect(mpArena->itemsBoundingRect());

  mCenter = sceneRect().center();
  centerOn(mCenter);
}


/*!
  \brief Reimplementation of the mouse press event handler
  \param event The given event

  Memorizes the mouse position for further processing in mouseMoveEvent().
*/ 
void
QCS2DDisplay::mousePressEvent(QMouseEvent *event)
{
  lastpos = event->pos();
}


/*!
  \brief Reimplementation of the mouse move event handler
  \param event The given event

  Handles translations when a mouse move occurs with left mouse button pressed.
  Memorizes the last mouse pointer position for further calls to mouseMoveEvent().
*/
void
QCS2DDisplay::mouseMoveEvent(QMouseEvent *event)
{
  if (event->buttons() & Qt::LeftButton)
    {
      int dx = event->x() - lastpos.x();
      int dy = event->y() - lastpos.y();

      mCenter -= QPointF(dx/mScale,dy/mScale);
      centerOn(mCenter);

      lastpos = event->pos();

      //QPointF upperleft = mapToScene(QPoint(0,0));
      //      std::cout << "Upper left corner has coordinates:  ( " << (int)upperleft.x() << ", " << (int)upperleft.y() << " )\n";
    }
}


/*!
  \brief Reimplementation of the mouse wheel event handler
  \param event The given event

  Translates wheel events into zooms.
*/
void
QCS2DDisplay::wheelEvent(QWheelEvent *event)
{
  qreal dist = ((qreal)event->delta())/1000.;

  if (dist <= -1.)
    dist += 1e-8;

  qreal nscale = 1/(1+dist);

  scale(nscale,nscale);

  mScale *= nscale;
  mCenter = mapToScene(QRect(0,0,width(),height()).center());
		       //sceneRect().center();
}


/*!
  \brief Reimplementation of the mouse double click event handler
  \param event The given event

  Recenters the viewport when zoomed.
*/
void
QCS2DDisplay::mouseDoubleClickEvent(QMouseEvent * event)
{
  centerOn(event->x(),event->y());
  mCenter = QPointF(0,0);
}


/*!
  \brief Reimplementation of the key press event handler
  \param event The given event

  Pressing Escape will rescale the scene to initial layout.
*/
void
QCS2DDisplay::keyPressEvent(QKeyEvent *event)
{
  if (event->key() == Qt::Key_Escape )
    {
      scale(1/mScale,1/mScale);
      mScale = 1.;
    }
  else
    {
      QGraphicsView::keyPressEvent(event);
    }
}


void
QCS2DDisplay::resizeEvent(QResizeEvent *event)
{
  // give the same size to the scene that his widget has
  if (scene())
    scene()->setSceneRect(QRect(QPoint(0, 0), event->size()));
  QGraphicsView::resizeEvent(event);
  mpQVTKViewport->GetRenderWindow()->SetSize(event->size().width(), event->size().height());
  mpQVTKViewport->GetRenderWindow()->Render();
}
