
#ifndef QCS_2D_DISPLAY_H
#define QCS_2D_DISPLAY_H

#include<QGraphicsView>

class QGraphicsScene;
class QCS2DObject;
class QGLContext;

class QVTKWidget2;
class vtkRenderer;

/*!
  \brief A display for 2D shapes
*/
class QCS2DDisplay : public QGraphicsView
{
  Q_OBJECT

 public:
  QCS2DDisplay(QWidget *parent=0);

  ~QCS2DDisplay();

  //void show();
  void addObject(QCS2DObject *);

  void openImage(QString fileName);

  void setRenderer(vtkRenderer* renderer) { mpRenderer = renderer; };
  vtkRenderer* getRenderer() const { return mpRenderer; }

  QGraphicsScene* getScene() const { return mpArena; }

  QGLContext* getGLContext() const { return mpGLContext; }

 protected:
  void mousePressEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void mouseDoubleClickEvent(QMouseEvent *event);
  void wheelEvent(QWheelEvent *event);
  void keyPressEvent(QKeyEvent *event);
  void resizeEvent(QResizeEvent *event);

  QVTKWidget2* mpQVTKViewport;
  vtkRenderer* mpRenderer;
  QGraphicsScene* mpArena;
  QGLContext* mpGLContext;

  QPoint lastpos;
  QPointF mCenter;
  qreal mScale;
};

#endif // QCS_2D_DISPLAY_H
