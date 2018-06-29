#ifndef QCS_GL_DISPLAY_H
#define QCS_GL_DISPLAY_H

class CSGLArena;

// Link to Backend
#include "../Core.h"


#if defined(CS_BUILD_IMAGEPROCESSING)
#  define  GLWIDGET QVTKWidget2
#  include "QVTKWidget2.h"
#else
#  define  GLWIDGET QGLWidget
#  include <QGLWidget>
#endif


/*!
  \brief A display for GLObjects
*/
class QCSGLDisplay : public GLWIDGET
{
    Q_OBJECT

public:
    QCSGLDisplay(QWidget *parent =0);
    ~QCSGLDisplay();

    /*!
      \brief Set the GLArena to be displayed
      \param arena The new GLArena
      This function will NOT delete the old arena, because that may be used by
      another model, or might be hidden.
    */
    void setArena(CSGLArena *arena) { mpObjects = arena; };

    void removeModel() { setArena(mpDefaultArena); /* deregister display? */};


    enum Color {
        ColorBlack =0,
        ColorDarkGray,
        ColorLightGray,
        ColorMauve,
        ColorWhite
    };

    QColor mBackgroundColor;
    void setBackgroundColor( QCSGLDisplay::Color, bool init=false );


    /*!
      \brief Set the Model the data of which is to be displayed
      \param model The new model
    */
    void setModel( CSModel * model );

    /*!
      \brief Get the GLArena displayed
      \returns The GLArena that is displayed
    */
    CSGLArena *arena() const
    {
        return mpObjects;
    };

public slots:
    void update()
    { emit updated(); GLWIDGET::update(); };

signals:
    // mainly for dependent displays, i.e. which should have the same arena:
    void updated();
    void arenaChanged();

protected:
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);

    void setXRotation(int angle);
    void setYRotation(int angle);
    void setDistance(double distance);

    void setXTranslation(double xSkip);
    void setYTranslation(double ySkip);


private:
    CSGLArena *mpObjects;

    // empty fall-back arena
    CSGLArena *mpDefaultArena; 

    int mXRot;
    int mYRot;
    int mZRot;
    double mDistance;

    QPoint lastpos;
    QPointF mOrigin;
};


#endif // QCS_GL_DISPLAY_H
