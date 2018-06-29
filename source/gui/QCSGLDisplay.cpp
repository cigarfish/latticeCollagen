
#include <QMouseEvent>
#include <QWheelEvent>
#include <QVector3D>
#include <iostream>

#include "QCSGLDisplay.h"

#include "CSGLArena.h"

#include "../model/Model/CSModel.h"

#define FULL_CIRCLE 360*16

#if defined(__BUILD_MAC__)
# include <GLUT/glut.h>
#else
# include <GL/glut.h>
#endif


/*!
  \brief Normalize an angle
  \param angle The angle to be normalized

  keep the angle between 0 and 2*PI.
*/

static void
qNormalizeAngle(int & angle)
{
    while ( angle < 0 )
        angle += FULL_CIRCLE;
    while (angle > FULL_CIRCLE)
        angle -= FULL_CIRCLE;

}

/*!
  \brief Constructor
  \param parent The widget's parent widget

  Initializing the list of objects with a GLArena.
  To date the display is populated with two spheres and one cube (which covers one sphere).
*/
QCSGLDisplay::QCSGLDisplay(QWidget *parent) :
    GLWIDGET(QGLFormat(QGL::SampleBuffers),parent),
    mOrigin(0,0)
{
    mXRot = 0;
    mYRot = 0;
    mZRot = 0;
    mDistance = 10;

    mpDefaultArena = new CSGLArena();
    mpObjects = mpDefaultArena;
    setBackgroundColor(QCSGLDisplay::ColorMauve, true);
}


QCSGLDisplay::~QCSGLDisplay()
{}



void
QCSGLDisplay::setModel( CSModel * model )
{
  CSGLArena *arena = model->arena();
  arena->registerDisplay(this);
  setArena(arena);
  emit arenaChanged();
}


/*!
  \brief Callback for mouse button press events
  \param event The given event

  Memorize position for processing move events.
*/
void
QCSGLDisplay::mousePressEvent(QMouseEvent *event)
{
    lastpos = event->pos();
}


/*!
  \brief Callback for mouse move events (with pressed button)
  \param event The given event

  1. Calculate the move distance.\n
  2. If\n
     a) the left button is pressed, set the rotation according to move distance;\n
     b) the right button is pressed, set the translation according to move distance.\n
  3. Memorize position for processing further mouse move events.
 */
void
QCSGLDisplay::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - lastpos.x();
    int dy = event->y() - lastpos.y();

    if (event->buttons() & Qt::LeftButton )
    {
        /*rotate*/
        setXRotation(mXRot + 8*dy);
        setYRotation(mYRot + 8*dx);
    }
    else if (event->buttons() & Qt::RightButton )
    {
        setXTranslation(mOrigin.x() - (float)dx/100.);
        setYTranslation(mOrigin.y() + (float)dy/100.);
    }

    lastpos = event->pos();
}


/*!
  \brief Callback for mouse wheel events
  \param event The given event

  Set the viewer's distance according to wheel move
 */
void
QCSGLDisplay::wheelEvent(QWheelEvent * event)
{
    double distance = mDistance + (double)( event->delta() )/120.;
    if(distance < 0)
        distance = 0;
    setDistance(distance);
    update();
}


void
QCSGLDisplay::setBackgroundColor( QCSGLDisplay::Color color, bool init )
{
    switch ( color )
    {
    case ColorBlack:
        mBackgroundColor = QColor( Qt::black );
        break;
    case ColorDarkGray:
        mBackgroundColor = QColor( Qt::darkGray );
        break;
    case ColorLightGray:
        mBackgroundColor = QColor( Qt::lightGray );
        break;
    default:
    case ColorMauve:
        mBackgroundColor = QColor::fromRgb(210,210,245,255);
        break;
    case ColorWhite:
        mBackgroundColor = QColor( Qt::white );
        break;
    }

    if ( ! init )
        qglClearColor(mBackgroundColor);
}


/*!
  \brief Initialize GL rendering options
 */
void
QCSGLDisplay::initializeGL()
{
    //glutInitDisplayMode(GLUT_DEPTH | GLUT_RGB);
    // qglClearColor(QColor::fromCmykF(0.39,0.39,0.39,0.39).dark());
    //mBackgroundColor = QColor::fromRgb(210,210,245,0);
    qglClearColor(mBackgroundColor);
    // enable everything at the beginning
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    // the following is hardware dependent, supplied in glext.h
#ifdef GL_MULTISAMPLE
    glEnable(GL_MULTISAMPLE);
#endif
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);

    static GLfloat lightPosition[4] = { 0.5, 5.0, 7.0, 0.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

    glEnable(GL_LIGHT0);

    resizeGL(width(),height());
    //gluPerspective( 30, width()/(float)height(), 0.001, 1000 );
    setAutoBufferSwap(true);
}


/*!
  \brief Callback for paint events

  This is called, whenever the GLArena is redrawn.
  Setting Scene's view angle and distance.  Redrawing all objects.
 */
void
QCSGLDisplay::paintGL()
{
    glClear( GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );
    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();

    gluPerspective(30, width()/(float)height(),0.001, 200);
    gluLookAt(0,0,mDistance,mOrigin.x(),mOrigin.y(),0,0,1,0);

    glMatrixMode(GL_MODELVIEW);

    glPushMatrix();

    glRotatef(mXRot/16., 1., 0., 0.);
    glRotatef(mYRot/16., 0., 1., 0.);
    glRotatef(mZRot/16., 0., 0., 1.);
    //glMatrixMode(GL_MODELVIEW);
    mpObjects->draw();

    glPopMatrix();
}

/*!
  \brief Callback for resize events
  \param width The viewport's new width
  \param height The viewport's new height

  Handles how the viewport changes when the scene's window is resized.
 */
void
QCSGLDisplay::resizeGL(int width, int height)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, width, height);
}


/*!
  \brief Set the x rotation variable
  \param angle The angle of the scene's rotation around the x axis through the origin
 */
void
QCSGLDisplay::setXRotation(int angle)
{
    qNormalizeAngle(angle);
    if ( angle != mXRot )
    {
        mXRot = angle;
        updateGL();
    }
}


/*!
  \brief Set the y rotation variable
  \param angle The angle of the scene's rotation around the y axis through the origin
 */
void
QCSGLDisplay::setYRotation(int angle)
{
    qNormalizeAngle(angle);
    if ( angle != mYRot )
    {
        mYRot = angle;
        updateGL();
    }
}


/*!
  \brief Set the PoV's distance from the origin
  \param distance The Distance of the PoV from the scene's origin
 */
void
QCSGLDisplay::setDistance(double distance)
{
    if (distance != mDistance )
        //    if (distance > 5)
    {
        mDistance = distance;
        updateGL();
    }
}


/*!
  \brief Set the x variable of the transpose transformation
  \param xSkip The distance the scene has to be tranposed in x direction
 */
void
QCSGLDisplay::setXTranslation(double xSkip)
{
    if ( xSkip != mOrigin.x() )
    {
        mOrigin.setX(xSkip);
        updateGL();
    }
}


/*!
  \brief Set the x variable of the transpose transformation
  \param ySkip The distance the scene has to be tranposed in y direction
 */
void
QCSGLDisplay::setYTranslation(double ySkip)
{
    if ( ySkip != mOrigin.y() )
    {
        mOrigin.setY(ySkip);
        updateGL();
    }
}
