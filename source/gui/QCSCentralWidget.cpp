
#include <QVBoxLayout>
#include <QButtonGroup>
#include <QCheckBox>
#include <QMouseEvent>
#include <QPushButton>

//#include "QGLDMainWindow.h"
#include "QCSGLDisplay.h"
#include "QCSCentralWidget.h"
#include "QCS3DTabWidget.h"
#include "etcUIDemo.h"
#include "tabODEs/ODEsTab.h"
#include "tabAbout/About.h"
#include "tabMonolayer/monolayer.h"
#include "tabVascularization/Vascularization.h"
#include "tabComplexCells/tabComplexCells.h"

#if defined( CS_BUILD_IMAGEPROCESSING )
  #include "tabImageProcessing/ImageProcessing.h"
  #include "tabToolsWithItkVtk/toolsWithItkVtk.h"
#endif


/*!
  \brief Constructor

  Setting up tabs and their widgets.
*/
QCSCentralWidget::QCSCentralWidget(QWidget *parent)
    :QTabWidget(parent)
{
    QSizePolicy spolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    setSizePolicy(spolicy);

#ifndef CS_TI_QUANT_ONLY

    QWidget * pageView3D = new QCS3DTabWidget(this);

    /*
     * The 2D View Tab
     */
//   QWidget * pageView2D = new QCS2DTabWidget(this);

    /*
     * The Etc Tab
     *  Featuring a collection of interactive Widgets like buttons, text input fields,
     *  and pull-down menues
     */
    QWidget * pageEtc = new etcUIDemo(this);

    // A tab that lets Stefan test some things
    QWidget * pageViewStefansPlayground = new monolayer(this);

#endif // CS_TI_QUANT_ONLY

#if defined( CS_BUILD_IMAGEPROCESSING )
    // A tab that lets Adrian test some things
    QWidget * pageViewAdriansPlayground = new ImageProcessing(this);
#  if !defined( CS_TI_QUANT_ONLY )
    // A tab that for things with VTK/ITK
    QWidget * pageViewJohannesPlaygroundTwo = new toolsWithItkVtk(this);
#  endif
#endif

    //A tab to display version and contact infos
    QWidget * pageAbout = new About(this);

#ifndef CS_TI_QUANT_ONLY
    
    // A tab that lets Geraldine test some things
    QWidget * pageViewGeraldinePlayground = new ODEsTab(this);

	QWidget * pageViewNickPlayground = new Vascularization(this);
	QWidget * pageViewComplexCells = new complexCells(this);

    // Adding all the tabs to the Widget:
    addTab(pageView3D, "3D");
//    addTab(pageView2D, "2D");
    addTab(pageEtc, "etc.");

    addTab(pageViewGeraldinePlayground, "ODEs (Geraldine)");

    addTab(pageViewStefansPlayground, "CellsSpherical (Stefan)");
#endif // CS_TI_QUANT_ONLY

#if defined( CS_BUILD_IMAGEPROCESSING )
    addTab(pageViewAdriansPlayground, "ImageProcessing");
#  if !defined( CS_TI_QUANT_ONLY )
    addTab(pageViewJohannesPlaygroundTwo, "Tools with Itk/Vtk");
#  endif
#endif

#ifndef CS_TI_QUANT_ONLY
    addTab(pageViewComplexCells, "Complex Cells");
    addTab(pageViewNickPlayground, "Vascularization");
#endif // CS_TI_QUANT_ONLY

    addTab(pageAbout, "About");
}
 

QCSCentralWidget::~QCSCentralWidget()
{}
