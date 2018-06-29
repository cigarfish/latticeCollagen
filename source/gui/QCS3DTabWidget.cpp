
#pragma region Includes

#include "QCS3DTabWidget.h"
#include "QCSGLDisplay.h"
#include "QCSMainWindow.h"
#include "CSGLArena.h"

#include <QGridLayout>
#include <QButtonGroup>
#include <QCheckBox>
#include <QPushButton>
#include <QShortcut>
#include <QKeySequence>

// Backend
#include "QCSSimulationThread.h"
#include "QDebugStream.h"

#include "../tabMonolayer/monolayer.h"
#include "../model/Model/ModelCellsSpherical/ModelCellsSpherical.h"

#pragma endregion

/*!
  \brief Conststructor

  Adding PushButton to create and display a GLDisplay.
  Adding a group of QCheckBoxes to color objects.
*/
QCS3DTabWidget::QCS3DTabWidget(QWidget *parent)
    : QWidget(parent),
      mColor(QColor("Green"))
{
    setupUi(this);

    connect( newDisplayButton, SIGNAL(clicked()), this, SLOT(createNewDisplay()) );

    #pragma region Keyboard shortcuts

    QShortcut * shortcut = new QShortcut(QKeySequence(QObject::tr("s")), this );
    connect( shortcut, SIGNAL(activated()), this, SLOT(startSim100()));
    shortcut = new QShortcut(QKeySequence(Qt::SHIFT + Qt::Key_S), this );
    connect( shortcut, SIGNAL(activated()), this, SLOT(startSim1()));
    // shortcut = new QShortcut(QKeySequence(Qt::Key_R), this );
    // connect( shortcut, SIGNAL(activated()), this, SLOT(resetSim()));

    connect( mpDisplay, SIGNAL(arenaChanged()), this, SLOT(arenaUpdate()) );

    #pragma endregion

    mpDisplay->setObjectName("Main Display");
}


QCS3DTabWidget::~QCS3DTabWidget()
{
}

#pragma region Methods for Keyboard shortcuts

void
QCS3DTabWidget::startSim1()
{
	// No function without model init
	if (!core->models["ModelCellsSpherical"]) return;

    (* (ModelCellsSpherical *)core->models["ModelCellsSpherical"]).Simulate(1);

    #pragma region Cell staining

    // core->modelCellsSpherical->UpdateCellsStaining(comboBoxCellStaining->currentIndex());

    #pragma endregion


    mpDisplay->update();
}

void
QCS3DTabWidget::startSim100()
{
	// No function without model init
	if (!core->models["ModelCellsSpherical"]) return;

    (* (ModelCellsSpherical *)core->models["ModelCellsSpherical"]).Simulate(100);

    #pragma region Cell staining

    // core->modelCellsSpherical->UpdateCellsStaining(comboBoxCellStaining->currentIndex());

    #pragma endregion

    mpDisplay->update();
}

#pragma endregion

/*!
  \brief Create a new GL display
 */
void
QCS3DTabWidget::createNewDisplay()
{
    QCSGLDisplay * newDisp = new QCSGLDisplay();
    newDisp->show();
    mpDisplay->arena()->registerDisplay( newDisp );
    newDisp->setArena( mpDisplay->arena() );
    connect( mpDisplay, SIGNAL(updated()), newDisp, SLOT(update()) );

    mOtherDisplays.push_back(newDisp);
}


/*!
  \brief Update dependent displays' arenas
 */
void
QCS3DTabWidget::arenaUpdate()
{
  CSGLArena * newArena = mpDisplay->arena();

  std::vector<QCSGLDisplay *>::iterator dispIt;
  for ( dispIt = mOtherDisplays.begin(); dispIt != mOtherDisplays.end(); ++dispIt )
    {
      newArena->registerDisplay( *dispIt );
      (*dispIt)->setArena( newArena );
    }
}
