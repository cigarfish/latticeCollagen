///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Vascularization.cpp                                                  //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes@neitsch.de>                               //
//    Created:  2012-12-20                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "Vascularization.h"

#include <QVector3D>

#include "../GLTools/CSGLSphere.h"
#include "../GLTools/CSGLVascularization.h"
#include "../QCSGLDisplay.h"
#include "../CSGLArena.h"
#include "../QCSSimulationThread.h"
#include "../../model/BasicDatatypes/Color.h"
#include "../../model/Model/ModelVascularization/CSVesselGraph.h"
#include "../../model/Model/ModelVascularization/ModelVascularization.h"
#include "../../model/Model/ModelVascularization/Tumor.hpp"

#if defined( CS_BUILD_IMAGEPROCESSING )
#  include "vtkImageActor.h"
#  include "vtkGraphLayoutView.h"
#  include "vtkImageViewer2.h"
#  include "vtkGenericOpenGLRenderWindow.h"
#  include "vtkRenderer.h"
#  include "vtkRenderWindowInteractor.h"
#endif


// Constructor
Vascularization::Vascularization(QWidget *parent) : QWidget(parent),
  vascularization(NULL),           // will be initialized in init()
  mpSimulationThread(NULL), // will be initialized in init()
  mReset(false){

  // Setting up the Qt Designer code
  // Connecting SIGNAL valueChanged() of the horizontal slider
  setupUi(this);

  connect( initButton, SIGNAL(clicked(bool)), this, SLOT(init()) );
  connect( loadButton, SIGNAL(clicked(bool)), this, SLOT(load()) );
  connect( playButton, SIGNAL(clicked(bool)), this, SLOT(play()) );
  connect( showButton, SIGNAL(clicked(bool)), this, SLOT(createNewDisplay()) );
}
Vascularization::~Vascularization()
{}


void Vascularization::createNewDisplay()
{
  if( vascularization == NULL) return;

    QCSGLDisplay * newDisp = new QCSGLDisplay();
    newDisp->setModel(vascularization);
    newDisp->show();
    connect( mpSimulationThread, SIGNAL(updateStep(int)), newDisp, SLOT(update()) );
    mDisplays.push_back( newDisp );
}

void Vascularization::play()
{
  if( vascularization == NULL) return;

this->vascularization->Run();

  
}

void Vascularization::init(){
  
  if( vascularization != NULL) return;

  vascularization = new Model_Vascularization();

  createNewDisplay();
  vascularization->enableSimulation = true;

  if ( ! mpSimulationThread ) {
    mpSimulationThread = new QCSSimulationThread();
    connect( mpSimulationThread, SIGNAL( finished() ), this, SLOT(stateInitialized()) );
    connect( mpSimulationThread, SIGNAL( updateStep(int) ), this, SLOT(updateProgress(int)) );
  }

  mpSimulationThread->setModel(vascularization);
}

void Vascularization::load(){

  if( vascularization == NULL) return;


  this->vascularization->SetupSimulation();



}
