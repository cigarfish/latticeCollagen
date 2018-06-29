///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  complexCells.h                                                       //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2013-10-24                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

//#include <iostream>

#include "tabComplexCells.h"
#include "../QCSGLDisplay.h"// Backend
#include "../QCSParameterDelegate.h"// Backend
#include "../QCSParameterModel.h"// Backend
#include "../QCSSimulationThread.h"// Backend
#include "../QDebugStream.h"// Backend
#include "../../model/Model/ModelCellsTriangulated/Model3D.h"
#include "../../tools/parameters/CSParameterContext.h"// Backend
#include "../../tools/parameters/CSParameterContextTemporary.h"// Backend




complexCells::complexCells(QWidget *parent)
    : QWidget(parent),
      mpModel(NULL),
      mIsReset(true),
      mStateParametersChanged(false),
      mThreadPaused(false),
      mpParameterModel(NULL),
      mpParameterContext(NULL)
{

  // Setting up the Qt Designer code
  setupUi(this);

  setObjectName("tabComplexCells");

  mName = Model3D::xmlType;

  // Parameters:
  // add a Context for this kind of model
  // at the begining we have no model, therefore we create a temporary tree
  // which has to be deleted after the context has been replaced by the
  // model's parameter tree:
  CSParameterContextTemporary * dummyContext = new CSParameterContextTemporary( mName );
  mpParameterContext = dynamic_cast<CSParameterContext *>( dummyContext );
  Model3D::DefaultParameters( mpParameterContext );

  mpParameterModel = new QCSParameterModel( mpParameterContext );
  parameterTreeView->setModel( mpParameterModel );
  parameterTreeView->setItemDelegate( new QCSParameterDelegate() );

  for ( int i =0; i < 3; ++i )
    parameterTreeView->resizeColumnToContents(i);


  // Parameterization

  // Simulation control
  connect( buttonStartSimulation, SIGNAL(clicked(bool)),   this, SLOT(startSimulationClicked()) );
  connect( buttonResetModel, SIGNAL(clicked(bool)),   this, SLOT(resetSimulation()) );
  connect( buttonAbortSimulation, SIGNAL(clicked(bool)), this, SLOT(buttonAbortSimulationButtonClicked()) );

  // Reset the ParameterTreeView with default values
  connect( pushButtonResetParametersToDefaults, SIGNAL( clicked() ), this, SLOT( ResetParametersToDefaults() ) );

  connect( buttonNewDisplay, SIGNAL(clicked(bool)), this, SLOT(createNewDisplay()) );

  GUIStatePristine();

}
complexCells::~complexCells()
{
  delete mpModel;
}


/*
 * Simulation control
 */

/*
 * region Init / Reset model
 */
void complexCells::initFromModel( CSModel *newModel )
{

  if ( !newModel ) return;

  if ( mpModel ){
    // deregister model from core
    core->models.erase( core->models.find(mpModel->name) );
    mpParameterContext = NULL;
    delete mpModel;
  }

  mpModel = dynamic_cast<Model3D *>( newModel );

  if (!mpModel){
    std::cerr << "Incompatible Model passed to tab complexCells\n";
    return;
  }

  // register model in core
  mName = newModel->name;

  core->models[mName] = mpModel;



  ConnectModelToGUI();


  // initializing the GUI's parameterTreeView with the model's parameter context.
  //
  // store the old parameterContext's pointer
  // (for later 'deep' deletion, i.e. incl. the data pointers)
  CSParameterContextTemporary * dummyContext = static_cast<CSParameterContextTemporary *>(mpParameterContext);

  QCSParameterModel *dummyModel = mpParameterModel;



  mpParameterContext = mpModel->GetParameters();
  mpParameterModel   = new QCSParameterModel(mpParameterContext);
  parameterTreeView->setModel( mpParameterModel );

  connect( mpParameterModel, SIGNAL( dataChanged(const QModelIndex &, const QModelIndex &) ), this, SLOT( ParametersChanged() ) );

  if ( dummyContext )
    delete dummyContext;

  if ( dummyModel )
    delete dummyModel;

  mIsReset = true;




  /*
   * Enable buttons
   */
  GUIStateInitialized();



}

void complexCells::resetSimulation()
{

  if ( ! mpModel ){

    mpModel = new Model3D();

    if ( mvDisplays.size() > 0 )
        {
            // add model to all existing displays:
            std::vector<QCSGLDisplay *>::iterator dispIt;
            for ( dispIt = mvDisplays.begin(); dispIt != mvDisplays.end(); ++dispIt )
                {
                    (*dispIt)->setModel(mpModel);
                    (*dispIt)->update();
                }
        }
    else
        createNewDisplay();

    ConnectModelToGUI();

    // ToDo:  prevent overwriting of existing models of the same name:
    core->models[mName] = mpModel;

    this->PushParametersToModel();

    mpModel->SetupSimulation();

  }else{
    if ( mpModel->SimulationThread()->isRunning() ){
      emit abortThread();
      mpModel->SimulationThread()->wait();
    }
  }

  mIsReset = true;
  mpModel->Reset();

  progressBar->setValue(0);



  /*
   * Enable buttons
   */
    GUIStateInitialized();


}



/*
 * region Start simulation (Button)
 */

void complexCells::startSimulationClicked()
{

  // No function without model init
  if (!mpModel) return;

  // simulation started => pause/resume when clicked
  if ( mpModel->SimulationThread()->isRunning() )
    pauseResumeSimulation();
  else{ // start the simulation:{

	  /*
	   * Prepare simulation control / observation in model
	   */

    mpModel->time_simulation = doubleSpinBoxSimulateUntil->value();
    //mpModel->enableObservation = checkBoxEnableObservation->isChecked();
    //mpModel->observeEveryDays  = doubleSpinBoxObserveEvery->value();


    /*
     * Prepare progress bar
     */

    progressBar->setMaximum(1000);
    progressBar->setValue(0);


    /*
     * Prepare and start thread that runs simulation
     */

    // initialize model parameters from GUI
    // PushAllParametersToBiolink();

    mpModel->SetupSimulation();

    mpModel->Run();

    // Enable/disable buttons
    GUIStateRunning();


    mIsReset = false;
  }
}



/*
 *  region Pause / Resume simulation (Button)
 */
void complexCells::pauseResumeSimulation()
{

  mThreadPaused = !mThreadPaused;
  emit pauseResumeThread( mThreadPaused );

  if ( mThreadPaused )
    GUIStatePaused();
  else
    GUIStateRunning();

}



/*
 * Abort simulation (Button)
 */
void complexCells::buttonAbortSimulationButtonClicked()
{

  // No function without model init
  if (!mpModel) return;

  emit abortThread();

  GUIStatePristine();

  // to avoid race conditions with abortThread() we should not use:
  // mpModel->enableSimulation = false;
  mThreadPaused = false;

}


/*
 * Method called after thread ended
 */
void complexCells::doAfterThreadIsFinished()
{

  //Set progress bar to 100%
  progressBar->setValue(1000);

  // Enable/disable buttons

  GUIStateFinished();
  mIsReset = true;

}



/*
 * region GUI Update (during thread update)
 */

void complexCells::threadUpdate()
{

  // No update without model init
  if (!mpModel) return;

  //Update simulation progress spinners in monolayer tab
  spinBoxSimulationProgressN->setValue(mpModel->cells.size());
  doubleSpinBoxSimulationProgressTime->setValue(mpModel->time);


  //Update progress bar
  progressBar->setValue((int)(1000.0 * mpModel->GetSimulationProgress()));

}








/*
 * Parameters
 */

/*
 * Set / Load
 */

void complexCells::ResetParametersToDefaults()
{

  Model3D::DefaultParameters( mpParameterContext );

  parameterTreeView->update();

  ParametersChanged();

}


void complexCells::ParametersChanged()
{

  mStateParametersChanged = true;

}


// only has to be done once after the Model is created
// but the parameters are not yet synchronised to the right pointers!
void complexCells::PushParametersToModel()
{

  // initializing the model's parameters with the possibly changed values in the GUI:
  mpModel->InitParameters( mpParameterContext );

  // destroying the temporary parameter tree
  CSParameterContextTemporary * dummyContext = static_cast<CSParameterContextTemporary *>(mpParameterContext);

  if ( dummyContext )
    delete dummyContext;

  if ( mpParameterModel )
    delete mpParameterModel;

  // updating the pointer to the parameter tree
  mpParameterContext = mpModel->GetParameters();

  // updating the GUI's data structures
  mpParameterModel = new QCSParameterModel( mpParameterContext );

   parameterTreeView->setModel( mpParameterModel );

  connect( mpParameterModel, SIGNAL( dataChanged(const QModelIndex &, const QModelIndex &) ),this, SLOT( ParametersChanged() ) );

}

void complexCells::createNewDisplay()
{

  if( mpModel != NULL){

    QCSGLDisplay * newDisp = new QCSGLDisplay();

    newDisp->setModel(mpModel);
    newDisp->show();

    connect( mpModel->SimulationThread() , SIGNAL(updateStep(int)), newDisp, SLOT(update()) );

    mvDisplays.push_back( newDisp );

  }

}

void complexCells::ConnectModelToGUI()
{

  if ( !mpModel )
    return;

  QCSSimulationThread * simulationThread = mpModel->SimulationThread();

  connect( simulationThread, SIGNAL(updateStep(int)), this, SLOT(threadUpdate()) );
  connect( this, SIGNAL(pauseResumeThread(bool)), simulationThread, SLOT(pause(bool)) );
  connect( this, SIGNAL(abortThread()), simulationThread, SLOT(reset()) );
  connect( simulationThread, SIGNAL(finishedNormally()), this, SLOT(doAfterThreadIsFinished()) ); // To enable bottons when thread is stopped

  QCSGLDisplay * display = dynamic_cast<QCSGLDisplay *>( parent()->findChild<QCSGLDisplay *>("Main Display") );
  if ( display ){
    display->setModel( dynamic_cast<CSModel *>( mpModel ) );
    connect( simulationThread, SIGNAL(updateStep(int)), display, SLOT(update()) );
    display->update();
  }

}

void complexCells::GUIStatePristine()
{
  buttonResetModel->setEnabled(true);

  buttonStartSimulation->setText("Start simulation");
  buttonStartSimulation->setEnabled(false);

  buttonAbortSimulation->setEnabled(false);


  parameterTreeView->setEnabled( true );

}
void complexCells::GUIStateInitialized()
{
  GUIStatePristine();

  buttonStartSimulation->setEnabled(true);

}
void complexCells::GUIStateRunning()
{
  buttonResetModel->setEnabled(false);

  buttonStartSimulation->setText("Pause simulation");
  buttonStartSimulation->setEnabled(true);

  buttonAbortSimulation->setEnabled(true);

  parameterTreeView->setEnabled(false);

}
void complexCells::GUIStatePaused()
{

  // Though we should only get into this state from GUIStateRunning()...
  GUIStateRunning();

  buttonStartSimulation->setText("Resume simulation");

  parameterTreeView->setEnabled(true);

}
void complexCells::GUIStateFinished()
{
  GUIStateInitialized();
}

