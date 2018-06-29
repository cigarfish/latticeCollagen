
#pragma region Includes

#include "monolayer.h"
#include "../../model/Model/ModelCellsSpherical/ModelCellsSpherical.h"
// Backend
#include "../../tools/parameters/CSParameterContext.h"
#include "../../tools/parameters/CSParameterContextTemporary.h"
#include "../QCSParameterModel.h"
#include "../QCSParameterDelegate.h"
#include "../QCSGLDisplay.h"
#include "../QCSSimulationThread.h"
#include "../QDebugStream.h"

#include "../../tools/discvor/discvor.h"

#include <iostream>

#pragma endregion

#pragma region Init

monolayer::monolayer(QWidget *parent)
    : QWidget(parent),
      mpModel(NULL),
      mIsReset(true),
      mStateParametersChanged(false),
      mThreadPaused(false)
{
    // Setting up the Qt Designer code
    setupUi(this);

    setObjectName("tabMonolayer");

    mName = ModelCellsSpherical::xmlType;

    int i = 2;
    // choose a name that is not yet registered
    //  will be useful when having multiple tabs of the same kind.
    if ( core->models.find( mName ) != core->models.end() )
        while ( true )
        {
            std::stringstream a;
            a << mName << " " << i;
            if ( core->models.find( a.str() ) == core->models.end() )
            {
                mName = a.str();
                break;
            }

            ++i;
        }

    // set the Name field:
    textSimulationName->setText( QString(mName.c_str()) );

    // Parameters:
    // add a Context for this kind of model
    // at the begining we have no model, therefore we create a temporary tree
    // which has to be deleted after the context has been replaced by the
    // model's parameter tree:
    CSParameterContextTemporary * dummyContext =
        new CSParameterContextTemporary( mName );
    mpParameterContext = dynamic_cast<CSParameterContext *>( dummyContext );
    ModelCellsSpherical::DefaultParameters( mpParameterContext );

    mpParameterModel = new QCSParameterModel( mpParameterContext );
    parameterTreeView->setModel( mpParameterModel );
    parameterTreeView->setItemDelegate( new QCSParameterDelegate() );

    for ( int i =0; i < 3; ++i )
        parameterTreeView->resizeColumnToContents(i);

    #pragma region Connects

    // Parameterization

    //  Following is taken care of by a call to PushAllParametersToBiolink() in startSimulationClicked()
    connect( pushButtonApplyParametersToAllExistingCells, SIGNAL(clicked()), this, SLOT(ApplyParametersToAllExistingCells()) );
    connect( textSimulationName, SIGNAL(editingFinished()), this, SLOT(simulationNameChanged()) );
    // connect( saveMXFButton, SIGNAL(clicked(bool)), this, SLOT(saveAsMXF()) );

    // Simulation control
    connect( buttonStartSimulation, SIGNAL(clicked(bool)),   this, SLOT(startSimulationClicked()) );
    connect( buttonResetModel, SIGNAL(clicked(bool)),   this, SLOT(resetSimulation()) );
    connect( buttonAbortSimulation, SIGNAL(clicked(bool)), this, SLOT(buttonAbortSimulationButtonClicked()) );
    connect( buttonWritePovray, SIGNAL(clicked(bool)), this, SLOT(observeCellPopulationPovButtonClicked()) );
    connect( observeCellPopulationSnapshotButton, SIGNAL(clicked(bool)), this, SLOT(observeCellPopulationSnapshotButtonClicked()) );
    //    connect( cellsVisible, SIGNAL(clicked()), this, SLOT(setCellsVisible()) );

    connect( comboBoxScreenBackgroundColor, SIGNAL(currentIndexChanged(int)), this, SLOT( setBackgroundColor(int) ) );

    // Reset the ParameterTreeView with default values
    connect( pushButtonResetParametersToDefaults, SIGNAL( clicked() ), this, SLOT( ResetParametersToDefaults() ) );

    // Visualisation
    connect( comboBoxCellStaining,  SIGNAL(currentIndexChanged(int)),   this, SLOT(updatedCellPopulationStaining(int)) );

    // Voxelization
    connect( pushButtonVoxelize, SIGNAL(clicked()), this, SLOT(voxelize()) );
    // Quantification
    connect( quantificationRunTestAnalysisButton,  SIGNAL(clicked()),   this, SLOT(runTestAnalysis()) );

    #pragma endregion


    // Monolayer console links
    mpDuplexer = new QDebugStream(core->tools->output->consoleModelCellsSpherical_Simulation, monolayerConsole_Simulation);
    mpDuplexer = new QDebugStream(core->tools->output->consoleModelCellsSpherical_Quantification, monolayerConsole_Quantification);

    GUIStatePristine();
}

#pragma endregion

#pragma region Destroy

monolayer::~monolayer()
{
    delete mpDuplexer;
    delete mpModel;
}

#pragma endregion

#pragma region Simulation control

#pragma region Init / Reset model

void
monolayer::initFromModel( CSModel *newModel )
{
    if ( !newModel )
        return;

    if ( mpModel )
    {
        // deregister model from core
        core->models.erase( core->models.find(mpModel->name) );
        mpParameterContext = NULL;
        delete mpModel;
    }

    mpModel = dynamic_cast<ModelCellsSpherical *>( newModel );

    if (!mpModel)
    {
        std::cerr << "Incompatible Model passed to tab monolayer\n";
        return;
    }

    // register model in core
    mName = newModel->name;

    // find a name that is not taken
    int i = 1;
    if ( core->models.find( mName ) != core->models.end() )
        while ( true )
        {
            std::stringstream a;
            a << mName << " " << i;
            if ( core->models.find( a.str() ) == core->models.end() )
            {
                mName = a.str();
                mpModel->SetName(mName);
                break;
            }

            ++i;
        }

    core->models[mName] = mpModel;
    textSimulationName->setText( QString(mName.c_str()) );


    ConnectModelToGUI();


    // initializing the GUI's parameterTreeView with the model's parameter context.
    //
    // store the old parameterContext's pointer
    // (for later 'deep' deletion, i.e. incl. the data pointers)
    CSParameterContextTemporary * dummyContext =
        static_cast<CSParameterContextTemporary *>(mpParameterContext);

    QCSParameterModel *dummyModel = mpParameterModel;

    mpParameterContext = mpModel->GetParameters();
    mpParameterModel   = new QCSParameterModel(mpParameterContext);
    parameterTreeView->setModel( mpParameterModel );

    connect( mpParameterModel, SIGNAL( dataChanged(const QModelIndex &, const QModelIndex &) ),
             this, SLOT( ParametersChanged() ) );

    if ( dummyContext )
        delete dummyContext;

    if ( dummyModel )
        delete dummyModel;

    // update visualization:
    mpModel->UpdateCellsStaining(comboBoxCellStaining->currentIndex());

    mIsReset = true;

    #pragma region Enable buttons

    GUIStateInitialized();

    #pragma endregion
}


void
monolayer::resetSimulation()
{
    if ( ! mpModel )
    {
        mpModel = new ModelCellsSpherical();

        ConnectModelToGUI();

        // ToDo:  prevent overwriting of existing models of the same name:
        core->models[mName] = mpModel;

        // read out the name of the model
        simulationNameChanged();

        PushParametersToModel();
    }
    else
    {
        if ( mpModel->SimulationThread()->isRunning() )
        {
            emit abortThread();
            mpModel->SimulationThread()->wait();
        }
    }

    mIsReset = true;
    mpModel->SwitchSubmodel(comboBoxSubModel->currentIndex());
    mpModel->SetScenario((ModelCellsSpherical::Scenario)comboBoxScenario->currentIndex());
    mpModel->Reset();

    mpModel->simulateUntilDays = doubleSpinBoxSimulateUntil->value();
    mpModel->enableObservation = checkBoxEnableObservation->isChecked();
    mpModel->observeEveryDays  = doubleSpinBoxObserveEvery->value();

    progressBar->setValue(0);

    #pragma region Enable buttons

    GUIStateInitialized();

    #pragma endregion
}

#pragma endregion

#pragma region Start simulation (Button)

void monolayer::startSimulationClicked()
{
    // No function without model init
    if (!mpModel) return;

    // simulation started => pause/resume when clicked
    if ( mpModel->SimulationThread()->isRunning() )
        pauseResumeSimulation();
    else // start the simulation:
    {
        #pragma region Prepare simulation control / observation in model

        mpModel->simulateUntilDays = doubleSpinBoxSimulateUntil->value();
        mpModel->enableObservation = checkBoxEnableObservation->isChecked();
        mpModel->observeEveryDays  = doubleSpinBoxObserveEvery->value();

        #pragma endregion

        #pragma region Prepare progress bar

        progressBar->setMaximum(1000);
        progressBar->setValue(0);

        #pragma endregion

        #pragma region Write info to console in GUI

        QString s;

        if (checkBoxEnableObservation->isChecked())
        {
            monolayerConsole_Simulation->append("Started simulation (" + textSimulationName->displayText() + ") from " + s.number(mpModel->biolink->getTimeInDays(mpModel->time)) + " until " + s.number(doubleSpinBoxSimulateUntil->value()) + " days (observe every " + s.number(doubleSpinBoxObserveEvery->value()) + " days).");
        }
        else
        {
            monolayerConsole_Simulation->append("Started simulation (" + textSimulationName->displayText() + ") from " + s.number(mpModel->biolink->getTimeInDays(mpModel->time)) + " until " + s.number(doubleSpinBoxSimulateUntil->value()) + " days without observations.");
        }
        #pragma endregion

        #pragma region Prepare and start thread that runs simulation

        if ( mIsReset )
            mpModel->UpdateParametersFromBioLink();

        // initialize model parameters from GUI
        // PushAllParametersToBiolink();

        doubleSpinBoxSimulationProgressTime->setMinimum( mpModel->biolink->getTimeInDays(mpModel->time) );

        mpModel->SetupSimulation();

        mpModel->Run();

        #pragma endregion

        #pragma region Enable/disable buttons

        GUIStateRunning();

        #pragma endregion

        mIsReset = false;
    }
}

#pragma endregion

#pragma region Pause / Resume simulation (Button)

void
monolayer::pauseResumeSimulation()
{
    mThreadPaused = !mThreadPaused;
    emit pauseResumeThread( mThreadPaused );

    if ( mThreadPaused )
        GUIStatePaused();
    else
        GUIStateRunning();
}

#pragma endregion

#pragma region Abort simulation (Button)

void
monolayer::buttonAbortSimulationButtonClicked()
{
    // No function without model init
    if (!mpModel) return;

    emit abortThread();

    GUIStatePristine();

    // to avoid race conditions with abortThread() we should not use:
    // mpModel->enableSimulation = false;
    mThreadPaused = false;
}

#pragma endregion

#pragma region Method called after thread ended

void
monolayer::doAfterThreadIsFinished()
{
    #pragma region Set progress bar to 100%

    progressBar->setValue(1000);

    #pragma endregion

    #pragma region Enable/disable buttons

    GUIStateFinished();

    mIsReset = true;

    #pragma endregion
}

#pragma endregion

#pragma region GUI Update (during thread update)

void
monolayer::threadUpdate()
{
    // No update without model init
    if (!mpModel) return;

    #pragma region Update simulation progress spinners in monolayer tab

    spinBoxSimulationProgressN->setValue(mpModel->cells.size());
    doubleSpinBoxSimulationProgressTime->setValue(mpModel->biolink->getTimeInDays(mpModel->time));

    #pragma endregion

    #pragma region Update progress bar

    progressBar->setValue((int)(1000.0 * mpModel->GetSimulationProgress()));

    #pragma endregion

    #pragma region Cell staining

    mpModel->UpdateCellsStaining(comboBoxCellStaining->currentIndex());

    #pragma endregion

}

#pragma endregion

#pragma endregion

#pragma region Debug: Obtain cell population snapshot

void
monolayer::observeCellPopulationSnapshotButtonClicked()
{
    // No function without model init
    if (!mpModel) return;

    mpModel->observe->ObserveCellPopulationSnapshot(mpModel->cells);
}


void
monolayer::observeCellPopulationPovButtonClicked()
{
    // No function without model init
    if (!mpModel) return;

    mpModel->observe->WritePOV(0);
}

#pragma endregion

#pragma region Visualization

#pragma region Method called upon cell staining update
void
monolayer::updatedCellPopulationStaining(int index)
{
    // No function without model init
    if (!mpModel) return;

    // Within a running simulation this is automatically updated
    if (mpModel->enableSimulation == false)
    {
        mpModel->UpdateCellsStaining(index);
    }
}
#pragma endregion

#pragma endregion

#pragma region Parameters

#pragma region Simulation name

void monolayer::simulationNameChanged()
{
    // No function without model init
    if ( !mpModel ) return;

    mpModel->SetName(textSimulationName->displayText().toStdString());

    monolayerConsole_Simulation->append("Changed simulation name to: " + textSimulationName->displayText());
}

#pragma endregion

#pragma region Set / Load

void monolayer::ApplyParametersToAllExistingCells()
{
    // No function without model init
    if (!mpModel) return;

    // GUI -> Biolink -> Dimensionless parameters in Model (Used in simulation)
    mpModel->UpdateParametersFromBioLink();

    // Dimensionless parameters in Model-> Biolink
    mpModel->UpdateParametersForAllCells();

    pushButtonApplyParametersToAllExistingCells->setEnabled(false);

    mStateParametersChanged = false;
}


void monolayer::ResetParametersToDefaults()
{
    ModelCellsSpherical::DefaultParameters( mpParameterContext );

    parameterTreeView->update();

    ParametersChanged();
}


void monolayer::ParametersChanged()
{
    if ( mpModel && mpModel->cells.size() )
        pushButtonApplyParametersToAllExistingCells->setEnabled(true);

    mStateParametersChanged = true;
}


// only has to be done once after the Model is created
// but the parameters are not yet synchronised to the right pointers!
void monolayer::PushParametersToModel()
{
    // initializing the model's parameters with the possibly changed values in the GUI:
    mpModel->InitParameters( mpParameterContext );

    // destroying the temporary parameter tree
    CSParameterContextTemporary * dummyContext =
        static_cast<CSParameterContextTemporary *>(mpParameterContext);

    if ( dummyContext )
        delete dummyContext;

    if ( mpParameterModel )
        delete mpParameterModel;

    // updating the pointer to the parameter tree
    mpParameterContext = mpModel->GetParameters();

    // updating the GUI's data structures
    mpParameterModel = new QCSParameterModel( mpParameterContext );

    parameterTreeView->setModel( mpParameterModel );

    connect( mpParameterModel, SIGNAL( dataChanged(const QModelIndex &, const QModelIndex &) ),
             this, SLOT( ParametersChanged() ) );

    mpModel->UpdateParametersFromBioLink();

}

// void
// monolayer::setCellsVisible(){

//   if( !mpModel )
//     return;

//   mpModel->setVisible( cellsVisible->isChecked() );

// }
// void monolayer::saveAsMXF(){

//   this->mpModel->writeXML2();

// }

void
monolayer::setBackgroundColor(int colorChoice)
{
    QCSGLDisplay * display = dynamic_cast<QCSGLDisplay *>( parent()->findChild<QCSGLDisplay *>("Main Display") );
    if ( !display )
        monolayerConsole_Simulation->append( "Main display not found\n" );
    else
    {
        display->setBackgroundColor( (QCSGLDisplay::Color)colorChoice );
        display->update();
    }
}


void
monolayer::ConnectModelToGUI()
{
    if ( !mpModel )
        return;

    QCSSimulationThread * simulationThread = mpModel->SimulationThread();

    connect( simulationThread, SIGNAL(updateStep(int)), this, SLOT(threadUpdate()) );
    connect( this, SIGNAL(pauseResumeThread(bool)), simulationThread, SLOT(pause(bool)) );
    connect( this, SIGNAL(abortThread()), simulationThread, SLOT(reset()) );
    connect( simulationThread, SIGNAL(finishedNormally()), this, SLOT(doAfterThreadIsFinished()) ); // To enable bottons when thread is stopped

    QCSGLDisplay * display = dynamic_cast<QCSGLDisplay *>( parent()->findChild<QCSGLDisplay *>("Main Display") );
    if ( !display )
        monolayerConsole_Simulation->append( "Main display not found\n" );
    else
    {
        display->setModel( dynamic_cast<CSModel *>( mpModel ) );
        connect( simulationThread, SIGNAL(updateStep(int)), display, SLOT(update()) );
        display->update();
    }
}


void
monolayer::GUIStatePristine()
{
    buttonResetModel->setEnabled(true);

    buttonStartSimulation->setText("Start simulation");
    buttonStartSimulation->setEnabled(false);

    buttonAbortSimulation->setEnabled(false);

    comboBoxSubModel->setEnabled(true);
    comboBoxScenario->setEnabled(true);

    parameterTreeView->setEnabled( true );

    pushButtonApplyParametersToAllExistingCells->setEnabled(false);
}


void
monolayer::GUIStateInitialized()
{
    GUIStatePristine();

    buttonStartSimulation->setEnabled(true);
}


void
monolayer::GUIStateRunning()
{
    buttonResetModel->setEnabled(false);

    buttonStartSimulation->setText("Pause simulation");
    buttonStartSimulation->setEnabled(true);

    buttonAbortSimulation->setEnabled(true);

    comboBoxSubModel->setEnabled(false);
    comboBoxScenario->setEnabled(false);

    parameterTreeView->setEnabled(false);
    pushButtonApplyParametersToAllExistingCells->setEnabled(false);
}


void
monolayer::GUIStatePaused()
{
    // Though we should only get into this state from GUIStateRunning()...
    GUIStateRunning();

    buttonStartSimulation->setText("Resume simulation");

    parameterTreeView->setEnabled(true);
    if ( mStateParametersChanged )
        pushButtonApplyParametersToAllExistingCells->setEnabled(true);
}


void
monolayer::GUIStateFinished()
{
    GUIStateInitialized();

    comboBoxSubModel->setEnabled( false );

    if ( mStateParametersChanged )
        pushButtonApplyParametersToAllExistingCells->setEnabled(true);
}

#pragma endregion

#pragma endregion

#pragma region Quantification

#pragma region Discrete voronoi analysis

void monolayer::runTestAnalysis()
{
    #pragma region Checks

    if ( !mpModel )
    {
        monolayerConsole_Quantification->append( "Error:  No initialized model." );
        return;
    }

    if ( !mpModel->mpGraphBloodVesselNetwork )
    {
        monolayerConsole_Quantification->append("Error, no blood vessel network found.\nPlease specify the path to an MXF file under point \"Path to Blood Vessel Network\" in the \"Parameters\" tab.\n");
        return;
    }

    #pragma endregion

	// Option A: Analyze cells only
	// DiscreteVoronoi dv(&mpModel->cells, &core->tools->output->consoleModelCellsSpherical_Quantification, quantificationDVA_ProgressBar);

	// Option B: Analyze cells and blood vessel network
	DiscreteVoronoi dv(&mpModel->cells, mpModel->mpGraphBloodVesselNetwork, &core->tools->output->consoleModelCellsSpherical_Quantification, quantificationDVA_ProgressBar);

	// Run analysis (Preliminarily locking GUI)
    dv.Run(discvor_VoxelSize->value(), discvor_CellCutOffRadius->value(),discvor_EnableDebugOutput->isChecked());
}

#pragma endregion

#pragma endregion

void monolayer::voxelize()
{
    if (mpModel)
    {
        mpModel->voxelize();
        // mpModel->UpdateCellsStaining( 6 );
    }
}
