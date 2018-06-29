//
//  ODEsTab.cpp
//  CSGUI
//
//  Created by celliere on 12/12/12.
//
//




#pragma region Includes

#include "ODEsTab.h"
#include "../../model/Model/ModelCellsSpherical/ModelCellsSpherical.h"
#include "../../model/Model/ModelCellsSphericalODE/ModelCellsSphericalODE.h"
// Backend
#include "../../tools/parameters/CSParameterContext.h"
#include "../../tools/parameters/CSParameterContextTemporary.h"
#include "../QCSParameterModel.h"
#include "../QCSParameterDelegate.h"
#include "../QCSGLDisplay.h"
#include "../QCSSimulationThread.h"
#include "../QDebugStream.h"
#include <QFileDialog>
#include <QTableWidgetItem>

// for some reason these three headers had to be included one by one, although
// they should have been pulled in by sbmlsolver/odeSolver.h in ODEsTab.h
#include <sbmlsolver/modelSimplify.h>
#include <sbmlsolver/odeConstruct.h>
#include <sbmlsolver/processAST.h>

#include <iostream>

#pragma endregion

#pragma region Init

int MyPredicate(const ASTNode_t *node); //prototype

ODEsTab::ODEsTab(QWidget *parent)
: QWidget(parent),
mpModel(NULL),
sbmlWithLink(NULL),
mIsReset(true),
mStateParametersChanged(false),
mThreadPaused(false)
{
    // Setting up the Qt Designer code
    setupUi(this);
    setObjectName("tabODEs");
    
    mName = ModelCellsSpherical::xmlType;
    
    int i = 2;
    // choose a name that is not yet registered
    //  will be useful when having multiple tabs of the same kind.
    if ( core->models.find( mName ) != core->models.end() )
        while ( true )
        {
            std::stringstream a;
            a << mName << " " << i;
            if ( core->models.find( a.str() ) == core->models.end())
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
    CSParameterContextTemporary * dummyContext = new CSParameterContextTemporary( mName );
    mpParameterContext = dynamic_cast<CSParameterContext *>( dummyContext );
    ModelCellsSpherical::DefaultParameters( mpParameterContext );
    
    mpParameterModel = new QCSParameterModel( mpParameterContext );
    parameterTreeView->setModel( mpParameterModel );
    parameterTreeView->setItemDelegate( new QCSParameterDelegate() );
    
    for ( int i =0; i < 3; ++i )
        parameterTreeView->resizeColumnToContents(i);
    
    
#pragma region Connects
    
    // Parameterization
    // connect( comboBoxSubModel,  SIGNAL(currentIndexChanged(int)),   this, SLOT(subModelIndexChanged(int)) );
    
    //  Following is taken care of by a call to PushAllParametersToBiolink() in startSimulationClicked()
    connect( pushButtonApplyParametersToAllExistingCells, SIGNAL(clicked()), this, SLOT(ApplyParametersToAllExistingCells()) );
    connect( textSimulationName, SIGNAL(editingFinished()), this, SLOT(simulationNameChanged()) );
    
    // Simulation control
    connect( buttonStartSimulation, SIGNAL(clicked(bool)),   this, SLOT(startSimulationClicked()) );
    connect( buttonResetModel, SIGNAL(clicked(bool)),   this, SLOT(resetSimulation()) );
    connect( buttonAbortSimulation, SIGNAL(clicked(bool)), this, SLOT(buttonAbortSimulationButtonClicked()) );
    connect( buttonWritePovray, SIGNAL(clicked(bool)), this, SLOT(observeCellPopulationPovButtonClicked()) );
    connect( observeCellPopulationSnapshotButton, SIGNAL(clicked(bool)), this, SLOT(observeCellPopulationSnapshotButtonClicked()) );
    //connect( cellsVisible, SIGNAL(clicked()), this, SLOT(setCellsVisible()) );
    
    // Reset the ParameterTreeView with default values
    //connect( pushButtonResetParametersToDefaults, SIGNAL( clicked() ), this, SLOT( ResetParametersToDefaults() ) );
    

    // Visualisation
    connect( comboBoxCellStaining,  SIGNAL(currentIndexChanged(int)),   this, SLOT(updatedCellPopulationStaining(int)) );
    
    // ODEs
    connect(browse_sbml, SIGNAL(clicked()), this, SLOT(openFileDialog()));
    connect(import_sbml, SIGNAL(clicked()), this, SLOT(importSBML()));
    connect(check_useODEs, SIGNAL(stateChanged(int)), this, SLOT(enableSBMLImport(int)));
    connect(tableWidget_param, SIGNAL(cellChanged(int,int)), this, SLOT(changeParameterValue(int,int)));
    connect(tableWidget_eq, SIGNAL(cellChanged(int,int)), this, SLOT(changeInitialValue(int,int)));
        
#pragma endregion
    
    // Monolayer console link
    mpDuplexer = new QDebugStream(core->tools->output->consoleModelCellsSpherical_Simulation, monolayerConsole);
    
    GUIStatePristine();
}


ODEsTab::~ODEsTab()
{
    delete mpDuplexer;
    delete mpModel;
}

void
ODEsTab::initFromModel( CSModel *newModel )
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
ODEsTab::resetSimulation()
{
    if ( ! mpModel )
    {
        if (check_useODEs->isChecked()){
            
            if ( ! sbmlWithLink){
                monolayerConsole->append( "ERROR : you checked the box to use ODEs but have not yet imported an SBML model. The model could not be initialized/reseted.");
                return;
            }else {
                mpModel = new ModelCellsSphericalODE();
                
                odeWithLink = ODEModel_create(sbmlWithLink);
                
                double Error = lineEdit_Error->text().toDouble();
                double RError = lineEdit_RError->text().toDouble();
                int Method = comboBox_method->currentIndex();
                int useCompiled;
                if (checkBox_useCompiled->isChecked()) {
                    useCompiled = 1;
                }else{
                    useCompiled = 0;
                }
                std::cout << Error << endl;
                std::cout << RError << endl;
                static_cast<ModelCellsSphericalODE *>(mpModel)->SetCVODESettings(Error, RError, Method, useCompiled);
                static_cast<ModelCellsSphericalODE *>(mpModel)->SetSBMLModel(sbmlWithLink);
            }
        }else{
            mpModel = new ModelCellsSpherical();
        }
        
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
    mpModel->Reset(0);
    //static_cast<ModelCellsSphericalODE *>(mpModel)->Reset(1);
    
    progressBar->setValue(0);
    GUIStateInitialized();
}

void ODEsTab::startSimulationClicked()
{
    // No function without model init
    if (!mpModel) return;
    
    // simulation started => pause/resume when clicked
    if ( mpModel->SimulationThread()->isRunning()  )
        pauseResumeSimulation();
    else // start the simulation:
    {
        mpModel->simulateUntilDays = doubleSpinBoxSimulateUntil->value();
        mpModel->enableObservation = checkBoxEnableObservation->isChecked();
        mpModel->observeEveryDays  = doubleSpinBoxObserveEvery->value();

        /*if (check_useODEs->isChecked())
        {
            static_cast<ModelCellsSphericalODE *>(mpModel)->SetupSimulation();
        }else{*/
            //mpModel->SetupSimulation();
        //}
        
        progressBar->setMaximum(1000);
        progressBar->setValue(0);
        
        QString s;
        
        if (checkBoxEnableObservation->isChecked())
        {
            monolayerConsole->append("Started simulation (" + textSimulationName->displayText() + ") from " + s.number(mpModel->biolink->getTimeInDays(mpModel->time)) + " until " + s.number(doubleSpinBoxSimulateUntil->value()) + " days (observe every " + s.number(doubleSpinBoxObserveEvery->value()) + " days).");
        }
        else
        {
            monolayerConsole->append("Started simulation (" + textSimulationName->displayText() + ") from " + s.number(mpModel->biolink->getTimeInDays(mpModel->time)) + " until " + s.number(doubleSpinBoxSimulateUntil->value()) + " days without observations.");
        }
        
        if ( mIsReset )
            mpModel->UpdateParametersFromBioLink();
        
        mpModel->SetupSimulation();
        mpModel->Run();
        GUIStateRunning();
        mIsReset = false;
    }
}

void
ODEsTab::pauseResumeSimulation()
{
    mThreadPaused = !mThreadPaused;
    emit pauseResumeThread( mThreadPaused );
    
    if ( mThreadPaused )
        GUIStatePaused();
    else
        GUIStateRunning();
}

void
ODEsTab::buttonAbortSimulationButtonClicked()
{
    // No function without model init
    if (!mpModel) return;
    
    emit abortThread();
    
    GUIStatePristine();
    
    // to avoid race conditions with abortThread() we should not use:
    // mpModel->enableSimulation = false;
    mThreadPaused = false;
}

void
ODEsTab::doAfterThreadIsFinished()
{
    progressBar->setValue(1000);
    GUIStateFinished();
    mIsReset = true;
}

void
ODEsTab::threadUpdate()
{
    // No update without model init
    if (!mpModel) return;
    
    spinBoxSimulationProgressN->setValue(mpModel->cells.size());
    doubleSpinBoxSimulationProgressTime->setValue(mpModel->biolink->getTimeInDays(mpModel->time));
    progressBar->setValue((int)(1000.0 * mpModel->GetSimulationProgress()));

    if (check_useODEs->isChecked())
    {
        static_cast<ModelCellsSphericalODE *>(mpModel)->UpdateCellsStaining(comboBoxCellStaining->currentIndex());
    }else{
        mpModel->UpdateCellsStaining(comboBoxCellStaining->currentIndex());
    }    
}

void
ODEsTab::observeCellPopulationSnapshotButtonClicked()
{
    // No function without model init
    if (!mpModel) return;
    
    mpModel->observe->ObserveCellPopulationSnapshot(mpModel->cells);
}


void
ODEsTab::observeCellPopulationPovButtonClicked()
{
    // No function without model init
    if (!mpModel) return;
    
    mpModel->observe->WritePOV(0);
}

void
ODEsTab::updatedCellPopulationStaining(int index)
{
	// No function without model init
	if (!mpModel) return;
    
    // Within a running simulation this is automatically updated
    if (mpModel->enableSimulation == false)
    {
        mpModel->UpdateCellsStaining(index);
    }
}

void ODEsTab::simulationNameChanged()
{
    // No function without model init
    if ( !mpModel ) return;
    
    mpModel->SetName(textSimulationName->displayText().toStdString());
    monolayerConsole->append("Changed simulation name to: " + textSimulationName->displayText());
}


void ODEsTab::ApplyParametersToAllExistingCells()
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

void ODEsTab::ResetParametersToDefaults()
{
    ModelCellsSpherical::DefaultParameters( mpParameterContext );
    
    parameterTreeView->update();
    
    ParametersChanged();
}


void ODEsTab::ParametersChanged()
{
    if ( mpModel && mpModel->cells.size() )
        pushButtonApplyParametersToAllExistingCells->setEnabled(true);
    
    mStateParametersChanged = true;
}

// only has to be done once after the Model is created
// but the parameters are not yet synchronised to the right pointers!
void ODEsTab::PushParametersToModel()
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

void
ODEsTab::setCellsVisible(){
    
    if( !mpModel )
        return;
    
    //mpModel->setVisible( cellsVisible->isChecked() );
    
}

void
ODEsTab::ConnectModelToGUI()
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
        monolayerConsole->append( "Main display not found\n" );
    else
    {
        display->setModel( dynamic_cast<CSModel *>( mpModel ) );
        connect( simulationThread, SIGNAL(updateStep(int)), display, SLOT(update()) );
        display->update();
    }
}


void
ODEsTab::GUIStatePristine()
{
    buttonResetModel->setEnabled(true);
    
    buttonStartSimulation->setText("Start simulation");
    buttonStartSimulation->setEnabled(false);
    
    buttonAbortSimulation->setEnabled(false);
    
    comboBoxSubModel->setEnabled( true );
    
    parameterTreeView->setEnabled( true );
    
    pushButtonApplyParametersToAllExistingCells->setEnabled(false);
}


void
ODEsTab::GUIStateInitialized()
{
    GUIStatePristine();
    
    buttonStartSimulation->setEnabled(true);
}


void
ODEsTab::GUIStateRunning()
{
    buttonResetModel->setEnabled(false);
    
    buttonStartSimulation->setText("Pause simulation");
    buttonStartSimulation->setEnabled(true);
    
    buttonAbortSimulation->setEnabled(true);
    
    comboBoxSubModel->setEnabled(false);
    
    parameterTreeView->setEnabled(false);
    pushButtonApplyParametersToAllExistingCells->setEnabled(false);
}


void
ODEsTab::GUIStatePaused()
{
    // Though we should only get into this state from GUIStateRunning()...
    GUIStateRunning();
    
    buttonStartSimulation->setText("Resume simulation");
    
    parameterTreeView->setEnabled(true);
    if ( mStateParametersChanged )
        pushButtonApplyParametersToAllExistingCells->setEnabled(true);
}


void
ODEsTab::GUIStateFinished()
{
    GUIStateInitialized();
    
    comboBoxSubModel->setEnabled( false );
    
    if ( mStateParametersChanged )
        pushButtonApplyParametersToAllExistingCells->setEnabled(true);
}



void ODEsTab::openFileDialog()
{
    QString sbmlFile = QFileDialog::getOpenFileName(this,"Open an SBML file", "~", "XML files (*.xml)");
    editSBMLPath->setText(sbmlFile);
}

/* Use the SOSlib and the libSBML libraries to import the SBML file et create a Model. The main Model's features (ODE equation, species initial amounts, 
 global parameters values) are then displayed in the GUI. If the user make changes in these equations, they will NOT be applied to the Model.*/

void ODEsTab::importSBML()
{
    
    groupBox_12->setEnabled(true);
    groupBox_13->setEnabled(true);
    tableWidget_param->setEnabled(true);
    
    tableWidget_param->clearContents();
    tableWidget_eq->clearContents();
    
    QString sbmlFile = editSBMLPath->text(); // takes the text in the QLineEdit
    // convert QString to char*, needed as input of readSBML :
    const char *ch;
    ch=sbmlFile.toStdString().c_str();
    
    SBMLReader_t *sr = SBMLReader_create();
    SBMLDocument_t *d = SBMLReader_readSBML(sr, ch);
    SBMLReader_free(sr);
    sbml = SBMLDocument_getModel(d); // creates Model (class of libSBML)
    sbmlWithLink = SBMLDocument_getModel(d);
    
    odesbml = ODEModel_create(sbml); // creates ODEModel (class of SOSlib where all the ODEs are written in RateRules)

    /*for (int i=0; i<ODEModel_getNumValues(odesbml); i++)
    {std::cout << odesbml->names[i] << "  " << odesbml->values[i] << endl;}
    std::cout << odesbml->neq << "  " << odesbml->nass << "  " << odesbml->ninitAss << "  " << odesbml->nconst << endl;*/
    
    
    //****** Display of the equation/initial value Table :
    
    tableWidget_eq->setColumnCount(2);
    QStringList cheaders;
    cheaders.append("Reactions");
    cheaders.append("Initial conc."); //amount are transformed into concentrations except if compartment size is 0 or species hasOnlySubstanceUnits
    tableWidget_eq->setHorizontalHeaderLabels (cheaders);
    tableWidget_eq->setRowCount(sbml->getNumSpecies());
    
    
    // display the dc/dt= and the equations using the odeModel
    for ( int i=0; i<odesbml->neq; i++ )
    {
        // Replace the ODEs in odemodel, written as for ex. R0-R1 (only name of reaction appears) by explicit equation. Global param as displayed with their name, local parameters are replaced by their value (to avoid duplicated local parameters).
        ASTNode_t *ode = copyAST(odesbml->ode[i]);
        for ( int j=odesbml->nass-1; j>=0; j-- )
            AST_replaceNameByFormula(ode,odesbml->names[odesbml->neq + j], odesbml->assignment[j]);
        
        // display equation rhs
        char* formula = SBML_formulaToString(ode); // ASTNode to char*
        QString *qformula = new QString(formula);  // char* to QString
        QTableWidgetItem *newItem = new QTableWidgetItem(*qformula);
        newItem->setFlags(newItem->flags() ^ Qt::ItemIsEditable);
        tableWidget_eq->blockSignals(true);
        tableWidget_eq->setItem(i, 0, newItem);
        tableWidget_eq->blockSignals(false);
        
        // display dc/dt=
        variableIndex_t * vi = ODEModel_getOdeVariableIndex(odesbml, i);
        const char * protName = ODEModel_getVariableName(odesbml, vi);
        VariableIndex_free(vi);
        QString *qProtName = new QString("d%1/dt=");
        QString qName = *qProtName;
        QString qName2 = qName.arg(protName); // replace %1 by protein name
        QTableWidgetItem *newHeader = new QTableWidgetItem(qName2);
        tableWidget_eq->setVerticalHeaderItem (i , newHeader); // put dc/dt= as header
        
    }
    
    for (unsigned int i=0; i<(sbml->getNumSpecies()-odesbml->neq); i++ )
    {
        
        // display equation rhs
        QTableWidgetItem *newItem = new QTableWidgetItem("0");
        newItem->setFlags(newItem->flags() ^ Qt::ItemIsEditable);
        tableWidget_eq->blockSignals(true);
        tableWidget_eq->setItem(i+odesbml->neq, 0, newItem);
        tableWidget_eq->blockSignals(false);
        
        // display dc/dt=
        const std::string id = odesbml->names[odesbml->neq+odesbml->nass+sbml->getNumCompartments()+i];
        QString qId = QString::fromStdString(id);
        QString *qProtName = new QString("d%1/dt=");
        QString qName = *qProtName;
        QString qName2 = qName.arg(qId); // replace %1 by protein name
        QTableWidgetItem *newHeader = new QTableWidgetItem(qName2);
        tableWidget_eq->setVerticalHeaderItem (i+odesbml->neq, newHeader); // put dc/dt= as header
        
    }
    
    // display the initial amounts using the odesbml model. InitialAmount is converted in concentration by SOSlib
    
    for (int i=0; i<odesbml->neq; i++)
    {
        QString sIniVal;
        const double val = odesbml->values[i];
        sIniVal = sIniVal.setNum(val);
        QTableWidgetItem *newIniValue = new QTableWidgetItem(sIniVal);
        tableWidget_eq->blockSignals(true);
        tableWidget_eq->setItem(i,1,newIniValue);
        tableWidget_eq->blockSignals(false);
    }
    
    for (unsigned int i=0; i<(sbml->getNumSpecies()-odesbml->neq); i++)
    {
        // get initial concentrations or amounts for each species 
        QString sIniVal;
        const double val = odesbml->values[odesbml->neq+odesbml->nass+sbml->getNumCompartments()+i];
        sIniVal = sIniVal.setNum(val);
        QTableWidgetItem *newIniValue = new QTableWidgetItem(sIniVal);
        tableWidget_eq->blockSignals(true);
        tableWidget_eq->setItem(i+odesbml->neq,1,newIniValue);
        tableWidget_eq->blockSignals(false);
        
        
        

    }
    

    tableWidget_eq->resizeColumnsToContents(); // Adjust the column width.
    tableWidget_eq->setColumnWidth( 0, 1050 ); //fix the width (can only be done after resizeColumnToContents)
    tableWidget_eq->setColumnWidth( 1, 100 );
    
    // ******** display the parameter table
    tableWidget_param->setRowCount(2);
    QStringList rheaders;
    rheaders.append("Global Parameter");
    rheaders.append("Value");
    tableWidget_param->setVerticalHeaderLabels (rheaders);
    tableWidget_param->horizontalHeader()->hide();
    tableWidget_param->setColumnCount(sbml->getNumParameters()+sbml->getNumCompartments());

    // for each global parameter
    for (unsigned int i=0; i<sbml->getNumParameters(); i++)
    {
        //Parameter *p = sbml->getParameter(i);
        //const double val = p->getValue();
        const double val = odesbml->values[ODEModel_getNumValues(odesbml)-sbml->getNumParameters()+i];
        QString sVal;
        sVal= sVal.setNum(val); //transform double into QString needed for QTableWidgetItem
        //const std::string id = p->getId();
        const std::string id = odesbml->names[ODEModel_getNumValues(odesbml)-sbml->getNumParameters()+i];
        QString qId = QString::fromStdString(id); // from sting to QString
        QTableWidgetItem *newValue = new QTableWidgetItem(sVal);
        QTableWidgetItem *newID = new QTableWidgetItem(qId);
        newID->setFlags(newID->flags() ^ Qt::ItemIsEditable);
        tableWidget_param->blockSignals(true);
        tableWidget_param->setItem(0, i, newID);
        tableWidget_param->setItem(1, i, newValue);
        tableWidget_param->blockSignals(false);
    }
    
    for (unsigned int i=0; i< sbml->getNumCompartments(); i++)
    {
        const double val = odesbml->values[odesbml->neq+odesbml->nass+i];
        QString sVal;
        sVal= sVal.setNum(val); //transform double into QString needed for QTableWidgetItem
        const std::string id = odesbml->names[odesbml->neq+odesbml->nass+i];
        QString qId = QString::fromStdString(id); // from sting to QString
        QTableWidgetItem *newValue = new QTableWidgetItem(sVal);
        QTableWidgetItem *newID = new QTableWidgetItem(qId);
        newID->setFlags(newID->flags() ^ Qt::ItemIsEditable);
        tableWidget_param->blockSignals(true);
        int nPara = sbml->getNumParameters();
        tableWidget_param->setItem(0, i+nPara, newID);
        tableWidget_param->setItem(1, i+nPara, newValue);
        tableWidget_param->blockSignals(false);
    }
    tableWidget_param->resizeColumnsToContents();
    
    
}

void ODEsTab::enableSBMLImport(int state)
{
    if (state == 0)
    {
    groupBox_11->setEnabled(false);
    groupBox_12->setEnabled(false);
    groupBox_13->setEnabled(false);
    ResetParametersToDefaults();
    parameterTreeView->setModel( new QCSParameterModel( mpParameterContext ) );
        
    }else if (state == 2)
    {
        groupBox_11->setEnabled(true);
        ResetParametersToDefaults();
        parameterTreeView->setModel( new QCSParameterModel( mpParameterContext ) );
        
        QTableWidgetItem *item = tableWidget_eq->item(0,0);
        if (item)
        {
            groupBox_12->setEnabled(true);
            groupBox_13->setEnabled(true);
        }
    }
    
}

void ODEsTab::changeParameterValue(int row ,int column)
{
    if (row == 1) //names cannot be changed, only values
    {
        double val = tableWidget_param->item(row,column)->text().toDouble();
        if ((unsigned int)column < sbml->getNumParameters()) // parameter is changed
        {
            odesbml->values[ODEModel_getNumValues(odesbml)-sbml->getNumParameters()+column] = val;
            const char *id = odesbml->names[ODEModel_getNumValues(odesbml)-sbml->getNumParameters()+column];
            const char *rid = "";
            Model_setValue(sbmlWithLink, id, rid, val);
        }else // compartment size is changed
        {
            odesbml->values[odesbml->neq+odesbml->nass+column] = val;
            const char *id = odesbml->names[odesbml->neq+odesbml->nass+column];
            const char *rid = "";
            Model_setValue(sbmlWithLink, id, rid, val);
        }
    }
}






void ODEsTab::changeInitialValue(int row, int column)
{
    if (column == 1) //only initial value can be changed, not teh equations
    {
        double val = tableWidget_eq->item(row,column)->text().toDouble();
        char *id;
        Species_t *s;
        if (row < odesbml->neq) // initial value of unconstant species is changed
        {
            odesbml->values[row] = val;
            id = odesbml->names[row];
        }else // initial value of constant species is changed
        {
            odesbml->values[odesbml->neq+odesbml->nass+sbml->getNumCompartments()+row] = val;
            id = odesbml->names[odesbml->neq+odesbml->nass+sbml->getNumCompartments()+row];
        }
        if ( (s = Model_getSpeciesById(sbmlWithLink, id)) != NULL )
        {
            if ( Species_isSetInitialAmount(s) ) // val is a concentration, but we need the amount (because it is what Model_setValue does) so multiply by the size of the compartment
            {
                const std::string comp = s->getCompartment();
                val = val * sbmlWithLink->getCompartment(comp)->getSize();
            }
        }
            
        const char *rid = "";
        Model_setValue(sbmlWithLink, id, rid, val);
    }
    
}

int MyPredicate(const ASTNode_t *node)
{
    return 1;
}



/* Checks if the nodes of the ASTNode (i.e. tree of nodes) that are names (and not values) really correspond to the id of a species, parameter, compartment, reaction of functiondefinition. The output string contains the names that do not correspond to an id.
 */
std::string ODEsTab::checkNamesInASTNode(ASTNode_t *node)
{
    int (*pointeurSurFonction)(const ASTNode_t *);  /*declare pointer of a fonction*/
    pointeurSurFonction = MyPredicate;              // initialize the pointer to the function MyPredicate
    List *list = node->getListOfNodes(MyPredicate);  // list of all nodes of the node because predicate return 1(true) for all nodes
    std::string listOfWrongNames = "";
    
    for (unsigned int i = 0; i<list->getSize(); i++) // going through all nodes
    {
        ASTNode_t *nodeInList = (ASTNode_t*)list->get(i);
        if (nodeInList->isName()) // if node is a name (and not a value), check if id correspond
        {
            if (sbmlWithLink->getSpecies(nodeInList->getName())==NULL && sbmlWithLink->getParameter(nodeInList->getName())==NULL && sbmlWithLink->getCompartment(nodeInList->getName())==NULL && sbmlWithLink->getReaction(nodeInList->getName())==NULL && sbmlWithLink->getFunctionDefinition(nodeInList->getName())==NULL) //name correspond to no id
            {
                listOfWrongNames += nodeInList->getName();
                listOfWrongNames += ", ";
            }
        }
    }
    return listOfWrongNames;

}
    
/* Takes an ASTNode, and if this ASTNode is null (usually because there was an error during parsing) return the string given as input
 */
std::string ODEsTab::checkIfASTNodeIsNULL(ASTNode_t *node, std::string paramName)
{
    if (node == NULL)
    {
        return paramName;
    }else
    {
        return "";
    }
}

