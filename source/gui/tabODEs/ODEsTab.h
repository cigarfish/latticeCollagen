//
//  ODEsTab.h
//  CSGUI
//
//  Created by celliere on 12/12/12.
//
//

#ifndef __CSGUI__ODEsTab__
#define __CSGUI__ODEsTab__

#include <iostream>

#include "ui_ODEsTab.h"
#include "sbmlsolver/odeSolver.h"

#include <string>

class ModelCellsSpherical;
class CSParameterContext;
class QCSParameterModel;
class QCSSimulationThread;
class QDebugStream;
class CSModel;

class ODEsTab : public QWidget, private Ui::ODEsTab
{
    Q_OBJECT
    
public:
    ODEsTab(QWidget *parent=0);
    ~ODEsTab();
    
    void initFromModel( CSModel * );
    
    private slots:
	// Parameterization
    void ResetParametersToDefaults();
    void PushParametersToModel();
	void ApplyParametersToAllExistingCells();
    void ParametersChanged();
    
	void simulationNameChanged();
    
    // Simulation
    void startSimulationClicked();
    void pauseResumeSimulation();
    void resetSimulation();
    void threadUpdate();
	void buttonAbortSimulationButtonClicked();
	void doAfterThreadIsFinished();
	void observeCellPopulationSnapshotButtonClicked();
	void observeCellPopulationPovButtonClicked();
    
    // Visualisation
    void updatedCellPopulationStaining(int index);
    void setCellsVisible();
    
    // ODEs
    void openFileDialog();
    void importSBML();
    void enableSBMLImport(int state);
    void changeParameterValue(int row ,int column);
    //void changeInitialValueType(int index);
    void changeInitialValue(int row, int column);
    //void unenableSBMLImport();
    std::string checkNamesInASTNode(ASTNode_t *node);
    std::string checkIfASTNodeIsNULL(ASTNode_t *node, std::string paramName);
    
private:
    void ConnectModelToGUI();
    
    void GUIStatePristine();
    void GUIStateInitialized();
    void GUIStateRunning();
    void GUIStatePaused();
    void GUIStateFinished();
    
signals:
    void abortThread();
    void pauseResumeThread(bool);
    
private:
    std::string mName;
    ModelCellsSpherical *mpModel;
    
    bool mIsReset;
    bool mStateParametersChanged;
    
    CSParameterContext *mpParameterContext;
    QCSParameterModel * mpParameterModel;
    
    QDebugStream * mpDuplexer;
    bool mThreadPaused;
    
    Model_t *sbml;
    Model_t *sbmlWithLink;
    odeModel_t *odesbml;
    odeModel_t *odeWithLink;
};

#endif /* defined(__CSGUI__ODEsTab__) */
