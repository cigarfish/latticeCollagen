
#ifndef MONOLAYER_H
#define MONOLAYER_H

#include "ui_monolayer.h"

#include <string>

class ModelCellsSpherical;
class CSParameterContext;
class QCSParameterModel;
class QCSSimulationThread;
class QDebugStream;
class CSModel;


class monolayer : public QWidget, private Ui::monolayer
{
    Q_OBJECT

    friend class monolayerDemo;

public:
    monolayer(QWidget *parent=0);
    ~monolayer();

    void initFromModel( CSModel * );


private slots:
	// Parameterization
    void ResetParametersToDefaults();
    void PushParametersToModel();
	void ApplyParametersToAllExistingCells();
    void ParametersChanged();
    //void saveAsMXF();

	void simulationNameChanged();

    // Simulation
    void startSimulationClicked();
    void pauseResumeSimulation();
    void resetSimulation();

	void runTestAnalysis();

    void threadUpdate();

	void buttonAbortSimulationButtonClicked();
	void doAfterThreadIsFinished();

	void observeCellPopulationSnapshotButtonClicked();
	void observeCellPopulationPovButtonClicked();

    // Visualisation
    void updatedCellPopulationStaining(int index);

    void setBackgroundColor( int index );

    void voxelize();
//    void setCellsVisible();

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
};


#endif
