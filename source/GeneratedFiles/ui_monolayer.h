/********************************************************************************
** Form generated from reading UI file 'monolayer.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MONOLAYER_H
#define UI_MONOLAYER_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QTreeView>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_monolayer
{
public:
    QGridLayout *gridLayout;
    QTabWidget *tabWidget;
    QWidget *tab;
    QGridLayout *gridLayout_3;
    QGroupBox *groupBoxParameters;
    QGridLayout *gridLayout_14;
    QTreeView *parameterTreeView;
    QPushButton *pushButtonResetParametersToDefaults;
    QPushButton *pushButtonApplyParametersToAllExistingCells;
    QGridLayout *gridLayout_5;
    QGroupBox *groupBoxSubmodel;
    QComboBox *comboBoxSubModel;
    QGroupBox *groupBoxSimulation;
    QGridLayout *gridLayout_11;
    QLabel *label_28;
    QLineEdit *textSimulationName;
    QGroupBox *groupBoxScenario;
    QComboBox *comboBoxScenario;
    QSpacerItem *horizontalSpacer_6;
    QWidget *tab_2;
    QGridLayout *gridLayout_4;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout_6;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_23;
    QSpinBox *spinBoxSimulationProgressN;
    QLabel *label_22;
    QSpacerItem *horizontalSpacer_2;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_7;
    QDoubleSpinBox *doubleSpinBoxSimulationProgressTime;
    QLabel *label_21;
    QProgressBar *progressBar;
    QGroupBox *groupBox_3;
    QGridLayout *gridLayout_7;
    QPushButton *buttonResetModel;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QDoubleSpinBox *doubleSpinBoxSimulateUntil;
    QLabel *label_2;
    QCheckBox *checkBoxEnableObservation;
    QSpacerItem *verticalSpacer;
    QHBoxLayout *horizontalLayout_2;
    QLabel *labelOutputObservablesEvery;
    QDoubleSpinBox *doubleSpinBoxObserveEvery;
    QLabel *labelOutputObservablesEveryUnit;
    QPushButton *buttonStartSimulation;
    QPushButton *buttonAbortSimulation;
    QPushButton *observeCellPopulationSnapshotButton;
    QCheckBox *checkBoxEnablePOVOutput;
    QPushButton *buttonWritePovray;
    QGroupBox *groupBox_4;
    QHBoxLayout *horizontalLayout_5;
    QTextBrowser *monolayerConsole_Simulation;
    QWidget *tab_3;
    QGridLayout *gridLayout_8;
    QGroupBox *groupBoxVisualizationScreen;
    QGridLayout *gridLayout_10;
    QComboBox *comboBoxScreenBackgroundColor;
    QLabel *label_27;
    QSpacerItem *horizontalSpacer;
    QCheckBox *checkBoxScreenShowGrid;
    QGroupBox *groupBoxPovRayOptions;
    QGridLayout *gridLayout_12;
    QComboBox *comboBox;
    QLabel *label_29;
    QSpacerItem *horizontalSpacer_4;
    QGroupBox *groupBoxVisualizationParameter;
    QTreeView *VisualizationTreeView;
    QPushButton *saveMXFButton;
    QGroupBox *groupBoxVisualizationElements;
    QGridLayout *gridLayout_9;
    QLabel *label_24;
    QSpacerItem *horizontalSpacer_3;
    QComboBox *comboBoxCellStaining;
    QPushButton *pushButtonVoxelize;
    QWidget *tab_4;
    QGroupBox *groupBox;
    QToolButton *quantificationRunTestAnalysisButton;
    QWidget *formLayoutWidget;
    QFormLayout *formLayout;
    QDoubleSpinBox *discvor_VoxelSize;
    QDoubleSpinBox *discvor_CellCutOffRadius;
    QLabel *label_3;
    QLabel *label_4;
    QCheckBox *discvor_EnableDebugOutput;
    QProgressBar *quantificationDVA_ProgressBar;
    QGroupBox *groupBox_5;
    QHBoxLayout *horizontalLayout_6;
    QTextBrowser *monolayerConsole_Quantification;

    void setupUi(QWidget *monolayer)
    {
        if (monolayer->objectName().isEmpty())
            monolayer->setObjectName(QStringLiteral("monolayer"));
        monolayer->resize(1045, 727);
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(monolayer->sizePolicy().hasHeightForWidth());
        monolayer->setSizePolicy(sizePolicy);
        gridLayout = new QGridLayout(monolayer);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        tabWidget = new QTabWidget(monolayer);
        tabWidget->setObjectName(QStringLiteral("tabWidget"));
        tabWidget->setMinimumSize(QSize(300, 0));
        tab = new QWidget();
        tab->setObjectName(QStringLiteral("tab"));
        gridLayout_3 = new QGridLayout(tab);
        gridLayout_3->setObjectName(QStringLiteral("gridLayout_3"));
        groupBoxParameters = new QGroupBox(tab);
        groupBoxParameters->setObjectName(QStringLiteral("groupBoxParameters"));
        gridLayout_14 = new QGridLayout(groupBoxParameters);
        gridLayout_14->setObjectName(QStringLiteral("gridLayout_14"));
        parameterTreeView = new QTreeView(groupBoxParameters);
        parameterTreeView->setObjectName(QStringLiteral("parameterTreeView"));

        gridLayout_14->addWidget(parameterTreeView, 0, 0, 1, 2);

        pushButtonResetParametersToDefaults = new QPushButton(groupBoxParameters);
        pushButtonResetParametersToDefaults->setObjectName(QStringLiteral("pushButtonResetParametersToDefaults"));

        gridLayout_14->addWidget(pushButtonResetParametersToDefaults, 1, 0, 1, 1);

        pushButtonApplyParametersToAllExistingCells = new QPushButton(groupBoxParameters);
        pushButtonApplyParametersToAllExistingCells->setObjectName(QStringLiteral("pushButtonApplyParametersToAllExistingCells"));
        pushButtonApplyParametersToAllExistingCells->setEnabled(false);

        gridLayout_14->addWidget(pushButtonApplyParametersToAllExistingCells, 1, 1, 1, 1);


        gridLayout_3->addWidget(groupBoxParameters, 3, 0, 1, 3);

        gridLayout_5 = new QGridLayout();
        gridLayout_5->setObjectName(QStringLiteral("gridLayout_5"));
        groupBoxSubmodel = new QGroupBox(tab);
        groupBoxSubmodel->setObjectName(QStringLiteral("groupBoxSubmodel"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(groupBoxSubmodel->sizePolicy().hasHeightForWidth());
        groupBoxSubmodel->setSizePolicy(sizePolicy1);
        groupBoxSubmodel->setMinimumSize(QSize(0, 55));
        comboBoxSubModel = new QComboBox(groupBoxSubmodel);
        comboBoxSubModel->setObjectName(QStringLiteral("comboBoxSubModel"));
        comboBoxSubModel->setEnabled(true);
        comboBoxSubModel->setGeometry(QRect(10, 20, 261, 22));

        gridLayout_5->addWidget(groupBoxSubmodel, 2, 0, 1, 1);

        groupBoxSimulation = new QGroupBox(tab);
        groupBoxSimulation->setObjectName(QStringLiteral("groupBoxSimulation"));
        QSizePolicy sizePolicy2(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(groupBoxSimulation->sizePolicy().hasHeightForWidth());
        groupBoxSimulation->setSizePolicy(sizePolicy2);
        groupBoxSimulation->setMinimumSize(QSize(280, 0));
        gridLayout_11 = new QGridLayout(groupBoxSimulation);
        gridLayout_11->setObjectName(QStringLiteral("gridLayout_11"));
        label_28 = new QLabel(groupBoxSimulation);
        label_28->setObjectName(QStringLiteral("label_28"));

        gridLayout_11->addWidget(label_28, 1, 0, 1, 1);

        textSimulationName = new QLineEdit(groupBoxSimulation);
        textSimulationName->setObjectName(QStringLiteral("textSimulationName"));
        textSimulationName->setEnabled(true);

        gridLayout_11->addWidget(textSimulationName, 1, 1, 1, 1);


        gridLayout_5->addWidget(groupBoxSimulation, 1, 0, 1, 1);

        groupBoxScenario = new QGroupBox(tab);
        groupBoxScenario->setObjectName(QStringLiteral("groupBoxScenario"));
        comboBoxScenario = new QComboBox(groupBoxScenario);
        comboBoxScenario->setObjectName(QStringLiteral("comboBoxScenario"));
        comboBoxScenario->setGeometry(QRect(10, 20, 261, 24));

        gridLayout_5->addWidget(groupBoxScenario, 2, 1, 1, 1);

        horizontalSpacer_6 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_6, 1, 1, 1, 1);


        gridLayout_3->addLayout(gridLayout_5, 1, 0, 1, 1);

        tabWidget->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QStringLiteral("tab_2"));
        gridLayout_4 = new QGridLayout(tab_2);
        gridLayout_4->setObjectName(QStringLiteral("gridLayout_4"));
        groupBox_2 = new QGroupBox(tab_2);
        groupBox_2->setObjectName(QStringLiteral("groupBox_2"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Maximum);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(groupBox_2->sizePolicy().hasHeightForWidth());
        groupBox_2->setSizePolicy(sizePolicy3);
        gridLayout_6 = new QGridLayout(groupBox_2);
        gridLayout_6->setObjectName(QStringLiteral("gridLayout_6"));
        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        label_23 = new QLabel(groupBox_2);
        label_23->setObjectName(QStringLiteral("label_23"));

        horizontalLayout_4->addWidget(label_23);

        spinBoxSimulationProgressN = new QSpinBox(groupBox_2);
        spinBoxSimulationProgressN->setObjectName(QStringLiteral("spinBoxSimulationProgressN"));
        spinBoxSimulationProgressN->setMinimumSize(QSize(70, 0));
        spinBoxSimulationProgressN->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinBoxSimulationProgressN->setKeyboardTracking(false);
        spinBoxSimulationProgressN->setMaximum(1000000);

        horizontalLayout_4->addWidget(spinBoxSimulationProgressN);

        label_22 = new QLabel(groupBox_2);
        label_22->setObjectName(QStringLiteral("label_22"));

        horizontalLayout_4->addWidget(label_22);


        gridLayout_6->addLayout(horizontalLayout_4, 0, 1, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_6->addItem(horizontalSpacer_2, 0, 3, 1, 1);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        label_7 = new QLabel(groupBox_2);
        label_7->setObjectName(QStringLiteral("label_7"));

        horizontalLayout_3->addWidget(label_7);

        doubleSpinBoxSimulationProgressTime = new QDoubleSpinBox(groupBox_2);
        doubleSpinBoxSimulationProgressTime->setObjectName(QStringLiteral("doubleSpinBoxSimulationProgressTime"));
        doubleSpinBoxSimulationProgressTime->setMinimumSize(QSize(70, 0));
        doubleSpinBoxSimulationProgressTime->setButtonSymbols(QAbstractSpinBox::NoButtons);
        doubleSpinBoxSimulationProgressTime->setKeyboardTracking(false);
        doubleSpinBoxSimulationProgressTime->setDecimals(3);
        doubleSpinBoxSimulationProgressTime->setMaximum(10000);

        horizontalLayout_3->addWidget(doubleSpinBoxSimulationProgressTime);

        label_21 = new QLabel(groupBox_2);
        label_21->setObjectName(QStringLiteral("label_21"));

        horizontalLayout_3->addWidget(label_21);


        gridLayout_6->addLayout(horizontalLayout_3, 0, 0, 1, 1);

        progressBar = new QProgressBar(groupBox_2);
        progressBar->setObjectName(QStringLiteral("progressBar"));
        QSizePolicy sizePolicy4(QSizePolicy::Expanding, QSizePolicy::MinimumExpanding);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(progressBar->sizePolicy().hasHeightForWidth());
        progressBar->setSizePolicy(sizePolicy4);
        progressBar->setValue(50);

        gridLayout_6->addWidget(progressBar, 2, 0, 1, 4);


        gridLayout_4->addWidget(groupBox_2, 1, 0, 1, 2);

        groupBox_3 = new QGroupBox(tab_2);
        groupBox_3->setObjectName(QStringLiteral("groupBox_3"));
        QSizePolicy sizePolicy5(QSizePolicy::Preferred, QSizePolicy::MinimumExpanding);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(groupBox_3->sizePolicy().hasHeightForWidth());
        groupBox_3->setSizePolicy(sizePolicy5);
        gridLayout_7 = new QGridLayout(groupBox_3);
        gridLayout_7->setObjectName(QStringLiteral("gridLayout_7"));
        buttonResetModel = new QPushButton(groupBox_3);
        buttonResetModel->setObjectName(QStringLiteral("buttonResetModel"));
        QSizePolicy sizePolicy6(QSizePolicy::Minimum, QSizePolicy::Maximum);
        sizePolicy6.setHorizontalStretch(0);
        sizePolicy6.setVerticalStretch(0);
        sizePolicy6.setHeightForWidth(buttonResetModel->sizePolicy().hasHeightForWidth());
        buttonResetModel->setSizePolicy(sizePolicy6);

        gridLayout_7->addWidget(buttonResetModel, 0, 0, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label = new QLabel(groupBox_3);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout->addWidget(label);

        doubleSpinBoxSimulateUntil = new QDoubleSpinBox(groupBox_3);
        doubleSpinBoxSimulateUntil->setObjectName(QStringLiteral("doubleSpinBoxSimulateUntil"));
        doubleSpinBoxSimulateUntil->setMaximum(10000);
        doubleSpinBoxSimulateUntil->setValue(8);

        horizontalLayout->addWidget(doubleSpinBoxSimulateUntil);

        label_2 = new QLabel(groupBox_3);
        label_2->setObjectName(QStringLiteral("label_2"));

        horizontalLayout->addWidget(label_2);


        gridLayout_7->addLayout(horizontalLayout, 1, 0, 1, 1);

        checkBoxEnableObservation = new QCheckBox(groupBox_3);
        checkBoxEnableObservation->setObjectName(QStringLiteral("checkBoxEnableObservation"));
        checkBoxEnableObservation->setChecked(true);

        gridLayout_7->addWidget(checkBoxEnableObservation, 2, 0, 1, 1);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_7->addItem(verticalSpacer, 7, 0, 1, 1);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        labelOutputObservablesEvery = new QLabel(groupBox_3);
        labelOutputObservablesEvery->setObjectName(QStringLiteral("labelOutputObservablesEvery"));

        horizontalLayout_2->addWidget(labelOutputObservablesEvery);

        doubleSpinBoxObserveEvery = new QDoubleSpinBox(groupBox_3);
        doubleSpinBoxObserveEvery->setObjectName(QStringLiteral("doubleSpinBoxObserveEvery"));
        doubleSpinBoxObserveEvery->setDecimals(2);
        doubleSpinBoxObserveEvery->setMinimum(0);
        doubleSpinBoxObserveEvery->setMaximum(100);
        doubleSpinBoxObserveEvery->setSingleStep(0.1);
        doubleSpinBoxObserveEvery->setValue(0.1);

        horizontalLayout_2->addWidget(doubleSpinBoxObserveEvery);

        labelOutputObservablesEveryUnit = new QLabel(groupBox_3);
        labelOutputObservablesEveryUnit->setObjectName(QStringLiteral("labelOutputObservablesEveryUnit"));

        horizontalLayout_2->addWidget(labelOutputObservablesEveryUnit);


        gridLayout_7->addLayout(horizontalLayout_2, 4, 0, 1, 1);

        buttonStartSimulation = new QPushButton(groupBox_3);
        buttonStartSimulation->setObjectName(QStringLiteral("buttonStartSimulation"));
        buttonStartSimulation->setEnabled(false);

        gridLayout_7->addWidget(buttonStartSimulation, 5, 0, 1, 1);

        buttonAbortSimulation = new QPushButton(groupBox_3);
        buttonAbortSimulation->setObjectName(QStringLiteral("buttonAbortSimulation"));
        buttonAbortSimulation->setEnabled(false);

        gridLayout_7->addWidget(buttonAbortSimulation, 6, 0, 1, 1);

        observeCellPopulationSnapshotButton = new QPushButton(groupBox_3);
        observeCellPopulationSnapshotButton->setObjectName(QStringLiteral("observeCellPopulationSnapshotButton"));
        QFont font;
        font.setItalic(false);
        observeCellPopulationSnapshotButton->setFont(font);

        gridLayout_7->addWidget(observeCellPopulationSnapshotButton, 8, 0, 1, 1);

        checkBoxEnablePOVOutput = new QCheckBox(groupBox_3);
        checkBoxEnablePOVOutput->setObjectName(QStringLiteral("checkBoxEnablePOVOutput"));
#ifndef QT_NO_TOOLTIP
        checkBoxEnablePOVOutput->setToolTip(QStringLiteral(""));
#endif // QT_NO_TOOLTIP

        gridLayout_7->addWidget(checkBoxEnablePOVOutput, 3, 0, 1, 1);

        buttonWritePovray = new QPushButton(groupBox_3);
        buttonWritePovray->setObjectName(QStringLiteral("buttonWritePovray"));

        gridLayout_7->addWidget(buttonWritePovray, 9, 0, 1, 1);


        gridLayout_4->addWidget(groupBox_3, 0, 0, 1, 1);

        groupBox_4 = new QGroupBox(tab_2);
        groupBox_4->setObjectName(QStringLiteral("groupBox_4"));
        QSizePolicy sizePolicy7(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
        sizePolicy7.setHorizontalStretch(0);
        sizePolicy7.setVerticalStretch(0);
        sizePolicy7.setHeightForWidth(groupBox_4->sizePolicy().hasHeightForWidth());
        groupBox_4->setSizePolicy(sizePolicy7);
        horizontalLayout_5 = new QHBoxLayout(groupBox_4);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(2, 2, 2, 2);
        monolayerConsole_Simulation = new QTextBrowser(groupBox_4);
        monolayerConsole_Simulation->setObjectName(QStringLiteral("monolayerConsole_Simulation"));
        sizePolicy7.setHeightForWidth(monolayerConsole_Simulation->sizePolicy().hasHeightForWidth());
        monolayerConsole_Simulation->setSizePolicy(sizePolicy7);
        monolayerConsole_Simulation->setReadOnly(false);

        horizontalLayout_5->addWidget(monolayerConsole_Simulation);


        gridLayout_4->addWidget(groupBox_4, 0, 1, 1, 1);

        tabWidget->addTab(tab_2, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QStringLiteral("tab_3"));
        gridLayout_8 = new QGridLayout(tab_3);
        gridLayout_8->setObjectName(QStringLiteral("gridLayout_8"));
        groupBoxVisualizationScreen = new QGroupBox(tab_3);
        groupBoxVisualizationScreen->setObjectName(QStringLiteral("groupBoxVisualizationScreen"));
        QSizePolicy sizePolicy8(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy8.setHorizontalStretch(0);
        sizePolicy8.setVerticalStretch(0);
        sizePolicy8.setHeightForWidth(groupBoxVisualizationScreen->sizePolicy().hasHeightForWidth());
        groupBoxVisualizationScreen->setSizePolicy(sizePolicy8);
        groupBoxVisualizationScreen->setMinimumSize(QSize(280, 0));
        gridLayout_10 = new QGridLayout(groupBoxVisualizationScreen);
        gridLayout_10->setObjectName(QStringLiteral("gridLayout_10"));
        comboBoxScreenBackgroundColor = new QComboBox(groupBoxVisualizationScreen);
        comboBoxScreenBackgroundColor->setObjectName(QStringLiteral("comboBoxScreenBackgroundColor"));

        gridLayout_10->addWidget(comboBoxScreenBackgroundColor, 1, 1, 1, 1);

        label_27 = new QLabel(groupBoxVisualizationScreen);
        label_27->setObjectName(QStringLiteral("label_27"));

        gridLayout_10->addWidget(label_27, 1, 0, 1, 1);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_10->addItem(horizontalSpacer, 1, 2, 1, 1);

        checkBoxScreenShowGrid = new QCheckBox(groupBoxVisualizationScreen);
        checkBoxScreenShowGrid->setObjectName(QStringLiteral("checkBoxScreenShowGrid"));
        checkBoxScreenShowGrid->setEnabled(false);

        gridLayout_10->addWidget(checkBoxScreenShowGrid, 2, 0, 1, 1);


        gridLayout_8->addWidget(groupBoxVisualizationScreen, 0, 0, 2, 1);

        groupBoxPovRayOptions = new QGroupBox(tab_3);
        groupBoxPovRayOptions->setObjectName(QStringLiteral("groupBoxPovRayOptions"));
        gridLayout_12 = new QGridLayout(groupBoxPovRayOptions);
        gridLayout_12->setObjectName(QStringLiteral("gridLayout_12"));
        comboBox = new QComboBox(groupBoxPovRayOptions);
        comboBox->setObjectName(QStringLiteral("comboBox"));
        QSizePolicy sizePolicy9(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy9.setHorizontalStretch(0);
        sizePolicy9.setVerticalStretch(0);
        sizePolicy9.setHeightForWidth(comboBox->sizePolicy().hasHeightForWidth());
        comboBox->setSizePolicy(sizePolicy9);
        comboBox->setMinimumSize(QSize(220, 0));

        gridLayout_12->addWidget(comboBox, 0, 1, 1, 1);

        label_29 = new QLabel(groupBoxPovRayOptions);
        label_29->setObjectName(QStringLiteral("label_29"));

        gridLayout_12->addWidget(label_29, 0, 0, 1, 1);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_12->addItem(horizontalSpacer_4, 0, 2, 1, 1);


        gridLayout_8->addWidget(groupBoxPovRayOptions, 6, 0, 1, 1);

        groupBoxVisualizationParameter = new QGroupBox(tab_3);
        groupBoxVisualizationParameter->setObjectName(QStringLiteral("groupBoxVisualizationParameter"));
        VisualizationTreeView = new QTreeView(groupBoxVisualizationParameter);
        VisualizationTreeView->setObjectName(QStringLiteral("VisualizationTreeView"));
        VisualizationTreeView->setGeometry(QRect(5, 15, 991, 101));

        gridLayout_8->addWidget(groupBoxVisualizationParameter, 5, 0, 1, 1);

        saveMXFButton = new QPushButton(tab_3);
        saveMXFButton->setObjectName(QStringLiteral("saveMXFButton"));

        gridLayout_8->addWidget(saveMXFButton, 7, 0, 1, 1);

        groupBoxVisualizationElements = new QGroupBox(tab_3);
        groupBoxVisualizationElements->setObjectName(QStringLiteral("groupBoxVisualizationElements"));
        gridLayout_9 = new QGridLayout(groupBoxVisualizationElements);
        gridLayout_9->setObjectName(QStringLiteral("gridLayout_9"));
        label_24 = new QLabel(groupBoxVisualizationElements);
        label_24->setObjectName(QStringLiteral("label_24"));

        gridLayout_9->addWidget(label_24, 2, 0, 1, 1);

        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_9->addItem(horizontalSpacer_3, 2, 3, 1, 1);

        comboBoxCellStaining = new QComboBox(groupBoxVisualizationElements);
        comboBoxCellStaining->setObjectName(QStringLiteral("comboBoxCellStaining"));
        sizePolicy1.setHeightForWidth(comboBoxCellStaining->sizePolicy().hasHeightForWidth());
        comboBoxCellStaining->setSizePolicy(sizePolicy1);

        gridLayout_9->addWidget(comboBoxCellStaining, 2, 2, 1, 1);

        pushButtonVoxelize = new QPushButton(groupBoxVisualizationElements);
        pushButtonVoxelize->setObjectName(QStringLiteral("pushButtonVoxelize"));

        gridLayout_9->addWidget(pushButtonVoxelize, 3, 0, 1, 1);


        gridLayout_8->addWidget(groupBoxVisualizationElements, 4, 0, 1, 1);

        tabWidget->addTab(tab_3, QString());
        tab_4 = new QWidget();
        tab_4->setObjectName(QStringLiteral("tab_4"));
        groupBox = new QGroupBox(tab_4);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        groupBox->setGeometry(QRect(10, 15, 196, 176));
        quantificationRunTestAnalysisButton = new QToolButton(groupBox);
        quantificationRunTestAnalysisButton->setObjectName(QStringLiteral("quantificationRunTestAnalysisButton"));
        quantificationRunTestAnalysisButton->setGeometry(QRect(10, 120, 121, 19));
        formLayoutWidget = new QWidget(groupBox);
        formLayoutWidget->setObjectName(QStringLiteral("formLayoutWidget"));
        formLayoutWidget->setGeometry(QRect(10, 25, 163, 66));
        formLayout = new QFormLayout(formLayoutWidget);
        formLayout->setObjectName(QStringLiteral("formLayout"));
        formLayout->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
        formLayout->setContentsMargins(0, 0, 0, 0);
        discvor_VoxelSize = new QDoubleSpinBox(formLayoutWidget);
        discvor_VoxelSize->setObjectName(QStringLiteral("discvor_VoxelSize"));
        discvor_VoxelSize->setDecimals(2);
        discvor_VoxelSize->setMinimum(0.01);
        discvor_VoxelSize->setMaximum(1);
        discvor_VoxelSize->setSingleStep(0.01);
        discvor_VoxelSize->setValue(0.2);

        formLayout->setWidget(0, QFormLayout::FieldRole, discvor_VoxelSize);

        discvor_CellCutOffRadius = new QDoubleSpinBox(formLayoutWidget);
        discvor_CellCutOffRadius->setObjectName(QStringLiteral("discvor_CellCutOffRadius"));
        discvor_CellCutOffRadius->setDecimals(2);
        discvor_CellCutOffRadius->setMinimum(0.5);
        discvor_CellCutOffRadius->setMaximum(2);
        discvor_CellCutOffRadius->setSingleStep(0.1);
        discvor_CellCutOffRadius->setValue(0.7);

        formLayout->setWidget(1, QFormLayout::FieldRole, discvor_CellCutOffRadius);

        label_3 = new QLabel(formLayoutWidget);
        label_3->setObjectName(QStringLiteral("label_3"));

        formLayout->setWidget(0, QFormLayout::LabelRole, label_3);

        label_4 = new QLabel(formLayoutWidget);
        label_4->setObjectName(QStringLiteral("label_4"));

        formLayout->setWidget(1, QFormLayout::LabelRole, label_4);

        discvor_EnableDebugOutput = new QCheckBox(groupBox);
        discvor_EnableDebugOutput->setObjectName(QStringLiteral("discvor_EnableDebugOutput"));
        discvor_EnableDebugOutput->setGeometry(QRect(10, 95, 88, 17));
        discvor_EnableDebugOutput->setChecked(true);
        quantificationDVA_ProgressBar = new QProgressBar(groupBox);
        quantificationDVA_ProgressBar->setObjectName(QStringLiteral("quantificationDVA_ProgressBar"));
        quantificationDVA_ProgressBar->setGeometry(QRect(10, 150, 181, 16));
        sizePolicy4.setHeightForWidth(quantificationDVA_ProgressBar->sizePolicy().hasHeightForWidth());
        quantificationDVA_ProgressBar->setSizePolicy(sizePolicy4);
        quantificationDVA_ProgressBar->setValue(50);
        groupBox_5 = new QGroupBox(tab_4);
        groupBox_5->setObjectName(QStringLiteral("groupBox_5"));
        groupBox_5->setGeometry(QRect(220, 15, 791, 211));
        sizePolicy7.setHeightForWidth(groupBox_5->sizePolicy().hasHeightForWidth());
        groupBox_5->setSizePolicy(sizePolicy7);
        horizontalLayout_6 = new QHBoxLayout(groupBox_5);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        horizontalLayout_6->setContentsMargins(2, 2, 2, 2);
        monolayerConsole_Quantification = new QTextBrowser(groupBox_5);
        monolayerConsole_Quantification->setObjectName(QStringLiteral("monolayerConsole_Quantification"));
        sizePolicy7.setHeightForWidth(monolayerConsole_Quantification->sizePolicy().hasHeightForWidth());
        monolayerConsole_Quantification->setSizePolicy(sizePolicy7);
        monolayerConsole_Quantification->setReadOnly(false);

        horizontalLayout_6->addWidget(monolayerConsole_Quantification);

        tabWidget->addTab(tab_4, QString());

        gridLayout->addWidget(tabWidget, 0, 0, 1, 1);


        retranslateUi(monolayer);

        tabWidget->setCurrentIndex(2);
        comboBoxScenario->setCurrentIndex(1);
        comboBoxScreenBackgroundColor->setCurrentIndex(3);
        comboBoxCellStaining->setCurrentIndex(3);


        QMetaObject::connectSlotsByName(monolayer);
    } // setupUi

    void retranslateUi(QWidget *monolayer)
    {
        monolayer->setWindowTitle(QApplication::translate("monolayer", "Form", Q_NULLPTR));
        groupBoxParameters->setTitle(QApplication::translate("monolayer", "Parameters", Q_NULLPTR));
        pushButtonResetParametersToDefaults->setText(QApplication::translate("monolayer", "Reset to default values", Q_NULLPTR));
        pushButtonApplyParametersToAllExistingCells->setText(QApplication::translate("monolayer", "Apply parameters to all existing cells", Q_NULLPTR));
        groupBoxSubmodel->setTitle(QApplication::translate("monolayer", "Submodel", Q_NULLPTR));
        comboBoxSubModel->clear();
        comboBoxSubModel->insertItems(0, QStringList()
         << QApplication::translate("monolayer", "Strict monolayer", Q_NULLPTR)
         << QApplication::translate("monolayer", "Tumor spheroid", Q_NULLPTR)
        );
        groupBoxSimulation->setTitle(QApplication::translate("monolayer", "Simulation", Q_NULLPTR));
        label_28->setText(QApplication::translate("monolayer", "Name:", Q_NULLPTR));
        textSimulationName->setText(QApplication::translate("monolayer", "default", Q_NULLPTR));
        groupBoxScenario->setTitle(QApplication::translate("monolayer", "Scenario", Q_NULLPTR));
        comboBoxScenario->clear();
        comboBoxScenario->insertItems(0, QStringList()
         << QApplication::translate("monolayer", "None (File Open... )", Q_NULLPTR)
         << QApplication::translate("monolayer", "Single Cell Origin", Q_NULLPTR)
         << QApplication::translate("monolayer", "Population in Embedding Medium", Q_NULLPTR)
        );
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("monolayer", "Parameters", Q_NULLPTR));
        groupBox_2->setTitle(QApplication::translate("monolayer", "Simulation progress", Q_NULLPTR));
        label_23->setText(QApplication::translate("monolayer", "N:", Q_NULLPTR));
        label_22->setText(QApplication::translate("monolayer", "cells", Q_NULLPTR));
        label_7->setText(QApplication::translate("monolayer", "Time:", Q_NULLPTR));
        label_21->setText(QApplication::translate("monolayer", "days", Q_NULLPTR));
        groupBox_3->setTitle(QApplication::translate("monolayer", "Simulation control", Q_NULLPTR));
        buttonResetModel->setText(QApplication::translate("monolayer", "Init / Reset model (N=1, t=0)", Q_NULLPTR));
        label->setText(QApplication::translate("monolayer", "Simulate until", Q_NULLPTR));
        doubleSpinBoxSimulateUntil->setPrefix(QString());
        label_2->setText(QApplication::translate("monolayer", "days", Q_NULLPTR));
        checkBoxEnableObservation->setText(QApplication::translate("monolayer", "Enable observation", Q_NULLPTR));
        labelOutputObservablesEvery->setText(QApplication::translate("monolayer", "Output observables every", Q_NULLPTR));
        labelOutputObservablesEveryUnit->setText(QApplication::translate("monolayer", "days", Q_NULLPTR));
        buttonStartSimulation->setText(QApplication::translate("monolayer", "Start simulation", Q_NULLPTR));
        buttonAbortSimulation->setText(QApplication::translate("monolayer", "Abort simulation", Q_NULLPTR));
        observeCellPopulationSnapshotButton->setText(QApplication::translate("monolayer", "Write cell population snapshot", Q_NULLPTR));
        checkBoxEnablePOVOutput->setText(QApplication::translate("monolayer", "Enable POV output", Q_NULLPTR));
        buttonWritePovray->setText(QApplication::translate("monolayer", "Write povray", Q_NULLPTR));
        groupBox_4->setTitle(QApplication::translate("monolayer", "Console", Q_NULLPTR));
        monolayerConsole_Simulation->setHtml(QApplication::translate("monolayer", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Noto Sans'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:'MS Shell Dlg 2'; font-size:8pt;\"><br /></p></body></html>", Q_NULLPTR));
        tabWidget->setTabText(tabWidget->indexOf(tab_2), QApplication::translate("monolayer", "Simulation", Q_NULLPTR));
        groupBoxVisualizationScreen->setTitle(QApplication::translate("monolayer", "Screen", Q_NULLPTR));
        comboBoxScreenBackgroundColor->clear();
        comboBoxScreenBackgroundColor->insertItems(0, QStringList()
         << QApplication::translate("monolayer", "black", Q_NULLPTR)
         << QApplication::translate("monolayer", "dark grey", Q_NULLPTR)
         << QApplication::translate("monolayer", "light grey", Q_NULLPTR)
         << QApplication::translate("monolayer", "mauve", Q_NULLPTR)
         << QApplication::translate("monolayer", "white", Q_NULLPTR)
        );
        label_27->setText(QApplication::translate("monolayer", "Background color: ", Q_NULLPTR));
        checkBoxScreenShowGrid->setText(QApplication::translate("monolayer", "Show grid", Q_NULLPTR));
        groupBoxPovRayOptions->setTitle(QApplication::translate("monolayer", "PovRay options", Q_NULLPTR));
        comboBox->clear();
        comboBox->insertItems(0, QStringList()
         << QApplication::translate("monolayer", "Camera settings of current 3D view", Q_NULLPTR)
        );
        label_29->setText(QApplication::translate("monolayer", "View perspective:", Q_NULLPTR));
        groupBoxVisualizationParameter->setTitle(QApplication::translate("monolayer", "Visualization Parameter", Q_NULLPTR));
        saveMXFButton->setText(QApplication::translate("monolayer", "save as mxf format", Q_NULLPTR));
        groupBoxVisualizationElements->setTitle(QApplication::translate("monolayer", "Model elements", Q_NULLPTR));
        label_24->setText(QApplication::translate("monolayer", "Cells are colored", Q_NULLPTR));
        comboBoxCellStaining->clear();
        comboBoxCellStaining->insertItems(0, QStringList()
         << QApplication::translate("monolayer", "all white", Q_NULLPTR)
         << QApplication::translate("monolayer", "by proliferation status (white = proliferating, grey = quiescent)", Q_NULLPTR)
         << QApplication::translate("monolayer", "by absolute force in last simulation step (hue)", Q_NULLPTR)
         << QApplication::translate("monolayer", "by volume (red: smallest just after division, green: largest just before division)", Q_NULLPTR)
         << QApplication::translate("monolayer", "Lobule layer", Q_NULLPTR)
         << QApplication::translate("monolayer", "Lesion edge", Q_NULLPTR)
        );
        pushButtonVoxelize->setText(QApplication::translate("monolayer", "Voxelize", Q_NULLPTR));
        tabWidget->setTabText(tabWidget->indexOf(tab_3), QApplication::translate("monolayer", "Visualization", Q_NULLPTR));
        groupBox->setTitle(QApplication::translate("monolayer", "Discrete Voronoi Analysis", Q_NULLPTR));
        quantificationRunTestAnalysisButton->setText(QApplication::translate("monolayer", "Run test analysis", Q_NULLPTR));
#ifndef QT_NO_ACCESSIBILITY
        discvor_VoxelSize->setAccessibleName(QString());
#endif // QT_NO_ACCESSIBILITY
        discvor_VoxelSize->setPrefix(QString());
        label_3->setText(QApplication::translate("monolayer", "Voxel size", Q_NULLPTR));
        label_4->setText(QApplication::translate("monolayer", "Cell cutoff radius", Q_NULLPTR));
        discvor_EnableDebugOutput->setText(QApplication::translate("monolayer", "DebugOutput", Q_NULLPTR));
        groupBox_5->setTitle(QApplication::translate("monolayer", "Console", Q_NULLPTR));
        monolayerConsole_Quantification->setHtml(QApplication::translate("monolayer", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Noto Sans'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:'MS Shell Dlg 2'; font-size:8pt;\"><br /></p></body></html>", Q_NULLPTR));
        tabWidget->setTabText(tabWidget->indexOf(tab_4), QApplication::translate("monolayer", "Quantification", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class monolayer: public Ui_monolayer {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MONOLAYER_H
