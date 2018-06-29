/********************************************************************************
** Form generated from reading UI file 'ODEsTab.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ODESTAB_H
#define UI_ODESTAB_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFrame>
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
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QTreeView>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ODEsTab
{
public:
    QGridLayout *gridLayout;
    QTabWidget *tabODEs;
    QWidget *tab;
    QGridLayout *gridLayout_3;
    QGroupBox *groupBox_7;
    QGridLayout *gridLayout_11;
    QLabel *label_28;
    QLineEdit *textSimulationName;
    QSpacerItem *verticalSpacer_3;
    QGroupBox *groupBox_6;
    QComboBox *comboBoxSubModel;
    QGroupBox *groupBoxParameters;
    QGridLayout *gridLayout_5;
    QGridLayout *gridLayout_2;
    QLabel *label_5;
    QDoubleSpinBox *doubleSpinBoxCellDiameter;
    QLabel *label_6;
    QLabel *label_16;
    QDoubleSpinBox *doubleSpinBoxCycleTime;
    QLabel *label_15;
    QLabel *label_19;
    QDoubleSpinBox *doubleSpinBoxYoungModulus;
    QLabel *label_20;
    QLabel *label_8;
    QLabel *label_9;
    QDoubleSpinBox *doubleSpinBoxDiffusionConstant;
    QLabel *label_10;
    QLabel *label_17;
    QDoubleSpinBox *doubleSpinBoxCellCellAdhesion;
    QLabel *label_18;
    QLabel *label_11;
    QDoubleSpinBox *doubleSpinBoxCellCellGamma;
    QLabel *label_12;
    QLabel *label_14;
    QDoubleSpinBox *doubleSpinBoxCellECMGamma;
    QLabel *label_13;
    QDoubleSpinBox *doubleSpinBoxCycleTimeSD;
    QLabel *label_26;
    QLabel *label_25;
    QDoubleSpinBox *doubleSpinBoxPoissonNumber;
    QPushButton *pushButtonApplyParametersToAllExistingCells;
    QGroupBox *groupBox_9;
    QGridLayout *gridLayout_13;
    QDoubleSpinBox *doubleSpinBox_3;
    QLabel *label_33;
    QDoubleSpinBox *doubleSpinBox;
    QDoubleSpinBox *doubleSpinBox_2;
    QLabel *label_32;
    QLabel *label_31;
    QLabel *label_34;
    QCheckBox *checkBoxEnableEvolution;
    QGroupBox *groupBox_10;
    QVBoxLayout *verticalLayout_3;
    QTreeView *parameterTreeView;
    QHBoxLayout *horizontalLayout_7;
    QLabel *label_35;
    QLineEdit *lineEdit_ApoTrigger;
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
    QLabel *label_3;
    QDoubleSpinBox *doubleSpinBoxObserveEvery;
    QLabel *label_4;
    QPushButton *buttonStartSimulation;
    QPushButton *buttonAbortSimulation;
    QPushButton *observeCellPopulationSnapshotButton;
    QCheckBox *checkBox;
    QPushButton *buttonWritePovray;
    QGroupBox *groupBox_4;
    QHBoxLayout *horizontalLayout_5;
    QTextBrowser *monolayerConsole;
    QWidget *tab_3;
    QGridLayout *gridLayout_8;
    QSpacerItem *verticalSpacer_6;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_10;
    QComboBox *comboBoxScreenBackgroundColor;
    QLabel *label_27;
    QSpacerItem *horizontalSpacer;
    QCheckBox *checkBoxScreenShowGrid;
    QGroupBox *groupBox_5;
    QGridLayout *gridLayout_9;
    QComboBox *comboBoxCellStaining;
    QLabel *label_24;
    QSpacerItem *horizontalSpacer_3;
    QGroupBox *groupBox_8;
    QGridLayout *gridLayout_12;
    QComboBox *comboBox;
    QLabel *label_29;
    QSpacerItem *horizontalSpacer_4;
    QWidget *tab_4;
    QVBoxLayout *verticalLayout;
    QCheckBox *check_useODEs;
    QHBoxLayout *horizontalLayout_8;
    QHBoxLayout *horizontalLayout_6;
    QGroupBox *groupBox_11;
    QGridLayout *gridLayout_17;
    QGridLayout *gridLayout_15;
    QPushButton *browse_sbml;
    QLineEdit *editSBMLPath;
    QPushButton *import_sbml;
    QSpacerItem *horizontalSpacer_5;
    QGroupBox *groupBox_12;
    QVBoxLayout *verticalLayout_2;
    QTableWidget *tableWidget_eq;
    QTableWidget *tableWidget_param;
    QLabel *label_30;
    QGroupBox *groupBox_13;
    QHBoxLayout *horizontalLayout_9;
    QComboBox *comboBox_method;
    QFrame *line;
    QLabel *label_36;
    QLineEdit *lineEdit_Error;
    QFrame *line_2;
    QLabel *label_37;
    QLineEdit *lineEdit_RError;
    QFrame *line_3;
    QCheckBox *checkBox_useCompiled;
    QSpacerItem *horizontalSpacer_9;

    void setupUi(QWidget *ODEsTab)
    {
        if (ODEsTab->objectName().isEmpty())
            ODEsTab->setObjectName(QStringLiteral("ODEsTab"));
        ODEsTab->resize(1423, 746);
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(ODEsTab->sizePolicy().hasHeightForWidth());
        ODEsTab->setSizePolicy(sizePolicy);
        gridLayout = new QGridLayout(ODEsTab);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        tabODEs = new QTabWidget(ODEsTab);
        tabODEs->setObjectName(QStringLiteral("tabODEs"));
        tabODEs->setMinimumSize(QSize(300, 0));
        tab = new QWidget();
        tab->setObjectName(QStringLiteral("tab"));
        gridLayout_3 = new QGridLayout(tab);
        gridLayout_3->setObjectName(QStringLiteral("gridLayout_3"));
        groupBox_7 = new QGroupBox(tab);
        groupBox_7->setObjectName(QStringLiteral("groupBox_7"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(groupBox_7->sizePolicy().hasHeightForWidth());
        groupBox_7->setSizePolicy(sizePolicy1);
        groupBox_7->setMinimumSize(QSize(280, 0));
        gridLayout_11 = new QGridLayout(groupBox_7);
        gridLayout_11->setObjectName(QStringLiteral("gridLayout_11"));
        label_28 = new QLabel(groupBox_7);
        label_28->setObjectName(QStringLiteral("label_28"));

        gridLayout_11->addWidget(label_28, 1, 0, 1, 1);

        textSimulationName = new QLineEdit(groupBox_7);
        textSimulationName->setObjectName(QStringLiteral("textSimulationName"));
        textSimulationName->setEnabled(true);

        gridLayout_11->addWidget(textSimulationName, 1, 1, 1, 1);


        gridLayout_3->addWidget(groupBox_7, 0, 0, 1, 1);

        verticalSpacer_3 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_3->addItem(verticalSpacer_3, 4, 0, 1, 1);

        groupBox_6 = new QGroupBox(tab);
        groupBox_6->setObjectName(QStringLiteral("groupBox_6"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(groupBox_6->sizePolicy().hasHeightForWidth());
        groupBox_6->setSizePolicy(sizePolicy2);
        groupBox_6->setMinimumSize(QSize(0, 55));
        comboBoxSubModel = new QComboBox(groupBox_6);
        comboBoxSubModel->setObjectName(QStringLiteral("comboBoxSubModel"));
        comboBoxSubModel->setEnabled(true);
        comboBoxSubModel->setGeometry(QRect(10, 20, 261, 22));

        gridLayout_3->addWidget(groupBox_6, 1, 0, 2, 1);

        groupBoxParameters = new QGroupBox(tab);
        groupBoxParameters->setObjectName(QStringLiteral("groupBoxParameters"));
        groupBoxParameters->setEnabled(true);
        gridLayout_5 = new QGridLayout(groupBoxParameters);
        gridLayout_5->setObjectName(QStringLiteral("gridLayout_5"));
        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        label_5 = new QLabel(groupBoxParameters);
        label_5->setObjectName(QStringLiteral("label_5"));

        gridLayout_2->addWidget(label_5, 0, 1, 1, 1);

        doubleSpinBoxCellDiameter = new QDoubleSpinBox(groupBoxParameters);
        doubleSpinBoxCellDiameter->setObjectName(QStringLiteral("doubleSpinBoxCellDiameter"));
        doubleSpinBoxCellDiameter->setEnabled(false);
        doubleSpinBoxCellDiameter->setDecimals(1);
        doubleSpinBoxCellDiameter->setMinimum(1);
        doubleSpinBoxCellDiameter->setMaximum(100);
        doubleSpinBoxCellDiameter->setValue(1);

        gridLayout_2->addWidget(doubleSpinBoxCellDiameter, 0, 2, 1, 1);

        label_6 = new QLabel(groupBoxParameters);
        label_6->setObjectName(QStringLiteral("label_6"));

        gridLayout_2->addWidget(label_6, 0, 3, 1, 1);

        label_16 = new QLabel(groupBoxParameters);
        label_16->setObjectName(QStringLiteral("label_16"));

        gridLayout_2->addWidget(label_16, 1, 1, 1, 1);

        doubleSpinBoxCycleTime = new QDoubleSpinBox(groupBoxParameters);
        doubleSpinBoxCycleTime->setObjectName(QStringLiteral("doubleSpinBoxCycleTime"));
        doubleSpinBoxCycleTime->setEnabled(false);
        doubleSpinBoxCycleTime->setDecimals(1);
        doubleSpinBoxCycleTime->setMinimum(1);
        doubleSpinBoxCycleTime->setMaximum(200);

        gridLayout_2->addWidget(doubleSpinBoxCycleTime, 1, 2, 1, 1);

        label_15 = new QLabel(groupBoxParameters);
        label_15->setObjectName(QStringLiteral("label_15"));

        gridLayout_2->addWidget(label_15, 1, 3, 1, 1);

        label_19 = new QLabel(groupBoxParameters);
        label_19->setObjectName(QStringLiteral("label_19"));

        gridLayout_2->addWidget(label_19, 3, 1, 1, 1);

        doubleSpinBoxYoungModulus = new QDoubleSpinBox(groupBoxParameters);
        doubleSpinBoxYoungModulus->setObjectName(QStringLiteral("doubleSpinBoxYoungModulus"));
        doubleSpinBoxYoungModulus->setEnabled(false);
        doubleSpinBoxYoungModulus->setDecimals(1);
        doubleSpinBoxYoungModulus->setMinimum(10);
        doubleSpinBoxYoungModulus->setMaximum(20000);
        doubleSpinBoxYoungModulus->setSingleStep(10);

        gridLayout_2->addWidget(doubleSpinBoxYoungModulus, 3, 2, 1, 1);

        label_20 = new QLabel(groupBoxParameters);
        label_20->setObjectName(QStringLiteral("label_20"));

        gridLayout_2->addWidget(label_20, 3, 3, 1, 1);

        label_8 = new QLabel(groupBoxParameters);
        label_8->setObjectName(QStringLiteral("label_8"));

        gridLayout_2->addWidget(label_8, 4, 1, 1, 1);

        label_9 = new QLabel(groupBoxParameters);
        label_9->setObjectName(QStringLiteral("label_9"));

        gridLayout_2->addWidget(label_9, 5, 1, 1, 1);

        doubleSpinBoxDiffusionConstant = new QDoubleSpinBox(groupBoxParameters);
        doubleSpinBoxDiffusionConstant->setObjectName(QStringLiteral("doubleSpinBoxDiffusionConstant"));
        doubleSpinBoxDiffusionConstant->setEnabled(false);
        doubleSpinBoxDiffusionConstant->setMaximum(10000);

        gridLayout_2->addWidget(doubleSpinBoxDiffusionConstant, 5, 2, 1, 1);

        label_10 = new QLabel(groupBoxParameters);
        label_10->setObjectName(QStringLiteral("label_10"));

        gridLayout_2->addWidget(label_10, 5, 3, 1, 1);

        label_17 = new QLabel(groupBoxParameters);
        label_17->setObjectName(QStringLiteral("label_17"));

        gridLayout_2->addWidget(label_17, 6, 1, 1, 1);

        doubleSpinBoxCellCellAdhesion = new QDoubleSpinBox(groupBoxParameters);
        doubleSpinBoxCellCellAdhesion->setObjectName(QStringLiteral("doubleSpinBoxCellCellAdhesion"));
        doubleSpinBoxCellCellAdhesion->setEnabled(false);
        doubleSpinBoxCellCellAdhesion->setMaximum(100);

        gridLayout_2->addWidget(doubleSpinBoxCellCellAdhesion, 6, 2, 1, 1);

        label_18 = new QLabel(groupBoxParameters);
        label_18->setObjectName(QStringLiteral("label_18"));

        gridLayout_2->addWidget(label_18, 6, 3, 1, 1);

        label_11 = new QLabel(groupBoxParameters);
        label_11->setObjectName(QStringLiteral("label_11"));

        gridLayout_2->addWidget(label_11, 7, 1, 1, 1);

        doubleSpinBoxCellCellGamma = new QDoubleSpinBox(groupBoxParameters);
        doubleSpinBoxCellCellGamma->setObjectName(QStringLiteral("doubleSpinBoxCellCellGamma"));
        doubleSpinBoxCellCellGamma->setEnabled(false);
        doubleSpinBoxCellCellGamma->setMinimum(0.1);
        doubleSpinBoxCellCellGamma->setMaximum(2000);
        doubleSpinBoxCellCellGamma->setSingleStep(0.1);

        gridLayout_2->addWidget(doubleSpinBoxCellCellGamma, 7, 2, 1, 1);

        label_12 = new QLabel(groupBoxParameters);
        label_12->setObjectName(QStringLiteral("label_12"));

        gridLayout_2->addWidget(label_12, 7, 3, 1, 1);

        label_14 = new QLabel(groupBoxParameters);
        label_14->setObjectName(QStringLiteral("label_14"));

        gridLayout_2->addWidget(label_14, 8, 1, 1, 1);

        doubleSpinBoxCellECMGamma = new QDoubleSpinBox(groupBoxParameters);
        doubleSpinBoxCellECMGamma->setObjectName(QStringLiteral("doubleSpinBoxCellECMGamma"));
        doubleSpinBoxCellECMGamma->setEnabled(false);
        doubleSpinBoxCellECMGamma->setMinimum(0.1);
        doubleSpinBoxCellECMGamma->setMaximum(2000);
        doubleSpinBoxCellECMGamma->setSingleStep(0.1);

        gridLayout_2->addWidget(doubleSpinBoxCellECMGamma, 8, 2, 1, 1);

        label_13 = new QLabel(groupBoxParameters);
        label_13->setObjectName(QStringLiteral("label_13"));

        gridLayout_2->addWidget(label_13, 8, 3, 1, 1);

        doubleSpinBoxCycleTimeSD = new QDoubleSpinBox(groupBoxParameters);
        doubleSpinBoxCycleTimeSD->setObjectName(QStringLiteral("doubleSpinBoxCycleTimeSD"));
        doubleSpinBoxCycleTimeSD->setEnabled(false);
        doubleSpinBoxCycleTimeSD->setMinimum(1);
        doubleSpinBoxCycleTimeSD->setMaximum(100);

        gridLayout_2->addWidget(doubleSpinBoxCycleTimeSD, 2, 2, 1, 1);

        label_26 = new QLabel(groupBoxParameters);
        label_26->setObjectName(QStringLiteral("label_26"));

        gridLayout_2->addWidget(label_26, 2, 3, 1, 1);

        label_25 = new QLabel(groupBoxParameters);
        label_25->setObjectName(QStringLiteral("label_25"));

        gridLayout_2->addWidget(label_25, 2, 1, 1, 1);

        doubleSpinBoxPoissonNumber = new QDoubleSpinBox(groupBoxParameters);
        doubleSpinBoxPoissonNumber->setObjectName(QStringLiteral("doubleSpinBoxPoissonNumber"));
        doubleSpinBoxPoissonNumber->setEnabled(false);
        doubleSpinBoxPoissonNumber->setMinimum(0.01);
        doubleSpinBoxPoissonNumber->setMaximum(0.5);
        doubleSpinBoxPoissonNumber->setSingleStep(0.01);

        gridLayout_2->addWidget(doubleSpinBoxPoissonNumber, 4, 2, 1, 1);


        gridLayout_5->addLayout(gridLayout_2, 1, 1, 1, 1);

        pushButtonApplyParametersToAllExistingCells = new QPushButton(groupBoxParameters);
        pushButtonApplyParametersToAllExistingCells->setObjectName(QStringLiteral("pushButtonApplyParametersToAllExistingCells"));
        pushButtonApplyParametersToAllExistingCells->setEnabled(false);

        gridLayout_5->addWidget(pushButtonApplyParametersToAllExistingCells, 2, 1, 1, 1);


        gridLayout_3->addWidget(groupBoxParameters, 3, 0, 1, 1);

        groupBox_9 = new QGroupBox(tab);
        groupBox_9->setObjectName(QStringLiteral("groupBox_9"));
        gridLayout_13 = new QGridLayout(groupBox_9);
        gridLayout_13->setObjectName(QStringLiteral("gridLayout_13"));
        doubleSpinBox_3 = new QDoubleSpinBox(groupBox_9);
        doubleSpinBox_3->setObjectName(QStringLiteral("doubleSpinBox_3"));

        gridLayout_13->addWidget(doubleSpinBox_3, 2, 6, 1, 1);

        label_33 = new QLabel(groupBox_9);
        label_33->setObjectName(QStringLiteral("label_33"));

        gridLayout_13->addWidget(label_33, 2, 5, 1, 1);

        doubleSpinBox = new QDoubleSpinBox(groupBox_9);
        doubleSpinBox->setObjectName(QStringLiteral("doubleSpinBox"));

        gridLayout_13->addWidget(doubleSpinBox, 2, 2, 1, 1);

        doubleSpinBox_2 = new QDoubleSpinBox(groupBox_9);
        doubleSpinBox_2->setObjectName(QStringLiteral("doubleSpinBox_2"));

        gridLayout_13->addWidget(doubleSpinBox_2, 2, 4, 1, 1);

        label_32 = new QLabel(groupBox_9);
        label_32->setObjectName(QStringLiteral("label_32"));

        gridLayout_13->addWidget(label_32, 2, 3, 1, 1);

        label_31 = new QLabel(groupBox_9);
        label_31->setObjectName(QStringLiteral("label_31"));

        gridLayout_13->addWidget(label_31, 2, 1, 1, 1);

        label_34 = new QLabel(groupBox_9);
        label_34->setObjectName(QStringLiteral("label_34"));

        gridLayout_13->addWidget(label_34, 2, 7, 1, 1);

        checkBoxEnableEvolution = new QCheckBox(groupBox_9);
        checkBoxEnableEvolution->setObjectName(QStringLiteral("checkBoxEnableEvolution"));

        gridLayout_13->addWidget(checkBoxEnableEvolution, 2, 0, 1, 1);


        gridLayout_3->addWidget(groupBox_9, 1, 2, 1, 1);

        groupBox_10 = new QGroupBox(tab);
        groupBox_10->setObjectName(QStringLiteral("groupBox_10"));
        verticalLayout_3 = new QVBoxLayout(groupBox_10);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        parameterTreeView = new QTreeView(groupBox_10);
        parameterTreeView->setObjectName(QStringLiteral("parameterTreeView"));

        verticalLayout_3->addWidget(parameterTreeView);

        horizontalLayout_7 = new QHBoxLayout();
#ifndef Q_OS_MAC
        horizontalLayout_7->setSpacing(-1);
#endif
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        horizontalLayout_7->setSizeConstraint(QLayout::SetMinimumSize);
        label_35 = new QLabel(groupBox_10);
        label_35->setObjectName(QStringLiteral("label_35"));
        label_35->setEnabled(false);

        horizontalLayout_7->addWidget(label_35);

        lineEdit_ApoTrigger = new QLineEdit(groupBox_10);
        lineEdit_ApoTrigger->setObjectName(QStringLiteral("lineEdit_ApoTrigger"));
        lineEdit_ApoTrigger->setEnabled(false);

        horizontalLayout_7->addWidget(lineEdit_ApoTrigger);


        verticalLayout_3->addLayout(horizontalLayout_7);


        gridLayout_3->addWidget(groupBox_10, 3, 2, 1, 1);

        tabODEs->addTab(tab, QString());
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
        doubleSpinBoxSimulationProgressTime->setMaximum(1000);

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
        doubleSpinBoxSimulateUntil->setMaximum(1000);
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
        label_3 = new QLabel(groupBox_3);
        label_3->setObjectName(QStringLiteral("label_3"));

        horizontalLayout_2->addWidget(label_3);

        doubleSpinBoxObserveEvery = new QDoubleSpinBox(groupBox_3);
        doubleSpinBoxObserveEvery->setObjectName(QStringLiteral("doubleSpinBoxObserveEvery"));
        doubleSpinBoxObserveEvery->setDecimals(2);
        doubleSpinBoxObserveEvery->setMinimum(0);
        doubleSpinBoxObserveEvery->setMaximum(100);
        doubleSpinBoxObserveEvery->setSingleStep(0.1);
        doubleSpinBoxObserveEvery->setValue(0.1);

        horizontalLayout_2->addWidget(doubleSpinBoxObserveEvery);

        label_4 = new QLabel(groupBox_3);
        label_4->setObjectName(QStringLiteral("label_4"));

        horizontalLayout_2->addWidget(label_4);


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

        checkBox = new QCheckBox(groupBox_3);
        checkBox->setObjectName(QStringLiteral("checkBox"));
#ifndef QT_NO_TOOLTIP
        checkBox->setToolTip(QStringLiteral(""));
#endif // QT_NO_TOOLTIP

        gridLayout_7->addWidget(checkBox, 3, 0, 1, 1);

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
        horizontalLayout_5->setContentsMargins(2, 2, 2, 2);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        monolayerConsole = new QTextBrowser(groupBox_4);
        monolayerConsole->setObjectName(QStringLiteral("monolayerConsole"));
        sizePolicy7.setHeightForWidth(monolayerConsole->sizePolicy().hasHeightForWidth());
        monolayerConsole->setSizePolicy(sizePolicy7);
        monolayerConsole->setReadOnly(false);

        horizontalLayout_5->addWidget(monolayerConsole);


        gridLayout_4->addWidget(groupBox_4, 0, 1, 1, 1);

        tabODEs->addTab(tab_2, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QStringLiteral("tab_3"));
        gridLayout_8 = new QGridLayout(tab_3);
        gridLayout_8->setObjectName(QStringLiteral("gridLayout_8"));
        verticalSpacer_6 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_8->addItem(verticalSpacer_6, 5, 0, 1, 1);

        groupBox = new QGroupBox(tab_3);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        QSizePolicy sizePolicy8(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy8.setHorizontalStretch(0);
        sizePolicy8.setVerticalStretch(0);
        sizePolicy8.setHeightForWidth(groupBox->sizePolicy().hasHeightForWidth());
        groupBox->setSizePolicy(sizePolicy8);
        groupBox->setMinimumSize(QSize(280, 0));
        gridLayout_10 = new QGridLayout(groupBox);
        gridLayout_10->setObjectName(QStringLiteral("gridLayout_10"));
        comboBoxScreenBackgroundColor = new QComboBox(groupBox);
        comboBoxScreenBackgroundColor->setObjectName(QStringLiteral("comboBoxScreenBackgroundColor"));

        gridLayout_10->addWidget(comboBoxScreenBackgroundColor, 1, 1, 1, 1);

        label_27 = new QLabel(groupBox);
        label_27->setObjectName(QStringLiteral("label_27"));

        gridLayout_10->addWidget(label_27, 1, 0, 1, 1);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_10->addItem(horizontalSpacer, 1, 2, 1, 1);

        checkBoxScreenShowGrid = new QCheckBox(groupBox);
        checkBoxScreenShowGrid->setObjectName(QStringLiteral("checkBoxScreenShowGrid"));
        checkBoxScreenShowGrid->setEnabled(false);

        gridLayout_10->addWidget(checkBoxScreenShowGrid, 2, 0, 1, 1);


        gridLayout_8->addWidget(groupBox, 0, 0, 2, 1);

        groupBox_5 = new QGroupBox(tab_3);
        groupBox_5->setObjectName(QStringLiteral("groupBox_5"));
        gridLayout_9 = new QGridLayout(groupBox_5);
        gridLayout_9->setObjectName(QStringLiteral("gridLayout_9"));
        comboBoxCellStaining = new QComboBox(groupBox_5);
        comboBoxCellStaining->setObjectName(QStringLiteral("comboBoxCellStaining"));
        sizePolicy2.setHeightForWidth(comboBoxCellStaining->sizePolicy().hasHeightForWidth());
        comboBoxCellStaining->setSizePolicy(sizePolicy2);

        gridLayout_9->addWidget(comboBoxCellStaining, 2, 2, 1, 1);

        label_24 = new QLabel(groupBox_5);
        label_24->setObjectName(QStringLiteral("label_24"));

        gridLayout_9->addWidget(label_24, 2, 0, 1, 1);

        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_9->addItem(horizontalSpacer_3, 2, 3, 1, 1);


        gridLayout_8->addWidget(groupBox_5, 3, 0, 1, 1);

        groupBox_8 = new QGroupBox(tab_3);
        groupBox_8->setObjectName(QStringLiteral("groupBox_8"));
        gridLayout_12 = new QGridLayout(groupBox_8);
        gridLayout_12->setObjectName(QStringLiteral("gridLayout_12"));
        comboBox = new QComboBox(groupBox_8);
        comboBox->setObjectName(QStringLiteral("comboBox"));
        QSizePolicy sizePolicy9(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy9.setHorizontalStretch(0);
        sizePolicy9.setVerticalStretch(0);
        sizePolicy9.setHeightForWidth(comboBox->sizePolicy().hasHeightForWidth());
        comboBox->setSizePolicy(sizePolicy9);
        comboBox->setMinimumSize(QSize(220, 0));

        gridLayout_12->addWidget(comboBox, 0, 1, 1, 1);

        label_29 = new QLabel(groupBox_8);
        label_29->setObjectName(QStringLiteral("label_29"));

        gridLayout_12->addWidget(label_29, 0, 0, 1, 1);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_12->addItem(horizontalSpacer_4, 0, 2, 1, 1);


        gridLayout_8->addWidget(groupBox_8, 4, 0, 1, 1);

        tabODEs->addTab(tab_3, QString());
        tab_4 = new QWidget();
        tab_4->setObjectName(QStringLiteral("tab_4"));
        verticalLayout = new QVBoxLayout(tab_4);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        check_useODEs = new QCheckBox(tab_4);
        check_useODEs->setObjectName(QStringLiteral("check_useODEs"));

        verticalLayout->addWidget(check_useODEs);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));

        verticalLayout->addLayout(horizontalLayout_8);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        groupBox_11 = new QGroupBox(tab_4);
        groupBox_11->setObjectName(QStringLiteral("groupBox_11"));
        groupBox_11->setEnabled(false);
        gridLayout_17 = new QGridLayout(groupBox_11);
        gridLayout_17->setObjectName(QStringLiteral("gridLayout_17"));
        gridLayout_15 = new QGridLayout();
        gridLayout_15->setObjectName(QStringLiteral("gridLayout_15"));
        browse_sbml = new QPushButton(groupBox_11);
        browse_sbml->setObjectName(QStringLiteral("browse_sbml"));

        gridLayout_15->addWidget(browse_sbml, 0, 1, 1, 1);

        editSBMLPath = new QLineEdit(groupBox_11);
        editSBMLPath->setObjectName(QStringLiteral("editSBMLPath"));

        gridLayout_15->addWidget(editSBMLPath, 0, 0, 1, 1);

        import_sbml = new QPushButton(groupBox_11);
        import_sbml->setObjectName(QStringLiteral("import_sbml"));

        gridLayout_15->addWidget(import_sbml, 1, 1, 1, 1);


        gridLayout_17->addLayout(gridLayout_15, 0, 0, 1, 1);

        horizontalSpacer_5 = new QSpacerItem(700, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_17->addItem(horizontalSpacer_5, 0, 1, 1, 1);


        horizontalLayout_6->addWidget(groupBox_11);


        verticalLayout->addLayout(horizontalLayout_6);

        groupBox_12 = new QGroupBox(tab_4);
        groupBox_12->setObjectName(QStringLiteral("groupBox_12"));
        groupBox_12->setEnabled(false);
        verticalLayout_2 = new QVBoxLayout(groupBox_12);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        tableWidget_eq = new QTableWidget(groupBox_12);
        if (tableWidget_eq->columnCount() < 2)
            tableWidget_eq->setColumnCount(2);
        tableWidget_eq->setObjectName(QStringLiteral("tableWidget_eq"));
        tableWidget_eq->setColumnCount(2);

        verticalLayout_2->addWidget(tableWidget_eq);

        tableWidget_param = new QTableWidget(groupBox_12);
        tableWidget_param->setObjectName(QStringLiteral("tableWidget_param"));
        tableWidget_param->setEnabled(false);
        tableWidget_param->setMaximumSize(QSize(16777215, 60));

        verticalLayout_2->addWidget(tableWidget_param);

        label_30 = new QLabel(groupBox_12);
        label_30->setObjectName(QStringLiteral("label_30"));

        verticalLayout_2->addWidget(label_30);


        verticalLayout->addWidget(groupBox_12);

        groupBox_13 = new QGroupBox(tab_4);
        groupBox_13->setObjectName(QStringLiteral("groupBox_13"));
        groupBox_13->setEnabled(false);
        groupBox_13->setMinimumSize(QSize(0, 70));
        horizontalLayout_9 = new QHBoxLayout(groupBox_13);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        comboBox_method = new QComboBox(groupBox_13);
        comboBox_method->setObjectName(QStringLiteral("comboBox_method"));
        comboBox_method->setMaximumSize(QSize(150, 16777215));

        horizontalLayout_9->addWidget(comboBox_method);

        line = new QFrame(groupBox_13);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::VLine);
        line->setFrameShadow(QFrame::Sunken);

        horizontalLayout_9->addWidget(line);

        label_36 = new QLabel(groupBox_13);
        label_36->setObjectName(QStringLiteral("label_36"));

        horizontalLayout_9->addWidget(label_36);

        lineEdit_Error = new QLineEdit(groupBox_13);
        lineEdit_Error->setObjectName(QStringLiteral("lineEdit_Error"));
        lineEdit_Error->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_9->addWidget(lineEdit_Error);

        line_2 = new QFrame(groupBox_13);
        line_2->setObjectName(QStringLiteral("line_2"));
        line_2->setFrameShape(QFrame::VLine);
        line_2->setFrameShadow(QFrame::Sunken);

        horizontalLayout_9->addWidget(line_2);

        label_37 = new QLabel(groupBox_13);
        label_37->setObjectName(QStringLiteral("label_37"));

        horizontalLayout_9->addWidget(label_37);

        lineEdit_RError = new QLineEdit(groupBox_13);
        lineEdit_RError->setObjectName(QStringLiteral("lineEdit_RError"));
        lineEdit_RError->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_9->addWidget(lineEdit_RError);

        line_3 = new QFrame(groupBox_13);
        line_3->setObjectName(QStringLiteral("line_3"));
        line_3->setFrameShape(QFrame::VLine);
        line_3->setFrameShadow(QFrame::Sunken);

        horizontalLayout_9->addWidget(line_3);

        checkBox_useCompiled = new QCheckBox(groupBox_13);
        checkBox_useCompiled->setObjectName(QStringLiteral("checkBox_useCompiled"));

        horizontalLayout_9->addWidget(checkBox_useCompiled);

        horizontalSpacer_9 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_9->addItem(horizontalSpacer_9);


        verticalLayout->addWidget(groupBox_13);

        tabODEs->addTab(tab_4, QString());

        gridLayout->addWidget(tabODEs, 0, 1, 1, 1);


        retranslateUi(ODEsTab);

        tabODEs->setCurrentIndex(2);
        comboBoxCellStaining->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(ODEsTab);
    } // setupUi

    void retranslateUi(QWidget *ODEsTab)
    {
        ODEsTab->setWindowTitle(QApplication::translate("ODEsTab", "Form", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        tabODEs->setToolTip(QApplication::translate("ODEsTab", "<html><head/><body><p><br/></p></body></html>", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        groupBox_7->setTitle(QApplication::translate("ODEsTab", "Simulation", Q_NULLPTR));
        label_28->setText(QApplication::translate("ODEsTab", "Name:", Q_NULLPTR));
        textSimulationName->setText(QApplication::translate("ODEsTab", "default", Q_NULLPTR));
        groupBox_6->setTitle(QApplication::translate("ODEsTab", "Submodel", Q_NULLPTR));
        comboBoxSubModel->clear();
        comboBoxSubModel->insertItems(0, QStringList()
         << QApplication::translate("ODEsTab", "Strict monolayer", Q_NULLPTR)
         << QApplication::translate("ODEsTab", "Tumor spheroid", Q_NULLPTR)
        );
        groupBoxParameters->setTitle(QApplication::translate("ODEsTab", "Cell population (newly created cells)", Q_NULLPTR));
        label_5->setText(QApplication::translate("ODEsTab", "Cell diameter", Q_NULLPTR));
        doubleSpinBoxCellDiameter->setSpecialValueText(QString());
        doubleSpinBoxCellDiameter->setPrefix(QString());
        doubleSpinBoxCellDiameter->setSuffix(QString());
        label_6->setText(QApplication::translate("ODEsTab", "\302\265m", Q_NULLPTR));
        label_16->setText(QApplication::translate("ODEsTab", "Cycle time", Q_NULLPTR));
        doubleSpinBoxCycleTime->setSpecialValueText(QString());
        label_15->setText(QApplication::translate("ODEsTab", "hours", Q_NULLPTR));
        label_19->setText(QApplication::translate("ODEsTab", "Young modulus", Q_NULLPTR));
        doubleSpinBoxYoungModulus->setSpecialValueText(QString());
        label_20->setText(QApplication::translate("ODEsTab", "Pa", Q_NULLPTR));
        label_8->setText(QApplication::translate("ODEsTab", "Poisson number", Q_NULLPTR));
        label_9->setText(QApplication::translate("ODEsTab", "Diffusion constant", Q_NULLPTR));
        doubleSpinBoxDiffusionConstant->setSpecialValueText(QString());
        label_10->setText(QApplication::translate("ODEsTab", "10^-16 cm\302\262/s", Q_NULLPTR));
        label_17->setText(QApplication::translate("ODEsTab", "Cell-Cell adhesion", Q_NULLPTR));
        doubleSpinBoxCellCellAdhesion->setSpecialValueText(QString());
        label_18->setText(QApplication::translate("ODEsTab", "1e15/m\302\262", Q_NULLPTR));
        label_11->setText(QApplication::translate("ODEsTab", "Cell-Cell gamma", Q_NULLPTR));
        doubleSpinBoxCellCellGamma->setSpecialValueText(QString());
        label_12->setText(QApplication::translate("ODEsTab", "10^7 Ns/m\302\263", Q_NULLPTR));
        label_14->setText(QApplication::translate("ODEsTab", "Cell-ECM gamma", Q_NULLPTR));
        doubleSpinBoxCellECMGamma->setSpecialValueText(QString());
        label_13->setText(QApplication::translate("ODEsTab", "10^7 Ns/m\302\263", Q_NULLPTR));
        doubleSpinBoxCycleTimeSD->setSpecialValueText(QString());
        label_26->setText(QApplication::translate("ODEsTab", "hours", Q_NULLPTR));
        label_25->setText(QApplication::translate("ODEsTab", "Cycle time SD", Q_NULLPTR));
        doubleSpinBoxPoissonNumber->setSpecialValueText(QString());
        pushButtonApplyParametersToAllExistingCells->setText(QApplication::translate("ODEsTab", "Apply parameters to all existing cells", Q_NULLPTR));
        groupBox_9->setTitle(QApplication::translate("ODEsTab", "Evolution / Mutation at cell division", Q_NULLPTR));
        label_33->setText(QApplication::translate("ODEsTab", "Max:", Q_NULLPTR));
        label_32->setText(QApplication::translate("ODEsTab", "Min:", Q_NULLPTR));
        label_31->setText(QApplication::translate("ODEsTab", "DeltaMax:", Q_NULLPTR));
        label_34->setText(QApplication::translate("ODEsTab", "10^-16 cm\302\262/s", Q_NULLPTR));
        checkBoxEnableEvolution->setText(QApplication::translate("ODEsTab", "Cell diffusion constant:", Q_NULLPTR));
        groupBox_10->setTitle(QApplication::translate("ODEsTab", "Parameters", Q_NULLPTR));
        label_35->setText(QApplication::translate("ODEsTab", "Apoptosis triggered if ", Q_NULLPTR));
        lineEdit_ApoTrigger->setInputMask(QString());
        lineEdit_ApoTrigger->setPlaceholderText(QApplication::translate("ODEsTab", "[MAPK] > 3.5", Q_NULLPTR));
        tabODEs->setTabText(tabODEs->indexOf(tab), QApplication::translate("ODEsTab", "Parameters", Q_NULLPTR));
        groupBox_2->setTitle(QApplication::translate("ODEsTab", "Simulation progress", Q_NULLPTR));
        label_23->setText(QApplication::translate("ODEsTab", "N:", Q_NULLPTR));
        label_22->setText(QApplication::translate("ODEsTab", "cells", Q_NULLPTR));
        label_7->setText(QApplication::translate("ODEsTab", "Time:", Q_NULLPTR));
        label_21->setText(QApplication::translate("ODEsTab", "days", Q_NULLPTR));
        groupBox_3->setTitle(QApplication::translate("ODEsTab", "Simulation control", Q_NULLPTR));
        buttonResetModel->setText(QApplication::translate("ODEsTab", "Init / Reset model (N=1, t=0)", Q_NULLPTR));
        label->setText(QApplication::translate("ODEsTab", "Simulate until", Q_NULLPTR));
        doubleSpinBoxSimulateUntil->setPrefix(QString());
        label_2->setText(QApplication::translate("ODEsTab", "days", Q_NULLPTR));
        checkBoxEnableObservation->setText(QApplication::translate("ODEsTab", "Enable observation", Q_NULLPTR));
        label_3->setText(QApplication::translate("ODEsTab", "Output observables every", Q_NULLPTR));
        label_4->setText(QApplication::translate("ODEsTab", "days", Q_NULLPTR));
        buttonStartSimulation->setText(QApplication::translate("ODEsTab", "Start simulation", Q_NULLPTR));
        buttonAbortSimulation->setText(QApplication::translate("ODEsTab", "Abort simulation", Q_NULLPTR));
        observeCellPopulationSnapshotButton->setText(QApplication::translate("ODEsTab", "Write cell population snapshot", Q_NULLPTR));
        checkBox->setText(QApplication::translate("ODEsTab", "Enable POV output", Q_NULLPTR));
        buttonWritePovray->setText(QApplication::translate("ODEsTab", "Write povray", Q_NULLPTR));
        groupBox_4->setTitle(QApplication::translate("ODEsTab", "Console", Q_NULLPTR));
        monolayerConsole->setHtml(QApplication::translate("ODEsTab", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Lucida Grande'; font-size:13pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:'MS Shell Dlg 2'; font-size:8pt;\"><br /></p></body></html>", Q_NULLPTR));
        tabODEs->setTabText(tabODEs->indexOf(tab_2), QApplication::translate("ODEsTab", "Simulation", Q_NULLPTR));
        groupBox->setTitle(QApplication::translate("ODEsTab", "Screen", Q_NULLPTR));
        comboBoxScreenBackgroundColor->clear();
        comboBoxScreenBackgroundColor->insertItems(0, QStringList()
         << QApplication::translate("ODEsTab", "black", Q_NULLPTR)
         << QApplication::translate("ODEsTab", "dark grey", Q_NULLPTR)
         << QApplication::translate("ODEsTab", "light grey", Q_NULLPTR)
         << QApplication::translate("ODEsTab", "white", Q_NULLPTR)
        );
        label_27->setText(QApplication::translate("ODEsTab", "Background color: ", Q_NULLPTR));
        checkBoxScreenShowGrid->setText(QApplication::translate("ODEsTab", "Show grid", Q_NULLPTR));
        groupBox_5->setTitle(QApplication::translate("ODEsTab", "Model elements", Q_NULLPTR));
        comboBoxCellStaining->clear();
        comboBoxCellStaining->insertItems(0, QStringList()
         << QApplication::translate("ODEsTab", "all white", Q_NULLPTR)
         << QApplication::translate("ODEsTab", "by proliferation status (white = proliferating, grey = quiescent)", Q_NULLPTR)
         << QApplication::translate("ODEsTab", "by absolute force in last simulation step (hue)", Q_NULLPTR)
         << QApplication::translate("ODEsTab", "by volume (white: smallest just after division, red: largest just before division)", Q_NULLPTR)
         << QApplication::translate("ODEsTab", "by internal state", Q_NULLPTR)
        );
        label_24->setText(QApplication::translate("ODEsTab", "Cells are colored", Q_NULLPTR));
        groupBox_8->setTitle(QApplication::translate("ODEsTab", "PovRay options", Q_NULLPTR));
        comboBox->clear();
        comboBox->insertItems(0, QStringList()
         << QApplication::translate("ODEsTab", "Camera settings of current 3D view", Q_NULLPTR)
        );
        label_29->setText(QApplication::translate("ODEsTab", "View perspective:", Q_NULLPTR));
        tabODEs->setTabText(tabODEs->indexOf(tab_3), QApplication::translate("ODEsTab", "Visualization", Q_NULLPTR));
        check_useODEs->setText(QApplication::translate("ODEsTab", "Use ODEs for intracellular pathways", Q_NULLPTR));
        groupBox_11->setTitle(QApplication::translate("ODEsTab", "SBML Model import", Q_NULLPTR));
        browse_sbml->setText(QApplication::translate("ODEsTab", "Browse", Q_NULLPTR));
        editSBMLPath->setInputMask(QString());
        editSBMLPath->setPlaceholderText(QApplication::translate("ODEsTab", "Enter the path to your SBML model", Q_NULLPTR));
        import_sbml->setText(QApplication::translate("ODEsTab", "Import", Q_NULLPTR));
        groupBox_12->setTitle(QApplication::translate("ODEsTab", "SBML Model summary", Q_NULLPTR));
        label_30->setText(QApplication::translate("ODEsTab", "Go to the Parameters tab to link the parameters with the ODEs", Q_NULLPTR));
        groupBox_13->setTitle(QApplication::translate("ODEsTab", "Integration Settings", Q_NULLPTR));
        comboBox_method->clear();
        comboBox_method->insertItems(0, QStringList()
         << QApplication::translate("ODEsTab", "BDF (CVODE stiff)", Q_NULLPTR)
         << QApplication::translate("ODEsTab", "Adams-Moulton (CVODE non-stiff)", Q_NULLPTR)
        );
        label_36->setText(QApplication::translate("ODEsTab", "Absolute Error:", Q_NULLPTR));
        lineEdit_Error->setText(QApplication::translate("ODEsTab", "1.e-18", Q_NULLPTR));
        lineEdit_Error->setPlaceholderText(QApplication::translate("ODEsTab", "1.e-18", Q_NULLPTR));
        label_37->setText(QApplication::translate("ODEsTab", "Relative Error:", Q_NULLPTR));
        lineEdit_RError->setText(QApplication::translate("ODEsTab", "1.e-10", Q_NULLPTR));
        lineEdit_RError->setPlaceholderText(QApplication::translate("ODEsTab", "1.e-10", Q_NULLPTR));
        checkBox_useCompiled->setText(QApplication::translate("ODEsTab", "Use compiled ODEs", Q_NULLPTR));
        tabODEs->setTabText(tabODEs->indexOf(tab_4), QApplication::translate("ODEsTab", "ODEs", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class ODEsTab: public Ui_ODEsTab {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ODESTAB_H
