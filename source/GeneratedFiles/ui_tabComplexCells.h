/********************************************************************************
** Form generated from reading UI file 'tabComplexCells.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_TABCOMPLEXCELLS_H
#define UI_TABCOMPLEXCELLS_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTreeView>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_complexCells
{
public:
    QGridLayout *gridLayout;
    QTabWidget *tabWidget;
    QWidget *simulation;
    QGroupBox *groupBox_7;
    QGridLayout *gridLayout_11;
    QPushButton *buttonNewDisplay;
    QGroupBox *groupBox_3;
    QGridLayout *gridLayout_7;
    QPushButton *buttonResetModel;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_3;
    QDoubleSpinBox *doubleSpinBoxObserveEvery;
    QLabel *label_4;
    QPushButton *buttonStartSimulation;
    QPushButton *buttonAbortSimulation;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QDoubleSpinBox *doubleSpinBoxSimulateUntil;
    QLabel *label_2;
    QSpacerItem *verticalSpacer;
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
    QWidget *parameter;
    QGroupBox *groupBox_10;
    QGridLayout *gridLayout_14;
    QPushButton *pushButtonResetParametersToDefaults;
    QTreeView *parameterTreeView;

    void setupUi(QWidget *complexCells)
    {
        if (complexCells->objectName().isEmpty())
            complexCells->setObjectName(QStringLiteral("complexCells"));
        complexCells->resize(1425, 499);
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(complexCells->sizePolicy().hasHeightForWidth());
        complexCells->setSizePolicy(sizePolicy);
        gridLayout = new QGridLayout(complexCells);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        tabWidget = new QTabWidget(complexCells);
        tabWidget->setObjectName(QStringLiteral("tabWidget"));
        simulation = new QWidget();
        simulation->setObjectName(QStringLiteral("simulation"));
        groupBox_7 = new QGroupBox(simulation);
        groupBox_7->setObjectName(QStringLiteral("groupBox_7"));
        groupBox_7->setGeometry(QRect(10, 10, 280, 81));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(groupBox_7->sizePolicy().hasHeightForWidth());
        groupBox_7->setSizePolicy(sizePolicy1);
        groupBox_7->setMinimumSize(QSize(280, 0));
        gridLayout_11 = new QGridLayout(groupBox_7);
        gridLayout_11->setObjectName(QStringLiteral("gridLayout_11"));
        buttonNewDisplay = new QPushButton(groupBox_7);
        buttonNewDisplay->setObjectName(QStringLiteral("buttonNewDisplay"));

        gridLayout_11->addWidget(buttonNewDisplay, 1, 0, 1, 2);

        groupBox_3 = new QGroupBox(simulation);
        groupBox_3->setObjectName(QStringLiteral("groupBox_3"));
        groupBox_3->setGeometry(QRect(10, 110, 329, 221));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::MinimumExpanding);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(groupBox_3->sizePolicy().hasHeightForWidth());
        groupBox_3->setSizePolicy(sizePolicy2);
        gridLayout_7 = new QGridLayout(groupBox_3);
        gridLayout_7->setObjectName(QStringLiteral("gridLayout_7"));
        buttonResetModel = new QPushButton(groupBox_3);
        buttonResetModel->setObjectName(QStringLiteral("buttonResetModel"));
        QSizePolicy sizePolicy3(QSizePolicy::Minimum, QSizePolicy::Maximum);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(buttonResetModel->sizePolicy().hasHeightForWidth());
        buttonResetModel->setSizePolicy(sizePolicy3);

        gridLayout_7->addWidget(buttonResetModel, 0, 0, 1, 1);

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


        gridLayout_7->addLayout(horizontalLayout_2, 2, 0, 1, 1);

        buttonStartSimulation = new QPushButton(groupBox_3);
        buttonStartSimulation->setObjectName(QStringLiteral("buttonStartSimulation"));
        buttonStartSimulation->setEnabled(false);

        gridLayout_7->addWidget(buttonStartSimulation, 3, 0, 1, 1);

        buttonAbortSimulation = new QPushButton(groupBox_3);
        buttonAbortSimulation->setObjectName(QStringLiteral("buttonAbortSimulation"));
        buttonAbortSimulation->setEnabled(false);

        gridLayout_7->addWidget(buttonAbortSimulation, 4, 0, 1, 1);

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

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_7->addItem(verticalSpacer, 5, 0, 1, 1);

        groupBox_2 = new QGroupBox(simulation);
        groupBox_2->setObjectName(QStringLiteral("groupBox_2"));
        groupBox_2->setGeometry(QRect(10, 340, 995, 96));
        QSizePolicy sizePolicy4(QSizePolicy::Preferred, QSizePolicy::Maximum);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(groupBox_2->sizePolicy().hasHeightForWidth());
        groupBox_2->setSizePolicy(sizePolicy4);
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
        QSizePolicy sizePolicy5(QSizePolicy::Expanding, QSizePolicy::MinimumExpanding);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(progressBar->sizePolicy().hasHeightForWidth());
        progressBar->setSizePolicy(sizePolicy5);
        progressBar->setValue(50);

        gridLayout_6->addWidget(progressBar, 2, 0, 1, 4);

        tabWidget->addTab(simulation, QString());
        parameter = new QWidget();
        parameter->setObjectName(QStringLiteral("parameter"));
        groupBox_10 = new QGroupBox(parameter);
        groupBox_10->setObjectName(QStringLiteral("groupBox_10"));
        groupBox_10->setGeometry(QRect(0, 10, 995, 441));
        gridLayout_14 = new QGridLayout(groupBox_10);
        gridLayout_14->setObjectName(QStringLiteral("gridLayout_14"));
        pushButtonResetParametersToDefaults = new QPushButton(groupBox_10);
        pushButtonResetParametersToDefaults->setObjectName(QStringLiteral("pushButtonResetParametersToDefaults"));

        gridLayout_14->addWidget(pushButtonResetParametersToDefaults, 2, 0, 1, 1);

        parameterTreeView = new QTreeView(groupBox_10);
        parameterTreeView->setObjectName(QStringLiteral("parameterTreeView"));

        gridLayout_14->addWidget(parameterTreeView, 1, 0, 1, 1);

        tabWidget->addTab(parameter, QString());

        gridLayout->addWidget(tabWidget, 0, 0, 1, 1);


        retranslateUi(complexCells);

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(complexCells);
    } // setupUi

    void retranslateUi(QWidget *complexCells)
    {
        complexCells->setWindowTitle(QApplication::translate("complexCells", "Form", Q_NULLPTR));
        groupBox_7->setTitle(QApplication::translate("complexCells", "Simulation", Q_NULLPTR));
        buttonNewDisplay->setText(QApplication::translate("complexCells", "new Display", Q_NULLPTR));
        groupBox_3->setTitle(QApplication::translate("complexCells", "Simulation control", Q_NULLPTR));
        buttonResetModel->setText(QApplication::translate("complexCells", "Init / Reset model (N=1, t=0)", Q_NULLPTR));
        label_3->setText(QApplication::translate("complexCells", "Output observables every", Q_NULLPTR));
        label_4->setText(QApplication::translate("complexCells", "days", Q_NULLPTR));
        buttonStartSimulation->setText(QApplication::translate("complexCells", "Start simulation", Q_NULLPTR));
        buttonAbortSimulation->setText(QApplication::translate("complexCells", "Abort simulation", Q_NULLPTR));
        label->setText(QApplication::translate("complexCells", "Simulate until", Q_NULLPTR));
        doubleSpinBoxSimulateUntil->setPrefix(QString());
        label_2->setText(QApplication::translate("complexCells", "days", Q_NULLPTR));
        groupBox_2->setTitle(QApplication::translate("complexCells", "Simulation progress", Q_NULLPTR));
        label_23->setText(QApplication::translate("complexCells", "N:", Q_NULLPTR));
        label_22->setText(QApplication::translate("complexCells", "cells", Q_NULLPTR));
        label_7->setText(QApplication::translate("complexCells", "Time:", Q_NULLPTR));
        label_21->setText(QApplication::translate("complexCells", "days", Q_NULLPTR));
        tabWidget->setTabText(tabWidget->indexOf(simulation), QApplication::translate("complexCells", "Simulation", Q_NULLPTR));
        groupBox_10->setTitle(QApplication::translate("complexCells", "Parameters", Q_NULLPTR));
        pushButtonResetParametersToDefaults->setText(QApplication::translate("complexCells", "Reset to default values", Q_NULLPTR));
        tabWidget->setTabText(tabWidget->indexOf(parameter), QApplication::translate("complexCells", "Parameter", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class complexCells: public Ui_complexCells {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TABCOMPLEXCELLS_H
