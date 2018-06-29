/********************************************************************************
** Form generated from reading UI file 'etcUIDemo.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ETCUIDEMO_H
#define UI_ETCUIDEMO_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLCDNumber>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QTreeView>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_etcUIDemo
{
public:
    QGridLayout *gridLayout;
    QLabel *label_3;
    QLCDNumber *lcdNumber;
    QSlider *horizontalSlider;
    QProgressBar *progressBar;
    QDoubleSpinBox *doubleSpinBox;
    QComboBox *comboBoxSpin;
    QLineEdit *textInput;
    QLabel *label;
    QLineEdit *inputNumOnly;
    QLabel *label_2;
    QComboBox *comboBoxNumOnly;
    QRadioButton *radioButton;
    QGroupBox *groupBox;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QCheckBox *checkBeer;
    QCheckBox *checkCalpico;
    QCheckBox *checkCider;
    QCheckBox *checkClubMate;
    QCheckBox *checkCokeZero;
    QCheckBox *checkNone;
    QToolButton *helpButton;
    QPushButton *pushButton;
    QTreeView *exampleTreeView;

    void setupUi(QWidget *etcUIDemo)
    {
        if (etcUIDemo->objectName().isEmpty())
            etcUIDemo->setObjectName(QStringLiteral("etcUIDemo"));
        etcUIDemo->resize(702, 506);
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(etcUIDemo->sizePolicy().hasHeightForWidth());
        etcUIDemo->setSizePolicy(sizePolicy);
        gridLayout = new QGridLayout(etcUIDemo);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        label_3 = new QLabel(etcUIDemo);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_3, 0, 0, 1, 1);

        lcdNumber = new QLCDNumber(etcUIDemo);
        lcdNumber->setObjectName(QStringLiteral("lcdNumber"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Minimum);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(lcdNumber->sizePolicy().hasHeightForWidth());
        lcdNumber->setSizePolicy(sizePolicy1);

        gridLayout->addWidget(lcdNumber, 6, 0, 1, 1);

        horizontalSlider = new QSlider(etcUIDemo);
        horizontalSlider->setObjectName(QStringLiteral("horizontalSlider"));
        horizontalSlider->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(horizontalSlider, 6, 1, 1, 2);

        progressBar = new QProgressBar(etcUIDemo);
        progressBar->setObjectName(QStringLiteral("progressBar"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(progressBar->sizePolicy().hasHeightForWidth());
        progressBar->setSizePolicy(sizePolicy2);
        progressBar->setValue(24);

        gridLayout->addWidget(progressBar, 6, 3, 1, 1);

        doubleSpinBox = new QDoubleSpinBox(etcUIDemo);
        doubleSpinBox->setObjectName(QStringLiteral("doubleSpinBox"));
        sizePolicy2.setHeightForWidth(doubleSpinBox->sizePolicy().hasHeightForWidth());
        doubleSpinBox->setSizePolicy(sizePolicy2);

        gridLayout->addWidget(doubleSpinBox, 0, 1, 1, 1);

        comboBoxSpin = new QComboBox(etcUIDemo);
        comboBoxSpin->setObjectName(QStringLiteral("comboBoxSpin"));

        gridLayout->addWidget(comboBoxSpin, 0, 2, 1, 1);

        textInput = new QLineEdit(etcUIDemo);
        textInput->setObjectName(QStringLiteral("textInput"));
        sizePolicy2.setHeightForWidth(textInput->sizePolicy().hasHeightForWidth());
        textInput->setSizePolicy(sizePolicy2);

        gridLayout->addWidget(textInput, 1, 1, 1, 1);

        label = new QLabel(etcUIDemo);
        label->setObjectName(QStringLiteral("label"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy3);
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label, 1, 0, 1, 1);

        inputNumOnly = new QLineEdit(etcUIDemo);
        inputNumOnly->setObjectName(QStringLiteral("inputNumOnly"));
        sizePolicy2.setHeightForWidth(inputNumOnly->sizePolicy().hasHeightForWidth());
        inputNumOnly->setSizePolicy(sizePolicy2);

        gridLayout->addWidget(inputNumOnly, 2, 1, 1, 1);

        label_2 = new QLabel(etcUIDemo);
        label_2->setObjectName(QStringLiteral("label_2"));
        sizePolicy3.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy3);
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_2, 2, 0, 1, 1);

        comboBoxNumOnly = new QComboBox(etcUIDemo);
        comboBoxNumOnly->setObjectName(QStringLiteral("comboBoxNumOnly"));
        sizePolicy2.setHeightForWidth(comboBoxNumOnly->sizePolicy().hasHeightForWidth());
        comboBoxNumOnly->setSizePolicy(sizePolicy2);

        gridLayout->addWidget(comboBoxNumOnly, 2, 2, 1, 1);

        radioButton = new QRadioButton(etcUIDemo);
        radioButton->setObjectName(QStringLiteral("radioButton"));
        sizePolicy2.setHeightForWidth(radioButton->sizePolicy().hasHeightForWidth());
        radioButton->setSizePolicy(sizePolicy2);

        gridLayout->addWidget(radioButton, 3, 1, 1, 1);

        groupBox = new QGroupBox(etcUIDemo);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        sizePolicy3.setHeightForWidth(groupBox->sizePolicy().hasHeightForWidth());
        groupBox->setSizePolicy(sizePolicy3);
        groupBox->setMinimumSize(QSize(0, 0));
        verticalLayoutWidget = new QWidget(groupBox);
        verticalLayoutWidget->setObjectName(QStringLiteral("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(10, 30, 143, 161));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        checkBeer = new QCheckBox(verticalLayoutWidget);
        checkBeer->setObjectName(QStringLiteral("checkBeer"));

        verticalLayout->addWidget(checkBeer);

        checkCalpico = new QCheckBox(verticalLayoutWidget);
        checkCalpico->setObjectName(QStringLiteral("checkCalpico"));

        verticalLayout->addWidget(checkCalpico);

        checkCider = new QCheckBox(verticalLayoutWidget);
        checkCider->setObjectName(QStringLiteral("checkCider"));

        verticalLayout->addWidget(checkCider);

        checkClubMate = new QCheckBox(verticalLayoutWidget);
        checkClubMate->setObjectName(QStringLiteral("checkClubMate"));

        verticalLayout->addWidget(checkClubMate);

        checkCokeZero = new QCheckBox(verticalLayoutWidget);
        checkCokeZero->setObjectName(QStringLiteral("checkCokeZero"));

        verticalLayout->addWidget(checkCokeZero);

        checkNone = new QCheckBox(verticalLayoutWidget);
        checkNone->setObjectName(QStringLiteral("checkNone"));

        verticalLayout->addWidget(checkNone);


        gridLayout->addWidget(groupBox, 0, 3, 4, 1);

        helpButton = new QToolButton(etcUIDemo);
        helpButton->setObjectName(QStringLiteral("helpButton"));

        gridLayout->addWidget(helpButton, 7, 0, 1, 1);

        pushButton = new QPushButton(etcUIDemo);
        pushButton->setObjectName(QStringLiteral("pushButton"));

        gridLayout->addWidget(pushButton, 8, 0, 1, 1);

        exampleTreeView = new QTreeView(etcUIDemo);
        exampleTreeView->setObjectName(QStringLiteral("exampleTreeView"));

        gridLayout->addWidget(exampleTreeView, 7, 1, 2, 3);


        retranslateUi(etcUIDemo);
        QObject::connect(horizontalSlider, SIGNAL(valueChanged(int)), lcdNumber, SLOT(display(int)));
        QObject::connect(horizontalSlider, SIGNAL(valueChanged(int)), progressBar, SLOT(setValue(int)));

        QMetaObject::connectSlotsByName(etcUIDemo);
    } // setupUi

    void retranslateUi(QWidget *etcUIDemo)
    {
        etcUIDemo->setWindowTitle(QApplication::translate("etcUIDemo", "Form", Q_NULLPTR));
        label_3->setText(QApplication::translate("etcUIDemo", "Spinning Numbers", Q_NULLPTR));
        textInput->setText(QApplication::translate("etcUIDemo", "Input", Q_NULLPTR));
        label->setText(QApplication::translate("etcUIDemo", "Input", Q_NULLPTR));
        inputNumOnly->setText(QApplication::translate("etcUIDemo", "1.00e-7", Q_NULLPTR));
        label_2->setText(QApplication::translate("etcUIDemo", "Numbers only", Q_NULLPTR));
        radioButton->setText(QApplication::translate("etcUIDemo", "Radio?", Q_NULLPTR));
        groupBox->setTitle(QApplication::translate("etcUIDemo", "Favorite Beverage", Q_NULLPTR));
        checkBeer->setText(QApplication::translate("etcUIDemo", "Beer", Q_NULLPTR));
        checkCalpico->setText(QApplication::translate("etcUIDemo", "Calpico", Q_NULLPTR));
        checkCider->setText(QApplication::translate("etcUIDemo", "Cider", Q_NULLPTR));
        checkClubMate->setText(QApplication::translate("etcUIDemo", "Club Mate", Q_NULLPTR));
        checkCokeZero->setText(QApplication::translate("etcUIDemo", "Coca-Cola Zero", Q_NULLPTR));
        checkNone->setText(QApplication::translate("etcUIDemo", "None of these", Q_NULLPTR));
        helpButton->setText(QApplication::translate("etcUIDemo", "Help", Q_NULLPTR));
        pushButton->setText(QApplication::translate("etcUIDemo", "Push Me", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class etcUIDemo: public Ui_etcUIDemo {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ETCUIDEMO_H
