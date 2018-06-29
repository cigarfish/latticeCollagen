/********************************************************************************
** Form generated from reading UI file 'Vascularization.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_VASCULARIZATION_H
#define UI_VASCULARIZATION_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Vascularization
{
public:
    QGridLayout *gridLayout;
    QPushButton *playButton;
    QPushButton *initButton;
    QPushButton *loadButton;
    QPushButton *showButton;

    void setupUi(QWidget *Vascularization)
    {
        if (Vascularization->objectName().isEmpty())
            Vascularization->setObjectName(QStringLiteral("Vascularization"));
        Vascularization->resize(1425, 704);
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(Vascularization->sizePolicy().hasHeightForWidth());
        Vascularization->setSizePolicy(sizePolicy);
        gridLayout = new QGridLayout(Vascularization);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        playButton = new QPushButton(Vascularization);
        playButton->setObjectName(QStringLiteral("playButton"));

        gridLayout->addWidget(playButton, 4, 0, 1, 1);

        initButton = new QPushButton(Vascularization);
        initButton->setObjectName(QStringLiteral("initButton"));

        gridLayout->addWidget(initButton, 1, 0, 1, 1);

        loadButton = new QPushButton(Vascularization);
        loadButton->setObjectName(QStringLiteral("loadButton"));

        gridLayout->addWidget(loadButton, 2, 0, 1, 1);

        showButton = new QPushButton(Vascularization);
        showButton->setObjectName(QStringLiteral("showButton"));

        gridLayout->addWidget(showButton, 3, 0, 1, 1);


        retranslateUi(Vascularization);

        QMetaObject::connectSlotsByName(Vascularization);
    } // setupUi

    void retranslateUi(QWidget *Vascularization)
    {
        Vascularization->setWindowTitle(QApplication::translate("Vascularization", "Form", Q_NULLPTR));
        playButton->setText(QApplication::translate("Vascularization", "play", Q_NULLPTR));
        initButton->setText(QApplication::translate("Vascularization", "init", Q_NULLPTR));
        loadButton->setText(QApplication::translate("Vascularization", "load", Q_NULLPTR));
        showButton->setText(QApplication::translate("Vascularization", "show", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Vascularization: public Ui_Vascularization {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_VASCULARIZATION_H
