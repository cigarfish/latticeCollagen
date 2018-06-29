/********************************************************************************
** Form generated from reading UI file 'About.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ABOUT_H
#define UI_ABOUT_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_About
{
public:
    QGridLayout *gridLayout;
    QGridLayout *aboutLayout;
    QLabel *labelEmpty;
    QLabel *labelWebsite;
    QLabel *labelCopyright;
    QLabel *labelDevelopers;
    QLabel *labelRelease;
    QTextBrowser *licenseBrowser;

    void setupUi(QWidget *About)
    {
        if (About->objectName().isEmpty())
            About->setObjectName(QStringLiteral("About"));
        About->resize(1425, 704);
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(About->sizePolicy().hasHeightForWidth());
        About->setSizePolicy(sizePolicy);
        gridLayout = new QGridLayout(About);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        aboutLayout = new QGridLayout();
        aboutLayout->setObjectName(QStringLiteral("aboutLayout"));
        labelEmpty = new QLabel(About);
        labelEmpty->setObjectName(QStringLiteral("labelEmpty"));

        aboutLayout->addWidget(labelEmpty, 1, 0, 1, 1);

        labelWebsite = new QLabel(About);
        labelWebsite->setObjectName(QStringLiteral("labelWebsite"));

        aboutLayout->addWidget(labelWebsite, 2, 0, 1, 1);

        labelCopyright = new QLabel(About);
        labelCopyright->setObjectName(QStringLiteral("labelCopyright"));

        aboutLayout->addWidget(labelCopyright, 4, 0, 1, 1);

        labelDevelopers = new QLabel(About);
        labelDevelopers->setObjectName(QStringLiteral("labelDevelopers"));

        aboutLayout->addWidget(labelDevelopers, 3, 0, 1, 1);

        labelRelease = new QLabel(About);
        labelRelease->setObjectName(QStringLiteral("labelRelease"));
        QFont font;
        font.setPointSize(28);
        labelRelease->setFont(font);

        aboutLayout->addWidget(labelRelease, 0, 0, 1, 1);

        licenseBrowser = new QTextBrowser(About);
        licenseBrowser->setObjectName(QStringLiteral("licenseBrowser"));

        aboutLayout->addWidget(licenseBrowser, 5, 0, 1, 1);


        gridLayout->addLayout(aboutLayout, 0, 0, 1, 1);


        retranslateUi(About);

        QMetaObject::connectSlotsByName(About);
    } // setupUi

    void retranslateUi(QWidget *About)
    {
        About->setWindowTitle(QApplication::translate("About", "Form", Q_NULLPTR));
        labelEmpty->setText(QString());
        labelWebsite->setText(QApplication::translate("About", "http://www.msysbio.com", Q_NULLPTR));
        labelCopyright->setText(QApplication::translate("About", "\302\251 Multicellular Systems Biology Group, IZBI, University Leipzig and MAMBA Group, INRIA Paris-Rocquencourt.", Q_NULLPTR));
        labelDevelopers->setText(QString());
        labelRelease->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class About: public Ui_About {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ABOUT_H
