/********************************************************************************
** Form generated from reading UI file 'QCS3DTabWidget.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QCS3DTABWIDGET_H
#define UI_QCS3DTABWIDGET_H

#include <QCSGLDisplay.h>
#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_QCS3DTabWidget
{
public:
    QHBoxLayout *horizontalLayout;
    QFrame *frame;
    QVBoxLayout *verticalLayout;
    QPushButton *newDisplayButton;
    QCSGLDisplay *mpDisplay;

    void setupUi(QWidget *QCS3DTabWidget)
    {
        if (QCS3DTabWidget->objectName().isEmpty())
            QCS3DTabWidget->setObjectName(QStringLiteral("QCS3DTabWidget"));
        QCS3DTabWidget->resize(1037, 798);
        horizontalLayout = new QHBoxLayout(QCS3DTabWidget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        frame = new QFrame(QCS3DTabWidget);
        frame->setObjectName(QStringLiteral("frame"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(frame->sizePolicy().hasHeightForWidth());
        frame->setSizePolicy(sizePolicy);
        frame->setMaximumSize(QSize(110, 16777215));
        frame->setAutoFillBackground(true);
        frame->setFrameShape(QFrame::StyledPanel);
        frame->setFrameShadow(QFrame::Raised);
        verticalLayout = new QVBoxLayout(frame);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        newDisplayButton = new QPushButton(frame);
        newDisplayButton->setObjectName(QStringLiteral("newDisplayButton"));
        QSizePolicy sizePolicy1(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(newDisplayButton->sizePolicy().hasHeightForWidth());
        newDisplayButton->setSizePolicy(sizePolicy1);
        newDisplayButton->setMinimumSize(QSize(100, 0));

        verticalLayout->addWidget(newDisplayButton);


        horizontalLayout->addWidget(frame);

        mpDisplay = new QCSGLDisplay(QCS3DTabWidget);
        mpDisplay->setObjectName(QStringLiteral("mpDisplay"));
        QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(mpDisplay->sizePolicy().hasHeightForWidth());
        mpDisplay->setSizePolicy(sizePolicy2);
        frame->raise();

        horizontalLayout->addWidget(mpDisplay);


        retranslateUi(QCS3DTabWidget);

        QMetaObject::connectSlotsByName(QCS3DTabWidget);
    } // setupUi

    void retranslateUi(QWidget *QCS3DTabWidget)
    {
        QCS3DTabWidget->setWindowTitle(QApplication::translate("QCS3DTabWidget", "Form", Q_NULLPTR));
        newDisplayButton->setText(QApplication::translate("QCS3DTabWidget", "New Display", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class QCS3DTabWidget: public Ui_QCS3DTabWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QCS3DTABWIDGET_H
