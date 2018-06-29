
#ifndef QCS_3DTAB_WIDGET_H
#define QCS_3DTAB_WIDGET_H

#include "ui_QCS3DTabWidget.h"

#include <vector>

/*!
  \brief Widget for the 3D display handling

  Features a GLDisplay, creates a new GLDisplay when clicking a PushButton.
  Contains controls for changing color of the display's objects.
 */
class QCS3DTabWidget : public QWidget, private Ui::QCS3DTabWidget
{
    Q_OBJECT

public:
    QCS3DTabWidget(QWidget * parent=0);
    ~QCS3DTabWidget();

signals:
    void changeFGColor(QColor &);

private slots:
    void createNewDisplay();
    void arenaUpdate();

    //void initSim();
    void startSim1();
    //void startSim10();
    void startSim100();
    //void resetSim();


private:
    QColor mColor;
    std::vector<QCSGLDisplay *> mOtherDisplays;
};


#endif // QCS_3DTAB_WIDGET_H
