
#ifndef QCS_2DTAB_WIDGET_H
#define QCS_2DTAB_WIDGET_H

#include "ui_QCS2DTabWidget.h"


/*!
  \brief Widget for the 2D display handling

  Create a 2Display on clicking a PushButton.
 */
class QCS2DTabWidget : public QWidget, private Ui::QCS2DTabWidget
{
  Q_OBJECT

 public:
  QCS2DTabWidget(QWidget *parent=0);
  ~QCS2DTabWidget();

 private slots:
  void createNewDisplay();
  void loadImageSlot();
};

#endif // QCS_2DTAB_WIDGET_H
