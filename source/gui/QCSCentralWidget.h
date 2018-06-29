#ifndef QCS_CENTRALWIDGET_H
#define QCS_CENTRALWIDGET_H

#include <QTabWidget>

/*!
  \brief TabWidget as the central CellSys widget
 */
class QCSCentralWidget : public QTabWidget
{
 public:
  QCSCentralWidget(QWidget *parent);

  ~QCSCentralWidget();
};


#endif // QCS_CENTRALWIDGET_H
