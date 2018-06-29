#ifndef QCS_MAIN_WINDOW_H
#define QCS_MAIN_WINDOW_H

#include <QMainWindow>

//#include "ui_QCSMainWindow.h"

class QCSGLDisplay;
class QCSCentralWidget;
class QCloseEvent;

/*!
  \brief The main window container

  Mainly used for constructing the MenuBar
*/
class QCSMainWindow : public QMainWindow//, private Ui::MainWindow
{
  Q_OBJECT

  friend class QCSCentralWidget;

 public:
  QCSMainWindow(QWidget *parent =0, Qt::WindowFlags flags =0);

  ~QCSMainWindow();

  QCSGLDisplay * getDisplay();


 protected:
  void setDisplay( QCSGLDisplay * );
  void closeEvent( QCloseEvent * );

 protected slots:
  void fileOpenDialog();
  void saveToFile();
  void saveToOtherFile();
  void messageNotImplementedYet();

 protected:
  QCSGLDisplay * mpActiveDisplay;
  QList<QCSGLDisplay *> mDisplays;

  QString mDataDirName;
  bool mSaved;
};

#endif // QCS_MAIN_WINDOW_H
