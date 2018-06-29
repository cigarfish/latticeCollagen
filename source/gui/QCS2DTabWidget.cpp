
#include "QCS2DTabWidget.h"

#include <QGridLayout>
#include <QPushButton>
#include <QFileDialog>

#include "QCS2DDisplay.h"


/*!
  \brief Constructor

  Setting up the widget's elements.
*/
QCS2DTabWidget::QCS2DTabWidget(QWidget *parent)
  : QWidget(parent)
{
  setupUi(this);

  connect(new2DViewButton,SIGNAL(clicked()),this,SLOT(createNewDisplay()));
  connect(loadImageButton,SIGNAL(clicked()),this,SLOT(loadImageSlot()));
}


QCS2DTabWidget::~QCS2DTabWidget()
{}


/*!
  \brief Slot to create a new 2D display
 */
void
QCS2DTabWidget::createNewDisplay()
{
  QCS2DDisplay * display = new QCS2DDisplay();
  display->resize(800,600);
  display->show();
  // TODO: put display in some container, like std::vector
}


/*!
  \brief Slot to handle image loading.
*/
void
QCS2DTabWidget::loadImageSlot()
{
  QString fileName = QFileDialog::getOpenFileName(this,tr("Open Image"),"",tr("Images (*.png *.tif *.jpg *.jpeg)"));

  if ( !fileName.isEmpty() )
    mpDisplay->openImage(fileName);
}
