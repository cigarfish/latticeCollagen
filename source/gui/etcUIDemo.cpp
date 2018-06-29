
#include "etcUIDemo.h"

#include <QMessageBox>
#include <iostream>

#include "../tools/parameters/CSParameterContext.h"
#include "../tools/parameters/CSParameterChoice.h"


const QString etcUIDemo::bevString[] = {"Beer","Calpico","Cider","Club Mate","Coke Zero","None"};

/*!
  \brief Constructor
  \param parent Parent Widget
*/
etcUIDemo::etcUIDemo(QWidget *parent)
  : QWidget(parent),
    etcDisplay(NULL)
{
  // setting up the Qt Designer code
  // connecting SIGNAL valueChanged() of the horizontal slider
  // with SLOT display() of lcdNumber and SLOT setValue() of progressBar
  // is done in the ui code.
  setupUi(this);

  // setting a validator in order to only accept numbers (including dot, e, + and -);
  // inputNumOnly is of class QLineEdit.
  inputNumOnly->setValidator(new QDoubleValidator(NULL));

  // put the checkboxes in a buttongroup in order to make choice exclusive,
  // i.e. only one of the boxes in the buttongroup can be ticked.
  QButtonGroup * choiceGroup = new QButtonGroup(this);
  choiceGroup->addButton(checkBeer);
  choiceGroup->addButton(checkCalpico);
  choiceGroup->addButton(checkCider);
  choiceGroup->addButton(checkClubMate);
  choiceGroup->addButton(checkCokeZero);
  choiceGroup->addButton(checkNone);

  // default choice: beer!
  mBeverage = Beer;
  checkBeer->setChecked(true);

  // connecting to the slots (inline! s. header) that merely change the value of
  // mBeverage to the corresponding choice.
  connect(checkBeer,SIGNAL(clicked(bool)),this,SLOT(choiceBeer(bool)));
  connect(checkCalpico,SIGNAL(clicked(bool)),this,SLOT(choiceCalpico(bool)));
  connect(checkCider,SIGNAL(clicked(bool)),this,SLOT(choiceCider(bool)));
  connect(checkClubMate,SIGNAL(clicked(bool)),this,SLOT(choiceGingerBeer(bool)));
  connect(checkCokeZero,SIGNAL(clicked(bool)),this,SLOT(choiceWater(bool)));
  connect(checkNone,SIGNAL(clicked(bool)),this,SLOT(choiceNone(bool)));

  // setting up the possibility to change the value of the other Widgets that use numbers
  // inputNumOnly, doubleSpinBox, horizontalSlider
  comboBoxSpin->addItem("Change...");
  comboBoxSpin->addItem("Change Numbers Only Input");
  comboBoxSpin->addItem("Change Slider Value");
  comboBoxSpin->addItem("Change Both");

  // changeFromSpin() takes the value of doubleSpinBox to set the values of inputNumOnly,
  // horizontalSpinBox, or both, according to the index of the chosen item, which is given
  // as the int argument in signal and slot.
  connect(comboBoxSpin,SIGNAL(currentIndexChanged(int)),this,SLOT(changeFromSpin(int)));

  // this time the same thing for inputNumOnly, changing the respective other values.
  comboBoxNumOnly->addItem("Change...");
  comboBoxNumOnly->addItem("Change Spin Box Value");
  comboBoxNumOnly->addItem("Change Slider Value");
  comboBoxNumOnly->addItem("Change Both");

  connect(comboBoxNumOnly,SIGNAL(currentIndexChanged(int)),this,SLOT(changeFromNum(int)));

  // the radio button's toggled() signal is connected to an inline function (s. header)
  // that changes the value of mRadio.
  mRadio = false;
  connect(radioButton,SIGNAL(toggled(bool)),this,SLOT(toggleRadio(bool)));

  // putting progressbar to zero, the default seems not to be that low.
  progressBar->setValue(0);

  // help button connected to the function that will create the help dialog
  connect(helpButton,SIGNAL(clicked()),this,SLOT(showHelp()));

  // pushButton connected to a function that creates a message which
  // displays the values of the widgets (other than progressBar, horizontalSlider or lcdNumber)
  connect(pushButton,SIGNAL(clicked()),this,SLOT(showData()));


  // Example of a parameter tree and its visualization:
  // The root of the tree:
  CSParameterContext * contextRoot = new CSParameterContext("Root"); // name will not be shown

  // a choice between three items:
  const char *choiceStrings[] = { "First", "Second", "Third" };
  // create the choice container:
  CSParameterChoice * choice = new CSParameterChoice( choiceStrings, 3, 0 );
  // add it to the root context
  contextRoot->addParameter("Your choice", CSParameter::Choice, choice, "" );

  // first sub-context
  CSParameterContext * firstContext = contextRoot->addContext( "first" );
  // prepare the data memory for the parameter and initialize it:
  mParameterFile = "";
  firstContext->addParameter( "File Name Parameter", CSParameter::FileName,
                              &mParameterFile, "" );

  // second sub-context
  CSParameterContext * secondContext = contextRoot->addContext( "second" );
  // prepare the data memory for the parameter and initialize it:
  mParameterDouble = 0;
  // add the parameter to the second sub-context
  secondContext->addParameter( "Double Parameter", CSParameter::Double,
                               &mParameterDouble, "degrees" );

  // the third sub-context
  CSParameterContext * thirdContext = contextRoot->addContext( "third" );
  // prepare the data memory for the parameter and initialize it:
  mParameterBool = true;
  // add the parameter to the third sub-context
  thirdContext->addParameter( "Boolean", CSParameter::Bool, &mParameterBool, "");

  // This sets the choice item in control over the 'visibility' of the three sub-contexts
  // in the GUI, i.e. if the first item in the choice is chosen, only the first subcontext is
  // active, and the other two are disabled.
  choice->setControlSubContexts();
  contextRoot->setupGUI( exampleTreeView );
}


etcUIDemo::~etcUIDemo()
{}


/*!
  \brief Display the widgets' values in a pop-up window.
 */
void
etcUIDemo::showData()
{
  if (etcDisplay)
    delete etcDisplay;

  etcDisplay = new QWidget();
  QGridLayout * layout = new QGridLayout(etcDisplay);

  int row = 0;

  QLabel * spinLabel = new QLabel(QString("Spinning Numbers: "));
  layout->addWidget(spinLabel,row,0);
  QLabel * spinValue = new QLabel(QString::number(doubleSpinBox->value()));
  layout->addWidget(spinValue,row,1);

  ++row;

  QLabel * outLineInputLabel = new QLabel(QString("input:"));
  layout->addWidget(outLineInputLabel,row,0);
  QLabel * outLineInputValue = new QLabel(textInput->text());
  layout->addWidget(outLineInputValue,row,1);

  ++row;

  QLabel * outInputNumOnlyLabel = new QLabel(QString("Numbers only:"));
  layout->addWidget(outInputNumOnlyLabel);
  QLabel * outInputNumOnlyValue = new QLabel(inputNumOnly->text());
  layout->addWidget(outInputNumOnlyValue);

  ++row;

  QLabel * beverageLabel = new QLabel(QString("Favorite Beverage:"));
  layout->addWidget(beverageLabel,row,0);
  QLabel * beverageValue = new QLabel(bevString[mBeverage]);
  layout->addWidget(beverageValue,row,1);

  ++row;

  QLabel * radioLabel = new QLabel("Radio?");
  layout->addWidget(radioLabel,row,0);
  QLabel * radioValue = new QLabel(mRadio?"yes":"no");
  layout->addWidget(radioValue,row,1);

  etcDisplay->show();
}


/*!
  \brief Create and display a MessageBox with a help text.
 */
void
etcUIDemo::showHelp()
{
  QMessageBox * helpWidget = new QMessageBox();

  helpWidget->setText("This is a little zoo of Widgets some of which can influence others\n \
Feel free to play around with the comboBoxes.  The Push Me button will show a widget with the settings of the text fields and buttons");
  helpWidget->setIcon(QMessageBox::Information);

  // QMessageBoxes are blocking!!
  helpWidget->exec();
}


/*!
  \brief Take the value of QSpinBox doubleSpinBox and change other widgets' display value.
  \param index Index of chosen combo box item. 
*/
void
etcUIDemo::changeFromSpin(int index)
{
  double value = doubleSpinBox->value();

  switch((changeOthers)index)
    {
    case NoChoice:
      break;
    case OtherNumber:
      inputNumOnly->setText(QString::number(value));
      break;
    case Both:
      inputNumOnly->setText(QString::number(value));
    case Slider:
      value = (value > horizontalSlider->maximum()) ? horizontalSlider->maximum() : value;
      horizontalSlider->setValue(value);
      break;
    }
}


/*!
  \brief Take the value of QLineEdit inputNumOnly and change other widgets' display value.
  \param index Index of chosen combo box item. 
*/
void
etcUIDemo::changeFromNum(int index)
{
  bool ok = false;
  double value = inputNumOnly->text().toDouble(&ok);

  switch((changeOthers)index)
    {
    case NoChoice:
      break;
    case OtherNumber:
      doubleSpinBox->setValue(value);
      break;
    case Both:
      doubleSpinBox->setValue(value);
    case Slider:
      value = (value > horizontalSlider->maximum()) ? horizontalSlider->maximum() : value;
      horizontalSlider->setValue(value);
      break;
    }
}
