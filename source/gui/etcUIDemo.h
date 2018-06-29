
#ifndef ETC_UI_DEMO_H
#define ETC_UI_DEMO_H

#include "ui_etcUIDemo.h"

/*!
  \brief The widget for the example ETC tab.

  This widget features a little overview of the different Qt
  widget control elements.  There's SpinBoxes and LineEdits,
  ComboBoxes and CheckBoxes, a RadioButton, a Slider, an
  LcdNumber and a progressBar.
  The progressBar and the LcdNumber are responsive to the
  change of the Slider, which was done in the Designer file.
  The Rest of the widgets interact with signal/slot pairs
  connected in the constructor of etcUIDemo.
  Simple accessor slot functions are kept inline.
  An enum is used for ease of extendability of values of
  choice for favorite beverages.
*/
class etcUIDemo : public QWidget, private Ui::etcUIDemo
{
  Q_OBJECT

 public:
  etcUIDemo(QWidget *parent=0);
  ~etcUIDemo();

 private slots:
  void showData();
  void showHelp();
  void changeFromSpin(int index);
  void changeFromNum(int index);

  void choiceBeer(bool) {mBeverage = Beer;};
  void choiceCalpico(bool) {mBeverage = Calpico;};
  void choiceCider(bool) {mBeverage = Cider;};
  void choiceGingerBeer(bool) {mBeverage = ClubMate;};
  void choiceWater(bool) {mBeverage = CokeZero;};
  void choiceNone(bool) {mBeverage = None;};

  void toggleRadio(bool yn) {mRadio=yn;};


 private:
  enum beverage
  {
    Beer=0,
    Calpico,
    Cider,
    ClubMate,
    CokeZero,
    None
  };

  static const QString bevString[];

  enum changeOthers {NoChoice=0,OtherNumber,Slider,Both};

  std::string mParameterFile;
  double mParameterDouble;
  bool mParameterBool;

  bool mRadio;
  beverage mBeverage;

  QWidget * etcDisplay;
};


#endif // ETC_UI_DEMO_H
