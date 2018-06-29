///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSParameterRangeWidget.h                                            //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-05-22 15:34:02                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef Q_CS_PARAMETER_RANGE_WIDGET_H
#define Q_CS_PARAMETER_RANGE_WIDGET_H

#include <QtGui>

#include <iostream>


class RangeWidget : public QWidget
{
 protected:

  RangeWidget( QWidget * parent, CSParameter::DataType type )
    : QWidget( parent ),
      mDataType( type )
  {
    layout = new QHBoxLayout();

    mpStartBox = new QDoubleSpinBox( this );
    hyphen = new QLabel( tr(" - "), this );
    mpEndBox = new QDoubleSpinBox( this );

    mpEndBox->setValue(7);
    layout->addWidget(mpStartBox);
    layout->addWidget(hyphen);
    layout->addWidget(mpEndBox);

    this->setLayout( layout );

    show();

    std::cout << "in RangeWidget\n";
  };

 public:
  static RangeWidget * create( QWidget * parent, CSParameter::DataType type )
  {
    switch ( type )
      {
      case CSParameter::Int:
      case CSParameter::Long:
      case CSParameter::Float:
      case CSParameter::Double:
        return new RangeWidget(parent, type);

      default:
        return NULL;
      }
  };

  /* void show() */
  /* { */
  /*   mpStartBox->show(); */
  /*   hyphen->show(); */
  /*   mpEndBox->show(); */
  /*   QWidget::show(); */
  /*   std::cout << "in RangeWidget::show() \n"; */
  /* }; */


  void setStart( double start )
  {
    mStart = start;
    
    switch ( mDataType )
      {
      case CSParameter::Int:
      case CSParameter::Long:
        mpStartBox->setValue( (int) mStart );
        break;

      default:
        mpStartBox->setValue( mStart );
      }
  };

  double start() { return mStart; };

  void setEnd( double end )
  {
    mEnd = end;
    
    switch ( mDataType )
      {
      case CSParameter::Int:
      case CSParameter::Long:
        mpEndBox->setValue( (int) end );
        break;

      default:
        mpEndBox->setValue( mEnd );
      }
  };

  double end() { return mEnd; };

private:
  QHBoxLayout * layout;
  QDoubleSpinBox * mpStartBox;
  QDoubleSpinBox * mpEndBox;
  QLabel   * hyphen;
  
  CSParameter::DataType mDataType;
  double mStart;
  double mEnd;
};


#endif //Q_CS_PARAMETER_RANGE_WIDGET_H
