///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSParameterFileNameEditor.cpp                                       //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-08-06 18:24:32                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "QCSParameterFileNameEditor.h"

#include <QHBoxLayout>
#include <QToolButton>
#include <QLineEdit>
#include <QString>
#include <QFileDialog>
#include <QFileInfo>


QCSParameterFileNameEditor::QCSParameterFileNameEditor(CSParameter * parm, QWidget * parent )
  : QWidget(parent),
    mpParameter(parm)
{
  this->setAutoFillBackground(1);

  mpLine = new QLineEdit();

  QToolButton * button = new QToolButton();
  button->setText( tr("...") );

  QHBoxLayout * layout = new QHBoxLayout();
  layout->setContentsMargins(0,0,0,0);
  layout->addWidget( mpLine );
  layout->addWidget( button );

  CSParameter::DataType type;
  if (mpParameter)
      type = mpParameter->dataType();
  else
      type = CSParameter::FileName;

  if ( type == CSParameter::FileName )
  {
      connect( button, SIGNAL( clicked() ), this, SLOT( showFileDialog() ) );
      connect( mpLine, SIGNAL( textChanged(const QString&) ), this, SLOT( validateFileName(const QString&) ) );
  }
  else if ( type == CSParameter::DirName )
  {
      connect( button, SIGNAL( clicked() ), this, SLOT( showDirDialog() ) );
      connect( mpLine, SIGNAL( textChanged(const QString&) ), this, SLOT( validateDirName(const QString&) ) );
  }

  this->setLayout( layout );

  setFocusProxy( button );
}


void
QCSParameterFileNameEditor::showFileDialog()
{
  QFileDialog dial( this, tr( "Choose file..." ) );
  if (mpParameter)
      dial.setNameFilter(tr(mpParameter->getAttribute("FileTypes").c_str()));

  if ( mName.size() )
    {
      QString qName(mName.c_str());
      QFileInfo intel( qName );
      QDir directory = intel.absoluteDir();
      if ( directory.exists() )
        dial.setDirectory( directory );
    }

  if ( dial.exec() ) 
    setText( dial.selectedFiles().at(0) );

  emit doneEditing();
}


void
QCSParameterFileNameEditor::showDirDialog()
{
  QFileDialog dial( this, tr( "Choose directory..." ) );

  dial.setOptions( QFileDialog::ShowDirsOnly );
  dial.setFileMode( QFileDialog::Directory );

  QFileInfo finfo( mpLine->text() );
  if ( finfo.exists() && finfo.isDir() )
      dial.setDirectory( finfo.absoluteDir() );

  if ( dial.exec() )
    setText( dial.selectedFiles().at(0) );  
}


void
QCSParameterFileNameEditor::setText(std::string & name)
{
  mName = name;
  mpLine->setText( QString( name.c_str() ) );
}


void
QCSParameterFileNameEditor::setText(QString name)
{
  mName = name.toStdString();
  mpLine->setText( name );
}


QString
QCSParameterFileNameEditor::text()
{
  QString ret = mpLine->text();
  mName = ret.toStdString();
  return ret;
}


void
QCSParameterFileNameEditor::validateDirName( const QString & fileNameString )
{
  QPalette p = mpLine->palette();
  QFileInfo fileInfo( fileNameString );

  if ( fileInfo.exists() )
  {
    if ( fileInfo.isDir() )
    {
      p.setColor( QPalette::Base, Qt::white );
      mpLine->setPalette( p );
      return;
    }
  }

  p.setColor( QPalette::Base, QColor( 128, 0, 0 ) );
  mpLine->setPalette( p );
}


void
QCSParameterFileNameEditor::validateFileName( const QString & fileNameString )
{
  QPalette p = mpLine->palette();
  QFileInfo fileInfo( fileNameString );

  if ( fileInfo.exists() )
  {
    if ( fileInfo.isFile() )
    {
      p.setColor( QPalette::Base, Qt::white );
      mpLine->setPalette( p );
      return;
    }
  }

  p.setColor( QPalette::Base, QColor( 255, 64, 64 ) );
  mpLine->setPalette( p );
}
