///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSParameterFileNameEditor.h                                         //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-08-06 18:20:41                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef Q_CS_PARAMETER_FILENAME_EDITOR_H
#define Q_CS_PARAMETER_FILENAME_EDITOR_H

#include <QWidget>

#include <string>

#include "../tools/parameters/CSParameter.h"


class QLineEdit;

class QCSParameterFileNameEditor : public QWidget
{
  Q_OBJECT

 public:
  QCSParameterFileNameEditor(CSParameter *parm, QWidget * parent=0);

  QString text();
  void setText(std::string &);
  void setText(QString);

 public slots:
   void showFileDialog();
   void showDirDialog();

   void validateDirName(const QString & fileNameString);
   void validateFileName(const QString & fileNameString);

 signals:
   void doneEditing();

 private:
   std::string   mName;
   QLineEdit   * mpLine;
   CSParameter * mpParameter;
};

#endif // Q_CS_PARAMETER_FILENAME_EDITOR_H
