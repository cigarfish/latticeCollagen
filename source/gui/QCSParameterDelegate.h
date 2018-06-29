///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSParameterDelegate.h                                               //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-05-21 19:50:50                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef Q_CS_PARAMETER_DELEGATE_H
#define Q_CS_PARAMETER_DELEGATE_H

#include <QStyledItemDelegate>


class QCSParameterDelegate : public QStyledItemDelegate
{
    Q_OBJECT

 public:
  QCSParameterDelegate(QObject * parent = 0);

  QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                           const QModelIndex &index) const;

  bool editorEvent( QEvent *event, QAbstractItemModel * model,
                    const QStyleOptionViewItem &option,
                    const QModelIndex &index );

  void setEditorData( QWidget *editor, const QModelIndex &index ) const;

  void setModelData( QWidget *editor, QAbstractItemModel *model,
                     const QModelIndex &index ) const;

  void updateEditorGeometry( QWidget *editor,
                             const QStyleOptionViewItem &option,
                             const QModelIndex &index ) const;

  void paint( QPainter * painter, 
              const QStyleOptionViewItem & opt,
              const QModelIndex & index ) const;

  QSize sizeHint( const QStyleOptionViewItem & option,
                  const QModelIndex & index ) const;

 public slots:
  void closeFileNameEditor();
  void commitComboBoxData();
};

#endif // Q_CS_PARAMETER_DELEGATE_H
