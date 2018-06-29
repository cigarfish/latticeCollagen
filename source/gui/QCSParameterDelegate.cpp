///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  QCSParameterDelegate.cpp                                             //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-05-21 19:56:10                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#include "QCSParameterDelegate.h"

#include "../tools/parameters/CSParameter.h"
#include "../tools/parameters/CSParameterContext.h"

#include <limits>

#if QT_VERSION >= 0x050000
# include <QtWidgets>
#else
# include <QtGui>
#endif

#include "QCSParameterModel.h"
#include "QCSParameterRangeWidget.h"
#include "QCSParameterFileNameEditor.h"


QCSParameterDelegate::QCSParameterDelegate(QObject * parent)
: QStyledItemDelegate(parent)
{}


QWidget *
QCSParameterDelegate::createEditor( QWidget *parent,
                                    const QStyleOptionViewItem &/*option*/,
                                    const QModelIndex &index) const
{
  // this relies on the model not to make ParameterContexts editable
  CSParameter * parm = (CSParameter *) index.internalPointer();

  QWidget * editor;

  switch ( parm->dataType() )
    {
    case CSParameter::Bool:
      editor = NULL;
      // handled by the view's checkbox
      // editor = new QCheckBox(parent);
      // editor->installEventFilter( const_cast<QCSParameterDelegate *>(this) );
      break;

    case CSParameter::Int:
    case CSParameter::Long:
      editor = new QSpinBox(parent);
      ((QSpinBox *)editor)->setMinimum( std::numeric_limits<int>::min() );
      ((QSpinBox *)editor)->setMaximum( std::numeric_limits<int>::max() );
      break;

    case CSParameter::Float:
    case CSParameter::Double:
      {
        editor = new QLineEdit(parent);
        QDoubleValidator * val = new QDoubleValidator();
        ((QLineEdit *)editor)->setValidator( val );
      }
      break;

    case CSParameter::String:
      editor = new QLineEdit(parent);
      break;

    case CSParameter::DirName:
    case CSParameter::FileName:
      {
        editor = new QCSParameterFileNameEditor(parm, parent);
        editor->setFocus();
        connect( editor, SIGNAL(doneEditing()), this, SLOT(closeFileNameEditor()) );
        // QHBoxLayout *layout = new QHBoxLayout();
        // QLineEdit * line = new QLineEdit(editor);
        // layout->addWidget(line);
        // QToolButton * button = new QToolButton(editor);
        // layout->addWidget(button);
        // editor->setLayout(layout);
      }
      break;

    case CSParameter::RangeInt:
      editor = RangeWidget::create(parent,CSParameter::Int);
      break;

    case CSParameter::RangeLong:
      editor = RangeWidget::create(parent,CSParameter::Long);
      break;

    case CSParameter::RangeFloat:
      editor = RangeWidget::create(parent,CSParameter::Float);
      break;

    case CSParameter::RangeDouble:
      editor = RangeWidget::create(parent,CSParameter::Double);
      break;

    case CSParameter::Choice:
      editor = new QComboBox(parent);
      connect(editor, SIGNAL(currentIndexChanged(int)), this, SLOT(commitComboBoxData()) );
      break;
    }

  return editor;
}


QSize
QCSParameterDelegate::sizeHint( const QStyleOptionViewItem & option,
                                const QModelIndex & index ) const
{
    if ( !index.isValid() )
        return QSize();

    QFontMetrics fm( option.font );
    QString indexDataString = index.data( Qt::DisplayRole ).toString();
    int width = fm.width( index.data( Qt::DisplayRole ).toString() );
    int height = fm.height();

    return QSize( width+20, height );
}


bool
QCSParameterDelegate::editorEvent( QEvent *event, QAbstractItemModel * model,
                                   const QStyleOptionViewItem &option,
                                   const QModelIndex &index )
{
  Q_ASSERT(event);
  Q_ASSERT(model);
 
  // make sure that the item is checkable
  Qt::ItemFlags flags = model->flags(index);
  if (!(flags & Qt::ItemIsUserCheckable) || !(flags & Qt::ItemIsEnabled))
    return false;

  // make sure that we have a check state
  QVariant value = index.data(Qt::CheckStateRole);
  if (!value.isValid())
    return false;

  // make sure that we have the right event type
  switch ( event->type() )
    {
    case QEvent::MouseButtonRelease:
      {
        const int textMargin = QApplication::style()->pixelMetric(QStyle::PM_FocusFrameHMargin) + 1;
        QRect checkRect = QStyle::alignedRect( option.direction, option.displayAlignment,
                                               option.decorationSize,
                                               QRect( option.rect.x() + (2 * textMargin), option.rect.y(),
                                                      option.rect.width() - (2 * textMargin),
                                                      option.rect.height()));

        if ( !checkRect.contains(static_cast<QMouseEvent*>(event)->pos()) )
          return false;

      }
      break;

    case QEvent::KeyPress:
      {
        if ( static_cast<QKeyEvent*>(event)->key() != Qt::Key_Space
             && static_cast<QKeyEvent*>(event)->key() != Qt::Key_Select
             && static_cast<QKeyEvent*>(event)->key() != Qt::Key_Enter
             && static_cast<QKeyEvent*>(event)->key() != Qt::Key_Return )
          return false;
      }
      break;

    default:
      return false;
    }

  Qt::CheckState state = ( static_cast<Qt::CheckState>(value.toInt()) ==
                           Qt::Checked ? Qt::Unchecked : Qt::Checked );
  QVariant vstate = QVariant(state);
  return static_cast<QCSParameterModel *>(model)->setData(index, vstate, (int)Qt::CheckStateRole);
}


void
QCSParameterDelegate::setEditorData( QWidget *editor, const QModelIndex &index )
  const
{
  CSParameter * parm = (CSParameter *) index.internalPointer();
  void * dataPtr = parm->dataPointer();

  switch ( parm->dataType() )
    {
    case CSParameter::Bool:
      // {
      //   bool checked = * (bool *) dataPtr;
      //   static_cast<QCheckBox *>(editor)->setCheckState( checked ? Qt::Checked : Qt::Unchecked );
      // }
      break;

    case CSParameter::Int:
      {
        long data = * (int *) dataPtr;
        ((QSpinBox *)editor)->setValue( (int)data );
      }
      break;

    case CSParameter::Long:
      {
        long data = * (long *) dataPtr;
        ((QSpinBox *)editor)->setValue( data );
      }
      break;

    case CSParameter::Float:
      {
        double data = * (float *) dataPtr;
        ((QLineEdit *)editor)->setText( QString("%1").arg(data, 0, 'g', 5 ) );
      }
      break;

    case CSParameter::Double:
      {
        double data = * (double *) dataPtr;
        ((QLineEdit *)editor)->setText( QString("%1").arg(data, 0, 'g', 5 ) );
      }
      break;

    case CSParameter::String:
      {
        std::string data = * (std::string *) dataPtr;
        ((QLineEdit *)editor)->setText( tr(data.c_str()) );
      }
      break;

    case CSParameter::DirName:
    case CSParameter::FileName:
      {
        std::string data = * (std::string *) dataPtr;
        ((QCSParameterFileNameEditor *) editor)->setText( data );
      }
      break;

    case CSParameter::RangeInt:
      {
        std::pair<int,int> data = * (std::pair<int,int> *) dataPtr;
        RangeWidget * ed = (RangeWidget *) editor;
        ed->setStart( data.first );
        ed->setEnd( data.second );
//        std::cout << "setting range\n";
      }
      break;

    case CSParameter::RangeLong:
      {
        std::pair<long,long> data = * (std::pair<long,long> *) dataPtr;
        RangeWidget * ed = (RangeWidget *) editor;
        ed->setStart( data.first );
        ed->setEnd( data.second );
      }
      break;

    case CSParameter::RangeFloat:
      {
        std::pair<float,float> data = * (std::pair<float,float> *) dataPtr;
        RangeWidget * ed = (RangeWidget *) editor;
        ed->setStart( data.first );
        ed->setEnd( data.second );
      }
      break;

    case CSParameter::RangeDouble:
      {
        std::pair<double,double> data = * (std::pair<double,double> *) dataPtr;
        RangeWidget * ed = (RangeWidget *) editor;
        ed->setStart( data.first );
        ed->setEnd( data.second );
      }
      break;

    case CSParameter::Choice:
      {
        QStringList items;
        CSParameterChoice *ptr = static_cast<CSParameterChoice *>(dataPtr);
        std::vector<std::string> choices = ptr->choices();
        std::vector<std::string>::const_iterator it;
        for ( it=choices.begin(); it!=choices.end(); it++ )
          items += QString( it->c_str() );

        QComboBox * ed =
            static_cast<QComboBox *>(editor);
        ed->clear();
        ed->addItems( items );
        ed->setCurrentIndex( ptr->currentIndex() );
      }
      break;
    }
  
}


void
QCSParameterDelegate::setModelData( QWidget *editor, QAbstractItemModel *model,
                                    const QModelIndex &index )
  const
{
  CSParameter * dataPtr = (CSParameter *) index.internalPointer();

  QVariant data;

  switch ( dataPtr->dataType() )
    {
    case CSParameter::Bool:
      data = QVariant();
      break;

    case CSParameter::Int:
      data = QVariant( (int)((QSpinBox *) editor)->value() );
      break;

    case CSParameter::Long:
      data = QVariant( (int) ((QSpinBox *) editor)->value() );
      break;

    case CSParameter::Float:
      data = QVariant( ((QLineEdit *) editor)->text().toFloat() );
      break;

    case CSParameter::Double:
      data = QVariant( ((QLineEdit *) editor)->text().toDouble() );
      break;

    case CSParameter::String:
      data = QVariant( ((QLineEdit *)editor)->text() );
      break;

    case CSParameter::DirName:
    case CSParameter::FileName:
      data = QVariant( ((QCSParameterFileNameEditor *)editor)->text() );
      break;

    case CSParameter::RangeInt:
    case CSParameter::RangeLong:
    case CSParameter::RangeFloat:
    case CSParameter::RangeDouble:
      data = QVariant();
      break;

    case CSParameter::Choice:
      data = QVariant( static_cast<QComboBox *>(editor)->currentIndex() );
      break;

    default:
      data = QVariant();
    }

  static_cast<QCSParameterModel *>(model)->setData( index, data, (int)Qt::EditRole );
}


void
QCSParameterDelegate::updateEditorGeometry( QWidget *editor,
                                            const QStyleOptionViewItem &option,
                                            const QModelIndex &index )
  const
{
  if (index.isValid())
    editor->setGeometry(option.rect);
  else
    QStyledItemDelegate::updateEditorGeometry(editor, option, index);
}


void
QCSParameterDelegate::paint( QPainter *painter,
                             const QStyleOptionViewItem &opt,
                             const QModelIndex & index ) const
{
  QStyleOptionViewItem revisedOpt(opt);

  CSParameterTreeItem * ptr =
    static_cast<CSParameterTreeItem *>(index.internalPointer());

  Qt::Alignment alignment = Qt::AlignLeft;

  if ( ptr->type() == CSParameterTreeItem::Parameter )
    {
      CSParameter * parm = static_cast<CSParameter *>(ptr);

      if ( index.column() == 1 )
        {
          alignment = Qt::AlignRight;

          if ( parm->dataType() == CSParameter::FileName )
          {
            alignment = Qt::AlignLeft;
            QFileInfo finfo( QString( parm->dataString().c_str() ) );
            if ( finfo.exists() && finfo.isFile() )
              painter->fillRect( opt.rect, Qt::white );
            else
              painter->fillRect( opt.rect, QColor( 255, 64, 64) );
          }
          else if ( parm->dataType() == CSParameter::DirName )
          {
            alignment = Qt::AlignLeft;
            QFileInfo finfo( QString( parm->dataString().c_str() ) );
            if ( finfo.exists() && finfo.isDir() )
              painter->fillRect( opt.rect, Qt::white );
            else
              painter->fillRect( opt.rect, QColor( 255, 64, 64 ) );
          }
        }
    }

  revisedOpt.displayAlignment = alignment;

  QStyledItemDelegate::paint( painter, revisedOpt, index );
}


void
QCSParameterDelegate::closeFileNameEditor()
{
    QCSParameterFileNameEditor * editor = qobject_cast<QCSParameterFileNameEditor *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}


void
QCSParameterDelegate::commitComboBoxData()
{
    QComboBox * editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
}
