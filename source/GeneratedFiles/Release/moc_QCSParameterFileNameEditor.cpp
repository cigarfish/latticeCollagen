/****************************************************************************
** Meta object code from reading C++ file 'QCSParameterFileNameEditor.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../gui/QCSParameterFileNameEditor.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'QCSParameterFileNameEditor.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_QCSParameterFileNameEditor_t {
    QByteArrayData data[8];
    char stringdata0[117];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_QCSParameterFileNameEditor_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_QCSParameterFileNameEditor_t qt_meta_stringdata_QCSParameterFileNameEditor = {
    {
QT_MOC_LITERAL(0, 0, 26), // "QCSParameterFileNameEditor"
QT_MOC_LITERAL(1, 27, 11), // "doneEditing"
QT_MOC_LITERAL(2, 39, 0), // ""
QT_MOC_LITERAL(3, 40, 14), // "showFileDialog"
QT_MOC_LITERAL(4, 55, 13), // "showDirDialog"
QT_MOC_LITERAL(5, 69, 15), // "validateDirName"
QT_MOC_LITERAL(6, 85, 14), // "fileNameString"
QT_MOC_LITERAL(7, 100, 16) // "validateFileName"

    },
    "QCSParameterFileNameEditor\0doneEditing\0"
    "\0showFileDialog\0showDirDialog\0"
    "validateDirName\0fileNameString\0"
    "validateFileName"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_QCSParameterFileNameEditor[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   39,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       3,    0,   40,    2, 0x0a /* Public */,
       4,    0,   41,    2, 0x0a /* Public */,
       5,    1,   42,    2, 0x0a /* Public */,
       7,    1,   45,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    6,
    QMetaType::Void, QMetaType::QString,    6,

       0        // eod
};

void QCSParameterFileNameEditor::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        QCSParameterFileNameEditor *_t = static_cast<QCSParameterFileNameEditor *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->doneEditing(); break;
        case 1: _t->showFileDialog(); break;
        case 2: _t->showDirDialog(); break;
        case 3: _t->validateDirName((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: _t->validateFileName((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (QCSParameterFileNameEditor::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&QCSParameterFileNameEditor::doneEditing)) {
                *result = 0;
                return;
            }
        }
    }
}

const QMetaObject QCSParameterFileNameEditor::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_QCSParameterFileNameEditor.data,
      qt_meta_data_QCSParameterFileNameEditor,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *QCSParameterFileNameEditor::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *QCSParameterFileNameEditor::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_QCSParameterFileNameEditor.stringdata0))
        return static_cast<void*>(const_cast< QCSParameterFileNameEditor*>(this));
    return QWidget::qt_metacast(_clname);
}

int QCSParameterFileNameEditor::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 5;
    }
    return _id;
}

// SIGNAL 0
void QCSParameterFileNameEditor::doneEditing()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
