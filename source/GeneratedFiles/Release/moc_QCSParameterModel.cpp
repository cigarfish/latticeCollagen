/****************************************************************************
** Meta object code from reading C++ file 'QCSParameterModel.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../gui/QCSParameterModel.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'QCSParameterModel.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_QCSParameterModel_t {
    QByteArrayData data[5];
    char stringdata0[72];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_QCSParameterModel_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_QCSParameterModel_t qt_meta_stringdata_QCSParameterModel = {
    {
QT_MOC_LITERAL(0, 0, 17), // "QCSParameterModel"
QT_MOC_LITERAL(1, 18, 13), // "choiceChanged"
QT_MOC_LITERAL(2, 32, 0), // ""
QT_MOC_LITERAL(3, 33, 11), // "dataChanged"
QT_MOC_LITERAL(4, 45, 26) // "changeActivatedSubContexts"

    },
    "QCSParameterModel\0choiceChanged\0\0"
    "dataChanged\0changeActivatedSubContexts"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_QCSParameterModel[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   29,    2, 0x06 /* Public */,
       3,    1,   32,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       4,    1,   35,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void, QMetaType::QModelIndex,    2,
    QMetaType::Void, QMetaType::Int,    2,

 // slots: parameters
    QMetaType::Void, QMetaType::QModelIndex,    2,

       0        // eod
};

void QCSParameterModel::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        QCSParameterModel *_t = static_cast<QCSParameterModel *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->choiceChanged((*reinterpret_cast< const QModelIndex(*)>(_a[1]))); break;
        case 1: _t->dataChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->changeActivatedSubContexts((*reinterpret_cast< const QModelIndex(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (QCSParameterModel::*_t)(const QModelIndex & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&QCSParameterModel::choiceChanged)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (QCSParameterModel::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&QCSParameterModel::dataChanged)) {
                *result = 1;
                return;
            }
        }
    }
}

const QMetaObject QCSParameterModel::staticMetaObject = {
    { &QAbstractItemModel::staticMetaObject, qt_meta_stringdata_QCSParameterModel.data,
      qt_meta_data_QCSParameterModel,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *QCSParameterModel::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *QCSParameterModel::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_QCSParameterModel.stringdata0))
        return static_cast<void*>(const_cast< QCSParameterModel*>(this));
    return QAbstractItemModel::qt_metacast(_clname);
}

int QCSParameterModel::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QAbstractItemModel::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void QCSParameterModel::choiceChanged(const QModelIndex & _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void QCSParameterModel::dataChanged(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
