/****************************************************************************
** Meta object code from reading C++ file 'tabComplexCells.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../gui/tabComplexCells/tabComplexCells.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'tabComplexCells.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_complexCells_t {
    QByteArrayData data[14];
    char stringdata0[260];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_complexCells_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_complexCells_t qt_meta_stringdata_complexCells = {
    {
QT_MOC_LITERAL(0, 0, 12), // "complexCells"
QT_MOC_LITERAL(1, 13, 11), // "abortThread"
QT_MOC_LITERAL(2, 25, 0), // ""
QT_MOC_LITERAL(3, 26, 17), // "pauseResumeThread"
QT_MOC_LITERAL(4, 44, 25), // "ResetParametersToDefaults"
QT_MOC_LITERAL(5, 70, 21), // "PushParametersToModel"
QT_MOC_LITERAL(6, 92, 17), // "ParametersChanged"
QT_MOC_LITERAL(7, 110, 22), // "startSimulationClicked"
QT_MOC_LITERAL(8, 133, 21), // "pauseResumeSimulation"
QT_MOC_LITERAL(9, 155, 15), // "resetSimulation"
QT_MOC_LITERAL(10, 171, 34), // "buttonAbortSimulationButtonCl..."
QT_MOC_LITERAL(11, 206, 12), // "threadUpdate"
QT_MOC_LITERAL(12, 219, 23), // "doAfterThreadIsFinished"
QT_MOC_LITERAL(13, 243, 16) // "createNewDisplay"

    },
    "complexCells\0abortThread\0\0pauseResumeThread\0"
    "ResetParametersToDefaults\0"
    "PushParametersToModel\0ParametersChanged\0"
    "startSimulationClicked\0pauseResumeSimulation\0"
    "resetSimulation\0buttonAbortSimulationButtonClicked\0"
    "threadUpdate\0doAfterThreadIsFinished\0"
    "createNewDisplay"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_complexCells[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   74,    2, 0x06 /* Public */,
       3,    1,   75,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       4,    0,   78,    2, 0x08 /* Private */,
       5,    0,   79,    2, 0x08 /* Private */,
       6,    0,   80,    2, 0x08 /* Private */,
       7,    0,   81,    2, 0x08 /* Private */,
       8,    0,   82,    2, 0x08 /* Private */,
       9,    0,   83,    2, 0x08 /* Private */,
      10,    0,   84,    2, 0x08 /* Private */,
      11,    0,   85,    2, 0x08 /* Private */,
      12,    0,   86,    2, 0x08 /* Private */,
      13,    0,   87,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void complexCells::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        complexCells *_t = static_cast<complexCells *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->abortThread(); break;
        case 1: _t->pauseResumeThread((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 2: _t->ResetParametersToDefaults(); break;
        case 3: _t->PushParametersToModel(); break;
        case 4: _t->ParametersChanged(); break;
        case 5: _t->startSimulationClicked(); break;
        case 6: _t->pauseResumeSimulation(); break;
        case 7: _t->resetSimulation(); break;
        case 8: _t->buttonAbortSimulationButtonClicked(); break;
        case 9: _t->threadUpdate(); break;
        case 10: _t->doAfterThreadIsFinished(); break;
        case 11: _t->createNewDisplay(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (complexCells::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&complexCells::abortThread)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (complexCells::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&complexCells::pauseResumeThread)) {
                *result = 1;
                return;
            }
        }
    }
}

const QMetaObject complexCells::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_complexCells.data,
      qt_meta_data_complexCells,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *complexCells::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *complexCells::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_complexCells.stringdata0))
        return static_cast<void*>(const_cast< complexCells*>(this));
    return QWidget::qt_metacast(_clname);
}

int complexCells::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 12)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 12;
    }
    return _id;
}

// SIGNAL 0
void complexCells::abortThread()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}

// SIGNAL 1
void complexCells::pauseResumeThread(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
