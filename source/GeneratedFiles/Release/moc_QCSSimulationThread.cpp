/****************************************************************************
** Meta object code from reading C++ file 'QCSSimulationThread.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../gui/QCSSimulationThread.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'QCSSimulationThread.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_QCSSimulationThread_t {
    QByteArrayData data[9];
    char stringdata0[80];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_QCSSimulationThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_QCSSimulationThread_t qt_meta_stringdata_QCSSimulationThread = {
    {
QT_MOC_LITERAL(0, 0, 19), // "QCSSimulationThread"
QT_MOC_LITERAL(1, 20, 10), // "updateStep"
QT_MOC_LITERAL(2, 31, 0), // ""
QT_MOC_LITERAL(3, 32, 6), // "paused"
QT_MOC_LITERAL(4, 39, 16), // "finishedNormally"
QT_MOC_LITERAL(5, 56, 4), // "stop"
QT_MOC_LITERAL(6, 61, 5), // "pause"
QT_MOC_LITERAL(7, 67, 6), // "resume"
QT_MOC_LITERAL(8, 74, 5) // "reset"

    },
    "QCSSimulationThread\0updateStep\0\0paused\0"
    "finishedNormally\0stop\0pause\0resume\0"
    "reset"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_QCSSimulationThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   54,    2, 0x06 /* Public */,
       3,    0,   57,    2, 0x06 /* Public */,
       4,    0,   58,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       5,    0,   59,    2, 0x0a /* Public */,
       6,    1,   60,    2, 0x0a /* Public */,
       6,    0,   63,    2, 0x2a /* Public | MethodCloned */,
       7,    0,   64,    2, 0x0a /* Public */,
       8,    0,   65,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    6,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void QCSSimulationThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        QCSSimulationThread *_t = static_cast<QCSSimulationThread *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->updateStep((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->paused(); break;
        case 2: _t->finishedNormally(); break;
        case 3: _t->stop(); break;
        case 4: _t->pause((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 5: _t->pause(); break;
        case 6: _t->resume(); break;
        case 7: _t->reset(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (QCSSimulationThread::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&QCSSimulationThread::updateStep)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (QCSSimulationThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&QCSSimulationThread::paused)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (QCSSimulationThread::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&QCSSimulationThread::finishedNormally)) {
                *result = 2;
                return;
            }
        }
    }
}

const QMetaObject QCSSimulationThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_QCSSimulationThread.data,
      qt_meta_data_QCSSimulationThread,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *QCSSimulationThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *QCSSimulationThread::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_QCSSimulationThread.stringdata0))
        return static_cast<void*>(const_cast< QCSSimulationThread*>(this));
    return QThread::qt_metacast(_clname);
}

int QCSSimulationThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 8)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 8;
    }
    return _id;
}

// SIGNAL 0
void QCSSimulationThread::updateStep(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void QCSSimulationThread::paused()
{
    QMetaObject::activate(this, &staticMetaObject, 1, Q_NULLPTR);
}

// SIGNAL 2
void QCSSimulationThread::finishedNormally()
{
    QMetaObject::activate(this, &staticMetaObject, 2, Q_NULLPTR);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
