TEMPLATE = app

TARGET = TiSimDemo

ROOTDIR = ../../..
SRCDIR  = $$ROOTDIR/source
LIBDIR  = $$ROOTDIR/lib

LIBS += -L$$LIBDIR

contains(QMAKE_CXX, g++) | contains(QMAKE_CXX, clang++) {
  LIBS += -lCSGUI -lCSCore $$LIBDIR/libCSGUI.a
  macx:QMAKE_LFLAGS += -mmacosx-version-min=10.7 -stdlib=libc++
  PRE_TARGETDEPS += $$LIBDIR/libCSGUI.a $$LIBDIR/libCSCore.a
  LIBS += -lz
} else {
  contains( QMAKE_CXX, cl ) {
    LIBS += -L$(IntDir)
    INCLUDEPATH += $$SRCDIR/GeneratedFiles
  }
  LIBS += -lCSCore -lCSGUI
}

include($$SRCDIR/common.pri)


CONFIG += qt

DESTDIR = $$ROOTDIR/bin

DEPENDPATH += . \
              $$SRCDIR \
              $$SRCDIR/gui \
              $$SRCDIR/images \
              $$SRCDIR/tools \
              $$SRCDIR/gui/2DTools \
              $$SRCDIR/gui/GLTools \
              $$SRCDIR/gui/graphTools \
              $$SRCDIR/gui/tabImageProcessing \
              $$SRCDIR/gui/tabToolsWithItkVtk \
              $$SRCDIR/gui/tabMonolayer \
              $$SRCDIR/gui/tabComplexCells \
              $$SRCDIR/gui/tabVascularization \
              $$SRCDIR/images/pipelines \
              $$SRCDIR/images/tools \
              $$SRCDIR/model/BasicDatatypes \
              $$SRCDIR/model/Cell \
              $$SRCDIR/model/CellPopulation \
              $$SRCDIR/model/Model \
              $$SRCDIR/model/Observation \
              $$SRCDIR/tools/triangulation \
              $$SRCDIR/tools/Vascularization \
              $$SRCDIR/tools/new_triangulation \
              $$SRCDIR/images/filters/convertFilters \
              $$SRCDIR/images/filters/graphFilters \
              $$SRCDIR/images/filters/imageFilters

INCLUDEPATH += $$DEPENDPATH

HEADERS += $$SRCDIR/gui/TiSimDemoApp/TiSimDemoMainWindow.h \
           $$SRCDIR/gui/TiSimDemoApp/TiSimDemoCentralWidget.h \
           $$SRCDIR/gui/TiSimDemoApp/monolayerDemo.h

SOURCES += $$SRCDIR/gui/TiSimDemoApp/main.cpp \
           $$SRCDIR/gui/TiSimDemoApp/TiSimDemoMainWindow.cpp \
           $$SRCDIR/gui/TiSimDemoApp/TiSimDemoCentralWidget.cpp

FORMS += $$SRCDIR/gui/tabMonolayer/monolayer.ui
