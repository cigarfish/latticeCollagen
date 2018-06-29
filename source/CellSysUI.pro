######################################################################
# Automatically generated by qmake (2.01a) Fri Jan 13 13:02:39 2012
######################################################################

TEMPLATE = app

TARGET = CellSysUI

#!isEmpty( TI_QUANT_ONLY ) {
#  TARGET = TiQuant
#}


build_info.target = commit.h
build_info.depends = FORCE
build_info.commands = sh ./prepare_build_info.sh

!win32:QMAKE_EXTRA_TARGETS += build_info
!win32:PRE_TARGETDEPS += commit.h
!win32:INCLUDES += commit.h tools/dataIO/buildInformation.h


contains(QMAKE_CXX, g++) | contains(QMAKE_CXX, clang++) | macx-xcode {
  LIBS += -lCSGUI -lCSCore ../lib/libCSGUI.a
  macx:QMAKE_LFLAGS += -mmacosx-version-min=10.7 -stdlib=libc++
  PRE_TARGETDEPS += ../lib/libCSGUI.a ../lib/libCSCore.a
  LIBS += -lz
} else {
  contains( QMAKE_CXX, cl ) {
    LIBS += -L$(IntDir)
  }
  LIBS += -lCSCore -lCSGUI
}

macx:PRE_TARGETDEPS += ../lib/libCSGUI.a ../lib/libCSCore.a


include(common.pri)

contains(QMAKE_CXX,g++) {
  !macx:LIBS += -Wl,-Bstatic $$STATICLIBS -Wl,-Bdynamic -ldl
  macx:LIBS  += $$STATICLIBS
} else {
  LIBS += $$STATICLIBS
}

CONFIG += qt

DESTDIR = ../bin

DEPENDPATH += . \
              gui \
              images \
              tools \
              gui/2DTools \
              gui/GLTools \
              gui/graphTools \
              gui/tabImageProcessing \
              gui/tabToolsWithItkVtk \
              gui/tabMonolayer \
              gui/tabComplexCells \
              gui/tabVascularization \
              images/pipelines \
              images/tools \
              model/BasicDatatypes \
              model/Cell \
              model/CellPopulation \
              model/Model \
              model/Observation \
              tools/triangulation \
              tools/Vascularization \
              tools/new_triangulation \
              images/filters/convertFilters \
              images/filters/graphFilters \
              images/filters/imageFilters
INCLUDEPATH += . \
               gui \
               tools \
               model/Model \
               model/CellPopulation \
               model/Cell \
               model/BasicDatatypes \
               gui/2DTools \
               images \
               images/filters/imageFilters \
               gui/tabMonolayer \
               gui/tabImageProcessing \
               gui/tabToolsWithItkVtk \
               gui/tabVascularization \
               gui/tabComplexCells \
               gui/GLTools \
               gui/graphTools \
               images/pipelines \
               images/filters/convertFilters \
               images/filters/graphFilters \
               images/tools \
               model/Observation \
               tools/triangulation \
               tools/new_triangulation \
               tools/Vascularization

# Input
HEADERS += 

SOURCES += gui/main.cpp

