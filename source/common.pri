#######################################################################################
##                                                                                   ##
##  File Name:  common.pri                                                           ##
##                                                                                   ##
##     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               ##
##    Created:  2011-11-26 00:23:13                                                  ##
##                                                                                   ##
##  This file is part of the CellSys7 code.                                          ##
##                                                                                   ##
##  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. ##
##                                                                                   ##
#######################################################################################

# the install prefix will be:
INSTALLDIR = ..


# CellSys libraries:
LIBS += -L$$INSTALLDIR/lib


INCLUDEPATH += .


linux-*      { DEFINES += __BUILD_LINUX__ }
macx         { DEFINES += __BUILD_MAC__ }
macx-xcode   { DEFINES += __BUILD_XCODE__ }
win32        { DEFINES += __BUILD_WINDOWS__ }

# DEFINES += _GLIBCXX_USE_CXX11_ABI=0

QT += concurrent

unix {
  CONFIG += debug_and_release
  QMAKE_CXXFLAGS_RELEASE = -O3 
  QMAKE_CXXFLAGS_WARN_ON += -Wno-unknown-pragmas

  mac | macx {
    QMAKE_CXXFLAGS += -mmacosx-version-min=10.7 -stdlib=libc++ -std=c++14
    LIBS += -framework Cocoa -framework IOKit
  } else {
    QMAKE_LFLAGS_RELEASE   = -Wl,-O3
    QMAKE_CXXFLAGS += -std=c++14
  }
}

include (deployment.pri)

# include(sundials.pri)
include(hdf5.pri)

include(opengl.pri)

include(soslib.pri)

include(multicellds.pri)

# Creating a TiQuant-only binary
isEmpty(TI_QUANT_ONLY) {
  TI_QUANT_ONLY = $$(TI_QUANT_ONLY)
}

!isEmpty(TI_QUANT_ONLY) {
  contains(TI_QUANT_ONLY, "yes") || contains(TI_QUANT_ONLY, "y") || contains(TI_QUANT_ONLY, "1") {
    DEFINES += CS_TI_QUANT_ONLY
    USE_IMAGEPROCESSING = 1
  }
}

# Compiling and linking image processing;  default is not to.
isEmpty(USE_IMAGEPROCESSING) {
  USE_IMAGEPROCESSING = $$(USE_IMAGEPROCESSING)
}

!isEmpty(USE_IMAGEPROCESSING) {
  contains(USE_IMAGEPROCESSING, "yes") || contains(USE_IMAGEPROCESSING, "y") || contains(USE_IMAGEPROCESSING, "1") {
    include(imageProcessing.pri)
  }
}

# Include Python API.
isEmpty(USE_PYTHON) {
  USE_PYTHON = $$(USE_PYTHON)
}

!isEmpty(USE_PYTHON) {
  contains(USE_PYTHON, "yes") || contains(USE_PYTHON, "y") || contains(USE_PYTHON, "1") {
    include(python.pri)
  }
}
