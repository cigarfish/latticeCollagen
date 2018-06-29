
DEFINES += CS_BUILD_PYTHON

# only supported version so far
PYTHON_VERSION = 2.7


# look for environment variable definition:
isEmpty (PYTHON_PREFIX) {
  ENV_PYTHON_PREFIX = $$(PYTHON_PREFIX)

  # Fall-back on GNU/Linux and MacOSX:
  isEmpty(ENV_PYTHON_PREFIX) {  
    ENV_PYTHON_PREFIX = $$CSDEPENDENCIES
    isEmpty(ENV_PYTHON_PREFIX) {
      ENV_PYTHON_PREFIX = $$(CSDEPENDENCIES)
      isEmpty(ENV_PYTHON_PREFIX) {
        unix:ENV_PYTHON_PREFIX = /usr
        message( Warning:  Neither CSDEPENDENCIES nor PYTHON_PREFIX is set! )
        message( Trying PYTHON_PREFIX = $$ENV_PYTHON_PREFIX .)
     }
    }
  }

} else {
  ENV_PYTHON_PREFIX = $$PYTHON_PREFIX
}

win32 {
  PYTHON_LIBS = -lpython27
  
  PYTHON_INC_DIR = \$\(PYTHON_PREFIX\)\\include
  PYTHON_LIB_DIR = \$\(PYTHON_PREFIX\)\\libs
  
  isEmpty( PYTHON_INC_DIR ) {
    error( Could not find valid Python installation in supplied PYTHON_PREFIX = $$(PYTHON_PREFIX) )
  }

} else {
  PYTHON_LIBS = -lpython2.7

  PYTHON_INC_DIR = $$ENV_PYTHON_PREFIX/include/python2.7
  PYTHON_LIB_DIR = $$ENV_PYTHON_PREFIX/lib64/
  
  isEmpty( PYTHON_LIB_DIR ) {
    PYTHON_LIB_DIR = $$ENV_PYTHON_PREFIX/lib/python2.7
  
    !exists( PYTHON_LIB_DIR ) {
      error( Could not find valid Python installation in supplied PYTHON_PREFIX = $$(PYTHON_PREFIX) )
    }
  }
  message( Detected Python version 2.7 in $$ENV_PYTHON_PREFIX )
}

INCLUDEPATH += $$PYTHON_INC_DIR

LIBS += -L$$PYTHON_LIB_DIR $$PYTHON_LIBS
linux-g++:LIBS += -Wl,-rpath=$$PYTHON_LIB_DIR

