
SOSLIB_LIBS = -lODES -lsbml -lxml2


win32 {
      SOSLIB_INC_DIR = \$\(SOSLIB_PREFIX\)\\include
      SOSLIB_LIB_DIR = \$\(SOSLIB_PREFIX\)\\lib
      INCLUDEPATH += \$\(LIBSBML_PREFIX\)\\include
      SOSLIB_LIBS = -lODES -llibsbml -lexpat
      LIBS += -L\$\(LIBSBML_PREFIX\)\\lib
} else {

  isEmpty (SOSLIB_PREFIX) {
    ENV_SOSLIB_PREFIX = $$(SOSLIB_PREFIX)

    isEmpty (ENV_SOSLIB_PREFIX) {
      ENV_SOSLIB_PREFIX = $$CSDEPENDENCIES
      isEmpty (ENV_SOSLIB_PREFIX) {
        ENV_SOSLIB_PREFIX = $$(CSDEPENDENCIES)
        isEmpty(ENV_SOSLIB_PREFIX) { 
          ENV_SOSLIB_PREFIX = /usr/local
          message( Warning:  Neither CSDEPENDENCIES nor SOSLIB_PREFIX is set! )
          message(  Trying SOSLIB_PREFIX = "$$ENV_SOSLIB_PREFIX".)
        }
      }
    }

  } else {
    ENV_SOSLIB_PREFIX = $$SOSLIB_PREFIX
  }

  !isEmpty( ENV_SOSLIB_PREFIX ) {
    exists( $$ENV_SOSLIB_PREFIX/include/sbmlsolver ) {
       SOSLIB_INC_DIR = $$ENV_SOSLIB_PREFIX/include/
       SOSLIB_LIB_DIR = $$ENV_SOSLIB_PREFIX/lib
       message( SBML ODE solver library found in $$ENV_SOSLIB_PREFIX )
    } else {
      error( Could not find an installation of SOSLib in $$ENV_SOSLIB_PREFIX. )
    }
  } else {
    error( SOSLIB LIBRARIES:  This should never happen );
  }
}


INCLUDEPATH += $$SOSLIB_INC_DIR

LIBS += -L$$SOSLIB_LIB_DIR $$SOSLIB_LIBS

include(sundials.pri)
