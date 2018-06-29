#  qmake include file for dependencies on the sundials solver

SUNDIALS_LIBS = -lsundials_kinsol -lsundials_nvecserial -lsundials_cvodes -lsundials_ida


win32 {

      SUNDIALS_INC_DIR = \$\(SUNDIALS_PREFIX\)\\include
      SUNDIALS_LIB_DIR = \$\(SUNDIALS_PREFIX\)\\lib

} else {

  isEmpty (SUNDIALS_PREFIX) {
    ENV_SUNDIALS_PREFIX = $$(SUNDIALS_PREFIX)

    isEmpty (ENV_SUNDIALS_PREFIX) {
      ENV_SUNDIALS_PREFIX = $$CSDEPENDENCIES
      isEmpty (ENV_SUNDIALS_PREFIX) {
        ENV_SUNDIALS_PREFIX = $$(CSDEPENDENCIES)
        isEmpty(ENV_SUNDIALS_PREFIX) {
          ENV_SUNDIALS_PREFIX = /usr/local
          message( Warning:  Neither CSDEPENDENCIES nor SUNDIALS_PREFIX is set! )
          message(  Trying SUNDIALS_PREFIX = "$$ENV_SUNDIALS_PREFIX".)
        }
      }
    }

  } else {
    ENV_SUNDIALS_PREFIX = $$SUNDIALS_PREFIX
  }

  !isEmpty( ENV_SUNDIALS_PREFIX ) {
    exists( $$ENV_SUNDIALS_PREFIX/include/sundials ) {
       SUNDIALS_INC_DIR = $$ENV_SUNDIALS_PREFIX/include
       SUNDIALS_LIB_DIR = $$ENV_SUNDIALS_PREFIX/lib
       message( Sundials found in $$ENV_SUNDIALS_PREFIX )
    } else {
      error( Could not find an installation of Sundials in $$ENV_SUNDIALS_PREFIX. )
    }
  } else {
    error( SUNDIALS LIBRARIES:  This should never happen );
  }
}


INCLUDEPATH += $$SUNDIALS_INC_DIR

STATICLIBS += -L$$SUNDIALS_LIB_DIR $$SUNDIALS_LIBS