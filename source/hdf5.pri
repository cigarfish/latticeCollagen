#  qmake include file for dependencies on the sundials solver

HDF5_LIBS = -lhdf5_cpp -lhdf5


win32 {
      HDF5_LIBS = -llibhdf5_cpp -llibhdf5 -llibszip -llibzlib
      HDF5_INC_DIR = \$\(HDF5_PREFIX\)\\include \$\(HDF5_PREFIX\)\\include\\cpp
      HDF5_LIB_DIR = \$\(HDF5_PREFIX\)\\lib

} else {

  isEmpty (HDF5_PREFIX) {
    ENV_HDF5_PREFIX = $$(HDF5_PREFIX)

    isEmpty (ENV_HDF5_PREFIX) {
      ENV_HDF5_PREFIX = $$CSDEPENDENCIES
      isEmpty (ENV_HDF5_PREFIX) {
        ENV_HDF5_PREFIX = $$(CSDEPENDENCIES)
        isEmpty(ENV_HDF5_PREFIX) {
          ENV_HDF5_PREFIX = /usr/local
          message( Warning:  Neither CSDEPENDENCIES nor HDF5_PREFIX is set! )
          message(  Trying HDF5_PREFIX = "$$ENV_HDF5_PREFIX".)
        }
      }
    }

  } else {
    ENV_HDF5_PREFIX = $$HDF5_PREFIX
  }

  !isEmpty( ENV_HDF5_PREFIX ) {
    exists( $$ENV_HDF5_PREFIX/include/H5Cpp.h ) {
       HDF5_INC_DIR = $$ENV_HDF5_PREFIX/include
       HDF5_LIB_DIR = $$ENV_HDF5_PREFIX/lib
       message( HDF5 found in $$ENV_HDF5_PREFIX )
    } else {
      error( Could not find an installation of HDF5 in $$ENV_HDF5_PREFIX. )
    }
  } else {
    error( HDF5 LIBRARIES:  This should never happen );
  }
}


INCLUDEPATH += $$HDF5_INC_DIR

STATICLIBS += -L$$HDF5_LIB_DIR $$HDF5_LIBS
