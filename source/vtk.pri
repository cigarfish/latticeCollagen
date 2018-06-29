
VTK_LIBS = -lQVTK -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkIO -lvtkFiltering -lvtkCommon -lvtksys -lvtkpng -lvtktiff -lvtkjpeg -lvtkexpat -lvtkzlib -lvtkfreetype -lvtkftgl -lvtkViews -lvtkInfovis -lvtkWidgets -lvtkHybrid -lvtkVolumeRendering -lvtkverdict -lvtkDICOMParser -lvtkNetCDF -lvtkNetCDF_cxx -lvtkmetaio -lvtksqlite -lvtklibxml2 -lvtkalglib -lvtkexoIIc


# look for environment variable definition:
isEmpty (VTK_PREFIX) {
  ENV_VTK_PREFIX = $$(VTK_PREFIX)

  # Fall-back on GNU/Linux and MacOSX:
  isEmpty(ENV_VTK_PREFIX) {
    ENV_VTK_PREFIX = $$CSDEPENDENCIES
    isEmpty (ENV_VTK_PREFIX) {
      ENV_VTK_PREFIX = $$(CSDEPENDENCIES)
      isEmpty(ENV_VTK_PREFIX) {
        unix:ENV_VTK_PREFIX = /usr/local
        message( Warning:  Neither CSDEPENDENCIES nor VTK_PREFIX is set! )
        message( Trying VTK_PREFIX = $$ENV_VTK_PREFIX .)
      }
    }
  }
  
} else {
  ENV_VTK_PREFIX = $$VTK_PREFIX
}

# default version
VTK_VERSION = 5.8

win32 {

  INC_DIR = $$files( $$(VTK_PREFIX)/include/vtk-5.* )
  isEmpty( INC_DIR ) {
    error( Could not find valid VTK installation in supplied VTK_PREFIX = $$(VTK_PREFIX) )
  }
  VTK_VERSION = $$split(INC_DIR,"-")
  VTK_VERSION = $$last(VTK_VERSION)
  VTK_INC_DIR = \$\(VTK_PREFIX\)\\include\\vtk-$$VTK_VERSION
  VTK_LIB_DIR = \$\(VTK_PREFIX\)\\lib\\vtk-$$VTK_VERSION

} else {

  !isEmpty( ENV_VTK_PREFIX ) {
    exists( $$ENV_VTK_PREFIX/include/vtk-5.* ) {
      VTK_VERSION = $$files( $$ENV_VTK_PREFIX/include/vtk-5.* )
      VTK_VERSION = $$last( VTK_VERSION )
      VTK_VERSION = $$split( VTK_VERSION, "-" )
      VTK_VERSION = $$last( VTK_VERSION )
      message( Using VTK version $$VTK_VERSION in $$ENV_VTK_PREFIX )
    } else {
      error( Could not find valid VTK installation in supplied VTK_PREFIX = $$(VTK_PREFIX) )
    }

    !exists( $$ENV_VTK_PREFIX/lib/vtk-$$VTK_VERSION ) {
      error( supplied VTK_LIBRARY_DIR does not exist: $$(VTK_LIBRARY_DIR) )
    }

    VTK_INC_DIR = $$ENV_VTK_PREFIX/include/vtk-$$VTK_VERSION
    VTK_LIB_DIR = $$ENV_VTK_PREFIX/lib/vtk-$$VTK_VERSION
  } else {
    !exists( $$VTK_INC_DIR ) { error( Cannot find VTK include folder $$VTK_INC_DIR - please specify an VTK_PREFIX ) }
    !exists( $$VTK_LIB_DIR ) { error( Cannot find VTK library folder $$VTK_INC_DIR - please specify an VTK_PREFIX ) }
  }
}

INCLUDEPATH += $$VTK_INC_DIR

LIBS += -L$$VTK_LIB_DIR $$VTK_LIBS
linux-g++:LIBS += -Wl,-rpath=$$VTK_LIB_DIR

