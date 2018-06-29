
ITK_LIBS = -litksys -lITKCommon -litktiff -litkvnl -litkvnl_algo -litkvcl -lITKDICOMParser -lITKniftiio -litkzlib -litkopenjpeg -litkpng -lITKNrrdIO -lITKIONRRD -litkv3p_netlib -litkNetlibSlatec -lITKStatistics -litkv3p_lsqr -lITKMetaIO -lITKznz -lITKEXPAT -lITKWatersheds -lITKSpatialObjects -lITKPath -lITKMesh -lITKVNLInstantiation -litkgdcmjpeg8 -litkgdcmjpeg12 -litkgdcmjpeg16 -litkjpeg -lITKIOImageBase -lITKIOTIFF -lITKIOLSM -lITKVTK -lITKLabelMap

win32:ITK_SYSLIBS = -lwsock32 -lRpcrt4 -lSnmpapi



# look for environment variable definition:
isEmpty (ITK_PREFIX) {
  ENV_ITK_PREFIX = $$(ITK_PREFIX)

  # Fall-back on GNU/Linux and MacOSX:
  isEmpty(ENV_ITK_PREFIX) {
    ENV_ITK_PREFIX = $$CSDEPENDENCIES
    isEmpty(ENV_ITK_PREFIX) {
      ENV_ITK_PREFIX = $$(CSDEPENDENCIES)
      isEmpty(ENV_ITK_PREFIX) {
        unix:ENV_ITK_PREFIX = /usr/local
        message( Warning:  Neither CSDEPENDENCIES nor ITK_PREFIX is set! )
        message( Trying ITK_PREFIX = $$ENV_ITK_PREFIX .)
      }
    }
  }

} else {
  ENV_ITK_PREFIX = $$ITK_PREFIX
}


win32 {
  ITK_INC_DIR = $$files( $$(ITK_PREFIX)/include/ITK-4.* )
  isEmpty( ITK_INC_DIR ) {
    error( Could not find valid ITK-4 installation in supplied ITK_PREFIX = $$(ITK_PREFIX) )
  }
  ITK_VERSION = $$split( ITK_INC_DIR, "-" )
  ITK_VERSION = $$last( ITK_VERSION )
  ITK_LIB_DIR = \$\(ITK_PREFIX\)\\lib
  ITK_INC_DIR = \$\(ITK_PREFIX\)\\include\\ITK-$$ITK_VERSION

} else {

  !isEmpty( ENV_ITK_PREFIX ) {
    exists( $$ENV_ITK_PREFIX/include/ITK-4.* ) {
      ITK_VERSION = $$files( $$ENV_ITK_PREFIX/include/ITK-4.* )
      ITK_VERSION = $$last( ITK_VERSION )
      ITK_VERSION = $$split( ITK_VERSION, "-" )
      ITK_VERSION = $$last( ITK_VERSION )
      
      message( Detected ITK version $$ITK_VERSION in $$ENV_ITK_PREFIX )

      ITK_INC_DIR = $$ENV_ITK_PREFIX/include/ITK-$$ITK_VERSION
      ITK_LIB_DIR = $$ENV_ITK_PREFIX/lib

    } else {
      error( could not find ITK-4 installation in $$ENV_ITK_PREFIX )
    }
  } else {
    error("Unknown ITK version!")
  }
}



#!exists( $$ITK_INC_DIR ) { error( Cannot find ITK include folder $$ITK_INC_DIR - please specify an ITK_PREFIX ) }
#!exists( $$ITK_LIB_DIR ) { error( Cannot find ITK library folder $$ITK_INC_DIR - please specify an ITK_PREFIX ) }


message("Libraries have postfix versioning... appending.")
DUMMYVAR =
for(libflag,ITK_LIBS) {
  DUMMYVAR += $${libflag}-$$ITK_VERSION
}
ITK_LIBS = $$DUMMYVAR

INCLUDEPATH += $$ITK_INC_DIR

LIBS += -L$$ITK_LIB_DIR $$ITK_LIBS $$ITK_SYSLIBS
