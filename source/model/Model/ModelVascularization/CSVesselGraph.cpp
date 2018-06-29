#include "CSVesselGraph.h"

#include "../../gui/GLTools/CSGLSphere.h"
#include "../../gui/GLTools/CSGLVascularization.h"
CSVesselGraph::CSVesselGraph( const char *pFilename, Tumor *&tumor )
  : VesselGraph( pFilename, tumor),
    ModelElement_VesselGraph( 0,0,0 )
{
  mpGLObject =new  CSGLVascularization( this );// new CSGLSphere( this );
}
