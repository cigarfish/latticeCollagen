#include "ModelElement_VesselGraph.h"

ModelElement_VesselGraph::ModelElement_VesselGraph(double x, double y, double z)
    : ModelElement( x, y, z, ModelElement::TypeBarrierTriangle )
{}


BoundingBox *
ModelElement_VesselGraph::boundingBox()
{
  return &mBoundingBox;
}
