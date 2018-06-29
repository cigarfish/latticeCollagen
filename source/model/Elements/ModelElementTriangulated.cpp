#include "ModelElementTriangulated.h"

ModelElementTriangulatedCell::ModelElementTriangulatedCell(double x, double y, double z)
    : ModelElement( x, y, z, ModelElement::TypeTriangulated )
{
}


BoundingBox *
ModelElementTriangulatedCell::boundingBox()
{
  return &mBoundingBox;
}
