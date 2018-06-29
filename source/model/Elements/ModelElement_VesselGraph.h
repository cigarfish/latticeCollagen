#ifndef MODEL_ELEMENT_VASCULARIZATION_H
#define MODEL_ELEMENT_VASCULARIZATION_H

#include "ModelElement.h"

struct BoundingBox;

class ModelElement_VesselGraph : public ModelElement
{
public:
  
  ModelElement_VesselGraph(double x, double y, double z);
  
  BoundingBox * boundingBox();
};

#endif
