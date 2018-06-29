#ifndef MODEL_ELEMENT_TRIANGULATEDCELL_H
#define MODEL_ELEMENT_TRIANGULATEDCELL_H

#include "ModelElement.h"

struct BoundingBox;

class ModelElementTriangulatedCell : public ModelElement
{
public:
	ModelElementTriangulatedCell(double x, double y, double z);

    BoundingBox * boundingBox();
};

#endif
