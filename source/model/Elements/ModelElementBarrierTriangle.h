#ifndef MODEL_ELEMENT_BARRIER_TRIANGLE_H
#define MODEL_ELEMENT_BARRIER_TRIANGLE_H

#include "ModelElement.h"

#include <sstream>

namespace H5 { class CompType; };

struct BoundingBox;

class ModelElementBarrierTriangle : public ModelElement
{
public:
  
  ModelElementBarrierTriangle(double x=0., double y=0., double z=0.);

  ~ModelElementBarrierTriangle();

  void setPoint(double x, double y, double z, int index);
  void setBoundingBox(double epsilon);
  double getD();
  void setNormalVector();
  double getNormalVector( int i);
  double* getNormalVector();

  BoundingBox * boundingBox();

  static void HDF5DataFormat( H5::CompType & );
  static H5::CompType ParseHDF5DataFormat( H5::CompType & inputTypeDefinition,
                                           std::stringstream &,
                                           std::stringstream & warnings );
private:

  double mpPoints[3][3];
  double mpNormalVector[3];

  double mD;//distance (hessian normal form)

};

#endif
