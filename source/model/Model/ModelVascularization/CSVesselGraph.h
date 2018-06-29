#ifndef CS_VESSELGRAPH_H
#define CS_VESSELGRAPH_H


#include "VesselGraph.hpp"
#include "Tumor.hpp"

#include "../Elements/ModelElement_VesselGraph.h"

#include "VoronoiDiagram.hpp"

class CSVesselGraph : public VesselGraph, public ModelElement_VesselGraph
{
public:
  CSVesselGraph( const char *pFilename, Tumor *&tumor );

};

#endif
