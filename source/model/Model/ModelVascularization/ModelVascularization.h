#ifndef MODEL_Vascularization_H
#define MODEL_Vascularization_H


#include "../Model/CSModel.h"

#include <sstream>

#include "VesselGraph.hpp"


class QCSSimulationThread;
class CSVesselGraph;

class VoronoiDiagram;




class Model_Vascularization : public CSModel
{
	friend class Vascularization;
public:
	Model_Vascularization();

	void Reset();
//	void Simulate();
	void writeXML( QXmlStreamWriter * ) const {};
	void SimulateInThread();

    void writeHDF5( H5::H5File * /*outputFile*/ ) const {};

    void readModelData( H5::H5File * /* inputFile */,
                        std::stringstream & /* errors */,
                        std::stringstream & /* warnings */ ) {};

    void SetupSimulation();

public:
  CSVesselGraph *vg;
  VoronoiDiagram *vornoi;
  CONCENTRATION_T ***concentration;
  int ***labelCells;

};
#endif
