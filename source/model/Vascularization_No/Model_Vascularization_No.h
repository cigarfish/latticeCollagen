#ifndef MODEL_Vascularization_No_H
#define MODEL_Vascularization_No_H


#include "../Model/CSModel.h"

#include <sstream>

class QCSSimulationThread;
class CSVesselGraph;

class VoronoiDiagram;

#include "../Model/ModelVascularization/VesselGraph.hpp";


class Model_Vascularization_No : public CSModel
{
	friend class Vascularization;
public:
	Model_Vascularization_No();

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
