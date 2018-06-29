#ifndef MODEL_H
#define MODEL_H

#include "../../tools/Tools.h"
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include <map>
#include <unordered_map>
#include <sstream>

class CSGLArena;
class QCSSimulationThread;
class QXmlStreamWriter;
class QXmlStreamReader;
class BoundingBoxList;

using namespace std::chrono;

namespace H5 {
    class H5File;
}

//! Generic class from which all models are derived from
class CSModel
{
public:

	#pragma region Name

    // Name of simulation instance, population etc.
    // Used as prefix for simulation output files
    std::string name;

    void SetName(std::string newName);

    #pragma endregion

    #pragma region Time

    // Current timestep (dt) usually in internal units
    double timeStep;

    // Current simulation time
    double time;

    #pragma endregion

    // Spatial dimension
    int dimension;

    // Default constructor
    CSModel(int dimension);

    // Default destructor
    ~CSModel();

    // Reset model
    virtual void Reset( bool addFirstCell =true);

    // Should the simulation be continued (e.g. by thread)
    bool enableSimulation;

    // Exactly one simulation step (t -> t+1)
    virtual void Simulate();

    // Includes a number (or one) simulation steps plus observation etc.
    // Method is called by GUI from parallel thread
    virtual void SimulateInThread();

    // Do anything necessary for initialising a batch job:
    virtual void SetupSimulation() =0;

    // Run model simulation loop.
    // A batch job would want to block until the loop returns.  When called from
    // the GUI, the default can be used, i.e. model->Run();
    int Run( bool blocking =false );


    // routine to write XML code in order to save the Model
    virtual void writeXML( QXmlStreamWriter * ) const =0;

    // routine for writing the extensive data to an HDF5 file
    virtual void writeHDF5( H5::H5File */*outputFile*/ ) const =0;

    // method to fill the model with data from an HDF5 file.
    // the Model should be initialised with parameters, e.g. from XML,
    // the Model HAS TO HAVE ITS name SET.
    virtual void readModelData(H5::H5File * inputFile,
                               std::stringstream & errors,
                               std::stringstream & warnings ) =0;

    // The factory for creating models according to the XML information.
    // Every model to be created from XML has to have an entry in this method.
    static CSModel * ModelFromXML( QXmlStreamReader *,
                                 std::stringstream & errors,
                                 std::stringstream & warnings );

    // give access to the collection of CSGLObjects of the model:
    CSGLArena * arena() const
    {
        return mpArena;
    };

    QCSSimulationThread * SimulationThread() const
    { return mpSimulationThread; };

    Random mRandom;

	// Public interface for cell populations
	//void registerCellPopulation(const std::string name);
	//void removeCellPopulation(const std::string name);
	// register cell in a boundingBoxList
	//void registerCell(ModelElement * element, const std::string populationName);
	//void removeCell(ModelElement * element, const std::string populationName);

	// update Cell Populations
	//void updateCellPopulations();

    BoundingBoxList *mpBBList;

	// Cell populations
	// TODO make this a smart pointer!
	//std::unordered_map< std::string, vector< ModelElement* >* > mCellPopulations;

	//vector< ModelElement * >* getCellPopulation(const std::string name);
	//std::size_t               getCellPopulationSize(const std::string name);



protected:
    CSGLArena * mpArena;
    QCSSimulationThread * mpSimulationThread;
};

#endif
