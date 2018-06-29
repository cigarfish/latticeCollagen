#include <cstdio>
#include <fstream>

#include "CSModel.h"

#include "Models.h"

#include <QtCore>

#include "../../../gui/CSGLArena.h"
#include "../../../gui/QCSSimulationThread.h"


CSModel::CSModel(int dimension_)
    : timeStep(1.),
      mpBBList(NULL)
{
  dimension = dimension_;
  mpArena = new CSGLArena();
  mpSimulationThread = new QCSSimulationThread();
}

CSModel::~CSModel()
{
  delete mpArena;
}


CSModel *
CSModel::ModelFromXML( QXmlStreamReader *xml,
                     std::stringstream & errors,
                     std::stringstream & warnings )
{
    Q_ASSERT( xml->isStartElement() && xml->name() == "Model" );

    CSModel * model = NULL;

    std::string modelName =
        xml->attributes().value("name").toString().toStdString();

    std::string modelType = xml->attributes().value("type").toString().toStdString();

    unsigned int modelDims = xml->attributes().value("dimensions").toString().toUInt();

    ////////////////////////////////////////////////////////////////////////////
    // Every Model ought to have two lines in this if-else-if cascade:
    ////////////////////////////////////////////////////////////////////////////
    if ( modelType == "ModelCellsSpherical" )
        model = ModelCellsSpherical::createFromXML( xml, errors, warnings );
    else if ( modelType == "Model3D" )
        model = Model3D::createFromXML( xml, errors, warnings );
    else
        errors << "Model::ModelFromXML:  Unknown model type \""
                  << modelType << "\"" << std::endl;
    ////////////////////////////////////////////////////////////////////////////

    if (model)
    {
        model->SetName(modelName);
        model->dimension = modelDims;
        model->Reset(false);
    }

    return model;
}


void CSModel::Reset(bool) {}

void CSModel::Simulate() {}

void CSModel::SimulateInThread() {}

void CSModel::SetName(string newName)
{
  name = newName;

//	core->tools->output->logfile << " In SetName: " << name  << "\n";; 
}

int CSModel::Run( bool blocking )
{
    if ( blocking )
    {
        while (enableSimulation)
            SimulateInThread();
        return 0;
    }

    mpSimulationThread->setModel( this );

    mpSimulationThread->start();

    return 0;
}

/*
void CSModel::registerCellPopulation(const std::string name)
{
	auto it(mCellPopulations.find(name));

	if (it != mCellPopulations.end())
	{
		//console->warn("Cell Population {0:s} already registered!", name);
		cout << "	-> Cell Population " << name << " already registered!" << endl;
		return;
	}

	//console->info("Registering cell population {0:s}.", name);
	cout << "	-> Registering cell population " << name << "." << endl;
	mCellPopulations[name] = new std::vector<ModelElement *>;
}

void CSModel::removeCellPopulation(const std::string name)
{
	auto it(mCellPopulations.find(name));

	// Pointers are not owned by this vector. So don't call delete!
	if (it != mCellPopulations.end())
		it->second->clear();
}

void CSModel::registerCell(ModelElement * element, const std::string populationName)
{
	//b_cells->add(element);
	mpBBList->add(element);

	auto it(mCellPopulations.find(populationName));

	if (it == mCellPopulations.end())
		registerCellPopulation(populationName);

	mCellPopulations[populationName]->push_back(element);
}

void CSModel::removeCell(ModelElement * element, const std::string hint)
{
	//b_cells->remove(element);
	mpBBList->remove(element);

	// see if it is registered anywhere. Use the hint
	auto it(mCellPopulations.find(hint));
	if (it != mCellPopulations.end())
	{
		auto cellPop = it->second;
		cellPop->erase(std::remove(cellPop->begin(), cellPop->end(), element), cellPop->end());
	}
	else // if hint unsuccessful check em all
	{
		for (auto& kv : mCellPopulations)
		{
			auto cellPop = kv.second;
			cellPop->erase(std::remove(cellPop->begin(), cellPop->end(), element), cellPop->end());
		}
	}
}

std::vector<ModelElement*> *
CSModel::getCellPopulation(const std::string name)
{
	auto it(mCellPopulations.find(name));

	if (it != mCellPopulations.end())
		return it->second;

	//console->error("Cell Population {0:s} not found!", name);
	cout << "	-> Cell Population " << name << " not found!" << endl;
	return nullptr;
}

std::size_t CSModel::getCellPopulationSize(const std::string name)
{
	auto population = getCellPopulation(name);
	if (population)
		return population->size();
	else
		return 0;
}

void CSModel::updateCellPopulations()
{
	assert(false);
}
*/
