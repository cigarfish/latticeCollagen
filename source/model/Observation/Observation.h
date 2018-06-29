#ifndef OBSERVATION_H
#define OBSERVATION_H

#include "../../model/Cell/CellSpherical.h"
#include <vector>
#include <stdio.h>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>


class CSModel;
class ModelCellsSpherical;


class Observation
{
private:

	// Used by ObserveModelDescription and ObserveModelMeasures to distinguish models
	// Model: 0: Monolayer model with spherical cells
	int modelType;

    ModelCellsSpherical * mpModelCellsSpherical;

	long int timeZero;

#pragma region Collection of observation methods

// Prelim: Shouldnt std::vector<CellSpherical *> better be std::vector<Cell *> to make the methods more general?

// Mode:
// 0 ... Count all cells
// 1 ... Count only proliferating cells
int GetNumberOfCells(std::vector<CellSpherical *> cells, int mode);

// 0 ... Use all cells
double GetRadiusOfGyration(std::vector<CellSpherical *> cells, int mode);

// 0 ... Use all cells in 2D
double GetDiameterInMicrometerBasedOnRgyr(std::vector<CellSpherical *> cells, int mode);

void GetCenterOfMass(std::vector<CellSpherical *> cells, double *comX, double *comY, double *comZ);

// Measure:
// 0 ...Absolute force
void GetCellsMeasure(std::vector<CellSpherical *> cells, int measure, double *min, double *mean, double *max);

#pragma endregion

#pragma region POVRAY

public:
void WritePOV(int mode);

#pragma endregion

public:

    Observation( CSModel *observedModel, int type );

// Writes a description file for the observed simulation.
// File ending: *.desc
void ObserveModelDescription();

// Writes the measurements
// File ending: *.txt
void ObserveModelMeasures();

void ObserveCellPopulationSnapshot(std::vector<CellSpherical *> cells);

};

#endif

