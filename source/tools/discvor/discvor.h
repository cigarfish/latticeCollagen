#ifndef DISCVOR_H
#define DISCVOR_H

#include <vector>
#include <iostream>
#include <QTextBrowser>
#include <QProgressBar>
#include <list>

class CellSpherical;
class GraphSphere;
class ModelCellsSpherical;

// Implementation of discrete Voronoi limited by a network
// Currently the class is tested by calling Test()
class DiscreteVoronoi
{
public:

	// Init for cells and sinusoids
	DiscreteVoronoi(std::vector<CellSpherical *> * cells, GraphSphere * graph, std::ofstream * console, QProgressBar * progressBar_);

	// Init for cells
	DiscreteVoronoi(std::vector<CellSpherical *> * cells, std::ofstream * console, QProgressBar * progressBar_);

    void Run(float widthOfVoxel_, float cellCutOffRadius_, bool enableDebugOutput_); // For testing functionality

private:

    #pragma region Init

	void InitDefault(std::vector<CellSpherical *> * cells_, std::ofstream * console_, QProgressBar * progressBar_);

    #pragma endregion

    #pragma region Interface to (new) CellSys

	// Link to cells
	std::vector<CellSpherical *> *cells;

	// Link to graph
	GraphSphere *graph;
    
	// For debug and info output (e.g. in gui)
	std::ofstream *console;

	// For progress output in gui
	QProgressBar *progressBar;

	#pragma endregion

    #pragma region Core parameters

    // Defines resolution of discretization
    float widthOfVoxel;

    // Maximal radius a Voronoi-cell for cells can have
    float cellCutOffRadius;

	bool enableDebugOutput;

    #pragma endregion

    #pragma region Interface methods (to be replaced later with CellSys7 methods)

    int GetNumberOfCells();
    // Id's are assumed to run from 0 to GetNumberOfCells()-1
    float GetPositionXOfCell(int id);
    float GetPositionYOfCell(int id);
    float GetPositionZOfCell(int id);
    float GetRadiusOfCell(int id);

    unsigned int IsOverlappingWithNetwork(float x, float y, float z);

    #pragma endregion

    #pragma region Helper methods

    // Basic math
    float MathDot(float x1, float y1, float z1, float x2, float y2, float z2);
    float MathNorm(float vx, float vy, float vz);
    float MathDist(float x1, float y1, float z1, float x2, float y2, float z2);
    float MathDistQ(float x1, float y1, float z1, float x2, float y2, float z2);

    // Get Distance from Point to Line segment
    float GetDistanceFromPointToLineSegment(float px, float py, float pz, float s0x, float s0y, float s0z, float s1x, float s1y, float s1z); // Point P, Segment S)

    #pragma endregion

    // (1) Measure the extension of the volume of interest
    float minX, maxX, minY, maxY, minZ, maxZ;    // Parameters describing the rectangular discrete volume in which all calculations take place
    void InitMinMaxXYZ();

    // (2) Init voxels
    // This includes memory allocation
    int resX, resY, resZ; // Resolution of discrete volume used for further calculations
    int *voxels; // Data array
    void InitVoxels();
	// Voxels: 0 = Empty space or Yet unassigned
	//         1 ... GetNumberOfCells() = Cells
	//         GetNumberOfCells() + 1 ... GetNumberOfCells() + 1 + graph.size() = Edges of graph


    // (3) Determine which voxels belong to the network
	bool useNetwork; // If set false, only cells are analyzed (e.g. tumor spheroid without blood vessels setting)
	void DefineNetwork();

    // (4) Determine which voxels belong to the cells
    int GetClosestCellIndexOrZero(float x, float y, float z);   // Prelim: Significant speedup is possible with cell-position precomputation in int-space (check if all cells have different pos)
    void DefineCells();

    // (5) Refine voxels and compute statistics for each cell
    bool IsVoxelInContactWithNetwork(int x, int y, int z); // Checks whether a voxel is surrounded by ANY network node (6-neigborhood)
	bool IsVoxelInContactWithNetwork(int x, int y, int z, int edge); // Checks whether a voxel is surrounded by A SPECIFIC network node (6-neigborhood)
    bool IsVoxelInContactWithOtherCell(int otherThanIndex, int x, int y, int z); // Checks wheter a voxel is surrounded by ANY 'other' cell node (6-neigborhood)
    bool IsVoxelInner(int x, int y, int z); // Checks wheter a voxel is surrounded ONLY by nodes of the same type (int value) (6-neigborhood)
    void Quantify(const char* fileNamePerCell, const char* fileNamePerCellAndEdge);

    #pragma region Debug

    // Output povray representation of result
    void OutputPovray(const char* filename, int version);

    #pragma endregion
};

#endif
