
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "discvor.h"
#include "../../model/Cell/CellSpherical.h"
#include "../../model/BasicDatatypes/GraphSphere.h"
#include "../../model/Model/ModelCellsSpherical/ModelCellsSpherical.h"

#pragma region Init

// Init for cells (no sinusoids)
DiscreteVoronoi::DiscreteVoronoi(std::vector<CellSpherical *> * cells_, std::ofstream * console_, QProgressBar * progressBar_)
{
    InitDefault(cells_, console_, progressBar_);

    // Do not interface graph
    graph = NULL;
    useNetwork = false;
}

// Init for cells and sinusoids
DiscreteVoronoi::DiscreteVoronoi(std::vector<CellSpherical *> * cells_, GraphSphere * graph_, std::ofstream * console_, QProgressBar * progressBar_)
{
    InitDefault(cells_, console_, progressBar_);

    // Additionally interface graph
    graph = graph_;
    useNetwork = true;
}

// Init
void DiscreteVoronoi::InitDefault(std::vector<CellSpherical *> * cells_, std::ofstream * console_, QProgressBar * progressBar_)
{
    #pragma region Parameters

    // Defines resolution of discretization
    widthOfVoxel = 0.2;

    // Maximal radius a Voronoi-cell for cells can have
    cellCutOffRadius = 0.70;

    // Debug output
    enableDebugOutput = true;

    #pragma endregion

    #pragma region Interface

    // Set interface
    cells = cells_;

    useNetwork = false;

    console = console_;

    progressBar = progressBar_;

    #pragma endregion
}

void DiscreteVoronoi::Run(float widthOfVoxel_, float cellCutOffRadius_, bool enableDebugOutput_)
{
    if (enableDebugOutput) (*console) << "DiscreteVoronoi::Test().\n";

    // Set parameters
    widthOfVoxel = widthOfVoxel_; // Defines resolution of discretization
    cellCutOffRadius = cellCutOffRadius_;// Maximal radius a Voronoi-cell for cells can have
    enableDebugOutput = enableDebugOutput_;// Debug output

    if (enableDebugOutput)
    {
        (*console) << "Use network: " << useNetwork << "\n";
        (*console) << "Width of voxels: " << widthOfVoxel << "\n";
        (*console) << "Cell cutoff radius: " << cellCutOffRadius << "\n";
        (*console) << "Enable debug output: " << enableDebugOutput << "\n";
    }

    // First tests
    if (cells==NULL) return;
    if ((*cells).size() < 1) return;
    if (useNetwork) if (graph==NULL) return;

    // (1) Measure the extension of the volume of interest
    InitMinMaxXYZ();

    // (2) Init voxels
    InitVoxels();

    // (3) Determine which voxels belong to the network
    if (useNetwork) DefineNetwork();

    // (4) Determine which voxels belong to the cells
    DefineCells();

    // (5) Compute statistics for each cell
    Quantify("resultsDVA_Cells.txt", "resultsDVA_CellsAndEdges.txt");

    // Debug output
    OutputPovray("out1.pov", 1);
    OutputPovray("out2.pov", 2);
    OutputPovray("out3.pov", 3);

    if (enableDebugOutput) (*console) << "\nTest complete.\n";
}

#pragma region (1) Measure the extension of the volume of interest

void DiscreteVoronoi::InitMinMaxXYZ()
{
    #pragma region Measure extension of cells

    if (GetNumberOfCells() < 1)  if (enableDebugOutput) (*console) << "Error: Not possible without cells.";

    // X
    minX = GetPositionXOfCell(0) - GetRadiusOfCell(0);
    maxX = GetPositionXOfCell(0) + GetRadiusOfCell(0);
    // Y
    minY = GetPositionYOfCell(0) - GetRadiusOfCell(0);
    maxY = GetPositionYOfCell(0) + GetRadiusOfCell(0);
    // Z
    minZ = GetPositionZOfCell(0) - GetRadiusOfCell(0);
    maxZ = GetPositionZOfCell(0) + GetRadiusOfCell(0);

    for (int i = 1; i<GetNumberOfCells(); i++)
    {
        // X
        if (GetPositionXOfCell(i) - GetRadiusOfCell(i) < minX) minX = GetPositionXOfCell(i) - GetRadiusOfCell(i);
        if (GetPositionXOfCell(i) + GetRadiusOfCell(i) > maxX) maxX = GetPositionXOfCell(i) + GetRadiusOfCell(i);
        // Y
        if (GetPositionYOfCell(i) - GetRadiusOfCell(i) < minY) minY = GetPositionYOfCell(i) - GetRadiusOfCell(i);
        if (GetPositionYOfCell(i) + GetRadiusOfCell(i) > maxY) maxY = GetPositionYOfCell(i) + GetRadiusOfCell(i);
        // Z
        if (GetPositionZOfCell(i) - GetRadiusOfCell(i) < minZ) minZ = GetPositionZOfCell(i) - GetRadiusOfCell(i);
        if (GetPositionZOfCell(i) + GetRadiusOfCell(i) > maxZ) maxZ = GetPositionZOfCell(i) + GetRadiusOfCell(i);
    }

    if (enableDebugOutput) (*console) << "Extension of " << GetNumberOfCells() << " cells is (" << minX << " " << maxX << "),(" << minY << " " << maxY << "),(" << minZ << " " << maxZ << ")\n";

    #pragma endregion
}

#pragma endregion

#pragma region (2) Init voxels

void DiscreteVoronoi::InitVoxels()
{
    if (enableDebugOutput) (*console) << "InitVoxels()\n";

    // Compute resolution
    resX = (int)ceil((maxX - minX) / widthOfVoxel);
    resY = (int)ceil((maxY - minY) / widthOfVoxel);
    resZ = (int)ceil((maxZ - minZ) / widthOfVoxel);

    voxels = (int *)malloc(resX * resY * resZ * sizeof(int));
    if (voxels == NULL)
    {
        if (enableDebugOutput) (*console) << "Error: Memory for " << resX * resY * resZ << " voxels (" << (float)(resX * resY * resZ * sizeof(int)) / (1024. * 1024.) << " MB) could not be allocated.\n";
    }
    else
    {
        if (enableDebugOutput) (*console) << "Memory for " << resX * resY * resZ << " voxels (" << (float)(resX * resY * resZ * sizeof(int)) / (1024. * 1024.) << " MB) allocated.\n";
    }

    // Init to 0
    for (int i=0; i<resX * resY * resZ; i++)
    {
        voxels[i] = 0;
    }
}

#pragma endregion

#pragma region (3) Determine which voxels belong to the network

void DiscreteVoronoi::DefineNetwork()
{
    int count = 0;
    int countOver = 0;
    int x;
    int currentOverlappingEdgePlus1OrZero;

    if (enableDebugOutput)
    {
        (*console) << "DefineNetwork() [" << resX << "]: ";
        (*console).flush();
    }

    #pragma omp parallel for reduction (+:count, countOver)

    for (x = 0; x<resX; x++)
    {
        for (int y= 0; y<resY; y++)
        {
            for (int z = 0; z<resZ; z++)
            {
                count++;

                currentOverlappingEdgePlus1OrZero = IsOverlappingWithNetwork(x * widthOfVoxel+ minX,y * widthOfVoxel+ minY,z* widthOfVoxel+ minZ);

                if (currentOverlappingEdgePlus1OrZero != 0)
                {
                    voxels[x + y * resX + z * resX * resY] = GetNumberOfCells() + currentOverlappingEdgePlus1OrZero; // Voxels: 0=Empty space or Yet unassigned, 1...GetNumberOfCells()=Cells, GetNumberOfCells()+1...GetNumberOfCells()+1+graph.size() = Edges of graph
                    countOver++;
                }
            }
        }
        if (enableDebugOutput)
        {
            (*console) << "*";
            (*console).flush();
        }
    }

    if (enableDebugOutput)
    {
        (*console) << "\nDefineNetwork() " << countOver << " of " << count << " (" << 100. * ((float)countOver / (float)count) << " percent) voxels are network.\n";
        (*console).flush();
    }
}

#pragma endregion

#pragma region (4) Determine which voxels belong to the cells

#pragma region Helper: GetClosestCellIndexOrZero

// Returns from 1 ... CellNumber
int DiscreteVoronoi::GetClosestCellIndexOrZero(float x, float y, float z)
{
    float cellCutOffRadiusQ = cellCutOffRadius * cellCutOffRadius;
    int indexOfClosestCellSoFar = 0; // None found
    float distanceQOfClosestCellSoFar = cellCutOffRadiusQ;

    float distanceQOfCurrent;

    for (int i = 0; i<GetNumberOfCells(); i++)
    {
        distanceQOfCurrent = MathDistQ(x,y,z, GetPositionXOfCell(i), GetPositionYOfCell(i), GetPositionZOfCell(i));

        if (distanceQOfCurrent <= distanceQOfClosestCellSoFar)
        {
            distanceQOfClosestCellSoFar = distanceQOfCurrent;
            indexOfClosestCellSoFar = i+1; // Index is returned from 1 ... n+1 // Voxels: 0=Empty space or Yet unassigned, 1...GetNumberOfCells()=Cells, GetNumberOfCells()+1...GetNumberOfCells()+1+graph.size() = Edges of graph
        }
    }

    return indexOfClosestCellSoFar;
}

#pragma endregion

void DiscreteVoronoi::DefineCells()
{
    int count = 0;
    int countOver = 0;
    int x;
    int firstNetWorkIndexInVoxels = GetNumberOfCells()+1;

    if (enableDebugOutput)
    {
        (*console) << "DefineCells() [" << resX << "]: ";
        (*console).flush();
    }

    #pragma omp parallel for reduction (+:count, countOver)
    for (x = 0; x<resX; x++)
    {
        for (int y= 0; y<resY; y++)
        {
            for (int z = 0; z<resZ; z++)
            {
                int index = x + y * resX + z * resX * resY;

                // Define cells only for non-network voxels
                if (voxels[index] < firstNetWorkIndexInVoxels)
                {
                    // Returns from 1 ... CellNumber, or 0 if empty
                    voxels[index] = GetClosestCellIndexOrZero(x * widthOfVoxel+ minX,y * widthOfVoxel+ minY,z* widthOfVoxel+ minZ);

                    // Count cell overlapped voxels
                    if (voxels[index] > 0) countOver++;
                }
                count++;
            }
        }
        if (enableDebugOutput)
        {
            (*console) << "*";
            (*console).flush();
        }
    }

    if (enableDebugOutput)
    {
        (*console) << "\nDefineCells() " << countOver << " of " << count << " (" << 100. * ((float)countOver / (float)count) << " percent) voxels are cells.\n";
    }
}

#pragma endregion

#pragma region (5) Quantify

#pragma region Helper IsVoxel...

bool DiscreteVoronoi::IsVoxelInContactWithNetwork(int x, int y, int z)
{
    int firstNetWorkIndexInVoxels = GetNumberOfCells()+1;

    if (x > 0)      if (voxels[(x-1) + y * resX + z * resX * resY] >= firstNetWorkIndexInVoxels) return true;
    if (x < resX-1) if (voxels[(x+1) + y * resX + z * resX * resY] >= firstNetWorkIndexInVoxels) return true;

    if (y > 0)      if (voxels[x + (y-1) * resX + z * resX * resY] >= firstNetWorkIndexInVoxels) return true;
    if (y < resY-1) if (voxels[x + (y+1) * resX + z * resX * resY] >= firstNetWorkIndexInVoxels) return true;

    if (z > 0)      if (voxels[x + y * resX + (z-1) * resX * resY] >= firstNetWorkIndexInVoxels) return true;
    if (z < resZ-1) if (voxels[x + y * resX + (z+1) * resX * resY] >= firstNetWorkIndexInVoxels) return true;

    return false;
}

bool DiscreteVoronoi::IsVoxelInContactWithNetwork(int x, int y, int z, int edge)
{
    int valueOfEdgeInVoxels = GetNumberOfCells()+1+edge;

    if (x > 0)      if (voxels[(x-1) + y * resX + z * resX * resY] == valueOfEdgeInVoxels) return true;
    if (x < resX-1) if (voxels[(x+1) + y * resX + z * resX * resY] == valueOfEdgeInVoxels) return true;

    if (y > 0)      if (voxels[x + (y-1) * resX + z * resX * resY] == valueOfEdgeInVoxels) return true;
    if (y < resY-1) if (voxels[x + (y+1) * resX + z * resX * resY] == valueOfEdgeInVoxels) return true;

    if (z > 0)      if (voxels[x + y * resX + (z-1) * resX * resY] == valueOfEdgeInVoxels) return true;
    if (z < resZ-1) if (voxels[x + y * resX + (z+1) * resX * resY] == valueOfEdgeInVoxels) return true;

    return false;
}

bool DiscreteVoronoi::IsVoxelInContactWithOtherCell(int otherThanIndex, int x, int y, int z)
{
    int firstNetWorkIndexInVoxels = GetNumberOfCells()+1;

    // Prelim: Implement more efficient but computing target index only once per line
    if (x > 0)      if (voxels[(x-1) + y * resX + z * resX * resY] > 0) if (voxels[(x-1) + y * resX + z * resX * resY] < firstNetWorkIndexInVoxels) if (voxels[(x-1) + y * resX + z * resX * resY] != otherThanIndex) return true;
    if (x < resX-1) if (voxels[(x+1) + y * resX + z * resX * resY] > 0) if (voxels[(x+1) + y * resX + z * resX * resY] < firstNetWorkIndexInVoxels) if (voxels[(x+1) + y * resX + z * resX * resY] != otherThanIndex) return true;

    if (y > 0)      if (voxels[x + (y-1) * resX + z * resX * resY] > 0) if (voxels[x + (y-1) * resX + z * resX * resY] < firstNetWorkIndexInVoxels) if (voxels[x + (y-1) * resX + z * resX * resY] != otherThanIndex) return true;
    if (y < resY-1) if (voxels[x + (y+1) * resX + z * resX * resY] > 0) if (voxels[x + (y+1) * resX + z * resX * resY] < firstNetWorkIndexInVoxels) if (voxels[x + (y+1) * resX + z * resX * resY] != otherThanIndex) return true;

    if (z > 0)      if (voxels[x + y * resX + (z-1) * resX * resY] > 0) if (voxels[x + y * resX + (z-1) * resX * resY] < firstNetWorkIndexInVoxels) if (voxels[x + y * resX + (z-1) * resX * resY] != otherThanIndex) return true;
    if (z < resZ-1) if (voxels[x + y * resX + (z+1) * resX * resY] > 0) if (voxels[x + y * resX + (z+1) * resX * resY] < firstNetWorkIndexInVoxels) if (voxels[x + y * resX + (z+1) * resX * resY] != otherThanIndex) return true;

    return false;
}

bool DiscreteVoronoi::IsVoxelInner(int x, int y, int z)
{
    int type = voxels[x + y * resX + z * resX * resY];

    if (x > 0)      if (voxels[(x-1) + y * resX + z * resX * resY] != type) return false;
    if (x < resX-1) if (voxels[(x+1) + y * resX + z * resX * resY] != type) return false;

    if (y > 0)      if (voxels[x + (y-1) * resX + z * resX * resY] != type) return false;
    if (y < resY-1) if (voxels[x + (y+1) * resX + z * resX * resY] != type) return false;

    if (z > 0)      if (voxels[x + y * resX + (z-1) * resX * resY] != type) return false;
    if (z < resZ-1) if (voxels[x + y * resX + (z+1) * resX * resY] != type) return false;

    return true;
}

#pragma endregion

void DiscreteVoronoi::Quantify(const char* fileNamePerCell, const char* fileNamePerCellAndEdge)
{
    if (enableDebugOutput)
    {
        (*console) << "Quantify() [" << GetNumberOfCells() << "]: ";
        (*console).flush();
    }

    #pragma region Quantifications per cell

    #pragma region Prepare file

    FILE *f = fopen(fileNamePerCell, "w");
    if (f == NULL)
    {
        if (enableDebugOutput) (*console) << "Could not open file in DiscreteVoronoi::Quantify()";
        return;
    }

    #pragma endregion

    // Header
    fprintf(f, "CellIndex\tVoxelsOfCell\tInnerVoxelsOfCell\tVoxelsBorderingToNetwork\tPercentOfSurfaceVoxelsBorderingToNetwork\tVoxelsBorderingToOtherCells\tPercentOfSurfaceVoxelsBorderingToOtherCells\n");

    // For all cells
    for (int c = 0; c < GetNumberOfCells(); c++)
    {
        int countVoxelsOfCell = 0;
        int countVoxelsBorderingToNetwork = 0;
        int countVoxelsBorderingToOtherCells = 0;
        int countVoxelsInner = 0;
        int x;

        #pragma omp parallel for reduction (+:countVoxelsOfCell, countVoxelsBorderingToNetwork, countVoxelsBorderingToOtherCells, countVoxelsInner)
        for (x = 0; x<resX; x++)
        {
            for (int y= 0; y<resY; y++)
            {
                for (int z = 0; z<resZ; z++)
                {
                    int index = x + y * resX + z * resX * resY;

                    // Quantify only for cell voxels
                    if (voxels[index] == c+1)
                    {
                        countVoxelsOfCell++;

                        if (IsVoxelInContactWithNetwork(x,y,z))
                        {
                            countVoxelsBorderingToNetwork++;
                        }
                        else if (IsVoxelInner(x,y,z))
                        {
                            countVoxelsInner++;
                        }
                        else if (IsVoxelInContactWithOtherCell(voxels[index], x,y,z))
                        {
                            countVoxelsBorderingToOtherCells++;
                        }
                    }
                }
            }
        }

        if (countVoxelsOfCell != 0)
        {
            fprintf(f, "%i\t%i\t%i\t%i\t%.3f\t%i\t%.3f\n",
                    c+1, // CellIndex
                    countVoxelsOfCell, // VoxelsOfCell
                    countVoxelsInner, // InnerVoxelsOfCell
                    countVoxelsBorderingToNetwork, // VoxelsBorderingToNetwork
                    100. * (float)countVoxelsBorderingToNetwork / (float)(countVoxelsOfCell - countVoxelsInner), // PercentOfSurfaceVoxelsBorderingToNetwork
                    countVoxelsBorderingToOtherCells, // VoxelsBorderingToOtherCells
                    100. * (float)countVoxelsBorderingToOtherCells / (float)(countVoxelsOfCell - countVoxelsInner)); // PercentOfSurfaceVoxelsBorderingToOtherCells
        }
    }

    #pragma region Close file

    fclose(f);

    #pragma endregion

    #pragma endregion

    if (useNetwork)
    {
        #pragma region Quantification of contact voxels per cell and edge

        #pragma region Prepare file

        f = fopen(fileNamePerCellAndEdge, "w");
        if (f == NULL)
        {
            if (enableDebugOutput) (*console) << "Could not open file in DiscreteVoronoi::Quantify()";
            return;
        }

        #pragma endregion

        // Header
        fprintf(f, "CellIndex\tEdgeIndex\tVoxelsBorderingToThisNetworkEdge\tPercentOfSurfaceVoxelsBorderingToThisNetworkEdge\n");

        // For all cells
        for (int c = 0; c < GetNumberOfCells(); c++)
        {
            for ( unsigned int e=0; e<graph->mvEdge.size(); e++ )
            {
                int countVoxelsOfCell = 0;
                int countVoxelsBorderingToNetwork = 0;
                int countVoxelsInner = 0;
                int x;

                #pragma omp parallel for reduction (+:countVoxelsOfCell, countVoxelsBorderingToNetwork, countVoxelsInner)
                for (x = 0; x<resX; x++)
                {
                    for (int y= 0; y<resY; y++)
                    {
                        for (int z = 0; z<resZ; z++)
                        {
                            int index = x + y * resX + z * resX * resY;

                            // Quantify only for cell voxels
                            if (voxels[index] == c+1)
                            {
                                countVoxelsOfCell++;

                                if (IsVoxelInContactWithNetwork(x,y,z,e))
                                {
                                    countVoxelsBorderingToNetwork++;
                                }
                                else if (IsVoxelInner(x,y,z))
                                {
                                    countVoxelsInner++;
                                }
                            }
                        }
                    }
                }

                if (countVoxelsOfCell != 0)
                {
                    fprintf(f, "%i\t%i\t%i\t%i\t%.1f\t%i\t%.3f\n",
                            c+1, // CellIndex (1..)
                            e+1, // EdgeIndex (1..)
                            countVoxelsBorderingToNetwork, // VoxelsBorderingToThisNetworkEdge
                            100. * (float)countVoxelsBorderingToNetwork / (float)(countVoxelsOfCell - countVoxelsInner)); // PercentOfSurfaceVoxelsBorderingToThisNetworkEdge
                }
            }
        }


        #pragma region Close file

        fclose(f);

        #pragma endregion

        #pragma endregion
    }

    if (enableDebugOutput)
    {
        (*console) << "\nRefineAndQuantify() complete\n";
        (*console).flush();
    }
}

#pragma endregion

#pragma region Debug output using povray

void DiscreteVoronoi::OutputPovray(const char* fileName, int version)
{
    if (enableDebugOutput)
    {
        (*console) << "OutputPovray()";
        (*console).flush();
    }

    #pragma region Prepare file

    FILE *f = fopen(fileName, "w");
    if (f == NULL)
    {
        if (enableDebugOutput) (*console) << "Could not open file in DiscreteVoronoi::OutputPovray()";
        return;
    }

    #pragma endregion

    #pragma region Povray header

    fprintf(f, "global_settings { max_trace_level 5 }\n");
    fprintf(f, "#macro rotation() 100 #end\n");
    fprintf(f, "#macro cut_rotation() 20 #end\n");
    fprintf(f, "#macro cut_depth() 5.000000 #end\n");
    fprintf(f, "#macro fc1() texture {pigment { color rgb  #end\n");
    fprintf(f, "#macro fc2() }} finish { ambient 0.25 diffuse 0.8 roughness 0.001} #end\n");
    fprintf(f, "\n");
    fprintf(f, "camera { location <23, 30, 2.8> look_at  <2, 5, 0>}\n");
    fprintf(f, "light_source { 1*x color rgb <0.7 0.5 0.5> area_light <8, 8, 0> <0, 8, 8>  4, 4 adaptive 0 jitter circular orient translate <40, 80, -40> }\n");
    fprintf(f, "light_source { 1*x color rgb 0.2 area_light <8, 8, 0> <0, 8, 8>  4, 4 adaptive 0 jitter circular orient translate < 0,  0, -50> }\n");
    fprintf(f, "light_source { 1*x color rgb <0.4 0.4 0.4> area_light <8, 8, 0> <0, 8, 8>  4, 4 adaptive 0 jitter circular orient translate <15, 23, 15> }\n");
    fprintf(f, "background { color rgb 0.2 }\n");
    fprintf(f, "\n");
    fprintf(f, "union { // Start of union\n");

    #pragma endregion

    #pragma region Version 1 - Per voxel output (Show cells and network)

    if (version == 1)
    {
        for (int x = 0; x<resX; x++)
        {
            for (int y= 0; y<resY; y++)
            {
                for (int z = 0; z<resZ; z++)
                {
                    // Network voxel?
                    if (voxels[x + y * resX + z * resX * resY] == 1)
                    {
                        fprintf(f, "box { <%.4f, %.4f, %.4f>,  <%.4f, %.4f, %.4f> fc1() <1.0 0.5 0.5> fc2() }\n",x * widthOfVoxel + minX,y * widthOfVoxel + minY,z* widthOfVoxel + minZ,(x+1) * widthOfVoxel + minX,(y+1) * widthOfVoxel + minY,(z+1)* widthOfVoxel + minZ);
                    }

                    // Cell voxel?
                    else if (voxels[x + y * resX + z * resX * resY] > 1)
                    {
                        fprintf(f, "box { <%.4f, %.4f, %.4f>,  <%.4f, %.4f, %.4f> fc1() <0.5 0.5 1.0> fc2() }\n",x * widthOfVoxel + minX,y * widthOfVoxel + minY,z* widthOfVoxel + minZ,(x+1) * widthOfVoxel + minX,(y+1) * widthOfVoxel + minY,(z+1)* widthOfVoxel + minZ);
                    }
                }
            }
        }
    }

    #pragma endregion

    #pragma region Version 2 - Colored per cells output

    else if (version == 2)
    {
        int countTo = GetNumberOfCells() / 10;

        if (enableDebugOutput) (*console) << " Cells [" << countTo << "]: ";

        // Cells
        for (int c = 0; c < GetNumberOfCells(); c++)
        {
            fprintf(f, "// Cell %i of %i\n",c+1, GetNumberOfCells());

            for (int x = 0; x<resX; x++)
            {
                for (int y= 0; y<resY; y++)
                {
                    for (int z = 0; z<resZ; z++)
                    {
                        int index = x + y * resX + z * resX * resY;

                        // Cells
                        if (voxels[index] == c+2)
                        {
                            if (IsVoxelInContactWithNetwork(x,y,z))
                            {
                                fprintf(f, "box { <%.4f, %.4f, %.4f>,  <%.4f, %.4f, %.4f> fc1() <1.0 0.3 1.0> fc2() }\n",x * widthOfVoxel + minX,y * widthOfVoxel + minY,z* widthOfVoxel + minZ,(x+1) * widthOfVoxel + minX,(y+1) * widthOfVoxel + minY,(z+1)* widthOfVoxel + minZ);
                            }
                            else if (IsVoxelInContactWithOtherCell(voxels[index], x,y,z))
                            {
                                fprintf(f, "box { <%.4f, %.4f, %.4f>,  <%.4f, %.4f, %.4f> fc1() <0.3 1.0 0.3> fc2() }\n",x * widthOfVoxel + minX,y * widthOfVoxel + minY,z* widthOfVoxel + minZ,(x+1) * widthOfVoxel + minX,(y+1) * widthOfVoxel + minY,(z+1)* widthOfVoxel + minZ);
                            }
                            else if (IsVoxelInner(x,y,z) == false)
                            {
                                fprintf(f, "box { <%.4f, %.4f, %.4f>,  <%.4f, %.4f, %.4f> fc1() <0.5 0.5 1.0> fc2() }\n",x * widthOfVoxel + minX,y * widthOfVoxel + minY,z* widthOfVoxel + minZ,(x+1) * widthOfVoxel + minX,(y+1) * widthOfVoxel + minY,(z+1)* widthOfVoxel + minZ);
                            }
                            // Do not draw inner
                        }
                    }
                }
            }

            if (enableDebugOutput) if (c % 10 == 0) (*console) << "*";

        }
        // Network
        fprintf(f, "\n\n// Network\n");
        for (int x = 0; x<resX; x++)
        {
            for (int y= 0; y<resY; y++)
            {
                for (int z = 0; z<resZ; z++)
                {
                    // Network voxel?
                    if (voxels[x + y * resX + z * resX * resY] == 1)
                    {
                        fprintf(f, "box { <%.4f, %.4f, %.4f>,  <%.4f, %.4f, %.4f> fc1() <1.0 0.5 0.5> fc2() }\n",x * widthOfVoxel + minX,y * widthOfVoxel + minY,z* widthOfVoxel + minZ,(x+1) * widthOfVoxel + minX,(y+1) * widthOfVoxel + minY,(z+1)* widthOfVoxel + minZ);
                    }
                }
            }
        }

    }

    #pragma endregion

    #pragma region Version 3 - Colored per cells output but only central slice

    else if (version == 3)
    {

        int countTo = GetNumberOfCells() / 10;
        int z = (int)(0.55f * (float)resZ);

        if (enableDebugOutput) (*console) << " Cells [" << countTo << "]: ";

        // Cells
        for (int c = 0; c < GetNumberOfCells(); c++)
        {
            fprintf(f, "// Cell %i of %i\n",c+1, GetNumberOfCells());

            for (int x = 0; x<resX; x++)
            {
                for (int y= 0; y<resY; y++)
                {
                    int index = x + y * resX + z * resX * resY;

                    // Cells
                    if (voxels[index] == c+2)
                    {
                        if (IsVoxelInContactWithNetwork(x,y,z))
                        {
                            fprintf(f, "box { <%.4f, %.4f, %.4f>,  <%.4f, %.4f, %.4f> fc1() <1.0 0.3 1.0> fc2() }\n",x * widthOfVoxel + minX,y * widthOfVoxel + minY,z* widthOfVoxel + minZ,(x+1) * widthOfVoxel + minX,(y+1) * widthOfVoxel + minY,(z+1)* widthOfVoxel + minZ);
                        }
                        else if (IsVoxelInContactWithOtherCell(voxels[index], x,y,z))
                        {
                            fprintf(f, "box { <%.4f, %.4f, %.4f>,  <%.4f, %.4f, %.4f> fc1() <0.3 1.0 0.3> fc2() }\n",x * widthOfVoxel + minX,y * widthOfVoxel + minY,z* widthOfVoxel + minZ,(x+1) * widthOfVoxel + minX,(y+1) * widthOfVoxel + minY,(z+1)* widthOfVoxel + minZ);
                        }
                        else if (IsVoxelInner(x,y,z) == false)
                        {
                            fprintf(f, "box { <%.4f, %.4f, %.4f>,  <%.4f, %.4f, %.4f> fc1() <0.5 0.5 1.0> fc2() }\n",x * widthOfVoxel + minX,y * widthOfVoxel + minY,z* widthOfVoxel + minZ,(x+1) * widthOfVoxel + minX,(y+1) * widthOfVoxel + minY,(z+1)* widthOfVoxel + minZ);
                        }
                        // Do not draw inner
                    }
                }
            }

            if (enableDebugOutput) if (c % 10 == 0) (*console) << "*";

        }
        // Network
        fprintf(f, "\n\n// Network\n");
        for (int x = 0; x<resX; x++)
        {
            for (int y= 0; y<resY; y++)
            {
                // Network voxel?
                if (voxels[x + y * resX + z * resX * resY] == 1)
                {
                    fprintf(f, "box { <%.4f, %.4f, %.4f>,  <%.4f, %.4f, %.4f> fc1() <1.0 0.5 0.5> fc2() }\n",x * widthOfVoxel + minX,y * widthOfVoxel + minY,z* widthOfVoxel + minZ,(x+1) * widthOfVoxel + minX,(y+1) * widthOfVoxel + minY,(z+1)* widthOfVoxel + minZ);
                }
            }
        }
    }

    #pragma endregion

    #pragma region Povray footer

    fprintf(f, "rotate <270, 360*clock+rotation(), 0>\n");
    fprintf(f, "} // End of union\n");

    #pragma endregion

    #pragma region Close file

    fclose(f);

    #pragma endregion

    if (enableDebugOutput) (*console) << " complete.\n";

}

#pragma endregion

#pragma region Helper

float DiscreteVoronoi::MathDot(float x1, float y1, float z1, float x2, float y2, float z2)
{
    return x1 * x2 + y1 * y2 + z1 * z2;
}

float DiscreteVoronoi::MathNorm(float vx, float vy, float vz)
{
    return sqrt(MathDot(vx, vy, vz, vx, vy, vz));
}

float DiscreteVoronoi::MathDist(float x1, float y1, float z1, float x2, float y2, float z2)
{
    return MathNorm(x1 - x2, y1-y2, z1-z2);
}

float DiscreteVoronoi::MathDistQ(float x1, float y1, float z1, float x2, float y2, float z2)
{
    return (x1 - x2) * (x1 - x2) +
           (y1 - y2) * (y1 - y2) +
           (z1 - z2) * (z1 - z2);
}

float DiscreteVoronoi::GetDistanceFromPointToLineSegment(float px, float py, float pz, float s0x, float s0y, float s0z, float s1x, float s1y, float s1z)
{
    float vx = s1x - s0x;
    float vy = s1y - s0y;
    float vz = s1z - s0z;

    float wx = px - s0x;
    float wy = py - s0y;
    float wz = pz - s0z;

    // Vector v = S.P1 - S.P0;
    // Vector w = P - S.P0;

    double c1 = MathDot(wx, wy, wz, vx, vy, vz);
    if ( c1 <= 0 )
        return MathDist(px, py, pz, s0x, s0y, s0z);

    double c2 = MathDot(vx, vy, vz,vx, vy, vz);
    if ( c2 <= c1 )
        return MathDist(px, py, pz, s1x, s1y, s1z);

    double b = c1 / c2;

    float pbx = s0x + b * vx;
    float pby = s0y + b * vy;
    float pbz = s0z + b * vz;
    //Point Pb = S.P0 + b * v;

    return MathDist(px, py, pz, pbx, pby, pbz);
}

#pragma endregion

#pragma region Interface methods

#pragma region Cells (e.g. Hepatocytes)

int DiscreteVoronoi::GetNumberOfCells()
{
    return (*cells).size();

    // This part is specific to the OLD CellSys and must be replaced by a corresponding method of the NEW CellSys
    // return zm.currentNumber;
}

float DiscreteVoronoi::GetPositionXOfCell(int id)
{
    return (*cells).at(id)->position.x;

    // This part is specific to the OLD CellSys and must be replaced by a corresponding method of the NEW CellSys
    // return (*zm.zelle[id-1]).x; // In old Cellsys array starts with 1
}

float DiscreteVoronoi::GetPositionYOfCell(int id)
{
    return (*cells).at(id)->position.y;

    // This part is specific to the OLD CellSys and must be replaced by a corresponding method of the NEW CellSys
    //return (*zm.zelle[id-1]).y; // In old Cellsys array starts with 1
}

float DiscreteVoronoi::GetPositionZOfCell(int id)
{
    return (*cells).at(id)->position.z;

    // This part is specific to the OLD CellSys and must be replaced by a corresponding method of the NEW CellSys
    // return (*zm.zelle[id-1]).z; // In old Cellsys array starts with 1
}

float DiscreteVoronoi::GetRadiusOfCell(int id)
{
    return (*cells).at(id)->mRadius;

    // This part is specific to the OLD CellSys and must be replaced by a corresponding method of the NEW CellSys
    // return (*zm.zelle[id-1]).radius; // In old Cellsys array starts with 1
}

#pragma endregion

#pragma region Networks (e.g. Sinusoids)

unsigned int DiscreteVoronoi::IsOverlappingWithNetwork(float x, float y, float z)
{
    for ( unsigned int i=0; i<graph->mvEdge.size(); ++i )
    {
        double avgSegmentRadius = (graph->mvEdge[i]->mRadius_start + graph->mvEdge[i]->mRadius_end)/2;
        if ( GetDistanceFromPointToLineSegment( x, y, z,
                                                graph->mvEdge[i]->mpStart->position.x,
                                                graph->mvEdge[i]->mpStart->position.y,
                                                graph->mvEdge[i]->mpStart->position.z,
                                                graph->mvEdge[i]->mpEnd->position.x,
                                                graph->mvEdge[i]->mpEnd->position.y,
                                                graph->mvEdge[i]->mpEnd->position.z) <= avgSegmentRadius )
            return i+1; // As 0 signals NO overlap
    }

    return 0;

    // This part is specific to the OLD CellSys and must be replaced by a corresponding method of the NEW CellSys
    /*
    for (int i = 1; i<=vasc.graph.currentNodeCount; i++)
    {
        for (int j = 1; j<=vasc.graph.nodes[i].numberOfLinksOutgoing; j++)
        {
            if (GetDistanceFromPointToLineSegment(x,y,z,
                                                  vasc.graph.nodes[i].x,
                                                  vasc.graph.nodes[i].y,
                                                  vasc.graph.nodes[i].z,
                                                  vasc.graph.nodes[vasc.graph.nodes[i].linksOutgoing[j].connectedId].x,
                                                  vasc.graph.nodes[vasc.graph.nodes[i].linksOutgoing[j].connectedId].y,
                                                  vasc.graph.nodes[vasc.graph.nodes[i].linksOutgoing[j].connectedId].z) <= vasc.graph.nodes[i].radius) return true;

        }
    }

    */
}



#pragma endregion

#pragma endregion
