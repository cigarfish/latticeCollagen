
#include "Observation.h"
#include "../../Core.h"

#include "../Model/ModelCellsSpherical/ModelCellsSpherical.h"

#include "../../../tools/Tools.h"



#pragma region Collection of observation methods

#pragma region Number of cells (all, subTypes, etc.)

int Observation::GetNumberOfCells(std::vector<CellSpherical *> cells, int mode)
{
    #pragma region 0: All cells (all subTypes)

    if (mode == 0)
    {
        return (int)cells.size();
    }
    #pragma endregion

    #pragma region 1: Only proliferating cells (all subTypes)

    else if (mode == 1)
    {
        int n=0;
        for ( unsigned int i=0; i<cells.size(); i++ )
        {
            if ( (cells.at(i)->cellcycleState & Cell::StateQuiescent) == false) n++;
        }
        return n;
    }
    #pragma endregion

    #pragma region 2: All cells (only subType 0)

    else if (mode == 2)
    {
        int n=0;
        for ( unsigned int i=0; i<cells.size(); i++ )
        {
            if (cells.at(i)->mSubType == 0) n++;
        }
        return n;
    }

    #pragma endregion

    #pragma region 3: All cells (only subType 1)

    else if (mode == 3)
    {
        int n=0;
        for ( unsigned int i=0; i<cells.size(); i++ )
        {
            if (cells.at(i)->mSubType == 1) n++;
        }
        return n;
    }

    #pragma endregion

    #pragma region 4: Only proliferating cells (only subType 0)

    else if (mode == 4)
    {
        int n=0;
        for ( unsigned int i=0; i<cells.size(); i++ )
        {
            if (cells.at(i)->mSubType == 0)
            {
                if ( (cells.at(i)->cellcycleState & Cell::StateQuiescent) == false) n++;
            }
        }
        return n;
    }

    #pragma endregion

    #pragma region 5: Only proliferating cells (only subType 1)

    else if (mode == 5)
    {
        int n=0;
        for ( unsigned int i=0; i<cells.size(); i++ )
        {
            if (cells.at(i)->mSubType == 1)
            {
                if ( (cells.at(i)->cellcycleState & Cell::StateQuiescent) == false) n++;
            }
        }
        return n;
    }

    #pragma endregion

    #pragma region Default: Error: Unknown mode

    else
    {
        // Prelim: Use output class to signal error here
        printf("ERROR: Unkown mode in Observation::GetNumberOfCells.\n");
        return 0;
    }

    #pragma endregion
}

#pragma endregion

#pragma region Radius of gyration (subType 0)

double Observation::GetRadiusOfGyration(std::vector<CellSpherical *> cells, int mode)
{
    #pragma region Mode 0: All cells in two dimensions

    if (mode == 0)
    {
        double r = 0;
        double comX, comY, comZ; // Center of mass
        int n = 0; // Number of cells of subType 0

        GetCenterOfMass(cells, &comX, &comY, &comZ);

        for ( unsigned int i=0; i<cells.size(); i++ )
        {
            if (cells.at(i)->mSubType == 0)
            {
                r += (cells.at(i)->position.x - comX) * (cells.at(i)->position.x - comX) +
                     (cells.at(i)->position.y - comY) * (cells.at(i)->position.y - comY);
				     (cells.at(i)->position.z - comZ) * (cells.at(i)->position.z - comZ);
                n++;
            }

        } // for all cells
        return sqrt(r / (double)n);
    }

    #pragma endregion

    #pragma region Default: Error: Unknown mode

    else
    {
        // Prelim: Use output class to signal error here
        printf("ERROR: Unkown mode in Observation::GetRadiusOfGyration.\n");
        return 0;
    }

    #pragma endregion
}

void Observation::GetCenterOfMass(std::vector<CellSpherical *> cells, double *comX, double *comY, double *comZ)
{
    *comX = 0;
    *comY = 0;
	*comZ = 0;

    int n= 0;

    for ( unsigned int i=0; i<cells.size(); i++ )
    {
        if (cells.at(i)->mSubType == 0)
        {
            *comX += cells.at(i)->position.x;
            *comY += cells.at(i)->position.y;
			*comZ += cells.at(i)->position.z;
            n++;
        }
    } // for all cells of subType 0

    *comX /= (double)n;
    *comY /= (double)n;
	*comZ /= (double)n;
}

#pragma endregion

#pragma endregion

Observation::Observation( CSModel * observedModel, int type )
{
    // t1m TODO:  make Observation to be able to handle other models, or make it
    //   a delegate that uses this implementation for handling
    //   ModelCellsSpherical, and implement Observation modules for each model.

    // t1m preliminary:  use the first model of type ModelCellsSpherical.
    mpModelCellsSpherical =
        dynamic_cast<ModelCellsSpherical *>( observedModel );

    // t1m not flexible, if this is the only point, where to set member model.
    //   The value may be changed in the GUI, while the same model is active.
    modelType = type;
}

#pragma region Output description file

// Writes a description file for the observed simulation.
void Observation::ObserveModelDescription()
{
    std::ofstream out;

    if (modelType == 0)
    {
        int column = 1;

        #pragma region Automatically create output directory if it does not already exist (Windows-only)

#ifdef WIN32

        LPCWSTR outputDir = L"output";
        CreateDirectory(outputDir, NULL);

#endif

        #pragma endregion

        // Construct filename
        std::string filename;
        filename.append("output/");
        filename.append(mpModelCellsSpherical->name);
        filename.append(".desc");

        core->tools->output->logfile << "Name in observation: " << mpModelCellsSpherical->name << "\n";

        out.open(filename.c_str());
        out.clear();

        out << "Column\tDescription\n\n";
        out << column << "\tTime in internal (dimensionless) units\n";
        column++;
        out << column << "\tTime in days\n";
        column++;
        out << column << "\tNumber of cells (subType 0)\n";
        column++;
        out << column << "\tNumber of cells (subType 1)\n";
        column++;
        out << column << "\tRadius of gyration of cell population in internal (dimensionless) units (subType 0)\n";
        column++;
        out << column << "\tCell population diameter in micrometer (based on radius of gyration) (subType 0)\n";
        column++;
        out << column << "\tCell population center of mass (x-coordinate)\n";
        column++;
        out << column << "\tCell population center of mass (y-coordinate)\n";
        column++;
		out << column << "\tCell population center of mass (z-coordinate)\n";
        column++;
        out << column << "\tLast absolute force (nN) Minimum of population\n";
        column++;
        out << column << "\tLast absolute force (nN) Average of population\n";
        column++;
        out << column << "\tLast absolute force (nN) Maximum of population\n";
        column++;
		out << column << "\tRealtime in seconds since last start of simulation\n";
        column++;
    }

    out.close();

	timeZero = time(NULL);
}

#pragma endregion

// Writes the measurements
void Observation::ObserveModelMeasures()
{
    std::ofstream out;

    // ModelType 0: ModelCellsSpherical
    if (modelType == 0)
    {
        core->tools->output->logfile << "Name in observation measures: " << mpModelCellsSpherical->name << "\n";

        #pragma region Automatically create output directory if it does not already exist (Windows-only)

#ifdef WIN32

        LPCWSTR outputDir = L"output";
        CreateDirectory(outputDir, NULL);

#endif

        #pragma endregion

        // Construct filename
        std::string filename;
        filename.append("output/");
        filename.append( mpModelCellsSpherical->name );
        filename.append(".txt");

        out.open(filename.c_str(), std::fstream::app);

        // Preparatory calculations
        double comx, comy, comz;
        GetCenterOfMass(mpModelCellsSpherical->cells, &comx, &comy, &comz);

        double minAF = 0;
        double meanAF = 0;
        double maxAF = 0;
        GetCellsMeasure(mpModelCellsSpherical->cells, 0, &minAF, &meanAF, &maxAF);

        out << mpModelCellsSpherical->time << "\t"
            << mpModelCellsSpherical->biolink->getTimeInDays(mpModelCellsSpherical->time) << "\t"
            << GetNumberOfCells(mpModelCellsSpherical->cells, 2) << "\t"
            << GetNumberOfCells(mpModelCellsSpherical->cells, 3) << "\t"
            << GetRadiusOfGyration(mpModelCellsSpherical->cells, 0) << "\t"
            << GetDiameterInMicrometerBasedOnRgyr(mpModelCellsSpherical->cells, 0) << "\t"
            << comx << "\t"
            << comy << "\t"
			<< comz << "\t"
            << minAF << "\t"
            << meanAF << "\t"
            << maxAF << "\t"
			<< (time(NULL) - timeZero) << "\t"

            << "\n";
    }

    out.close();
}

double Observation::GetDiameterInMicrometerBasedOnRgyr(std::vector<CellSpherical *> cells, int mode)
{
    if (mode == 0)
    {
        double rgyr_intern = GetRadiusOfGyration(cells, 0);
        return (double)(sqrt(2.) * ((ModelCellsSpherical *)mpModelCellsSpherical)->biolink->getLengthInMicrometers(1) * rgyr_intern * 2.);
    }
    return 0;
}


// Prelim: Writes to logfile
void Observation::ObserveCellPopulationSnapshot(std::vector<CellSpherical *> cells)
{
    core->tools->output->consoleModelCellsSpherical_Simulation << "\n\n" <<
            "Cell number\t" <<
            "Cell position x\t" <<
            "Cell position y\t" <<
            "Last friction coefficient\t" <<
            "Last contact area\t" <<
            "Last force absolute" <<
            "Young modulus" <<
            "Poisson ratio" <<
            "Current Radius" <<
            "Delta Radius" <<
            "\n";

    for ( unsigned int i=0; i<cells.size(); i++ )
    {
        core->tools->output->logfile <<
                                     i << "\t" <<
                                     cells.at(i)->position.x << "\t" <<
                                     cells.at(i)->position.y << "\t" <<
                                     cells.at(i)->frictionCoefficient << "\t" <<
                                     cells.at(i)->lastContactArea << "\t" <<
                                     cells.at(i)->lastForceAbsolute << "\t" <<
                                     cells.at(i)->youngModulus << "\t" <<
                                     cells.at(i)->poissonRatio << "\t" <<
                                     cells.at(i)->mRadius  << "\t" <<
                                     cells.at(i)->deltaRadius << "\t" <<
                                     "\n";
    }
}

// Measure:
// 0 ...Absolute force in nN
void Observation::GetCellsMeasure(std::vector<CellSpherical *> cells, int measure, double *min, double *mean, double *max)
{
    if (cells.size()==0) return;

    #pragma region 0: Last force absolute

    if (measure == 0)
    {
        *min = cells.at(0)->lastForceAbsolute;
        *mean = 0;
        *max = cells.at(0)->lastForceAbsolute;
        for ( unsigned int i=0; i<cells.size(); i++ )
        {
            if (cells.at(i)->lastForceAbsolute > *max) *max = mpModelCellsSpherical->biolink->getForceInNanoNewton(cells.at(i)->lastForceAbsolute);
            if (cells.at(i)->lastForceAbsolute < *min) *min = mpModelCellsSpherical->biolink->getForceInNanoNewton(cells.at(i)->lastForceAbsolute);
            *mean += mpModelCellsSpherical->biolink->getForceInNanoNewton(cells.at(i)->lastForceAbsolute);
        }

    }
    else
    {
        core->tools->output->PrintError("Unknown measure parameter in Observation::GetCellsMeasure");
    }

    #pragma endregion

    *mean /= (double)cells.size();
}

void Observation::WritePOV(int mode)
{
    #pragma region Automatically create output directory if it does not already exist (Windows-only)

#ifdef WIN32

    LPCWSTR outputDir = L"output";
    CreateDirectory(outputDir, NULL);

#endif

    #pragma endregion

    #pragma region Mode 0: Spherical cell model

    if (mode == 0)
    {
        #pragma region Construct filename

        std::string filename;
        filename.append("output/pov/");
        filename.append(mpModelCellsSpherical->name);
        // Prelim: Add number here
        filename.append(".pov");

        #pragma endregion

        #pragma region Open file stream

        std::fstream outfile;
        outfile.open(filename.c_str(), std::ios_base::out);
        if (outfile.bad())
        {
            core->tools->output->PrintError("File open failed in Observation::WritePOV()");
            return;
        }

        #pragma endregion

        core->tools->output->consoleModelCellsSpherical_Simulation << "Write povray snapshot in " << filename << ".\n";

        #pragma region Global settings and macros

        outfile << "global_settings { max_trace_level 5 }\n";

        outfile << "// Rotation macro for videos\n";
        outfile << "#macro \n";
        outfile << "rotation() 0\n";
        outfile << "#end\n";
        outfile << "// Defines the surface finish of cells\n";
        outfile << "#macro fin_cells()\n";
        outfile << "finish { ambient 0.7 diffuse 0.8 roughness 0.001}\n";
        outfile << "#end\n";

        #pragma endregion

        #pragma region Camera

        outfile << "camera { location <0, 0, 40> look_at  <0, 0, 0> }\n";

        #pragma endregion

        #pragma region Lights

        outfile << "light_source { 1*x color rgb 0.7 area_light <8, 8, 0> <0, 8, 8>  4, 4 adaptive 0 jitter circular orient translate <40, 80, -40> }\n";
        outfile << "light_source { 1*x color rgb 0.7 area_light <8, 8, 0> <0, 8, 8>  4, 4 adaptive 0 jitter circular orient translate <-40, -80, 40> }\n";

        #pragma endregion

        #pragma region Background color

        outfile << "background { color rgb 0.2 }\n";

        #pragma endregion

        #pragma region Cell population

        // Refresh colors
        mpModelCellsSpherical->UpdateCellsStaining();

        outfile << "union\n{\n";
        for (unsigned int i=0; i < mpModelCellsSpherical->cells.size(); i++)
        {
            outfile << "sphere {<" << mpModelCellsSpherical->cells.at(i)->position.x << "," <<
                    mpModelCellsSpherical->cells.at(i)->position.y << "," <<
                    mpModelCellsSpherical->cells.at(i)->position.z << ">, " <<
                    mpModelCellsSpherical->cells.at(i)->mRadius << " texture {pigment { color rgb <" <<
                    mpModelCellsSpherical->cells.at(i)->color.red << " " <<
                    mpModelCellsSpherical->cells.at(i)->color.green << " " <<
                    mpModelCellsSpherical->cells.at(i)->color.blue << "> }} fin_cells() }\n";
        }

        outfile << "rotate <0, 360*clock+rotation(), 0>\n";
        outfile << "}\n";

        #pragma endregion

        #pragma region Close file

        outfile.flush();
        outfile.close();

        #pragma endregion
    }

    #pragma endregion

    #pragma region Error

    else
    {
        core->tools->output->PrintError("Unkown mode in Observation::WritePOV()");
    }

    #pragma endregion


}
