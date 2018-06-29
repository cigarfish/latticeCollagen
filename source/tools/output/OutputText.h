#ifndef OUTPUTTEXT_H
#define OUTPUTTEXT_H

#include <string>
#include <sstream>
#include <fstream>

class OutputText
{
private:

    //! Level of verbosity (0: No output, 1: Only Errors, Warnings and System messages, 2 and 3 ... Levels for generic messages)
    int verboseLevel;

public:

#pragma region List of global debug switches used to control logfile output

	static const bool debugMonolayer = false; // Enable textual debug output when simulating a monolayer model

#pragma endregion

	//! Preliminary stream for logfile
	std::ofstream logfile;

	//! Stream that ends up on console of model wth spherical cells (Tab Simulation)
	std::ofstream consoleModelCellsSpherical_Simulation;

	//! Stream that ends up on console of model wth spherical cells (Tab Quantification)
	std::ofstream consoleModelCellsSpherical_Quantification;

    //! Method to generate generic output (Displayed only if verboseLevelOfMessage <= verboseLevel)
    void Print(std::string s, int verboseLevelOfMessage);

    //! Method to generate output when initializing something (auto verbose level 1)
    void PrintInit(std::string s);

    //! Method to generate output when initializing something (auto verbose level 1)
    void PrintSystem(std::string s);

    //! Method to generate output for a warning (auto verbose level 1)
    void PrintWarning(std::string s);

    //! Method to generate output for an error (auto verbose level 1)
    void PrintError(std::string s);

    OutputText();
};

#endif