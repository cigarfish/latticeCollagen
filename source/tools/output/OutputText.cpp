
#include "OutputText.h"
#include <iostream>
#include <fstream>

OutputText::OutputText()
{
    verboseLevel = 1;

	logfile.open("output/logfile.txt");
	logfile.clear();

	logfile << "OutputText::OutputText(). Init logfile.\n";
}

void OutputText::Print(std::string s, int verboseLevelOfMessage)
{
    if (verboseLevelOfMessage < 1)
    {
        PrintError("A message verbose level of smaller than 1 is not allowed.");
    }
    else if (verboseLevelOfMessage < 2)
    {
        PrintError("A message verbose level of 1 should only be used for Errors, Warnings, Init and Systems messages (Use the existing functions).");
    }
    else if (verboseLevelOfMessage <= verboseLevel)
    {
        std::cout << s;
    }
}

void OutputText::PrintInit(std::string s)
{
    if (1 <= verboseLevel)
    {
        std::cout << "Init: ";
    }
}

void OutputText::PrintSystem(std::string s)
{
    if (1 <= verboseLevel)
    {
        std::cout << "System: ";
    }
}

void OutputText::PrintError(std::string s)
{
    if (1 <= verboseLevel)
    {
        std::cout << "Error: " << s;
    }
}

void OutputText::PrintWarning(std::string s)
{
    if (1 <= verboseLevel)
    {
        std::cout << "Warning: " << s;
    }
}