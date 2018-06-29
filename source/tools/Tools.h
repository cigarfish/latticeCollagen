
#ifndef TOOLS_H
#define TOOLS_H

#include "random/Random.h"
#include "output/OutputText.h"
#include "colors/ColorTools.h"

//! Collection of all toolbox like classes
class Tools
{
public:
    // Konstruktor
    Tools();

    //! Random number generation
    Random *random;

    //! Output class
    OutputText *output;

	//! Color transformations etc.
	ColorTools *color;

	// Prelim: Create math tools
	double Clamp(double value, double min, double max);

};

#endif
