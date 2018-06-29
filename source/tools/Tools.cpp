
#include "Tools.h"

Tools::Tools()
{
	// Init random numbers
	random = new Random();

	// Init Output
	output = new OutputText();

	// Color transformations etc.
	color = new ColorTools();

}

double Tools::Clamp(double value, double min, double max)
{
if (value < min) return min;
if (value > max) return max;
return value;
}
