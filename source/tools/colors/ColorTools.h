#ifndef COLORTOOLS_H
#define COLORTOOLS_H

class ColorTools
{
public:

    ColorTools();

	// Color values set by methods below
	float r,g,b;

	// Creates a Red -> yellow -> green hue from a linear value (black = min)
	void CreateTrafficHue(float value, float min, float max);

	// Creates a black -> red -> yellow -> green -> cyan -> blue -> magenta -> white hue from a linear value (black = min)
    void CreateSuperTrafficHue(float value, float min, float max);

};

#endif