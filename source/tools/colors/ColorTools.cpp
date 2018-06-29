
#include "ColorTools.h"

ColorTools::ColorTools()
{}

#pragma region TrafficHue: Red -> yellow -> green (red = min)

void ColorTools::CreateTrafficHue(float value, float min, float max)
{
    if (min==max)
    {
        r = 1.;
        b = 1.;
        g = 1.;
        return;
    }
    // min = red (100), middle = yellow (110), max = green (010)
    float middle = ((max-min)/2)+min;
    float halfrange = (max-min)/2;
    if (value<=middle)
    {
        r = 1.0;
        g = (value-min)/halfrange;
        b = 0.0;
    }
    else
    {
        r = (max-value)/halfrange;
        g = 1.0;
        b = 0.0;
    }
    // catch overflows
    if (value>max)
    {
        r = 0.;
        g = 1.f;
        b = 0.; // temp
    }
    if (value<min)
    {
        r = 0.5;
        g = 0.;
        b = 0.;
    }
}

#pragma endregion

#pragma region SuperTrafficHue: Black -> red -> yellow -> green -> cyan -> blue -> magenta -> white hue from a linear value (black = min)

void ColorTools::CreateSuperTrafficHue(float value, float min, float max)
{
    if (min>max)
    {
        r = 0.5;
        b = 0.5;
        g = 0.5;
        return;
    }
    if (min==max)
    {
        r = 0.5;
        b = 0.5;
        g = 0.5;
        return;
    }

    //min - white/magenta - h1 - magentan/blue - h2 - blue/cyan - h3 - cyan/green - h4 - green/yellow - h5 - yellow/red - h6 - red/black - max
    float h1 = ((max-min)*(1./7.))+min;
    float h2 = ((max-min)*(2./7.))+min;
    float h3 = ((max-min)*(3./7.))+min;
    float h4 = ((max-min)*(4./7.))+min;
    float h5 = ((max-min)*(5./7.))+min;
    float h6 = ((max-min)*(6./7.))+min;
    float subrange = (max-min)/7.;
    //printf("min = [%f]  %f  %f  %f  [%f] (%f)\n", min, h1, h2, h3, max, subrange);
    if (value<=h1)
    {
        r = 1.0 ;
        g = 1.0 - (value-min)/subrange; // absteigend
        b = 1.0 ;
    }
    else if (value<=h2)
    {
        r = 1-(value-h1)/subrange;
        g = 0.0;
        b = 1.0;
    }
    else if (value<=h3)
    {
        r = 0.0;
        g = (value-h2)/subrange;
        b = 1.0;
    }
    else if (value<=h4)
    {
        r = 0.0;
        g = 1.0;
        b = 1.0 - (value-h3)/subrange;
    }
    else if (value<=h5)
    {
        r = (value-h4)/subrange;
        g = 1.0;
        b = 0.0;
    }
    else if (value<=h6)
    {
        r = 1.0;
        g = 1.0 - (value-h5)/subrange;
        b = 0.0;
    }
    else
    {
        r = 1.0 - (value-h6)/subrange;
        g = 0.0;
        b = 0.0;
    }

    // catch overflows
    if (value>max)
    {
        r = 0.;
        g = 0.;
        b = 0.;    // > max is black
    }
    if (value<min)
    {
        r = 1.;
        g = 1.;
        b = 1.;     // < min is white
    }
}

#pragma endregion

