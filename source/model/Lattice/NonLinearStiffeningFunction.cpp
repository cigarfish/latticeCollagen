
#include <cmath>

#include "NonLinearStiffeningFunction.h"

NonLinearStrainFunction::NonLinearStrainFunction()
    : StrainFunction()
{
    //k0.setXMLPath("k");
    //registerParameter(k0);

    //d0.setXMLPath("d0");
    //registerParameter(d0);

    //ls.setXMLPath("ls");
    //registerParameter(ls);

    //ds.setXMLPath("ds");
    //registerParameter(ds);
}

/*NonLinearStrainFunction::NonLinearStrainFunction(XMLNode& functionNode,
                   std::stringstream& errors,
                   std::stringstream& warnings)
    : NonLinearStrainFunction()
{
    load(functionNode, errors, warnings);
}*/

NonLinearStrainFunction::NonLinearStrainFunction(const NonLinearStrainFunction& other)
    : k0(other.k0),
    d0(other.d0),
    ls(other.ls),
    ds(other.ds)
{}

NonLinearStrainFunction::NonLinearStrainFunction(NonLinearStrainFunction&& other)
    : k0(other.k0),
    d0(other.d0),
    ls(other.ls),
    ds(other.ds)
{}

NonLinearStrainFunction::~NonLinearStrainFunction()
{}

double NonLinearStrainFunction::operator()(const double strain) const
{
    double ret;
    //if (strain < 0.) // disable for the moment
    //    ret = d0() * (std::exp(strain / d0()) - 1.); // TODO
    if (strain < ls)
        ret = strain;
    else // strain > ls()
        //ret = ls() + ds() * ExpMinusOne((strain - ls()) / ds());
        ret = ls + ds * (std::exp((strain - ls) / ds) - 1.);

    return ret * k0;
}

double NonLinearStrainFunction::stiffness()
{
    return k0;
}

void NonLinearStrainFunction::print(std::ostream& stream) const
{
    stream << "NonLinearStrainFunction("
        << "k0="  << k0
        << " ls=" << ls
        << " d0=" << d0
        << " ds=" << ds
        << ").";
}

bool NonLinearStrainFunction::ready() const
{
    return true;
}

double NonLinearStrainFunction::ExpMinusOne(const double x) const
{
    return x * (1. + 0.5 * x * (1. + 1. / 3. * x * (1. + 0.25 * x)));
}

