
#include "LinearFunction.h"

LinearStrainFunction::LinearStrainFunction()
    : StrainFunction()
{
    //k.setXMLPath("k");
    //registerParameter(k);
}

/*LinearStrainFunction::LinearStrainFunction(XMLNode& functionNode,
                   std::stringstream& errors,
                   std::stringstream& warnings)
    : LinearStrainFunction()
{
    load(functionNode, errors, warnings);
}*/

LinearStrainFunction::LinearStrainFunction(const LinearStrainFunction& other)
    : k(other.k)
{}

LinearStrainFunction::LinearStrainFunction(LinearStrainFunction&& other)
    : k(other.k)
{}

LinearStrainFunction::~LinearStrainFunction()
{}

double LinearStrainFunction::operator()(const double strain) const
{
    return k * strain;
}

double LinearStrainFunction::stiffness()
{
    return k;
}

void LinearStrainFunction::print(std::ostream& stream) const
{
    stream << "LinearStrainFunction(k=" << k << ").";
}

bool LinearStrainFunction::ready() const
{
    return true;
}


