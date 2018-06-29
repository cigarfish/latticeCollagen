#include "Data.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;
using namespace BasicDatatypes;

AverageData::AverageData(const string& name)
    : Data(name), avg(0.0), counter(0)
{}

AverageData::AverageData()
    : AverageData("")
{}

AverageData::AverageData(const AverageData& data)
    : Data(data), avg(data.avg), counter(data.counter), ids(data.ids)
{}

void AverageData::Reset()
{
    avg = 0.0;
    counter = 0;
    ids.clear();
}

void AverageData::Add(double value, pair<int, int> _ids)
{
    if ((find(ids.begin(), ids.end(), _ids) != ids.end()) ||
        find(ids.begin(), ids.end(), make_pair(_ids.second, _ids.first)) != ids.end())
    {
        cerr << "WARNING(" << getName() << "): Interaction(" <<
            _ids.first << ", " << _ids.second << ") already present." << endl;
        return;
    }

    Add(value);
    ids.push_back(_ids);
}

void AverageData::Add(double value)
{
    avg += value;
    counter++;
}

double AverageData::Average(void) const
{
    if (counter<=0.)
        return 0.;
    else
        return avg / counter;
}

double AverageData::operator()(void) const
{
    return Average();
}

StressData::StressData(double& v)
    : Data("Stress"), stress(), volume(v), empty(true)
{}

StressData::StressData(const StressData& data)
    : Data(data), stress(data.stress), volume(data.volume), empty(data.empty)
{}

void StressData::Reset()
{
    stress.set(0);
    ids.clear();
    empty = true;
}

void StressData::Add(const Vector3f& force, const Vector3f& distance)
{
    stress += outer(force, distance);
    empty = false;
}

void StressData::Add(const Vector3f& force, const Vector3f& distance,
                     pair<int, int> _ids)
{
    if ((find(ids.begin(), ids.end(), _ids) != ids.end()) ||
        find(ids.begin(), ids.end(), make_pair(_ids.second, _ids.first)) != ids.end())
    {
        cerr << "WARNING(" << getName() << "): Interaction(" <<
            _ids.first << ", " << _ids.second << ") already present." << endl;
        return;
    }

    ids.push_back(_ids);
    Add(force, distance);
}

double StressData::Pressure(void) const
{
    if (empty) return 0.0;
    //double pressure =  1. / (3. * volume) * trace(stress);
    double pressure = 1.;
    assert(std::isfinite(pressure));
    return pressure;
}

double StressData::operator()(void) const
{
    return Pressure();
}
