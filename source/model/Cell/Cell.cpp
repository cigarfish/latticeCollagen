
//#include "BasicDatatypes/Position.h"
#include "Cell.h"
#include "../../tools/model/CSModelTools.h"
#include "../../tools/random/Random.h"


Cell::Cell()
    : cellcycleState(0)
{
  this->mGeneration = 0;
}

void Cell::SetCycleTimeGaussClamped(double mean, double stddev)
{
	cycleTime = mpRandom->GetRandomGauss(mean, stddev);
	cycleTime = CSModelTools::ClampDouble(cycleTime, mean - 4 * stddev, mean + 4 * stddev); // Clamp values
}
