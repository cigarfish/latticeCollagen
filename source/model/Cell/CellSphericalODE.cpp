//
//  CellSphericalODE.cpp
//  CSCore
//
//  Created by celliere on 12/14/12.
//
//
#include "../../Core.h"
#include "CellSphericalODE.h"
#include "../Model/ModelCellsSphericalODE/ModelCellsSphericalODE.h"


CellSphericalODE::CellSphericalODE(double x, double y, double z) : CellSpherical(x,y,z)
{
    //! for the TRAIL example
    isDeathCommitted = false ;
}


int CellSphericalODE::runODEsandUpdate(double currentTime, double timeStep)
{
    //! TRAIL example : if cell committed to death, don't simulate
    if (isDeathCommitted) { return 1 ; }

    IntegratorInstance_setNextTimeStep(intInstance,currentTime+timeStep);
    //int canContinue = IntegratorInstance_integrateOneStep(intInstance);
    int canContinue = IntegratorInstance_integrateOneStepWithoutEventProcessing(intInstance);
    
    if (canContinue == 0)
    {
        core->tools->output->consoleModelCellsSpherical_Simulation << "ODE simulation failed";
        std::cout << "ODE simulation failed" << endl;
        return 0;
    }
    else
    {

        //! for the TRAIL example : look if cell commited to death
        if ( getcPARP () > 100000 )
        {
            isDeathCommitted = true ;
            timeOfDeathCommitment = currentTime + timeStep ;
        }

        //! treatments for the CCl4 example
        /*
        //retrieve concentration of the toxic compound
        concentrationOfEffector = IntegratorInstance_getVariableValue(intInstance,viEffector);
        //std::cout << CCl4 << " " << concentrationOfToxic << " " << IntegratorInstance_getVariableValue(intInstance,vis4) << " " << IntegratorInstance_getVariableValue(intInstance,viCYP2E1) << " " << IntegratorInstance_getVariableValue(intInstance,viCCl4inflow) << " " << IntegratorInstance_getVariableValue(intInstance,vikcat) << " " << IntegratorInstance_getVariableValue(intInstance,viKm) << endl;
        //std::cout << CCl4 << " " << concentrationOfToxic << " " << IntegratorInstance_getVariableValue(intInstance,viCCl4inflow) << endl;
        //set the CCl4inflow to the value calculated by flow simulation
        IntegratorInstance_setVariableValue(intInstance, viinflowFromOutside, inflowFromOutside);
        */

        return 1;
    }
}


//! definition of methods for the TRAIL example

void
CellSphericalODE::setTrail (double trail)
{
    IntegratorInstance_setVariableValue ( intInstance , viTrail , trail ) ;
}

double
CellSphericalODE::getcPARP ()
{
    return IntegratorInstance_getVariableValue ( intInstance , vicPARP ) ;
}



