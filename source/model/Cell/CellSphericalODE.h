//
//  CellSphericalODE.h
//  CSCore
//
//  Created by celliere on 12/14/12.
//
//

#ifndef __CSCore__CellSphericalODE__
#define __CSCore__CellSphericalODE__

#include <iostream>
#include "CellSpherical.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "sbmlsolver/odeSolver.h"

class ModelCellsSphericalODE;


class CellSphericalODE : public CellSpherical
{

    public :
    
    CellSphericalODE(double x, double y, double z);
    int runODEsandUpdate(double currentTime, double timeStep);
    
    // parameters that are defined only in the model and not for each cell in stefans model
    double diffusionCells;
    double IndivAdhesionDensity;
    double IndivGammaCells;
    double IndivGammaECM;
    
    /*variableIndex_t *viCellDiameter;
    variableIndex_t *viccTime;
    variableIndex_t *viccSD;
    variableIndex_t *viyoungmod;
    variableIndex_t *vipoisson;
    variableIndex_t *vidiffu;
    variableIndex_t *viccAdh;
    variableIndex_t *viccgamma;
    variableIndex_t *vicegamma;*/
    
    //! variables for the TRAIL example
    variableIndex_t * viTrail ;
    variableIndex_t * vicPARP ;
    bool isDeathCommitted ;
    double timeOfDeathCommitment ;

    //! methods for the TRAIL example
    void setTrail ( double trail ) ;
    double getcPARP () ;


    //! variables for the CCl4 example
    /*double concentrationOfEffector;
    double inflowFromOutside;
    double CCl4;
    variableIndex_t *viEffector;
    variableIndex_t *viinflowFromOutside;
    variableIndex_t *viCCl4;
    variableIndex_t *viCYP2E1;
    variableIndex_t *vis4;
    variableIndex_t *vikcat;
    variableIndex_t *viKm;*/
    
    integratorInstance_t *intInstance;
    


};

#endif /* defined(__CSCore__CellSphericalODE__) */
