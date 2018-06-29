//
//  ModelCellsSphericalODE.h
//  CSCore
//
//  Created by celliere on 12/14/12.
//
//

#ifndef __CSCore__ModelCellsSphericalODE__
#define __CSCore__ModelCellsSphericalODE__

#include <iostream>
#include <vector>
#include "../../Cell/CellSphericalODE.h"
#include "../ModelCellsSpherical/ModelCellsSpherical.h"

// have to include this here because of typedefs
#include "sbmlsolver/odeSolver.h"
#include "cvode/cvode.h"
#include "sbmlsolver/odeSolver.h"

class CSParameter;
class QXmlStreamWriter;
class QXmlStreamReader;

class BoundingBoxList;


class ModelCellsSphericalODE : public ModelCellsSpherical
{
    
public :
    
    ModelCellsSphericalODE();

    virtual void SetCVODESettings(double Error, double RError, int Method, int useCompiled);
    virtual void SetSBMLModel(Model_t * sbmlwithlink);
    virtual void AddCell(double x, double y, double z);
    virtual void AddCell(double x, double y, double z, CellSphericalODE *motherCell);
    virtual void DivideCell(int i);
    virtual void UpdateCellFriction( CellSphericalODE * );
    virtual inline void UpdateSingleCellEffects();
    virtual void UpdateCellsStaining(int mode);

    virtual void GrowAndDivide();
    
    std::vector<CellSphericalODE *> cells;
    
private :
    
    cvodeSettings_t *options; // option for the integration
    Model_t *sbmlModel; // the sbml model with links to the cellsys parameters
    odeModel_t *odeModel; //the ode model created from the sbml model

    // to allow different timestep for ODE relative to physics
    double timeStepODEs ;
    double durationSinceLastODE ;

    //! for the trail examples
    double cellNumberThresholdForTrailApplication ;
    double trailDose ;
    bool trailNotSetYet ;
    int numberOfClonesToSeed ;
    double CvOfProtLevelBetweenClones ;
    double radiusOfMonolayer ;
    void seedDifferentClones () ;

};


inline
void
ModelCellsSphericalODE::UpdateSingleCellEffects()
{
    std::vector<CellSphericalODE *>::iterator cell;

    bool doODEs = durationSinceLastODE >= timeStepODEs ;

    if ( (int) time % 300 == 0 )
    {
        cout << "time = " << time << " ( timestep = " << timeStep << " , cell number = " << cells.size () << " )" << endl ;
    }

    for ( cell = cells.begin(); cell != cells.end(); ++cell )
    {
        //! for TRAIL example: death committed cell dont move
        if (!(*cell)->isDeathCommitted )
        {
            UpdateCellFriction(*cell);
            UpdateForcesLangevin(*cell);
        }
        if ( doODEs )
        {
            (*cell)->runODEsandUpdate(time,durationSinceLastODE) ; // calculate the state of ODEs and parameters for the next timepoint
        }
    }

    if ( doODEs )
    {
        cout << " at this step, I did compute ODEs !" << endl ;
        durationSinceLastODE = 0. ;
    }
    else
    {
        durationSinceLastODE += timeStep ;
    }

};




#endif /* defined(__CSCore__ModelCellsSphericalODE__) */
