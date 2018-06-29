//
//  ModelCellsSphericalODE.cpp
//  CSCore
//
//  Created by celliere on 12/14/12.
//
//

#include "ModelCellsSphericalODE.h"
#include "../../tools/model/BoundingBoxList.h"
#include "../../model/Cell/CellSphericalODE.h"

#include "../Interactions/CSInteractionFrictionMatrix.h"


ModelCellsSphericalODE::ModelCellsSphericalODE() : ModelCellsSpherical()
{
    options = new cvodeSettings_t;
    sbmlModel = new Model_t();
    
}

void ModelCellsSphericalODE::SetCVODESettings(double Error, double RError, int Method, int useCompiled)
{
    //(double Time, int PrintStep, double Error, double RError, int Mxstep, int Method, int IterMethod, int UseJacobian, int Indefinitely, int HaltOnEvent, int HaltOnSteadyState, int StoreResults, int Sensitivity, int SensMethod)
    options = CvodeSettings_createWith (10, 10, Error, RError, 100000, Method, 0, 0, 1, 0, 0, 0, 0, 0);

    /*std::cout << CvodeSettings_getEndTime(options) << endl;
        std::cout << CvodeSettings_getTimeStep(options) << endl;
        std::cout << CvodeSettings_getMethod(options) << endl;
        std::cout << CvodeSettings_getMaxOrder(options) << endl;
        std::cout << CvodeSettings_getIterMethod(options) << endl;
        std::cout << CvodeSettings_getJacobian(options) << endl;
        std::cout << CvodeSettings_getTimeStep(options) << endl;*/

    CvodeSettings_setTStop(options, 1); // tells the solver to not integrate farer than the next output time, so that the setValue is not missed
    CvodeSettings_setCompileFunctions(options, useCompiled); // 1=use compiled functions, 0 = don't use compiled functions

    //! for TRAIL example
    trailDose = 1000 ;
    numberOfClonesToSeed = 40 ;
    cellNumberThresholdForTrailApplication = 10 * numberOfClonesToSeed ;
    CvOfProtLevelBetweenClones = 0.30 ;
    radiusOfMonolayer = 60. ;
    timeStepODEs = 120.0 ;
    trailNotSetYet = true ;

    cout << "timeStepODEs = " << timeStepODEs << endl ;
    
}

void ModelCellsSphericalODE::SetSBMLModel(Model_t * sbmlwithlink)
{
    if (sbmlModel)
    {
        delete sbmlModel;
    }
    sbmlModel = new Model_t;
    sbmlModel = sbmlwithlink;

    odeModel = ODEModel_create(sbmlModel);
    
}


void
ModelCellsSphericalODE::seedDifferentClones ()
{
    cout << "Starting to seed different clones..." << endl ;

    for (int ic = 1 ; ic < numberOfClonesToSeed ; ic++ )
    {
        // choose x, y , z to occupy a uniformly distributed disk of radius radiusMonolayer
        double x , y ;
        cells.at (0)->mpRandom->GetRandomUnitVector (&x,&y) ;
        double cellPosRad = radiusOfMonolayer * sqrt ( cells.at (0)->mpRandom->GetRandomUniform01 () ) ;
        x *= cellPosRad ;
        y *= cellPosRad ;
        cout << "position for the clone chosen." << endl ;

        // add the clone
        AddCell (x,y,0.) ;
        CellSphericalODE* addedCell = cells.back () ;
        CellSphericalODE* referenceClone = cells.at (0) ;
        cout << "Cell added" << endl ;


        // clone-to-clone variability in the level of native proteins

        // names of variable to change
        std::vector <std::string> nativeNames ;
        nativeNames.push_back ("R") ;
        nativeNames.push_back ("flip") ;
        nativeNames.push_back ("pC8") ;
        nativeNames.push_back ("BAR") ;
        nativeNames.push_back ("pC3") ;
        nativeNames.push_back ("pC6") ;
        nativeNames.push_back ("pC9") ;
        nativeNames.push_back ("Bid") ;
        nativeNames.push_back ("XIAP") ;
        nativeNames.push_back ("PARP") ;
        nativeNames.push_back ("Smacm") ;
        nativeNames.push_back ("Bcl2c") ;
        nativeNames.push_back ("Bax") ;
        nativeNames.push_back ("Bcl2") ;
        nativeNames.push_back ("M") ;
        nativeNames.push_back ("Apaf") ;
        nativeNames.push_back ("CytoCm") ;

        // iterate on native prot names
        for ( unsigned int nn = 0 ; nn < nativeNames.size () ; nn++ )
        {
            cout << "setting value in clone for prot : " << nativeNames.at (nn) << endl ;
            // get the mean from the ref clone
            variableIndex_t* vi = ODEModel_getVariableIndex ( odeModel , nativeNames.at (nn).c_str () ) ;
            double meanFromRefClone = IntegratorInstance_getVariableValue ( referenceClone->intInstance , vi ) ;
            // log normal sampling
            double sampledCloneValue = addedCell->mpRandom->GetRandomLogNormal ( meanFromRefClone , CvOfProtLevelBetweenClones ) ;
            // put sampled value into cell
            IntegratorInstance_setVariableValue ( addedCell->intInstance , vi , sampledCloneValue ) ;
            cout << "Value Set !" << endl ;
        }

    }
}

void ModelCellsSphericalODE::AddCell(double x, double y, double z)
{
    CellSphericalODE * c = new CellSphericalODE(x,y,z);

    c->mpRandom = &mRandom;

    //ODEs
    //std::cout << "Adding the first cell" << endl;
    //CvodeSettings_dump(c->optionscell);
    //CvodeSettings_setTime(c->optionscell,10,10);
    //CvodeSettings_setIndefinitely(c->optionscell, 1);
    //CvodeSettings_dump(c->optionscell);
    //std::cout << CvodeSettings_getEndTime(c->optionscell) << endl;
    //std::cout << CvodeSettings_getTimeStep(c->optionscell) << endl;
    c->intInstance = IntegratorInstance_create(odeModel, options);


    //! for the TRAIL example
    c->viTrail = ODEModel_getVariableIndex ( odeModel , "L" ) ;
    c->vicPARP = ODEModel_getVariableIndex ( odeModel , "CPARP" ) ;
    c->setTrail (0.) ;


    //! for the CCL4 example
    /* c->viinflowFromOutside = ODEModel_getVariableIndex(odeModel, "CCl4inflow");
    c->viEffector          = ODEModel_getVariableIndex(odeModel, "s2"); */
    //c->viCCl4       = ODEModel_getVariableIndex(odeModel, "s1");
    //c->vikcat              = ODEModel_getVariableIndex(odeModel, "kcat");
    /*c->viCYP2E1     = ODEModel_getVariableIndex(c->odesbmlcell, "s3");
    c->vis4         = ODEModel_getVariableIndex(c->odesbmlcell, "s4");
    c->vikcat       = ODEModel_getVariableIndex(c->odesbmlcell, "kcat");
    c->viKm         = ODEModel_getVariableIndex(c->odesbmlcell, "Km");*/
    //cell-to-cell variability in kcat
    /* IntegratorInstance_setVariableValue(c->intInstance,ODEModel_getVariableIndex(odeModel, "kcat"),c->mpRandom->GetRandomGauss(1,0.5));
    c->concentrationOfEffector=0;
    c->inflowFromOutside=1; */


    // Initially, cells are non-quiescent
    c->setState( Cell::StateQuiescent, false );

    c->initialRadius = defaultInitialCellRadius;
    c->mRadius = defaultInitialCellRadius;
    c->divisionRadius = defaultDivisionCellRadius;
    c->SetCycleTimeGaussClamped(defaultCellCycleTime, defaultCellCycleTimeStandardDeviation);
    c->UpdateDeltaRadius(timeStep);
    c->youngModulus = defaultYoungModulusCells;
    c->poissonRatio = defaultPoissonRatioCells;
    c->diffusionCells = defaultDiffusionConstantCells;
    c->IndivAdhesionDensity = adhesionDensity;
    c->IndivGammaCells = gammaCellsParallel;
    c->IndivGammaECM = gammaECM;

    // Add newly created cell to population
    cells.push_back(c);
    cells2->add(c);
    
    mpArena->addObject(c->GLObject());

    if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << "ModelMonolayer::AddCell(). Created a new model cell at (" << x << " " << y << ") with radius " << c->mRadius << " dividing at " << c->divisionRadius << " with delta-radius " << c->deltaRadius << "\n";

}

// add a mew cell that has the same parameters and ODE state as the mother cell
void ModelCellsSphericalODE::AddCell(double x, double y, double z, CellSphericalODE *motherCell)
{

    CellSphericalODE * c = new CellSphericalODE(x,y,z);

    c->mpRandom = &mRandom;

    //ODEs
    //std::cout << "Adding the second cell" << endl;
    //CvodeSettings_dump(c->optionscell);
    //CvodeSettings_setIndefinitely(c->optionscell, 1);
    c->intInstance = IntegratorInstance_create(odeModel, options);

    //    IntegratorInstance_dumpNames(c->intInstance);
    //    IntegratorInstance_dumpData(c->intInstance);
    //    IntegratorInstance_dumpData(motherCell->intInstance);

    // copy the state or mother cell into daughter. !!!! the time is not copied !!!! this is the reason for the two next lines
    IntegratorInstance_copyVariableState ( c->intInstance , motherCell->intInstance ) ;

    //    IntegratorInstance_dumpData(c->intInstance);
    //    exit (2) ;

    // this is necessary for setInitialTime to actually work:
    IntegratorInstance_setNextTimeStep(c->intInstance, time+timeStep);
    IntegratorInstance_setInitialTime(c->intInstance,time);
    
    //! for the TRAIL example
    c->viTrail = ODEModel_getVariableIndex ( odeModel , "L" ) ;
    c->vicPARP = ODEModel_getVariableIndex ( odeModel , "CPARP" ) ;

    cout << "TRAIL in mother = " << IntegratorInstance_getVariableValue ( motherCell->intInstance , motherCell->viTrail ) << endl ;
    cout << "TRAIL in daughter = " << IntegratorInstance_getVariableValue ( c->intInstance , c->viTrail ) << endl ;

    cout << "cPARP in mother = " << IntegratorInstance_getVariableValue ( motherCell->intInstance , motherCell->vicPARP ) << endl ;
    cout << "cPARP in daughter = " << IntegratorInstance_getVariableValue ( c->intInstance , c->vicPARP ) << endl ;

    cout << "C8 in mother = " << IntegratorInstance_getVariableValue ( motherCell->intInstance , ODEModel_getVariableIndex ( odeModel , "C8" ) ) << endl ;
    cout << "C8 in daughter = " << IntegratorInstance_getVariableValue ( c->intInstance , ODEModel_getVariableIndex ( odeModel , "C8" ) ) << endl ;


    //! for the CCL4 example
    /*
    c->viinflowFromOutside = ODEModel_getVariableIndex(odeModel, "CCl4inflow");
    c->viEffector          = ODEModel_getVariableIndex(odeModel, "s2");
    //c->viCCl4            = ODEModel_getVariableIndex(odeModel, "s1");
    //c->vikcat            = ODEModel_getVariableIndex(odeModel, "kcat");
    //cell-to-cell variability in kcat
    IntegratorInstance_setVariableValue(c->intInstance,ODEModel_getVariableIndex(odeModel, "kcat"),c->mpRandom->GetRandomGauss(1,0.5));
    IntegratorInstance_dumpData(c->intInstance);
    c->concentrationOfEffector= motherCell->concentrationOfEffector;
    c->inflowFromOutside=10;
    */


    c->initialRadius = motherCell->mRadius / 1.26;
    c->mRadius = motherCell->mRadius / 1.26;
    c->divisionRadius = motherCell->divisionRadius;
    c->SetCycleTimeGaussClamped(defaultCellCycleTime, defaultCellCycleTimeStandardDeviation);
    c->UpdateDeltaRadius(timeStep);
    c->youngModulus = motherCell->youngModulus;
    c->poissonRatio = motherCell->poissonRatio;
    c->diffusionCells = motherCell->diffusionCells ;
    c->IndivAdhesionDensity = motherCell->IndivAdhesionDensity;
    c->IndivGammaCells = motherCell->IndivGammaCells;
    c->IndivGammaECM = motherCell->IndivGammaECM;

        // Initially, cells are non-quiescent
    c->setState( Cell::StateQuiescent, false );

    // Add newly created cell to population
    cells.push_back(c);
    cells2->add(c);
    
    mpArena->addObject(c->GLObject());

    
    if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << "ModelMonolayer::AddCell(). Created a new model cell at (" << x << " " << y << ") with radius " << c->mRadius << " dividing at " << c->divisionRadius << " with delta-radius " << c->deltaRadius << "\n";

}


void ModelCellsSphericalODE::UpdateCellFriction( CellSphericalODE * cell )
{
    // Add cell-cell friction
    cell->frictionCoefficient += cell->IndivGammaCells * cell->surfaceContactArea;

    // Set surface area to the remaining surface area (not in contact with other cells)
    // Full surface of (spherical) cell is 4*pi*R^2
    cell->surfaceContactArea = 4 * M_PI * cell->mRadius * cell->mRadius - cell->surfaceContactArea;

    // Add cell-ecm friction
    if (cell->surfaceContactArea > 0) // If there is still some free surface available
    {
        cell->frictionCoefficient += cell->IndivGammaECM * cell->surfaceContactArea;
    }
}


void ModelCellsSphericalODE::DivideCell(int i)
{
    // If division radius was reached
    if ((*cells.at(i)).mRadius >= (*cells.at(i)).divisionRadius)
    {
#pragma region Debug output
        
        if (core->tools->output->debugMonolayer) (*core->tools->output).logfile << "Divide cell ";
        
#pragma endregion
        
#pragma region Preparations
        
        // Spatial shift of both daughter cells from mother cell
        Vector3f delta;
        
        // Obtain random direction of cell division
        if (is2D) mRandom.GetRandomUnitVector(&delta.x, &delta.y);
        else mRandom.GetRandomUnitVector(&delta.x, &delta.y, &delta.z);
        
        // Rescale to division distance
        delta.x *= 0.9 * cells.at(i)->divisionRadius / 1.26; //80% of the radius of the dauthers cells
        delta.y *= 0.9 * cells.at(i)->divisionRadius / 1.26;
        if (is2D == false) delta.z *= 0.9 * cells.at(i)->divisionRadius / 1.26;
        
#pragma endregion
        
#pragma region Create daughter cell 1
        
        if (is2D) AddCell(cells.at(i)->position.x + delta.x, cells.at(i)->position.y + delta.y, 0, cells.at(i));
        else  AddCell(cells.at(i)->position.x + delta.x, cells.at(i)->position.y + delta.y, cells.at(i)->position.z + delta.z, cells.at(i));
        
#pragma endregion
        
#pragma region Create daughter cell 2 (by modifying the existing cell)
        
        // Reset radius
        cells.at(i)->mRadius /= 1.26;
        cells.at(i)->initialRadius = cells.at(i)->mRadius;
        
        // Update position
        cells.at(i)->position.x -= delta.x;
        cells.at(i)->position.y -= delta.y;
        if (is2D == false) cells.at(i)->position.z -= delta.z;
        
        // Recaclulate cell cycle time
        //cells.at(i)->SetCycleTimeGaussClamped(biolink->scaleBiologyToInternal(IntegratorInstance_getVariableValue(cells.at(i)->intInstance, cells.at(i)->viccTime), BiologyLink::ScaleTime), biolink->scaleBiologyToInternal(IntegratorInstance_getVariableValue(cells.at(i)->intInstance, cells.at(i)->viccSD), BiologyLink::ScaleTime));
        cells.at(i)->SetCycleTimeGaussClamped(defaultCellCycleTime, defaultCellCycleTimeStandardDeviation);
        // Update corresponding growth step
        cells.at(i)->UpdateDeltaRadius(timeStep);
        
#pragma endregion
        
#pragma region Prelim: State transition (Growing -> Quiescent)
        
        if ((*cells.at(i)).lastForceAbsolute > 1000)
        {
            //cells.at(i)->cellcycleState[CELLCYCLE_STATE_QUIESCENT] = true;
            cells.at(i)->setState( Cell::StateQuiescent );
            // cells.at(cells.size()-1)->cellcycleState[CELLCYCLE_STATE_QUIESCENT] = true;
            cells.at(cells.size()-1)->setState( Cell::StateQuiescent );

        }
        
#pragma endregion
    }
}


// with the overriding of ModelCellsSpherical::cells by ModelCellsSphericalODE::cells
// it's now still necessary to duplicate this method
// Grows all cells (using the given time delta) and applies cell divisions
void ModelCellsSphericalODE::GrowAndDivide()
{

    //! for TRAIL example
    if ( time == 0. ) { seedDifferentClones () ; }
    bool doTrail = cells.size () > cellNumberThresholdForTrailApplication && trailNotSetYet ;
    if (doTrail) { trailNotSetYet = false ; }

    // For all cells
    for (unsigned int i = 0; i< cells.size(); i++)
    {

        //! for TRAIL example: just display cPARP level
        //        cout << "cell " << i << " : cPARP = " << cells[i]->getcPARP () << endl ;

        //! for TRAIL example : if cell number > X , apply TRAIL
        if ( doTrail )
        {
            cells[i]->setTrail (trailDose) ;
        }

        // If cell is not quiescent
        if ( cells[i]->getState( Cell::StateQuiescent ) == false )
        {
            //! additional condition from TRAIL example: should not be committed to death
            if ( !cells[i]->isDeathCommitted )

            {

#pragma region Debug output

                if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << "Cell " << i << " Grow(). Radius: " << (*cells.at(i)).mRadius << " -> ";

#pragma endregion

                // Grow the cell
                (*cells.at(i)).Grow();

                // Divide the cell if possible
                if ((*cells.at(i)).CanDivide())
                {
                    std::cout << "divide cell" << endl;
                    DivideCell(i);
                }

#pragma region Debug output

                if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << (*cells.at(i)).mRadius << "\n";

#pragma endregion

            }

        }
    }
}

void ModelCellsSphericalODE::UpdateCellsStaining(int mode)
{
    if (mode == 0) // all white
    {
        for (unsigned int i=0; i<cells.size(); i++)
        {
            cells.at(i)->SetColor(1, 1, 1);
        }
    }
    else if (mode == 1) // cell state quiescent or dividing
    {
        for (unsigned int i=0; i<cells.size(); i++)
        {
            if (cells.at(i)->getState(Cell::StateQuiescent))
            {
                cells.at(i)->SetColor(0.5, 0.5, 0.5);
            }
            else
            {
                cells.at(i)->SetColor(1, 1, 1);
            }
        }
    }
    else if (mode == 2) // last absolute force
    {
        for (unsigned int i=0; i<cells.size(); i++)
        {
            (*core->tools).color->CreateSuperTrafficHue(cells.at(i)->lastForceAbsolute, 0., 1000.);
            cells.at(i)->SetColor( (*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b );
        }
    }
    
    else if (mode == 3) // by volume
    {
        for (unsigned int i=0; i<cells.size(); i++)
        {
            // Prelim: Later: Actual volume calculation
            (*core->tools).color->CreateTrafficHue(cells.at(i)->mRadius, defaultInitialCellRadius, defaultDivisionCellRadius);
            cells.at(i)->SetColor( (*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b );
        }
    }
    
    else if (mode == 4) //internal concentration
    {
        for (unsigned int i=0; i<cells.size(); i++)
        {
            // Prelim: Later: Actual volume calculation

            //! for trail example
            double valueForColor = ( cells.at(i)->getcPARP() > 100000 ) ? 0. : 100000 - cells.at(i)->getcPARP() ;
            (*core->tools).color->CreateTrafficHue( valueForColor , 0. , 1. ) ;

            //! for CCL4 example
            //(*core->tools).color->CreateTrafficHue(cells.at(i)->concentrationOfToxic, 0, 10000);

            cells.at(i)->SetColor( (*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b );
        }
    }
    
#pragma endregion
    
    // Remember colormode
    lastColormodeCells= mode;
}
