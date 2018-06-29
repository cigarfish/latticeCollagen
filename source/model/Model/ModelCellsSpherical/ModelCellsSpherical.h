#ifndef MODEL_CELLS_SPHERICAL_H
#define MODEL_CELLS_SPHERICAL_H

#include "../../model/Model/CSModel.h"
#include "../../tools/Tools.h"
#include "../../gui/CSGLArena.h"
#include "../Cell/CellSpherical.h"
#include "../Cell/CellSphericalPolar.h"
#include "../Elements/ModelElementBarrierTriangle.h"
#include "../BiologyLink.h"
#include "../../Observation/Observation.h"
#include "../../../tools/random/Random.h"
#include "../../../tools/dataIO/vtp/CSVTPWriter.h"

#include "../../BasicDatatypes/GraphSphere.h"

// added by Jieling
#include "../../Elements/ModelElementECMSphere.h"
#include "model/Lattice/ModelLattice.h"

// #include "../../../gui/tabMonolayer/monolayer.h"
#include "../../../Core.h"

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <functional>

class CSParameterContext;
class QXmlStreamWriter;
class QXmlStreamReader;
class CSInteractionHertz;
class CSInteractionJKR;
class CSInteractionFrictionMatrix;
class ModelElementBarrierTriangle;
class BoundingBoxList;
class CSVTPWriter;
class LobuleLesionVoxelization;


//! Model class for 2D monolayer models
class ModelCellsSpherical : public CSModel
{
    // preliminary:  declare Cell classes as friends to access parameters
    // to be replaced by cell-local parameters.
    friend class CellSpherical;
    friend class CellSphericalPolar;

    friend class LobuleLesionVoxelization;
    friend class monolayer;

#pragma region Parameters

public:
    //! The model type.
    static const std::string xmlType;

    enum ContactModel { ContactModelHertz, ContactModelJKR };
    enum LobuleShape { None, Hexagonal, Quadric };

    ContactModel mContactModel;
    LobuleShape mLobuleShape;

    static const std::string mContactModelNames[];
    static const std::string mLobuleShapeNames[];

    //! The individual model instantiation name (in case of multiple instances):
    std::string xmlName;

    //! Change the name
    // ...and the parameter context's name!
    void SetName( std::string newName )
    {
        if ( newName != name )
        {
            core->models[newName] = this;
            if ( name.size() )
                core->models.erase( core->models.find(name) );
            xmlName = newName;
            CSModel::SetName(newName);
        }
    };

    enum Scenario {
        ScenarioReadFromData=0,
        ScenarioSingleCell,
        ScenarioEmbeddingMedium,
        ScenarioSpherePacking,
        ScenarioLobuleRegeneration
    };

    Scenario mScenario;

    static const std::string mScenarioNames[];

    void SetScenario(Scenario scen) {mScenario=scen;};

    //! Included default cell parameters with biological units and methods for unit conversions, is used to get dimensionless default parameters for model
    BiologyLink *biolink;

    // BioLink -> (default) parameters in Model class
    void UpdateParametersFromBioLink();

	// BioLink -> Dimensionless (default) parameters in Model class
	void UpdateDimensionlessParametersFromBioLink();

    // Push all dimensionless (default) parameters in Model class -> cell
    void UpdateParametersForCell( CellSpherical * );

	// Push all dimensionless (default) parameters in Model class -> All existing cells
	void UpdateParametersForAllCells();

#pragma endregion

    //! Model measurements
    Observation *observe;

    #pragma region Submodels

protected:
    bool is2D; // Can 2D optimizations be enabled (For strict monolayers)

    //border for lobule shape
    double mlobule_radius;
    double mlobule_height;

    bool mBloodVesselNetwork;
    std::string mBloodVesselNetworkPath;
    GraphSphere * mpGraphBloodVesselNetwork;

    bool mReadCells;

    //! Container used when a set of cells are to be killed.
    std::vector<CellSpherical*> mCellsToDelete;

    //! Parameter to change the radius of lesion after intoxication - Scenario Lobule Regeneration
    double mKillZoneRadius;

    double mKillCellsToDiePerInterval;
    //! Time interval for when to remove cells
    double mKillInterval;
    //! Next absolute time for when cells are removed, incremented by mKillInterval
    double mNextKillTime;
    //! The time period over which chosen cells will die. This period will be
    //! divided into mKillInterval large intervals.
    double mKillPeriod;
    //! Absolute time when the killing period will end.  After this time, all
    //! remaining cells will be removed.
    double mKillEnd;

    CSListContainer<ModelElementBarrierTriangle *> mBarriers;

    // t1m-debug for experimenting with scaling the S phase data to proliferation rates
    double mScaleProliferationRates;
    //!t1m-debug
    double mProliferationUpdateInterval;
    double nextObservablesUpdate;
    double nextProliferationRateUpdate;

    double nextLesionDiscretization;

    //! Additional active motion term toward the central vein (lobule model).
    //!  This force is only applied to cells at the edge of a lesion.
    double mCenterPullForce;

    double mParmCenterPullForce;

	// added by Jieling
	double mECMPullForce;
	int blockNumber;
	double blockSize;
	double sStrain;
	bool releaseStress; // true: change strain completely; false: change strain gradually

    //! Data structure to represent experimental data for proliferation patterns in CCl4 treated liver lobules.
    //  The data is taken
    //  - at specific time points
    //  - for specific layers (one of 14) within the lobule
    //  - gives rates in cells in S-phase per 24 hours
    // since the number of layers in the data is constant, the rates are saved in an array of 14 elements.
    static std::deque< std::pair<double /*time*/, std::vector<double /*rate by layer*/> > > mvProliferationData;

    std::deque< std::pair<double /*time*/, std::vector<double /*rate by layer*/> > > mvProliferationProfiles;

    double mpProliferationRatesPerLayer[14];
    double mpProbabilityToProliferateInLayer[14];

    LobuleLesionVoxelization * mpLesionVoxelization;

public:
    void SwitchSubmodel(int index);

    void CalcRadiusGyration();

    // Method to manually trigger voxelization from GUI
    void voxelize();

    #pragma endregion

    //! Default radius of the cells at t=0 and right after cell division
    double defaultInitialCellRadius;

protected:

    #pragma region Growth and Division

    //! Mean cell cycle time
    double defaultCellCycleTime;

    //! Standard deviation of cell cycle time (Gaussian distribution)
    double defaultCellCycleTimeStandardDeviation;

    //! Default radius beyond which the cells divide into two daughter cells
    double defaultDivisionCellRadius;

    //! Default distance of each of the two daughter cells center to the center of the mother cell
    double defaultDivisionDistance;

    //! Grows all cells and applies cell divisions (assuming one timestep)
    virtual void GrowAndDivide();

    //! Divide cell i
    // Modifies position and radis of cell i and creates an additional cell
    virtual void DivideCell(int i);

    //! Adds and initializes a new cell
    virtual void AddPolarCell(double x, double y, double z);
    virtual void AddCell(double x, double y, double z);


    //! Adds an existing cell
    virtual void AddCell( CellSpherical * newCell );

    //! Removes a cell from all containers.
    void RemoveCell( CellSpherical * );

    //Add Cells from MXF file
    void AddReadCells( std::string & filename );

    //Add shape of Lobule
    void AddBarrierTriangle( double mPoints[][3], bool visible=1, double epsilon=0.1, double r=1, double g=1, double b=0, double a=0.5);

    void AddBloodVesselNetwork( GraphSphere * );
    void AddBloodVesselNetwork( std::string & filename , int filetype=2 );

    void AddShape(double radius=10, double height=10, double r=1, double g=1, double b=0, double a=0.5, double epsilon=0.1, bool visible=1);
    void AddLobuleShape(double radius=10, double height=10, double r=1, double g=1, double b=0, double a=0.5, double epsilon=0.1, bool visible=1);
    void AddQuadricShape(double radius=10, double height=10, double r=1, double g=1, double b=0, double a=0.5, double epsilon=0.1, bool visible=1);


    #pragma endregion

    #pragma region Biophysics

    #pragma region Cell parameters (Dimensionless)

    double defaultYoungModulusSinusoids;      //! Default young modulus of sinusoids (dimensionless)
    double defaultPoissonRatioSinusoids;      //! Default poisson radtio of sinusoids

    double defaultYoungModulusCells;          //! Default young modulus of cells (dimensionless)
    double defaultPoissonRatioCells;          //! Default poisson radtio of cells
    double defaultDiffusionConstantCells;     //! Default diffusion constant of cells

	// added by Jieling
	double defaultYoungModulusECM;	// Default young modulus of ECM (dimensionless)
	double defaultPoissonRatioECM;	// Default poisson ratio of ECM

    //! Single bond energy (dimensionless)
    double singleBondEnergy;

    //! Default cell-cell adhesion density
    double adhesionDensity;

    //! Max angle of rotation for metropolis algorithm in case of Polar Cells
    double mMaxRotationAngleForMetropolis;

    //! Friction coefficient with other cells
    double gammaCells;

    double gammaCellsParallel;
    double gammaCellsPerpendicular;

    //! Fricition coefficient with ECM
    double gammaECM;
	// added by Jieling
	double gammaECMLattice;

    //pressure for quiescent
    double mQuiescentMin;
    double mQuiescentMax;

    // parameters for ScenarioEmbeddingmedium
    // Radius of edge of initial population in ScenarioEmbeddingmedium
    double mPopulationRadius;
    // distance between cells in initial configuration in Scenarioembeddingmedium
    double mPopulationInitialDistance;

    #pragma endregion

    #pragma region Force calculation

    //! Resets directional forces
    void InitForces();


    // Possible CSInteractions used in this model
    //  HertzForce interaction -> CSInteractionHertz
    //  ToDo:  JKR interaction
    //  FrictionMatrix from Hertz or JKR interactions -> CSInteractionFrictionMatrixHertz
    //  ...
    // These CSInteraction objects are not allocated in the constructor,
    //  but should be initialized in SetupSimulation depending on the options
    //  set by the user in the GUI or data bundle (xml).
    CSInteractionHertz          * mpInteractionHertz;
    CSInteractionJKR            * mpInteractionJKR;
    CSInteractionFrictionMatrix * mpInteractionFrictionMatrix;

    //! Loop over all cell-cell interactions.
    //! Calling CSInteractionHertz(.,.) and
    //! CSInteractionFrictionMatrixHertz(.,.)
    //! on every interacting Cell-Cell pair.
    virtual void UpdateInteractions();
//    void UpdateInteractionsLattice();
    void UpdateInteractionsCellsOnly();
    void UpdateInteractionsTest();
    void UpdateInteractionsTest2();
    void UpdateInteractionsOverlapTest();

    //! Loop over all cells.
    //! Calling UpdateCellFriction(.) and
    //! UpdateForcesLangevin(.) on every cell in
    //! std::vector<CellSpherical*> cells
    virtual void UpdateSingleCellEffects();
	// added by Jieling
	virtual void UpdateSingleECMEffects();


    //! Prepares cell friction based on accumulated contact areas
    virtual void UpdateCellFriction( CellSpherical * );
	// added by Jieling
	void UpdateLatticeNodeFriction( ModelElementLatticeNode * );

    //! Adds random Langevin-based forces
    virtual void UpdateForcesLangevin( CellSpherical * );

    virtual void AddDirectedMotion( CellSpherical * );

    //! Method to rescale the contribution of the Langevin force when adapting
    //! time step size.  Since the Langevin force scales with sqrt(1/timeStep).
    //! \param timeStepFactor  The factor by which the timeStep changed.
    void rescaleLangevinContrib( double timeStepFactor );

    #pragma endregion

    #pragma endregion


    #pragma region Simulation

    #pragma region Parameters for simulation control

public:

    // Remember start of simulation for progress bar
    double simulateFromDays;

    // Target time of current simulation (in days)
    double simulateUntilDays;

    // Start time w.r.t. t=0.
    double mStartTime;

    // Time between measurements (observations) in current simulation (in days)
    double observeEveryDays;

    bool enableObservation;

    // Next time of measurement / observation in current simulation (in days)
    double nextObservationTime;

    //! Flag for deciding to use CellSphericalPolar instead of CellSpherical.
    bool mUsePolarCells;

    //! Parameter to calculate the pole region angle of CellSphericalPolars
    //  from.  Unit is '%', i.e. a factor of 0.01 is to be used.  The angle of
    //  the polar region is calculated by
    //    mAngle = acos( 1 - .01*mParmAdhesiveSurface );
    double mParmAdhesiveSurface;

    //! Flag for deciding to use dumbbell or instantaneous division.
    bool mUseDumbbell;
    double mDefaultDumbbellPeriod;

    //! Flag for choosing dynamical step size adaptation per maximum allowed
    //! displacement for each time step.
    bool mUseDynamicalTimeSteps;
    double mDynamicTimeStepMin;//min and max Time step for dynamic timesteps
    double mDynamicTimeStepMax;

    //! Scaling factor for dynamical timeSteps.
    unsigned int mTimeScalingFactor;

    //! A variable to determine the maximum velocity in solveSystem() for time
    //! step adaptation.
    double mVelocityMaxSquared;

    //! A variable to hold the square of the maximum allowed displacement:
    //  (0.1 * defaultInitialCellRadius)^2.  Constant during simulation.
    double mMaximumDisplacementSquared;
    double mMaximumDisplacementCellRadius; //max allowed displacement

protected:

    #pragma endregion

    // One simulation step. t -> t+1
    virtual void Simulate();

    //! Solve the equation of motion
    virtual void solveSystem();

    virtual void FrictionMatrixTimesV( double *v_in, double *v_out );
    virtual void InitialResidualAndPreconditioner( double * v_in, double * residual_out, double *preconditioner_out );

    #pragma endregion

    //! Array for keeping track which cells have been handled in UpdateInteractions.
    //  This is only necessary for the JKR contact model, in order not to account for contacts twice.
    bool * mpJKRElementDone;
    unsigned long int mJKRElementDoneSize;

    //! Vector holding the velocities (for solving the equation of motion)
    //  size = 3*cells.size()
    double * mpVelocities;
    unsigned long int mProblemAllocationSize;

    // vectors for solving the overdamped equation of motion with friction
    // i.e. if mpInteractionFrictionMatrix != NULL;
    // declared as members, so that reallocation has to happen rarely
    double * mpTmpVectorP, * mpTmpVectorQ;
    double * mpResidualVector;
    double * mpPreconditioner;

    // storage for the directed forces
    double * mpDirectedForcesVector;

    CSParameterContext * mpParameters;
    CSParameterContext * mpParametersVisualization;

    virtual void RegisterParameters();

    struct observables {
        unsigned long numCells; // cells.size() might differ since dumbbell cells are 2 CellSphericals
        unsigned long numCellsPerLayer[14];
        unsigned long numSPhaseCells;
        unsigned long numSPhaseCellsPerLayer[14];
        unsigned long numProliferatingCells;
        unsigned long numProliferatingCellsPerLayer[14];
        unsigned long numCellsLesion;
        double radiusLesion;
        double areaDensity;
        double lesionArea;
        double areaDensityVoxelized;
        double areaDensityVoxelizedWithLesion;
    };
    observables mObservables;

    inline void ResetObservables();
    inline void UpdateObservables();
    inline void WriteoutObservables();

    #pragma region Scenarios
    void InitScenarioLobuleRegeneration( bool createInitialConditions );
    #pragma endregion



public:
    //! Main cell population
    std::vector<CellSpherical *> cells;

	// added by Jieling
	std::vector<ModelElementECMSphere*> HSCs; // HSC cell to deposit ECM
	//ECMGraph *mpECMNetwork; // stored ECM spheres and edges
	Lattice *mpCollagenNetwork; // stored Collagen nodes and springs
	double mECMNetworkTime = 172800; // 1 hours after injury
	void InitHSCsDistribution();
	void stressForStrainTest();
	void printStressStrain();
	bool stableStatus;
	bool unLoaded;
	int LoadedStep;
	int unLoadedStep;
	double getLatticeStrain();

    // boundary elements
    std::vector<ModelElementBarrierTriangle *> barrier;

    //! Alternative cell population:
    BoundingBoxList * cells2;

    // Default constructor
    ModelCellsSpherical();

    #pragma region Model interface (Reset, Simulation control)

    // Preliminary init
	virtual void Reset(bool loadScenario=true);

    // Returns the progress of the simulation (Normalized to 0...1)
    float GetSimulationProgress();

    // virtual method from class model (to deprecate InitSimulateInThread).
    // Sets parameters for simulation control
    virtual void SetupSimulation();

    // Method to do multiple simulation steps
    void Simulate(int numberOfSimulationSteps);

    // Method for simulation from GUI
    virtual void SimulateInThread();

    #pragma endregion

    #pragma region Cell population staining

	//! Stores the last color mode that was called from the gui
	int lastColormodeCells;

	//! Update cells colors using that last color mode
	virtual void UpdateCellsStaining();

	//! Update cell color using a new color mode
    virtual void UpdateCellsStaining(int mode);

    //! Update one cell's color using the last color mode;
    virtual void UpdateCellsStaining(CellSpherical *);

    void setVisible( bool visible);

    #pragma endregion

    #pragma region access to parameters

    virtual CSParameterContext * GetParameters( std::string contextName ="" );
    virtual void InitParameters( CSParameterContext * fromContext=0 );

    static void DefaultParameters( CSParameterContext * );

    #pragma region data I/O

    virtual void writeXML(QXmlStreamWriter *) const;

    static ModelCellsSpherical * createFromXML( QXmlStreamReader *xmlReader,
                                                std::stringstream & errors,
                                                std::stringstream & warnings );

    virtual void writeHDF5( H5::H5File * outputFile ) const;

    // method for reading the hdf5 data into an already allocated *Model
    // used after the createFromXML has created a certain Model and initialised
    // the parameters and, most importantly, its name.
    virtual void readModelData( H5::H5File * inputFile,
                        std::stringstream & errors,
                        std::stringstream & warnings );

    void writeXML2();

    //output class
    CSVTPWriter* mpParaview;
    std::string mOutputPath;
    std::string mOutputPrefix;

    size_t mOutputCounter;

    // preliminary vtp writing for ModelElementVesselSpheres:
    void writeVesselSphereVTP();

    #pragma endregion

    std::vector<double *> mFrictionMatrices;

    // The pointer to the array for all friction matrices created by the
    // cell-cell parallel/perpendicular friction calculation in
    // mpInteractionFrictionMatrix used by the conjugated gradient solver to
    // solve an over-damped equation of motion:
    // This should be 'connected' to the corresponding pointer in
    // mpInteractionFrictionMatrix:
    //   mpFrictionMatrices = mpInteractionFrictionMatrix->mpFrictionMatrices
    // after allocation of mpInteractionFrictionMatrix;
    double * mpFrictionMatrices;

    //global Radius of Gyration
    double mRadiusGyration;

    //adds cell in a sphere (closest sphere packing)
    void spherePacking(double radiusSphere);

};


inline
void
ModelCellsSpherical::UpdateSingleCellEffects()
{
    if ( mUsePolarCells )
    {
        std::function<unsigned long(unsigned long)>
            randFunc(std::bind(&Random::GetRandomUniform0x, &mRandom, std::placeholders::_1));
        std::random_shuffle(cells.begin(), cells.end(), randFunc );
    }

    for ( auto tmpCell: cells )
    {
        // CellSpherical* tmpCell = cells[existingCells[j]];
        tmpCell->lastForceAbsolute = tmpCell->accumulatedForceAbsolute;
        tmpCell->lastPressure = tmpCell->accumulatedPressure;

        UpdateCellFriction(tmpCell);
        if( this->defaultDiffusionConstantCells != 0)
            UpdateForcesLangevin(tmpCell);

        if ( (tmpCell)->mType == ModelElement::TypeCellSphericalPolar )
            static_cast<CellSphericalPolar *>(tmpCell)->OrientationMetropolis();

        if ( mScenario == ScenarioLobuleRegeneration )
        {
			if (time >= 0)
			{
				if ((tmpCell)->mLesionEdge)
				{
					Vector3f centerPull = (tmpCell)->position;
					centerPull[2] = 0.;
					double distance = centerPull.Normalize();
					centerPull.Multiply(-mCenterPullForce);
					(tmpCell)->directedForce.Add(centerPull);
				}
			}
        }

        // Applying forces on cells in dumbbell division phase onto the center of
        // mass of the two division partners.
        // This will be done twice unfortunately, the second time
        // c1Force==c2Force.  It will, though, remove the necessity to distribute
        // every single force contribution onto the division partners.
        if ( (tmpCell)->mpDivisionPartner)
        {
            Vector3f c1Force = (tmpCell)->directedForce;
            Vector3f c2Force = (tmpCell)->mpDivisionPartner->directedForce;

            (tmpCell)->directedForce = 0.5 * ( c1Force + c2Force );
            (tmpCell)->mpDivisionPartner->directedForce = (tmpCell)->directedForce;

            Vector3f c1Langevin = (tmpCell)->mLangevinForce;
            Vector3f c2Langevin = (tmpCell)->mpDivisionPartner->mLangevinForce;

            (tmpCell)->mLangevinForce = 0.5 * ( c1Langevin + c2Langevin );
            (tmpCell)->mpDivisionPartner->mLangevinForce = (tmpCell)->mLangevinForce;
        }
     //   i = i+1;
    }

};


inline
void
ModelCellsSpherical::ResetObservables()
{
    memset( &mObservables, 0, sizeof(observables) );
}


inline
void
ModelCellsSpherical::UpdateObservables()
{
    double backup = mObservables.lesionArea;
    ResetObservables();
    mObservables.lesionArea = backup;

    std::vector<CellSpherical *>::const_iterator cell;
    for ( auto cell: cells )
    {
        if ( mUseDumbbell )
        {
            if ( cell->mDaughterCell )
                continue;  // ignore dumbbell partner

            // cells.size() is accurate for non-dumbbell simulations
            ++mObservables.numCells;
        }

        if ( mScenario != ScenarioLobuleRegeneration &&
             !cell->getState(Cell::StateQuiescent) )
        {
            ++mObservables.numProliferatingCells;

            if ( cell->getState( Cell::StateSPhase ) )
                ++mObservables.numSPhaseCells;
        }
        else
        {
            unsigned int layer =
                (unsigned int) std::floor( cell->mLayerIndicator );
            if (layer > 13) layer = 13;

            ++mObservables.numCellsPerLayer[layer];

            if ( !cell->getState(Cell::StateQuiescent) )
            {
                ++mObservables.numProliferatingCells;
                ++mObservables.numProliferatingCellsPerLayer[layer];

                if ( cell->getState( Cell::StateSPhase ) )
                {
                    ++mObservables.numSPhaseCells;
                    ++mObservables.numSPhaseCellsPerLayer[layer];
                }
            }

            if ( cell->mLesionEdge )
            {
                Vector3f centerVector = cell->position;
                centerVector[2] = 0;
                mObservables.radiusLesion += centerVector.Norm() - cell->mRadius;
                ++mObservables.numCellsLesion;
            }
        }
    }

    if (!mUseDumbbell) mObservables.numCells = cells.size();

    if ( mScenario == ScenarioLobuleRegeneration )
    {
        if ( mObservables.numCellsLesion )
            mObservables.radiusLesion /= mObservables.numCellsLesion;
        else
            mObservables.radiusLesion = 0;

        // Area density:  1. divide by #levels in z direction
        mObservables.areaDensity
            = mObservables.numCells / (mlobule_height/(2*defaultInitialCellRadius));
        // Area density:  2. divide by area of the lobule in x-y-plane
        mObservables.areaDensity /= (mlobule_radius*mlobule_radius) * 2.598
            - M_PI * mObservables.radiusLesion*mObservables.radiusLesion;
        // Area is 6 * the area of an equilateral triangle of side length mlobule_radius:
        //  = 6 * mlobule_radius * h/2 = 3 *mlobule_radius * mlobule_radius * sin(Pi/3)
        //  = 3 * Sqrt(3/4) *mlobule_radius^2 = 2.59807621135 * mlobule_radius^2
        double lobuleRadiusInMeters = 1e-6 * biolink->getLengthInMicrometers(mlobule_radius);
        double lobuleHeightInMeters = 1e-6 * biolink->getLengthInMicrometers(mlobule_height);
        double initialCellRadiusInMeters = 1e-6 * biolink->getLengthInMicrometers(defaultInitialCellRadius);
        mObservables.areaDensityVoxelizedWithLesion = mObservables.numCells * 2 * initialCellRadiusInMeters / (lobuleHeightInMeters * lobuleRadiusInMeters*lobuleRadiusInMeters * 2.598);
        mObservables.areaDensityVoxelized = mObservables.numCells * 2 * initialCellRadiusInMeters / (lobuleHeightInMeters * ((lobuleRadiusInMeters*lobuleRadiusInMeters) * 2.598 - mObservables.lesionArea));
    }
}


inline
void
ModelCellsSpherical::WriteoutObservables()
{
    UpdateObservables();

    // ToDo:  open file
    if ( time == mStartTime )
        std::cout << "#"
                  << " Time"
                  << "\tnumCells"
                  << "\t"
                  << "\tnumCellsPerLayer:0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13"
                  << "\t"
                  << "\tnumSphaseCells(overall)"
                  << "\t"
                  << "\tnumSPhaseCellsPerLayer:0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13"
                  << "\t"
                  << "\tnumProliferatingCells(overall)"
                  << "\t"
                  << "\tnumProliferatingCellsPerLayer:0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13"
                  << "\t"
                  << "\tnumCellsLesion"
                  << "\tradiusLesion"
                  << "\t"
                  << "\tareaDensity"
                  << "\tlesionArea"
                  << "\tareaDensityVoxelized"
                  << "\tareaDensityVoxelizedWithLesion"
                  << std::endl;

    std::cout << time
              << "\t" << mObservables.numCells;
    for ( int i=0; i<14; ++i)
            std::cout << "\t" << mObservables.numCellsPerLayer[i];
    std::cout << "\t" << mObservables.numSPhaseCells
              << "\t";
    for ( int i=0; i<14; ++i)
        std::cout << "\t" << mObservables.numSPhaseCellsPerLayer[i];
    std::cout << "\t"
              << "\t" << mObservables.numProliferatingCells
              << "\t";
    for ( int i=0; i<14; ++i)
        std::cout << "\t" << mObservables.numProliferatingCellsPerLayer[i];
    std::cout << "\t"
              << "\t" << mObservables.numCellsLesion
              << "\t" << mObservables.radiusLesion
              << "\t"
              << "\t" << mObservables.areaDensity
              << "\t" << mObservables.lesionArea
              << "\t" << mObservables.areaDensityVoxelized
              << "\t" << mObservables.areaDensityVoxelizedWithLesion
              << std::endl;
}


inline
void
ModelCellsSpherical::rescaleLangevinContrib( double timeStepFactor )
{
    std::vector<CellSpherical *>::iterator cell;

    double scaleFactor = ( 1 - sqrt( 1/timeStepFactor ) );

    for ( cell = cells.begin(); cell != cells.end(); ++cell )
    {
        (*cell)->directedForce.x -= (*cell)->mLangevinForce.x * scaleFactor;
        (*cell)->directedForce.y -= (*cell)->mLangevinForce.y * scaleFactor;
        (*cell)->directedForce.z -= (*cell)->mLangevinForce.z * scaleFactor;
    }

	// added by Jieling
	for (int i = 0; i < mpCollagenNetwork->getNodes().size(); i++)
	{
		mpCollagenNetwork->getNode(i)->directedForce.x -= mpCollagenNetwork->getNode(i)->mLinearForce.x * scaleFactor;
		mpCollagenNetwork->getNode(i)->directedForce.y -= mpCollagenNetwork->getNode(i)->mLinearForce.y * scaleFactor;
		mpCollagenNetwork->getNode(i)->directedForce.z -= mpCollagenNetwork->getNode(i)->mLinearForce.z * scaleFactor;

		mpCollagenNetwork->getNode(i)->directedForce.x -= mpCollagenNetwork->getNode(i)->mRotationalForce.x * scaleFactor;
		mpCollagenNetwork->getNode(i)->directedForce.y -= mpCollagenNetwork->getNode(i)->mRotationalForce.y * scaleFactor;
		mpCollagenNetwork->getNode(i)->directedForce.z -= mpCollagenNetwork->getNode(i)->mRotationalForce.z * scaleFactor;
	}
}

#endif
