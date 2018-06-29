#include <cstdio>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cfloat>
#include <time.h>
#include <ctime>

#include <fstream>
#include <iostream>


#include "ModelCellsSpherical.h"
#include "tools/model/CSModelTools.h"
//#include "../../../tools/model/CSModelTools.h"
#include "../../Core.h"
#include "../../../tools/parameters/CSParameterContext.h"
#include "../../../tools/parameters/CSParameterContextTemporary.h"
#include "../../../tools/model/BoundingBoxList.h"
#include "../../Interactions/CSInteractionHertz.h"
#include "../../Interactions/CSInteractionJKR.h"
#include "../../Interactions/CSInteractionFrictionMatrix.h"
#include "../../../gui/QCSSimulationThread.h"
#include "../../Elements/ModelElementBarrierTriangle.h"
//#include "../../Elements/ModelElementVesselSphere.h"
#include "../../../tools/dataIO/vtp/CSVTPWriter.h"
#include "../../../tools/dataIO/xml/CSXMLWriter.h"

#include "../../BasicDatatypes/GraphSphere.h"

#include "LobuleLesionVoxelization.h"

// added by Jieling
#include "model/Lattice/LinearSpring.h"
#include "model/Lattice/RotationalSpring.h"

#include <QtCore>

#include <H5Cpp.h>
#include "../../../tools/dataIO/hdf5/CSHDF5Helpers.h"
#include "../../../tools/math/mathematics.h"

#include <locale.h>


const std::string ModelCellsSpherical::xmlType = "ModelCellsSpherical";

const std::string ModelCellsSpherical::mContactModelNames[] =
{
    "Hertz Model",
    "JKR Model",
    ""
};

const std::string ModelCellsSpherical::mLobuleShapeNames[] =
{
    "None",
    "Hexagonal",
    "Quadric",
    ""
};

const std::string ModelCellsSpherical::mScenarioNames[] =
{
    "Read from Data",
    "Single Cell",
    "Embedding Medium",
    "Sphere Packing",
    "Lobule Regeneration",
    ""
};


// Preliminary:  fixed CellSphericalPolar 'region angle'.
//   This value corresponds to an angle of 15 degrees around the polar axis
#define POLAR_REGION_ANGLE .2618


std::deque< std::pair<double, std::vector<double> > > ModelCellsSpherical::mvProliferationData = {
     //{ 0,
 //      { 5.555555, 1.666667, 2.128623, 2.528736, 2.5, 0., 1.111111, 0.,0.,0.,0.,0.,0.,0.} },
     //  { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.} },
     //{ .5,
     //  { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.} },
     //{ 1,
     //  { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.} },
     { 2,
       { 13.400002, 18.521447, 18.714073, 18.343575, 11.403219, 7.692308, 0.,0.,0.,0.,0.,0.,0.,0.} },
     { 3,
       { 35.180954, 26.545464, 21.775543, 13.474222, 19.712429, 8.571428, 8.475783, 0.,0.,0.,0.,0.,0.,0.} },
     { 4,
       { 5.749891, 5.193867, 7.450827, 7.892646, 3.926282, 4.6443, 2.561072, 0.,0.,0.,0.,0.,0.,0.} },
     { 8,
       { 0., 0., 0.806451, 0.625, 0., 0.480769, 0.,0.,0.,0.,0.,0.,0.,0.} },
     { 16,
       {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.} }
 };

/*    { 2,
      { 0.171386, 0.2368893, 0.239353, 0.2346143, 0.1458472, 0.09838462, 0., 0., 0., 0., 0., 0., 0., 0. } },
    { 3,
      { 0.4499644, 0.3395165, 0.2785092, 0.1723353, 0.252122, 0.1096286, 0.1084053, 0.2717875, 0.0799375, 0., 0., 0., 0., 0., 0. } },
    { 4,
      { 0.0735411, 0.06642956, 0.09529608, 0.1009469, 0.05021715, 0.0594006, 0.03275611, 0., 0., 0., 0., 0., 0., 0. } },
    { 5,
      { 0.05515583, 0.04982217, 0.07405069, 0.07770864, 0.03766286, 0.04608771, 0.02456708, 0., 0., 0., 0., 0., 0., 0. } },
    { 6,
      { 0.03677055, 0.03321478, 0.0528053, 0.05447035, 0.02510857, 0.03277482, 0.01637806, 0., 0., 0., 0., 0., 0., 0. } },
    { 7,
      { 0.01838528, 0.01660739, 0.0315599, 0.03123205, 0.01255429, 0.01946193, 0.008189028, 0., 0., 0., 0., 0., 0., 0. } },
    { 8,
      { 0., 0., 0.01031451, 0.00799375, 0., 0.006149036, 0., 0., 0.01065833, 0., 0., 0., 0., 0. } } };
*/
/*
 * Constructor
 */

ModelCellsSpherical::ModelCellsSpherical()
    : CSModel (3),
      cells2( NULL ),
      mScenario( ScenarioReadFromData ),
      mContactModel( ContactModelHertz ),
      mUseDynamicalTimeSteps(true),
      mTimeScalingFactor(2),
      mpInteractionHertz(NULL),
      mpInteractionJKR(NULL),
      mpInteractionFrictionMatrix(NULL),
      mpJKRElementDone(NULL),
      mJKRElementDoneSize(0),
      mpVelocities(NULL),
      mProblemAllocationSize(0),
      mpTmpVectorP(NULL),
      mpTmpVectorQ(NULL),
      mpResidualVector(NULL),
      mpPreconditioner(NULL),
      mpParameters(NULL),
      mpFrictionMatrices(NULL),
      mpGraphBloodVesselNetwork(NULL),
	  mpCollagenNetwork(NULL),
      mpParaview(NULL),
      mpLesionVoxelization(NULL)
{
    #pragma region Output

    // Console output
    core->tools->output->consoleModelCellsSpherical_Simulation<< "Init model.\n";

    if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << "ModelMonolayer::ModelMonolayer(). Create model\n";

    #pragma endregion

    // Unit conversions, Biological -> Dimensionless values for modeling etc.
    biolink = new BiologyLink(0); // Init for ModelCellsSpherical model

    // Model measurements
    observe = new Observation(this, 0);

    #pragma region Init Growth and Division parameters

    // register all parameters in our parameter context
    RegisterParameters();

    // Cycle time
    defaultCellCycleTime
    = biolink->scaleBiologyToInternal( biolink->cycletime_bio,
                                       BiologyLink::ScaleTime);

    // Cycle time deviation (Gaussian distribution)
    defaultCellCycleTimeStandardDeviation
    = biolink->scaleBiologyToInternal( biolink->cycletime_stddev_bio,
                                       BiologyLink::ScaleTime );

    // Set initial radius of the cells
    defaultInitialCellRadius = 0.5;

    // Set initial radius of the cells
    defaultDivisionCellRadius = 0.63;

    // Set initial radius of the cells
    defaultDivisionDistance = 0.4;

    #pragma endregion

    #pragma region Init Visualisation parameters

    lastColormodeCells = 0;

    #pragma endregion

    // Set model to strictly monolayer (Also calles Reset())
    SwitchSubmodel(0);


    this->mRadiusGyration = 0;

	stableStatus = false; // added by Jieling
	unLoaded = false;
	releaseStress = false; // change strain gradually
}

// end Constructor


/*
 *  region Setup
 */

void ModelCellsSpherical::Reset( bool loadScenario )
{
    #pragma region Init Biophysics parameters

    UpdateParametersFromBioLink();

    #pragma endregion

    #pragma region Output (Debug / Console)

    mOutputCounter = 0;

    if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << "ModelMonolayer::Reset(). Reset monolayer model\n";

    int dims = 3;

	is2D = false; // added by Jieling

    // Console output
    if (is2D)
    {
        core->tools->output->consoleModelCellsSpherical_Simulation << "Reset model (strictly monolayer).\n";
        dims = 2;
    }
    else core->tools->output->consoleModelCellsSpherical_Simulation << "Reset model (tumor spheroid).\n";

    #pragma endregion

    if ( cells2 && loadScenario )
    {
        delete cells2;
        cells2 = nullptr;
    }


    if (cells2 == nullptr)
        cells2 = new BoundingBoxList( is2D ? 2 : 3 );

	cells2->setDimensions(3); // added by Jieling

    mpBBList = cells2;

    #pragma region Reset interactions

    // the interactions will be allocated in SetupSimulation()
    if ( mpInteractionFrictionMatrix )
        delete mpInteractionFrictionMatrix;

    mpInteractionFrictionMatrix = NULL;

    if ( mpInteractionHertz )
        delete mpInteractionHertz;

    mpInteractionHertz = NULL;

    if ( mpInteractionJKR )
        delete mpInteractionJKR;

    mpInteractionJKR = NULL;


    // Tableaus for the conjugate gradient solver.
    // These will stay NULL, if friction was disabled.
    // They are (re-)allocated in chunks in InitForces()
    if ( mpPreconditioner )
    {
        free( mpPreconditioner );
        mpPreconditioner = NULL;
    }

    if ( mpResidualVector )
    {
        free( mpResidualVector );
        mpResidualVector = NULL;
    }

    if ( mpTmpVectorP )
    {
        free( mpTmpVectorP );
        mpTmpVectorP = NULL;
    }

    if ( mpTmpVectorQ )
    {
        free( mpTmpVectorQ );
        mpTmpVectorQ = NULL;
    }

    mProblemAllocationSize = 0;

    #pragma endregion

    #pragma resetting when started from scratch
    // reset only when starting from scratch - to prevent overwriting/deletion of read-in data:

    if ( loadScenario )
    {
        // The state of mRandom is read in from data, if loadScenario is false.
        mRandom.Init();

        // if we loaded from data bundle, time is set accordingly by CSXMLReader,
        // while mStartTime denotes the start time of the original simulation
        time = mStartTime; // Reset time

        // Prelim. Fixed time step // Normal: 0.1 - 1s
        //timeStep = 1;

        while (cells.size())
            RemoveCell(cells.back());

		// added by Jieling
		for (int i0 = 0; i0 < (int)HSCs.size(); i0++)
		{
			delete HSCs[i0];
			HSCs[i0] = NULL;
		}
		HSCs.clear();

        // Remove cells from arena
        mpArena->clear();

        #pragma endregion

        // handle requests to read in mxf to be like choosing ScenarioLobuleRegeneration:
        if (mReadCells && mBloodVesselNetwork)
            mScenario = ScenarioLobuleRegeneration;

		// added by Jieling
		if (this->mpCollagenNetwork != NULL)
			delete this->mpCollagenNetwork;
		this->mpCollagenNetwork = NULL;
		// using Andreas's Lattice to store collagen fibers
		mpCollagenNetwork = new Lattice();
		mpCollagenNetwork->setArena(mpArena);
		mpCollagenNetwork->setBoundingBoxList(cells2);

		blockSize = 3;
		blockNumber = 11; // at least 2
		double blockNumberHalf = (blockNumber - 1) * 0.5;
		for (int ix = 0; ix < blockNumber; ix++)
		{
			for (int iy = 0; iy < blockNumber; iy++)
			{
				for (int iz = 0; iz < blockNumber; iz++)
				{
					double x = (ix - blockNumberHalf) * blockSize;
					double y = (iy - blockNumberHalf) * blockSize;
					double z = (iz - blockNumberHalf) * blockSize;

					double sd = blockSize / 10;
					double dx = mRandom.GetRandomGauss(0., sd);
					double dy = mRandom.GetRandomGauss(0., sd);
					double dz = mRandom.GetRandomGauss(0., sd);
					double dd = std::sqrt(dx*dx + dy*dy + dz*dz);
					x += dx / dd * sd;
					y += dy / dd * sd;
					z += dz / dd * sd;
					
					ModelElementLatticeNode *CollagenNode = new ModelElementLatticeNode(x, y, z, 0.05);
					CollagenNode->mYoung = defaultYoungModulusECM;
					CollagenNode->mPoisson = defaultPoissonRatioECM;
					mpCollagenNetwork->addNode(CollagenNode);
					CollagenNode->latticeIndex = mpCollagenNetwork->getNodes().size();
					if (iy == 0)
						CollagenNode->bottom = true;
					if (iy == blockNumber - 1)
						CollagenNode->top = true;
					mpArena->addObject(CollagenNode->GLObject()); // for test
				}
			}
		}
		/**************
		    y 
		     | 1320~1330 (10, 0, 0)
		     o-----o
			 |     |
			 |     |
			 o-----o-->x
            0~10 (-10, 0, 0)
		**************/
		mpCollagenNetwork->createNetwork(blockSize * 1.8, 1); // each fibre has several springs
		// initialization
		stressForStrainTest();
		LoadedStep = 0;
		unLoadedStep = 0;
		// test for bar visualization
		//ModelElementLatticeNode * testNode = new ModelElementLatticeNode(5, 5, 0, 0.05);
		//mpCollagenNetwork->addNode(testNode);
		//LinearSpring *STest = static_cast<LinearSpring*>(mpCollagenNetwork->getSprings().at(1));
		//STest->changeN1(testNode);
    }
    #pragma endregion

    mContactModel
    = (ContactModel)((CSParameterChoice *)mpParameters->findParameter("Contact Model")->dataPointer())->currentIndex();

    enableSimulation = false;

	mUseDynamicalTimeSteps = false; // added by Jieling

    //output
    if( this->mpParaview)
      delete this->mpParaview;

    mOutputPrefix += (mOutputPrefix.size()) ? "_" : "";
    if (mOutputPath.size()==0)
        mOutputPath = ".";
    if (mOutputPath.back()!='/')
        mOutputPath += '/';
    this->mpParaview = new CSVTPWriter( this->mOutputPath + this->mOutputPrefix + "hepatocytes_" );
    this->mpParaview->setOutputRadius( true );
    this->mpParaview->setOutputQuiescence( true );
    this->mpParaview->setOutputPressure( true );

    if (mUsePolarCells)
        this->mpParaview->setOutputPolarVector(true);
    else
        this->mpParaview->setOutputPolarVector(false);

    this->mpParaview->setOutput_sphericalCells_1D_cut( 10 );

//    this->mpParaview->exec(this);//TODO no cell if using batch mode!


    //output xml file
    std::string xmlOutputFileName = this->mOutputPath + mOutputPrefix + "_parameter.xml";
     CSXMLWriter * writer = new CSXMLWriter(xmlOutputFileName.c_str());
    writer->exec(this);

    time_t now = std::time(0);
    char timestamp[22];
    strftime(timestamp, 22, "%d.%m.%Y \t %H:%M:%S", localtime(&now));

    std::ofstream F;
    std::string outputFileName = mOutputPath + mOutputPrefix + "_CellDivision.dat";
    F.open( outputFileName.c_str(), std::ios::out|std::ios::app );
    F << time << "\t"
      << std::time(0) << "\t"
      << this->cells.size() << "\t"
      << this->mRadiusGyration << "\t"
      << timestamp;
    F << std::endl;
    F.close();

}


void ModelCellsSpherical::SetupSimulation()
{
    // the maximum allowed displacement for dynamic step size adaptation:
    // set to 10% of a cell's initial radius, but squared.
    this->mMaximumDisplacementSquared =   this->mMaximumDisplacementCellRadius
                                        * this->mMaximumDisplacementCellRadius
                                        * this->defaultInitialCellRadius
                                        * this->defaultInitialCellRadius;

    // Remember simulation start for progress bar
    simulateFromDays = biolink->getTimeInDays(time);

    // Target simulation time has already been reached -> No initial observation
    if (simulateFromDays >= simulateUntilDays) return;

    mpSimulationThread->setUpdateInterval(50);

    mpSimulationThread->setUpdateInterval(50);

    enableSimulation = true;

    if (enableObservation)
    {
        core->tools->output->logfile << " In SetupSimulation: " << name  << "\n";;

        // Initial observation
        observe->ObserveModelDescription();
        // initial ObserveModelMeasures() now in SimulateInThread():
        //observe->ObserveModelMeasures();

        // Set next observation time
        nextObservationTime = simulateFromDays;
    }


    //
    // Initializing interactions:
    //
    mpInteractionFrictionMatrix = new CSInteractionFrictionMatrix();

    // Connecting model parameters to mpInteractionFrictionMatrix
    mpInteractionFrictionMatrix->mpGammaCellsParallel      = &gammaCellsParallel;
    mpInteractionFrictionMatrix->mpGammaCellsPerpendicular = &gammaCellsPerpendicular;


    if ( mContactModel == ContactModelHertz )
    {
        mpInteractionHertz = new CSInteractionHertz();

        // Connecting model parameters to mpInteractionHertz
        mpInteractionHertz->mIs2D = is2D;
        mpInteractionHertz->mpSingleBondEnergy = &singleBondEnergy;
        mpInteractionHertz->mpAdhesionDensity  = &adhesionDensity;

        // Connecting 'output' of mpInteractionHertz to 'input' of
        // mpInteractionFrictionMatrix from the Hertz force calculation:
        mpInteractionFrictionMatrix->mpDistance    = &mpInteractionHertz->mDistance;
        mpInteractionFrictionMatrix->mpContactArea = &mpInteractionHertz->mContactArea;
    }
    else
    {
        mpInteractionJKR = new CSInteractionJKR();

        // Connecting model parameters to mpInteractionHertz
        mpInteractionJKR->mIs2D = is2D;
        mpInteractionJKR->mpSingleBondEnergy = &singleBondEnergy;
        mpInteractionJKR->mpAdhesionDensity  = &adhesionDensity;

        // Connecting 'output' of mpInteractionHertz to 'input' of
        // mpInteractionFrictionMatrix from the JKR force calculation:
        mpInteractionFrictionMatrix->mpDistance    = &mpInteractionJKR->mDistance;
        mpInteractionFrictionMatrix->mpContactArea = &mpInteractionJKR->mContactArea;
    }


    UpdateParametersForAllCells();
}

// end Setup


/*
 * Simulation
 */

//! Calculates the next model state t + dt

void ModelCellsSpherical::Simulate()
{
	//if (time < 10) {
#pragma region Update forces

		InitForces(); // Set everything to 0
		//cout << " -> InitForces is done! " << endl;

		UpdateInteractions(); // Cell-cell interaction with Axis-Aligned Bounding Boxes.
		//cout << " -> UpdateInteractions is done! " << endl;

		UpdateSingleCellEffects();
		//cout << " -> updateSingleCellEffects is done! " << endl;

		// added by Jieling
		// add forces from the hepatocyte to drive the ECM sphere towards the CV along the vessel
		UpdateSingleECMEffects();
		//cout << " -> UpdateSingeECMEffects is done! " << endl;

		// Compute new position based on forces
		mpFrictionMatrices = mpInteractionFrictionMatrix->mpFrictionMatrices;

		solveSystem();
		//FO << "ModelCellsSpherical.cpp / Simulate()	after solveSystem() " << endl;
		if (mUseDynamicalTimeSteps)
		{
			double timeStepFactor = 1; // will be determined in the while(true) loop and applied to this.timeStep later

			bool shootOver = false;
			while (!shootOver)
			{
				solveSystem();

				double displacementSquared = mVelocityMaxSquared * (timeStep * timeStepFactor) * (timeStep * timeStepFactor);
				if (displacementSquared <= mMaximumDisplacementSquared) { // upscale

					if (mTimeScalingFactor * mTimeScalingFactor * displacementSquared > mMaximumDisplacementSquared) {
						break;
					}

					if (timeStepFactor < 1) { // prevent upscale-downscale-loops
						break;
					}

					if (timeStep * timeStepFactor >= mDynamicTimeStepMax) {
						break;
					}

					timeStepFactor *= mTimeScalingFactor;
				}
				else { // downscale

					if (timeStep * timeStepFactor <= mDynamicTimeStepMin) {
						break;
					}

					// the last timeStep scaling was up, we shot over, which means,
					// we need to scale down and re-solve the system
					shootOver = timeStepFactor > 1;

					timeStepFactor /= mTimeScalingFactor;
				}

				rescaleLangevinContrib(timeStepFactor);

				solveSystem();
			}

			timeStep *= timeStepFactor;
		}

		double maxDiff = 0.;
		for (unsigned int i = 0; i < cells2->size(); ++i)
		{
			ModelElement * element = cells2->element(i);

			if (element == NULL) {
				continue;
			}

			if (element->mStatic) {
				continue;
			}

			int j = 3 * i;

			// ToDo:  Code pertaining dumbbell should be handled by Cells.
			// ToDo:  Hierarchically handled bitmap for types, so that derived
			//        types will not be needed for this and similar tests.
			if (mUseDumbbell)
			{
				if (element->mType == ModelElement::TypeCellSpherical || element->mType == ModelElement::TypeCellSphericalPolar)
				{
					CellSpherical * cell = static_cast<CellSpherical *>(element);
					if (cell->mpDivisionPartner)
					{
						double *daughterVelocity = mpVelocities + 3 * element->mGlobalIndex;
						double *motherVelocity = mpVelocities + 3 * cell->mpDivisionPartner->mGlobalIndex;

						daughterVelocity[0] = 0.5 * (daughterVelocity[0] + motherVelocity[0]);
						daughterVelocity[1] = 0.5 * (daughterVelocity[1] + motherVelocity[1]);
						daughterVelocity[2] = 0.5 * (daughterVelocity[2] + motherVelocity[2]);

						motherVelocity[0] = daughterVelocity[0];
						motherVelocity[1] = daughterVelocity[1];
						motherVelocity[2] = daughterVelocity[2];
					}
				}
			}


			//std::cout << "	-> element " << element->mGlobalIndex << " velocity: "
			//	<< mpVelocities[j] << ", "
			//	<< mpVelocities[j + 1] << ", "
			//	<< mpVelocities[j + 2] << std::endl;
			if (element->mType == ModelElement::TypeCellLatticeNode)
			{
				if (((ModelElementLatticeNode*)element)->touched)// && !unLoaded) // the fixed nodes
				{
					if (!releaseStress) // when gradually releasing the stress, fix the nodes before unloading
						continue;
				}
			}

			element->position.x += mpVelocities[j] * timeStep;

			if (dimension != 1)
			{
				if (element->mType == ModelElement::TypeCellLatticeNode)
				{
					if (!unLoaded) // load the external force
					{
						if (((ModelElementLatticeNode*)element)->touched)
						{
							// top and bottom nodes
						}
						else
						{
							// other nodes
							element->position.y += mpVelocities[j + 1] * timeStep;
						}
					} 
					else // unload the external force
					{
						if (((ModelElementLatticeNode*)element)->touched)
						{
							// top and bottom nodes
						}
						else
						{
							// other nodes
							element->position.y += mpVelocities[j + 1] * timeStep;
						}
					}
				}
				else
				{
					element->position.y += mpVelocities[j + 1] * timeStep;
				}

				if (!is2D)
				{
					element->position.z += mpVelocities[j + 2] * timeStep;
				}
			}
			// added by Jieling
			if (element->mType == ModelElement::TypeCellLatticeNode)
			{
				double mVelocityDiff = (((ModelElementLatticeNode*)element)->mVelocity.x - mpVelocities[j]) *
					(((ModelElementLatticeNode*)element)->mVelocity.x - mpVelocities[j]) +
					(((ModelElementLatticeNode*)element)->mVelocity.y - mpVelocities[j + 1]) *
					(((ModelElementLatticeNode*)element)->mVelocity.y - mpVelocities[j + 1]) +
					(((ModelElementLatticeNode*)element)->mVelocity.z - mpVelocities[j + 2]) *
					(((ModelElementLatticeNode*)element)->mVelocity.z - mpVelocities[j + 2]);
				if (mVelocityDiff > maxDiff) maxDiff = mVelocityDiff;
				((ModelElementLatticeNode*)element)->mVelocity.x = mpVelocities[j];
				((ModelElementLatticeNode*)element)->mVelocity.y = mpVelocities[j + 1];
				((ModelElementLatticeNode*)element)->mVelocity.z = mpVelocities[j + 2];
			}

			std::cout << "	-> element " << element->mGlobalIndex << " new position: "
				<< element->position.x << ", "
				<< element->position.y << ", "
				<< element->position.z << std::endl;
		}
		if (maxDiff < 1e-12) stableStatus = true;
		else stableStatus = false;

		//Cell growth and division
		GrowAndDivide();//Disable for test cases
	//}

    #pragma region Update model time

    time += timeStep;

	//std::cout << "	-> simulation at time " << time << ", timeStep: " << timeStep << ", stable: " << stableStatus << std::endl;

    #pragma endregion
}


void ModelCellsSpherical::Simulate(int numberOfSimulationSteps)
{
    for (int i = 1; i<numberOfSimulationSteps; i++)
    {
        Simulate();
    }
}

// end Simulation
void ModelCellsSpherical::RemoveCell( CellSpherical * cell )
{
    std::vector<CellSpherical *>::iterator found
    = std::find( cells.begin(), cells.end(), cell );

    if ( found != cells.end() )
        cells.erase( found );
    else
        std::cerr << "Error in ModelCellsSpherical:  Attempted to remove unregistered cell from cells vector\n";

    cells2->remove( cell );

    mpArena->removeObject( cell->GLObject() );

    //delete cell;

    // Hack for being able to destroy CellSpherical/CellSphericalPolar, in case
    // they were read in by hdf5:
    if (cell->mType == ModelElement::TypeCellSphericalPolar)
        static_cast<CellSphericalPolar *>(cell)->~CellSphericalPolar();
    else
        cell->~CellSpherical();
}


#pragma region Grow and Divide

// Grows all cells (using the given time delta) and applies cell divisions
void ModelCellsSpherical::GrowAndDivide()
{
    bool triggerProliferation = false;

    // For all cells
    for (unsigned long i = 0; i<cells.size(); i++)
    {
        CellSpherical * cell = cells[i];

        // If cell is not quiescent
        if ( cell->getState( Cell::StateQuiescent ) == false )
        {
            // Divide the cell if possible
            if ( cell->CanDivide() )
            {
                CellSpherical * newCell =
                    static_cast<CellSpherical *>(cell->Divide(this));
                AddCell( newCell );

                time_t now = std::time(0);
                char timestamp[22];
                strftime(timestamp, 22, "%d.%m.%Y \t %H:%M:%S", localtime(&now));

                this->CalcRadiusGyration();

                std::ofstream F;
                std::string outputFileName = mOutputPath + mOutputPrefix + "_CellDivision.dat";
                F.open( outputFileName.c_str(), std::ios::out|std::ios::app );
                F << time << "\t"
                  << std::time(0) << "\t"
                  << this->cells.size() << "\t"
                  << this->mRadiusGyration << "\t"
                  << timestamp;
                F << std::endl;
                F.close();

            }
        }
    }

	if (mScenario == ScenarioLobuleRegeneration)
	{
		/*if ( time >= nextLesionDiscretization )
		{
			nextLesionDiscretization += 8640;
			mpLesionVoxelization->exec();
			mObservables.lesionArea = mpLesionVoxelization->lesionArea();
		}*/
	}

    if ( time >= nextObservablesUpdate )
    {
      WriteoutObservables();
      nextObservablesUpdate += mProliferationUpdateInterval;
    }
}


// Obsoleted by virtual ModelElement * Cell::Divide(Model*)!
void ModelCellsSpherical::DivideCell(int i)
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
        if( this->dimension == 1){
          delta.x = 1;
          delta.y = 0;
          delta.z = 0;
        }else{
          if (is2D) mRandom.GetRandomUnitVector(&delta.x, &delta.y);
          else mRandom.GetRandomUnitVector(&delta.x, &delta.y, &delta.z);
        }
        // Rescale to division distance
        delta.x *= defaultDivisionDistance;
        delta.y *= defaultDivisionDistance;
        if (is2D == false) delta.z *= defaultDivisionDistance;

        #pragma endregion

        #pragma region Create daughter cell 1

        switch ( cells.at(i)->mType )
        {

        case ModelElement::TypeCellSphericalPolar:
            AddPolarCell( cells.at(i)->position.x + delta.x,
                          cells.at(i)->position.y + delta.y,
                          (is2D) ? 0. : cells.at(i)->position.z + delta.z );
            break;

        default:
        case ModelElement::TypeCellSpherical:
            AddCell( cells.at(i)->position.x + delta.x,
                     cells.at(i)->position.y + delta.y,
                     (is2D) ? 0. : cells.at(i)->position.z + delta.z);

        }

        #pragma endregion


        #pragma region Create daughter cell 2 (by modifying the existing cell)

        // Reset radius
        cells.at(i)->mRadius = defaultInitialCellRadius;

        // Update position
        cells.at(i)->position.x -= delta.x;
        cells.at(i)->position.y -= delta.y;
        if (is2D == false) cells.at(i)->position.z -= delta.z;

        // Recaclulate cell cycle time
        cells.at(i)->SetCycleTimeGaussClamped(defaultCellCycleTime, defaultCellCycleTimeStandardDeviation);

        // Update corresponding growth step
        cells.at(i)->UpdateDeltaRadius(timeStep);

        #pragma endregion

        #pragma region Prelim: State transition (Growing -> Quiescent)

        if ( (*cells.at(i)).lastPressure > this->mQuiescentMax)
        {
            cells.at(i)->setState( Cell::StateQuiescent );
            cells.at(cells.size()-1)->setState( Cell::StateQuiescent );

            #pragma region Debug output

            // if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << "Made cell quiecent by absolute force " << (*cells.at(i)).lastForceAbsolut;

            #pragma endregion
        }

        #pragma endregion
    }
}


//! Adds and initializes a new cell
//! This method is used by the model to add cells at the given coordinates.
void ModelCellsSpherical::AddCell(double x, double y, double z)
{
    CellSpherical * newCell = new CellSpherical(x,y,z);

    // we do this in AddCell( CellSpherical ) again,
    // both are necessary: this one is for setting the cell's cycleTime,
    // the other is for cells generated from data
    newCell->mpRandom = &mRandom;

    #pragma region Size and growth specifics

    // Set current radius
    newCell->mRadius = defaultInitialCellRadius;

    // Set new cycle time
    newCell->SetCycleTimeGaussClamped(defaultCellCycleTime, defaultCellCycleTimeStandardDeviation);

    #pragma endregion

    #pragma region Cell cycle state

    // Initially, cells are non-quiescent
    newCell->setState(Cell::StateQuiescent, false);

    #pragma endregion

    newCell->mDetachTime = time;

    // Add newly created cell to population
    AddCell( newCell );
}


//! Adds and initializes a new cell
//! This method is used by the model to add cells at the given coordinates.
void ModelCellsSpherical::AddPolarCell(double x, double y, double z)
{
    CellSphericalPolar * newCell = new CellSphericalPolar(x,y,z);

    // we do this in AddCell( CellSpherical ) again,
    // both are necessary: this one is for setting the cell's cycleTime,
    // the other is for cells generated from data
    newCell->mpRandom = &mRandom;

    #pragma region Size and growth specifics

    // Set current radius
    newCell->mRadius = defaultInitialCellRadius;

    // Set new cycle time
    newCell->SetCycleTimeGaussClamped(defaultCellCycleTime, defaultCellCycleTimeStandardDeviation);

    #pragma endregion

    #pragma region Cell cycle state

    // Initially, cells are non-quiescent
    newCell->setState(Cell::StateQuiescent, false);

    #pragma endregion

    double lx, ly, lz;
    mRandom.GetRandomUnitVector(&lx,&ly,&lz);

    newCell->mPolarDirection.Set(lx,ly,lz);
    newCell->mNewPolarDirection.Set(lx,ly,lz);

    // Angle of the adhesive pole regions (spherical caps) calculated from a
    // surface area percentage (mParmAdhesivesurface).
    // Two caps with the area 2*Pi*R*h, where R is the cell's radius and h the
    // height of the spherical cap.  With h=R(1-cos(alpha)), the area percentage
    // is
    //   Q = 4 Pi R^2 (1-cos(alpha)) / ( 4 Pi R^2 ) = (1-cos(alpha)),
    // and therefore
    //   alpha = acos( 1 - Q ), where Q is mParmAdhesivesurface/100.
    newCell->SetPoleRegionAngle( std::acos( 1. - .01 * mParmAdhesiveSurface ) );

    newCell->setMaxRotationAngleForMetropolis(mMaxRotationAngleForMetropolis*M_PI/180);

    newCell->mDetachTime = time;

    // Add newly created cell to population
    AddCell( newCell );
}


//! Adds an already created cell.  Useful when loading data.
void ModelCellsSpherical::AddCell(CellSpherical * newCell)
{
    UpdateParametersForCell(newCell);
    UpdateCellsStaining(newCell);

    // we do this in AddCell( double, double, double ) too, both instances are
    // necessary: this one is for cells generated from data, the other is for
    // setting the cell's cycleTime,
    newCell->mpRandom = &mRandom;
    newCell->mpModel  = static_cast<CSModel *>(this);

    cells.push_back(newCell);
    cells2->add(newCell);

    mpArena->addObject(newCell->GLObject());

    #pragma region Debug output

    // if ((*core->tools->output).debugMonolayer) (*core->tools->output).logfile << "ModelMonolayer::AddCell(). Created a new model cell at (" << x << " " << y << ") with radius " << newCell->mRadius << " dividing at " << newCell->divisionRadius << " with delta-radius " << newCell->deltaRadius << "\n";

    #pragma endregion
}


//Add cells from mxf file
void ModelCellsSpherical::AddReadCells( std::string & filename )
{

    setlocale (LC_ALL,"C");

    //read mxf cell position
    std::fstream f;

    const int l = 512;
    char cstring[l];

    memset( cstring, 0, l );

    std::string s = filename.c_str();

    f.open(filename.c_str(), std::ios::in);

    if ( !f.is_open() )
    {
        std::cerr << "Warning:  Couldn't open file \""
                  << filename
                  << "\" for reading." << std::endl;
        return;
    }

    //search postion
    bool search = 1;
    while( !f.eof() && search )
    {
        if ( strncmp( cstring, "<CellList>", 10 ) == 0 )
            search = 0;
        f.getline(cstring, sizeof(cstring));
    }

    //read position
    search = 1;
    while ( !f.eof() && search )
    {
        if ( strncmp( cstring, "</CellList>", 11 ) == 0 )
            search = 0;
        else
        {
            char *pch;
            pch = strtok(cstring," \t");
            pch = strtok(NULL," ");
            pch = strtok(NULL,"\"");
            pch = strtok(NULL,"\"");

            double tmp_position[3];

            int index = 0;

            while( pch != NULL && index < 3)
            {

                tmp_position[index] = atof(pch);
                pch = strtok(NULL,"\"");
                pch = strtok(NULL,"\"");
                index++;
            }
            if ( mUsePolarCells )
                AddPolarCell(tmp_position[0],tmp_position[1],tmp_position[2]);
            else
                AddCell(tmp_position[0],tmp_position[1],tmp_position[2]);

            cells.back()->setQuality(10,10);
            cells.back()->setState( Cell::StateQuiescent );
            cells.back()->mProliferationControl = Cell::ProliferationControlCustom;

            f.getline(cstring, sizeof(cstring));
        }
    }

    f.close();
}



void ModelCellsSpherical::AddLobuleShape(double radius, double height, double, double, double, double, double, bool )
{
    //Add lobule shape
    if( this->is2D == false)
    {
        double d0[3][3] = { {0,0,-height*0.5},{-radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5},{0,-radius,-height*0.5} };
        AddBarrierTriangle( d0 ,0 );
        double u0[3][3] = { {0,0,height*0.5},{0,-radius,height*0.5},{-radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5} };
        AddBarrierTriangle( u0 ,0 );
        double d1[3][3] = { {0,0,-height*0.5},{-radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5},{-radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5} };
        AddBarrierTriangle( d1 ,0 );
        double u1[3][3] = { {0,0,height*0.5},{-radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5},{-radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5} };
        AddBarrierTriangle( u1 ,0 );
        double d2[3][3] = { {0,0,-height*0.5},{0,radius,-height*0.5},{-radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5} };
        AddBarrierTriangle( d2 ,0 );
        double u2[3][3] = { {0,0,height*0.5},{-radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5},{0,radius,height*0.5} };
        AddBarrierTriangle( u2 ,0 );
        double d3[3][3] = { {0,0,-height*0.5},{radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5},{0,radius,-height*0.5} };
        AddBarrierTriangle( d3 ,0 );
        double u3[3][3] = { {0,0,height*0.5},{0,radius,height*0.5},{radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5} };
        AddBarrierTriangle( u3 ,0 );
        double d4[3][3] = { {0,0,-height*0.5},{radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5},{radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5} };
        AddBarrierTriangle( d4 ,0 );
        double u4[3][3] = { {0,0,height*0.5},{radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5},{radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5} };
        AddBarrierTriangle( u4 ,0 );
        double d5[3][3] = { {0,0,-height*0.5},{0,-radius,-height*0.5} ,{radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5}};
        AddBarrierTriangle( d5 ,0 );
        double u5[3][3] = { {0,0,height*0.5},{radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5},{0,-radius,height*0.5} };
        AddBarrierTriangle( u5 ,0 );
    }

    double s0[3][3] = {{-radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5},{0,-radius,height*0.5},{0,-radius,-height*0.5} };
    AddBarrierTriangle( s0 );
    double s1[3][3] = {{-radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5},{-radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5},{0,-radius,height*0.5} };
    AddBarrierTriangle( s1 );
    double s2[3][3] = { {-radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5},{-radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5},{-radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5} };
    AddBarrierTriangle( s2 );
    double s3[3][3] = { {-radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5},{-radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5},{-radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5} };
    AddBarrierTriangle( s3 );
    double s4[3][3] = { {0,radius,-height*0.5},{-radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5},{-radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5} };
    AddBarrierTriangle( s4 );
    double s5[3][3] = { {0,radius,-height*0.5},{0,radius,height*0.5},{-radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5} };
    AddBarrierTriangle( s5 );
    double s6[3][3] = { {radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5},{0,radius,height*0.5},{0,radius,-height*0.5} };
    AddBarrierTriangle( s6 );
    double s7[3][3] = { {radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5},{radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5},{0,radius,height*0.5} };
    AddBarrierTriangle( s7 );
    double s8[3][3] = { {radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5},{radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5},{radius*std::sqrt(3.)*0.5,radius*0.5,-height*0.5} };
    AddBarrierTriangle( s8 );
    double s9[3][3] = { {radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5},{radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5},{radius*std::sqrt(3.)*0.5,radius*0.5,height*0.5} };
    AddBarrierTriangle( s9 );
    double s10[3][3] = { {0,-radius,-height*0.5},{radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5},{radius*std::sqrt(3.)*0.5,-radius*0.5,-height*0.5}};
    AddBarrierTriangle( s10 );
    double s11[3][3] = { {0,-radius,-height*0.5},{0,-radius,height*0.5},{radius*std::sqrt(3.)*0.5,-radius*0.5,height*0.5}};
    AddBarrierTriangle( s11 );




}
void ModelCellsSpherical::AddBarrierTriangle( double mPoints[][3], bool visible,  double epsilon, double r, double g, double b, double a )
{

    ModelElementBarrierTriangle *mebt = new ModelElementBarrierTriangle(0,0,0);
    for( int i = 0 ; i < 3 ; i++)
        mebt->setPoint(mPoints[i][0],mPoints[i][1],mPoints[i][2],i);
    mebt->setNormalVector();
    mebt->setBoundingBox(epsilon);
    mebt->SetColor(r,g,b,a);
    mebt->freeSurfaceArea = 1;
    cells2->add(mebt);
    mBarriers.push_back(mebt);
    if( visible )
        this->mpArena->addObject( mebt->GLObject() );


}


void ModelCellsSpherical::AddBloodVesselNetwork( GraphSphere * network )
{
	mpGraphBloodVesselNetwork = network;

    mpGraphBloodVesselNetwork->SampleGraph( 0.15 );

    mpGraphBloodVesselNetwork->setBoundingBoxList(cells2);

    mpGraphBloodVesselNetwork->defaultPoissonRatio = defaultPoissonRatioSinusoids;
    mpGraphBloodVesselNetwork->defaultYoungModulus = defaultYoungModulusSinusoids;

    mpGraphBloodVesselNetwork->updateParametersForAllSpheres();
    for( unsigned int i = 0 ; i< mpGraphBloodVesselNetwork->mvNode.size() ; i++)
        mpArena->addObject(mpGraphBloodVesselNetwork->mvNode[i]->GLObject());

    //set maximal values in each dimensions
    mpGraphBloodVesselNetwork->setBoundingBox();

    //connect nodes
    mpGraphBloodVesselNetwork->ConnectNodes();
}


void ModelCellsSpherical::AddBloodVesselNetwork( std::string & filename , int filetype )
{
    GraphSphere * network = new GraphSphere();

    network->setBoundingBoxList(cells2);

    network->read( filename.c_str(), filetype );

    AddBloodVesselNetwork(network);

	//FO << "	-> The number of vessel nodes right after AddBloodVessel: " << mpGraphBloodVesselNetwork->mvNode.size() << endl;

    //calc dimension of lobule
    double shift = mpGraphBloodVesselNetwork->calcDimensions(mlobule_radius, mlobule_height);

    //shift height to center
    mpGraphBloodVesselNetwork->shift( -shift , 2 );

    if( mScenario == ScenarioLobuleRegeneration )
    {
        for( unsigned int i = 0 ; i < cells.size() ; i++ )
            cells[i]->position.z -= shift;
    }
	//FO.close();
}

// added by Jieling
void ModelCellsSpherical::InitHSCsDistribution()
{
	// HSCs should be along the vessels
	if (mScenario == ScenarioLobuleRegeneration && mpGraphBloodVesselNetwork != NULL)
	{
		// VesselType = 1: CV; 2: PV; 3: Sinusoid; Add 1000 HSC, only along sinusoid
		int NumHSCs = 1000;
		int NumVesselSpheres = (int)mpGraphBloodVesselNetwork->mvNode.size();
		vector<int> sIndex;
		for (int i = 0; i < NumVesselSpheres; i++)
		{
			if (mpGraphBloodVesselNetwork->mvNode[i]->mVesselType == 3 && mpGraphBloodVesselNetwork->mvNode[i]->mRadius < 0.201)
			{
				// only sinusoid node with radius equals to 0.2 (not to consider nodes around the CV)
				sIndex.push_back(i);
			}
		}
		int NumSample = sIndex.size();
		for (int i = 0; i < NumHSCs; i++)
		{
			int randomIndex = sIndex[(int)(((double)rand() + 1.)/((RAND_MAX) + 1.) * NumSample)];
			double vx = mpGraphBloodVesselNetwork->mvNode[randomIndex]->position.x;
			double vy = mpGraphBloodVesselNetwork->mvNode[randomIndex]->position.y;
			double vz = mpGraphBloodVesselNetwork->mvNode[randomIndex]->position.z;
			double vr = mpGraphBloodVesselNetwork->mvNode[randomIndex]->mRadius * 1.0; // very close to the sinusoid
			int NumNeighbor = (int)mpGraphBloodVesselNetwork->mvNode[randomIndex]->mvpNeighbor.size();
			if (NumNeighbor == 0)
			{
				std::cerr << "Error: no neighbor vessels found for vessel cell " << randomIndex << endl;
				break;
			}
			int NeighborIndex = (int)(((double)rand() + 1.) / ((RAND_MAX)+1.) * NumNeighbor);
			double nvx = mpGraphBloodVesselNetwork->mvNode[randomIndex]->mvpNeighbor[NeighborIndex]->position.x;
			double nvy = mpGraphBloodVesselNetwork->mvNode[randomIndex]->mvpNeighbor[NeighborIndex]->position.y;
			double nvz = mpGraphBloodVesselNetwork->mvNode[randomIndex]->mvpNeighbor[NeighborIndex]->position.z;
			double nvd = sqrt((nvx - vx) * (nvx - vx) + (nvy - vy) * (nvy - vy) + (nvz - vz) * (nvz - vz));
			if (nvd == 0) break;
			/*****************************************************
			  o----------o
			  v    |    nv
				   o
				  HSC

			 Sampled HSC is along the direction from vessel cell to 
			 its neighbor vessel cell. We choose the middle position 
			 between vessel cell and its neighbor. The distance should
			 be more than the radius of the vessel cell. We need to 
			 calculate the rotation matrix to get the position by assuming
			 the direction of vessel cell to its neighbor is rotated from 
			 z-axis.

			 vector a: (0, 0, 1)
			 vector b: vessel cell -> its neighbor
			 vector v = a x b = (v1, v2, v3)
			 c = cos(a,b) = a * b
			 vx = |   0 -v3  v2 |
			      |  v3   0 -v1 |
				  | -v2  v1   0 |
			 Rotation matrix R that rotates a to b is
			 R = I + vx + vx^2 * 1/(1+c) 
			******************************************************/
			double R[3][3];
			double a[3]; a[0] = 0; a[1] = 0; a[2] = 1;
			double b[3]; b[0] = (nvx - vx) / nvd; b[1] = (nvy - vy) / nvd; b[2] = (nvz - vz) / nvd;
			double cos_a_b = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
			double v[3]; v[0] = a[1] * b[2] - a[2] * b[1]; v[1] = a[2] * b[0] - a[0] * b[2]; v[2] = a[0] * b[1] - a[1] * b[0];
			double v_x[3][3];
			v_x[0][0] = 0;    v_x[0][1] =-v[2]; v_x[0][2] = v[1];
			v_x[1][0] = v[2]; v_x[1][1] = 0;    v_x[1][2] =-v[0];
			v_x[2][0] =-v[1]; v_x[2][1] = v[0]; v_x[2][2] = 0;
			double v_x2[3][3];
			v_x2[0][0] = v_x[0][0] * v_x[0][0] + v_x[0][1] * v_x[1][0] + v_x[0][2] * v_x[2][0];
			v_x2[0][1] = v_x[0][0] * v_x[0][1] + v_x[0][1] * v_x[1][1] + v_x[0][2] * v_x[2][1];
			v_x2[0][2] = v_x[0][0] * v_x[0][2] + v_x[0][1] * v_x[1][2] + v_x[0][2] * v_x[2][2];
			v_x2[1][0] = v_x[1][0] * v_x[0][0] + v_x[1][1] * v_x[1][0] + v_x[1][2] * v_x[2][0];
			v_x2[1][1] = v_x[1][0] * v_x[0][1] + v_x[1][1] * v_x[1][1] + v_x[1][2] * v_x[2][1];
			v_x2[1][2] = v_x[1][0] * v_x[0][2] + v_x[1][1] * v_x[1][2] + v_x[1][2] * v_x[2][2];
			v_x2[2][0] = v_x[2][0] * v_x[0][0] + v_x[2][1] * v_x[1][0] + v_x[2][2] * v_x[2][0];
			v_x2[2][1] = v_x[2][0] * v_x[0][1] + v_x[2][1] * v_x[1][1] + v_x[2][2] * v_x[2][1];
			v_x2[2][2] = v_x[2][0] * v_x[0][2] + v_x[2][1] * v_x[1][2] + v_x[2][2] * v_x[2][2];
			double co_v_x2 = 1. / (1. + cos_a_b);
			R[0][0] = 1. + v_x[0][0] + v_x2[0][0] * co_v_x2;
			R[0][1] = 0. + v_x[0][1] + v_x2[0][1] * co_v_x2;
			R[0][2] = 0. + v_x[0][2] + v_x2[0][2] * co_v_x2;
			R[1][0] = 0. + v_x[1][0] + v_x2[1][0] * co_v_x2;
			R[1][1] = 1. + v_x[1][1] + v_x2[1][1] * co_v_x2;
			R[1][2] = 0. + v_x[1][2] + v_x2[1][2] * co_v_x2;
			R[2][0] = 0. + v_x[2][0] + v_x2[2][0] * co_v_x2;
			R[2][1] = 0. + v_x[2][1] + v_x2[2][1] * co_v_x2;
			R[2][2] = 1. + v_x[2][2] + v_x2[2][2] * co_v_x2;
			// sample one angle to locate the HSC
			double angleS = ((double)rand() + 1.) / (RAND_MAX + 1.) * M_PI * 2.; // 0 to 360 
			//angleS = 90; // test
			double oc[3]; 
			oc[0] = std::cos(angleS) * vr; 
			oc[1] = std::sin(angleS) * vr; 
			oc[2] = 0.; // the original coordinate of the sampled position
			double rc[3]; // rotate the sampled position with respect to the direction from vessel cell to its neighbor
			rc[0] = R[0][0] * oc[0] + R[0][1] * oc[1] + R[0][2] * oc[2] + (vx * 0.75 + nvx * 0.25);
			rc[1] = R[1][0] * oc[0] + R[1][1] * oc[1] + R[1][2] * oc[2] + (vy * 0.75 + nvy * 0.25);
			rc[2] = R[2][0] * oc[0] + R[2][1] * oc[1] + R[2][2] * oc[2] + (vz * 0.75 + nvz * 0.25);
			ModelElementECMSphere *HSCNew = new ModelElementECMSphere(rc[0],rc[1],rc[2]);
			HSCNew->mIndex = HSCs.size() + 1;
			HSCNew->setElementType(ModelElementECMSphere::HSC);
			// bind the vessel sphere to the HSC sphere
			HSCNew->vesselNeighbor = mpGraphBloodVesselNetwork->mvNode[randomIndex];
			HSCNew->vesselNeighbor->highlight = 2; // bound with HSC
			HSCs.push_back(HSCNew);
			//mpArena->addObject(HSCNew->GLObject()); // comment out after display test			
		}
		// test 
		for (int i = 0; i < HSCs.size(); i++)
		{
			//cout << "	-> HSC_" << i+1 << ": vessel " << HSCs[i]->vesselNeighbor->mRadius << ", radius: " << HSCs[i]->vesselNeighbor->mRadius << endl;
		}
	}
}

void ModelCellsSpherical::stressForStrainTest()
{
	// just for initialization
	sStrain = 0.3;
	double displace = blockSize * (blockNumber - 1) * sStrain;
	for (int ii = 0; ii < mpCollagenNetwork->getNodes().size(); ii++)
	{
		if (mpCollagenNetwork->getNode(ii)->top)
			mpCollagenNetwork->getNode(ii)->touched = true;
		if (mpCollagenNetwork->getNode(ii)->bottom)
			mpCollagenNetwork->getNode(ii)->touched = true;
	}
	/*for (int ix = 0; ix < blockNumber; ix++)
	{
		int iy = blockNumber - 1;
		double displaceY = displace / 2;
		for (int iz = 0; iz < blockNumber; iz++)
		{
			int nodeId = blockNumber * blockNumber * ix + blockNumber * iy + iz;
			mpCollagenNetwork->getNode(nodeId)->mStrainTestForce = 0.;
			mpCollagenNetwork->getNode(nodeId)->mStrainTestForceLoaded = 0.;
			mpCollagenNetwork->getNode(nodeId)->touched = true;
		}
		iy = 0;
		for (int iz = 0; iz < blockNumber; iz++)
		{
			int nodeId = blockNumber * blockNumber * ix + blockNumber * iy + iz;
			mpCollagenNetwork->getNode(nodeId)->mStrainTestForce = 0.;
			mpCollagenNetwork->getNode(nodeId)->mStrainTestForceLoaded = 0.;
			mpCollagenNetwork->getNode(nodeId)->touched = true;
		}
	}*/
}

void ModelCellsSpherical::printStressStrain()
{
	// test debug
	std::string outputFileName = "debugFile.txt";
	std::ofstream FO;
	FO.open(outputFileName.c_str(), std::ios::out | std::ios::app);
	//
	double topx = 0.;
	double topy = 0.;
	double topz = 0.;
	double bottomx = 0.;
	double bottomy = 0.;
	double bottomz = 0.;
	for (int ix = 0; ix < blockNumber; ix++)
	{
		int iy = blockNumber - 1; // for top nodes
		for (int iz = 0; iz < blockNumber; iz++)
		{
			int nodeId = blockNumber * blockNumber * ix + blockNumber * iy + iz;
			topx += mpCollagenNetwork->getNode(nodeId)->mStrainTestForce.x;
			topy += mpCollagenNetwork->getNode(nodeId)->mStrainTestForce.y;
			topz += mpCollagenNetwork->getNode(nodeId)->mStrainTestForce.z;
		}
		iy = 0; // for the bottom nodes
		for (int iz = 0; iz < blockNumber; iz++)
		{
			int nodeId = blockNumber * blockNumber * ix + blockNumber * iy + iz;
			bottomx += mpCollagenNetwork->getNode(nodeId)->mStrainTestForce.x;
			bottomy += mpCollagenNetwork->getNode(nodeId)->mStrainTestForce.y;
			bottomz += mpCollagenNetwork->getNode(nodeId)->mStrainTestForce.z;
		}
	}
	double topd = std::sqrt(topx*topx + topy*topy + topz*topz);
	double bottomd = std::sqrt(bottomx*bottomx + bottomy*bottomy + bottomz*bottomz);
	double forceT = (topd + bottomd) / 2;
	// get the strain
	double lStrain = getLatticeStrain();
	std::cout << "	-> ForceTop: " << topx << ", " << topy << ", " << topz << "; " << topd << std::endl;
	std::cout << "	-> ForceBottom: " << bottomx << ", " << bottomy << ", " << bottomz << "; " << bottomd << std::endl;
	std::cout << "	-> force: " << forceT << ", strain: " << lStrain << std::endl;
	FO << "Top," << topx << "," << topy << "," << topz << "," << topd
		<< ",Bottom," << bottomx << "," << bottomy << "," << bottomz << "," << bottomd
		<< ",Strain," << lStrain << ",Stable," << stableStatus << std::endl;
	FO.close();
}

double ModelCellsSpherical::getLatticeStrain()
{
	double aStrain = 0.;
	for (int ix = 0; ix < blockNumber; ix++)
	{
		// for top nodes
		int iy = blockNumber - 1;
		double aTopx = 0.;
		for (int iz = 0; iz < blockNumber; iz++)
		{
			int nodeId = blockNumber * blockNumber * ix + blockNumber * iy + iz;
			aTopx += mpCollagenNetwork->getNode(nodeId)->position.x;
		}
		aTopx /= blockNumber;
		// for bottom nodes
		iy = 0;
		double aBottomx = 0.;
		for (int iz = 0; iz < blockNumber; iz++)
		{
			int nodeId = blockNumber * blockNumber * ix + blockNumber * iy + iz;
			aBottomx += mpCollagenNetwork->getNode(nodeId)->position.x;
		}
		aBottomx /= blockNumber;
		//
		aStrain += (aTopx - aBottomx) / blockSize / (blockNumber - 1);
	}
	aStrain /= blockNumber;
	return aStrain;
}

void ModelCellsSpherical::AddShape(double radius, double height, double r, double g, double b, double a, double epsilon, bool visible)
{

    if( this->is2D == true)
        this->mlobule_height = 1;

    if( this->mLobuleShape == ModelCellsSpherical::Hexagonal )
        this->AddLobuleShape(this->mlobule_radius,this->mlobule_height,1,1,0,0.5,1);
    else if( this->mLobuleShape == ModelCellsSpherical::Quadric )
        this->AddQuadricShape(this->mlobule_radius,this->mlobule_height,1,1,0,0.5,1);
}
void ModelCellsSpherical::AddQuadricShape(double radius, double height, double r, double g, double b, double a, double epsilon, bool visible)
{

    //Add lobule shape
    if( this->is2D == false)
    {
        double d0[3][3] = { {-radius,-radius,-height*0.5},{radius,radius,-height*0.5},{-radius,radius,-height*0.5} };
        AddBarrierTriangle( d0, 0 );
        double d1[3][3] = { {-radius,-radius,-height*0.5},{radius,-radius,-height*0.5},{radius,radius,-height*0.5} };
        AddBarrierTriangle( d1, 0 );
        double u0[3][3] = { {-radius,-radius,height*0.5},{-radius,radius,height*0.5},{radius,radius,height*0.5} };
        AddBarrierTriangle( u0, 0 );
        double u1[3][3] = { {-radius,-radius,height*0.5},{radius,radius,height*0.5},{radius,-radius,height*0.5} };
        AddBarrierTriangle( u1, 0 );
    }

    double s0[3][3] = { {-radius,-radius,-height*0.5},{-radius,radius,height*0.5},{-radius,-radius,height*0.5} };
    AddBarrierTriangle( s0 );
    double s1[3][3] = { {-radius,radius,-height*0.5},{-radius,radius,height*0.5},{-radius,-radius,-height*0.5} };
    AddBarrierTriangle( s1 );
    double s2[3][3] = { {radius,-radius,-height*0.5},{radius,-radius,height*0.5},{radius,radius,height*0.5} };
    AddBarrierTriangle( s2 );
    double s3[3][3] = { {radius,radius,-height*0.5},{radius,-radius,-height*0.5},{radius,radius,height*0.5} };
    AddBarrierTriangle( s3 );
    double s4[3][3] = { {-radius,radius,-height*0.5},{radius,radius,height*0.5},{-radius,radius,height*0.5} };
    AddBarrierTriangle( s4 );
    double s5[3][3] = { {radius,radius,height*0.5},{-radius,radius,-height*0.5},{radius,radius,-height*0.5} };
    AddBarrierTriangle( s5 );
    double s6[3][3] = { {-radius,-radius,-height*0.5},{-radius,-radius,height*0.5},{radius,-radius,height*0.5} };
    AddBarrierTriangle( s6 );
    double s7[3][3] = { {radius,-radius,height*0.5},{radius,-radius,-height*0.5},{-radius,-radius,-height*0.5} };
    AddBarrierTriangle( s7 );
}


#pragma endregion



//Forces




void ModelCellsSpherical::InitForces()
{
    // For all cells
    for (unsigned int i = 0; i< cells2->EndIndex(); i++)
    {
		if (ModelElement *element = cells2->element(i)) // dead cells might remain in cells2 as NULL pointers!
		{
			element->Reset();
			// added by Jieling
			if (element->mType == ModelElement::TypeCellLatticeNode)
			{
				((ModelElementLatticeNode*)(element))->mLinearForce = 0;
				((ModelElementLatticeNode*)(element))->mRotationalForce = 0;
			}
		}
    }

    mpInteractionFrictionMatrix->Reset();

    std::size_t problemSize = 3*cells2->size();
    if ( problemSize >= mProblemAllocationSize )
    {
        // At least allocate  problemSize  elements, allow for 25% more, round up to the next 1024.
        std::size_t newAllocationSize = std::ceil( 1.25 * (problemSize >> 10) );
        newAllocationSize = ((newAllocationSize) ? newAllocationSize : 1) <<10;
        // newAllocationSize = newAllocationSize << 10; // allocation_chunk *= 1024;

        mpVelocities =
            (double *) realloc( mpVelocities,
                                newAllocationSize *sizeof(double));

        // set all new elements to nought
        memset( mpVelocities + mProblemAllocationSize, 0, (newAllocationSize-mProblemAllocationSize)*sizeof(double) );

        mProblemAllocationSize = newAllocationSize;
        mpPreconditioner =
            (double *) realloc( mpPreconditioner, mProblemAllocationSize*sizeof(double) );
        mpResidualVector =
            (double *) realloc( mpResidualVector, mProblemAllocationSize*sizeof(double) );

        mpTmpVectorP =
            (double *) realloc( mpTmpVectorP, mProblemAllocationSize*sizeof(double) );
        mpTmpVectorQ =
            (double *) realloc( mpTmpVectorQ, mProblemAllocationSize*sizeof(double) );

    }

    if (mpInteractionJKR)
    {
        if ( cells2->size() > mJKRElementDoneSize )
        {
            std::size_t difference = cells2->size() - mJKRElementDoneSize;
            difference = ((std::size_t) std::ceil( 1.25 * difference/1024 )) <<10;

            mJKRElementDoneSize += difference;
            // for JKR only
            mpJKRElementDone =
                (bool *) realloc( mpJKRElementDone, mJKRElementDoneSize*sizeof(bool) );
        }

        memset( mpJKRElementDone, 0, mJKRElementDoneSize*sizeof(bool) );
    }
}


void ModelCellsSpherical::UpdateInteractions()
{

	cells2->update();

    BoundingBoxList::iterator cellIterator = cells2->begin();

    while( cellIterator != cells2->end() )
    {
        ModelElement *object_1 = *cellIterator;

        CSListContainer< unsigned long > & interactingElements
            =  object_1->mIntersectingList;

        // prepare interactingElements for JKR
        size_t counter = 0;
        size_t switchToEstablishedContacts = interactingElements.size();

        if ( mpInteractionJKR )
        {
            // see if elements in the newly interactingElements are in the establishedContacts
            //  push_back into establishedContacts only if not.
            unsigned long * iterator;
            for ( iterator  = object_1->mEstablishedContacts.begin();
                  iterator != object_1->mEstablishedContacts.end();
                  ++iterator )
            {
                if ( !cells2->element(*iterator) )
                    continue;

                if ( !mpJKRElementDone[*iterator] )
                {
                    if ( interactingElements.find( *iterator ) == interactingElements.end() )
                        interactingElements.push_back( *iterator );
                }
            }

            // the first elements in interactingElements are from object_1->mEstablishedcontacts,
            // switch to true, when switchToEstablishedContacts is reached.
            mpInteractionJKR->mContactEstablished = false;
        }


        unsigned long * otherIterator = interactingElements.begin();

        while ( otherIterator != interactingElements.end() )
        {
            if ( *otherIterator < object_1->mGlobalIndex )
            {
                ++otherIterator;
                continue;
            }

            ModelElement * object_2 = cells2->element( *otherIterator );

            if ( mpInteractionHertz )
            {
                (*mpInteractionHertz)( object_1, object_2 );

                // if there's no contact area, there's no Hertz interaction,
                // next please!
                if ( mpInteractionHertz->mContactArea == 0 )
                {
                    ++otherIterator;
                    continue;
                }
            }
            else  // JKR!
            {
                if ( counter++ == switchToEstablishedContacts )
                    mpInteractionJKR->mContactEstablished = true;

                if ( mpJKRElementDone[*otherIterator] )
                {
                    ++otherIterator;
                    continue;
                }

                (*mpInteractionJKR)( object_1, object_2 );

                if ( mpInteractionJKR->mContactArea == 0 )
                {
                    ++otherIterator;
                    continue;
                }
            }

            // if ( useDetailedCellCellFriction )
            (*mpInteractionFrictionMatrix)( object_1, object_2 );

            ++otherIterator;
        }


        if (mpInteractionJKR)
            mpJKRElementDone[object_1->mGlobalIndex] = true;

        if ( object_1->mType != ModelElement::TypeBarrierTriangle )
        {
            ModelElementSphere *cell1 = static_cast<ModelElementSphere *>(object_1);
            cell1->freeSurfaceArea = 4 * M_PI * cell1->mRadius * cell1->mRadius - cell1->surfaceContactArea;
            if ( cell1->freeSurfaceArea <0 )
                cell1->freeSurfaceArea = 0;
        }

        ++cellIterator;
    }

}


void ModelCellsSpherical::UpdateCellFriction( CellSpherical * cell )
{
    // Add cell-cell friction
    cell->frictionCoefficient += gammaCellsPerpendicular * cell->surfaceContactArea;

    // Set surface area to the remaining surface area (not in contact with other cells)
    // Full surface of (spherical) cell is 4*pi*R^2
    cell->surfaceContactArea = 4 * M_PI * cell->mRadius * cell->mRadius - cell->surfaceContactArea;

    // Add cell-ecm friction
    if (cell->surfaceContactArea > 0) // If there is still some free surface available
    {
        cell->frictionCoefficient += gammaECM * cell->surfaceContactArea;
    }
}

// added by Jieling
void ModelCellsSpherical::UpdateLatticeNodeFriction(ModelElementLatticeNode * node)
{
	node->lastForceAbsolute = node->accumulatedForceAbsolute;
	node->lastPressure = node->accumulatedPressure;

	node->freeSurfaceArea = 4 * M_PI * node->mRadius * node->mRadius;
	// add latticeNode-ecm friction
	if (node->freeSurfaceArea > 0)
		node->frictionCoefficient += gammaECMLattice * node->freeSurfaceArea;
	else
		node->freeSurfaceArea = 0.;
}

// added by Jieling
inline 
void
ModelCellsSpherical::UpdateSingleECMEffects()
{
	if (mpCollagenNetwork != NULL)
	{
		for (int i = 0; i < mpCollagenNetwork->getNodes().size(); i++)
		{
			// add LatticeNode-ecm friction
			//UpdateLatticeNodeFriction(mpCollagenNetwork->getNode(i));
			//mpCollagenNetwork->getNode(i)->mStrainTestForce = 0.;
		}
		// store linear/rotational force for each lattice node
		for (int i = 0; i < mpCollagenNetwork->getSprings().size(); i++)
		{
			if (mpCollagenNetwork->getSprings().at(i)->mSpringType == LatticeSpring::TypeLinearSpring)
			{
				LinearSpring *lsI = (LinearSpring*)(mpCollagenNetwork->getSprings().at(i));
				lsI->strainTestForce();
			}
			else
			{
				RotationalSpring *rsI = (RotationalSpring*)(mpCollagenNetwork->getSprings().at(i));
				rsI->strainTestForce();
			}
		}
		if (!unLoaded) // gradually loading the strain
		{
			double displace = blockSize * (blockNumber - 1) * sStrain / 2;
			int stableStep = 200;
			int tracePoint = 10;
			int Remaining = LoadedStep % stableStep;
			int Fold = LoadedStep / stableStep;
			if (Remaining == 0 && Fold <= tracePoint)
			{
				double displaceY = displace / tracePoint;
				for (int ii = 0; ii < mpCollagenNetwork->getNodes().size(); ii++)
				{
					if (mpCollagenNetwork->getNode(ii)->top)
						mpCollagenNetwork->getNode(ii)->position.x -= displaceY;
					if (mpCollagenNetwork->getNode(ii)->bottom)
						mpCollagenNetwork->getNode(ii)->position.x += displaceY;
				}
				/*for (int ix = 0; ix < blockNumber; ix++)
				{
					int iy = blockNumber - 1; // top nodes
					for (int iz = 0; iz < blockNumber; iz++)
					{
						int nodeId = blockNumber * blockNumber * ix + blockNumber * iy + iz;
						mpCollagenNetwork->getNode(nodeId)->position.x += displaceY;
					}
					iy = 0;
					for (int iz = 0; iz < blockNumber; iz++)
					{
						int nodeId = blockNumber * blockNumber * ix + blockNumber * iy + iz;
						mpCollagenNetwork->getNode(nodeId)->position.x -= displaceY;
					}
				}*/
			}
			bool writeIt = false;
			if (Remaining == stableStep - 1 && Fold <= tracePoint) writeIt = true;
			if (writeIt) // printout the info
			{
				printStressStrain();
			}
			LoadedStep++;
			// check if the loading process is completed
			if (LoadedStep == stableStep * tracePoint)
			{
				unLoaded = true;
				// permanent elongation
				for (int i = 0; i < (int)mpCollagenNetwork->getSprings().size(); i++)
				{
					if (mpCollagenNetwork->getSprings().at(i)->mSpringType == LatticeSpring::TypeLinearSpring)
					{
						LinearSpring *lsI = (LinearSpring*)(mpCollagenNetwork->getSprings().at(i));
						lsI->elongation();
					}
				}
			}
			// update the unbinding and sliding
			//mpCollagenNetwork->updateFibreBound();
		}
		else if (unLoaded) // gradually unloading the strain
		{
			double displace = blockSize * (blockNumber - 1) * sStrain / 2;
			int stableStep = 200;
			int tracePoint = 10;
			int unRemaining = unLoadedStep % stableStep;
			int unFold = unLoadedStep / stableStep;
			if (unRemaining == 0 && unFold < tracePoint)
			{
				double displaceY = displace / tracePoint;
				for (int ii = 0; ii < mpCollagenNetwork->getNodes().size(); ii++)
				{
					if (mpCollagenNetwork->getNode(ii)->top)
					{
						if (!releaseStress) // gradually
							mpCollagenNetwork->getNode(ii)->position.x += displaceY;
					}
					if (mpCollagenNetwork->getNode(ii)->bottom)
					{
						if (!releaseStress) // gradually
							mpCollagenNetwork->getNode(ii)->position.x -= displaceY;
					}
				}
				/*for (int ix = 0; ix < blockNumber; ix++)
				{
					int iy = blockNumber - 1; // top nodes
					for (int iz = 0; iz < blockNumber; iz++)
					{
						int nodeId = blockNumber * blockNumber * ix + blockNumber * iy + iz;
						if (!releaseStress) // gradually
							mpCollagenNetwork->getNode(nodeId)->position.x -= displaceY;
					}
					iy = 0; // bottom nodes
					for (int iz = 0; iz < blockNumber; iz++)
					{
						int nodeId = blockNumber * blockNumber * ix + blockNumber * iy + iz;
						if (!releaseStress) // gradually
							mpCollagenNetwork->getNode(nodeId)->position.x += displaceY;
					}
				}*/
			}
			bool writeIt = false;
			if (releaseStress) writeIt = true; // completely release the stress, record every step
			else
			{
				if (unRemaining == stableStep - 1 && unFold <= tracePoint) writeIt = true; // release gradually, only record when stable is reached
			}
			if (writeIt) // printout the info
			{
				// sum up the shear force
				printStressStrain();
			}
			unLoadedStep++;
		}

		for (int i = 0; i < (int)mpCollagenNetwork->getSprings().size(); i++)
		{
			if (mpCollagenNetwork->getSprings().at(i)->mSpringType == LatticeSpring::TypeLinearSpring)
			{
				LinearSpring *lsI = (LinearSpring*)(mpCollagenNetwork->getSprings().at(i));
				lsI->update(1.); // Update the force of Linear Spring
			}
			else
			{
				RotationalSpring *rsI = (RotationalSpring*)(mpCollagenNetwork->getSprings().at(i));
				rsI->update(1.); // Update the force of Rotational Spring
			}
		}
		for (int i = 0; i < mpCollagenNetwork->getNodes().size(); i++)
		{
			// add LatticeNode-ecm friction
			UpdateLatticeNodeFriction(mpCollagenNetwork->getNode(i));
		}
	}
}

void ModelCellsSpherical::UpdateForcesLangevin( CellSpherical *cell )
{
    // Update sigma and langevin forces (after friction is completely known)
    double s;

    // For all cells
    // Earlier version: s^2 / timestep

    double    absoluteForce;

    if (is2D)
    {
        s = cell->frictionCoefficient * sqrt(2 * 2/*dimension*/ * defaultDiffusionConstantCells / this->timeStep);

        cell->mLangevinForce.x = mRandom.GetRandomGauss(0.0, s);
        cell->mLangevinForce.y = mRandom.GetRandomGauss(0.0, s);
        cell->mLangevinForce.z = 0;

        absoluteForce = sqrt(   cell->mLangevinForce.x*cell->mLangevinForce.x
                                + cell->mLangevinForce.y*cell->mLangevinForce.y );
    }
    else
    {
        s = cell->frictionCoefficient * sqrt(2 * dimension * defaultDiffusionConstantCells / this->timeStep);

        cell->mLangevinForce.x = mRandom.GetRandomGauss(0.0, s);
        cell->mLangevinForce.y = mRandom.GetRandomGauss(0.0, s);
        cell->mLangevinForce.z = mRandom.GetRandomGauss(0.0, s);

        absoluteForce = cell->mLangevinForce.Norm();

    }

    if (cell->mpDivisionPartner)
    {
        cell->directedForce.Add( cell->mLangevinForce.x/2,
                                 cell->mLangevinForce.y/2,
                                 cell->mLangevinForce.z/2 );
        cell->mpDivisionPartner->directedForce.Add( cell->mLangevinForce.x/2,
                cell->mLangevinForce.y/2,
                cell->mLangevinForce.z/2 );
    }
    else
        cell->directedForce.Add( cell->mLangevinForce.x,
                                 cell->mLangevinForce.y,
                                 cell->mLangevinForce.z );

    cell->lastForceAbsolute += absoluteForce;

    // Review?  how do the langevin forces contribute to pressure?  Omitted for now.
}

void ModelCellsSpherical::AddDirectedMotion( CellSpherical *cell )
{
    // ToDo:  change into a parameter
    double someFactor = 200;
    Vector3f directedForceVector = - cell->position;
    // force is driving to the z-axis (where the central vein should be):
    directedForceVector[2] = 0;
    // linear dependency on distance:
    double distance = directedForceVector.Normalize(); // |directedForcVector| = 1
    double directedMotionForce = someFactor * ( mlobule_radius - distance );
    directedForceVector *= directedMotionForce;

    cell->directedForce += directedForceVector;
}

//end Forces

//Solve equation of motion

// Solve the linear system with Conjugate Gradient method
// as described in Appendix B.2 of:
//  Schaller, G.:
//    On selected numerical approaches to cellular tissue,
//    PhD Thesis, 2005.
void
ModelCellsSpherical::solveSystem()
{
    int iteration = 0;
    double residualSquared;
    double maxResidualElementSquared; // the square of mpResidualVector's element with largest absolute value
    double rho;
    const unsigned int problemSize = 3*cells2->size();
    const double convergenceLimit = std::pow( DBL_EPSILON, 2/3. )/problemSize;

    // dummy variable to hold the square of the velocity vector:
    double velocitySquared;

    // t1m-debug
    // if (!time)
    //     PrintFrictionMatrix("matrix.txt");
    // !t1m-debug

    mVelocityMaxSquared = 0;

    // the initial residual:
    // q = A.v
    InitialResidualAndPreconditioner( mpVelocities, mpResidualVector, mpPreconditioner );

    residualSquared = 0;
    rho             = 0;
    for ( unsigned int i=0; i<problemSize; ++i )
    {
        residualSquared    += mpResidualVector[i]*mpResidualVector[i];
        mpTmpVectorP[i]     = mpPreconditioner[i]*mpResidualVector[i];
        rho                += mpTmpVectorP[i] * mpResidualVector[i];
    }

    if ( residualSquared < convergenceLimit )
        return;

    FrictionMatrixTimesV( mpTmpVectorP, mpTmpVectorQ );

    double qp = 0;
    for ( unsigned int i=0; i<problemSize; ++i )
    {
        qp += mpTmpVectorP[i] * mpTmpVectorQ[i];
    }

    double alpha = rho / qp;

    for ( unsigned int i=0; i<problemSize; i+=3 )
    {
        mpVelocities[i]     += alpha * mpTmpVectorP[i];
        mpVelocities[i+1]   += alpha * mpTmpVectorP[i+1];
        mpVelocities[i+2]   += alpha * mpTmpVectorP[i+2];

        velocitySquared = mpVelocities[i]   * mpVelocities[i]
                          + mpVelocities[i+1] * mpVelocities[i+1]
                          + mpVelocities[i+2] * mpVelocities[i+2];
        if  (velocitySquared > mVelocityMaxSquared)
        {
			// added by Jieling
			//std::cout << "		-> The first Large velocitySquared: " << velocitySquared << std::endl;
			int eIndex = (int)(i / 3);
			ModelElement * element = cells2->element(eIndex);
			if (element != NULL) {
				if (element->mType != ModelElement::TypeECMSphere) {
					mVelocityMaxSquared = velocitySquared;
				}
			}
        }

        mpResidualVector[i]   -= alpha * mpTmpVectorQ[i];
        mpResidualVector[i+1] -= alpha * mpTmpVectorQ[i+1];
        mpResidualVector[i+2] -= alpha * mpTmpVectorQ[i+2];
    }

    while ( true )
    {
        ++iteration;
        double lastRho = rho;

        residualSquared    = 0;
        maxResidualElementSquared = 0;
        rho = 0;
        for ( unsigned int i=0; i<problemSize; ++i )
        {
            double residualElementSquared  = mpResidualVector[i]*mpResidualVector[i];
            residualSquared               += residualElementSquared;
            if ( residualElementSquared > maxResidualElementSquared )
                maxResidualElementSquared = residualElementSquared;
            rho += mpPreconditioner[i] * residualElementSquared;
        }

        if ( maxResidualElementSquared < convergenceLimit )
            break;

        double beta = rho / lastRho;

        for ( unsigned int i=0; i<problemSize; ++i )
        {
            mpTmpVectorP[i] = mpPreconditioner[i]*mpResidualVector[i] + beta * mpTmpVectorP[i];
        }

        FrictionMatrixTimesV( mpTmpVectorP, mpTmpVectorQ );

        qp =0;
        for ( unsigned int i=0; i<problemSize; ++i )
        {
            qp += mpTmpVectorP[i] * mpTmpVectorQ[i];
        }

        alpha = rho / qp;

        for ( unsigned int i=0; i<problemSize; i+=3 )
        {
            mpVelocities[i]     += alpha * mpTmpVectorP[i];
            mpVelocities[i+1]   += alpha * mpTmpVectorP[i+1];
            mpVelocities[i+2]   += alpha * mpTmpVectorP[i+2];

            velocitySquared = mpVelocities[i]   * mpVelocities[i]
                              + mpVelocities[i+1] * mpVelocities[i+1]
                              + mpVelocities[i+2] * mpVelocities[i+2];
            if  (velocitySquared > mVelocityMaxSquared)
            {
                mVelocityMaxSquared = velocitySquared;
            }

            mpResidualVector[i]   -= alpha * mpTmpVectorQ[i];
            mpResidualVector[i+1] -= alpha * mpTmpVectorQ[i+1];
            mpResidualVector[i+2] -= alpha * mpTmpVectorQ[i+2];
        }
		// added by Jieling
		if (iteration >= 1000000)
		{
			cout << "	-> Endless while loop error!" << endl;
		}
    }
}

void
ModelCellsSpherical::FrictionMatrixTimesV( double *v_in, double *v_out )
{
    int problemVecSize = cells2->size();

    for ( int i = 0; i < problemVecSize; ++i )
    {
        int vindex0 = 3*i;
        int vindex1 = vindex0+1;
        int vindex2 = vindex0+2;

        ModelElement * element = cells2->element(i);

        if ( !element )
        {
            v_out[vindex0] = 0;
            v_out[vindex1] = 0;
            v_out[vindex2] = 0;

            continue;
        }

        // ecm friction on the diagonal of the friction matrix:
        if (element->mStatic)
        {
            v_out[vindex0] = 0;
            v_out[vindex1] = 0;
            v_out[vindex2] = 0;

            continue;
        }

		// added by Jieling
		if (element->mType == ModelElement::TypeVesselSphere)
		{
			v_out[vindex0] = 0;
			v_out[vindex1] = 0;
			v_out[vindex2] = 0;

			continue;
		}

        double ecmFrictionCoeff = gammaECM * element->freeSurfaceArea;

        v_out[vindex0] = ecmFrictionCoeff * v_in[vindex0];
        v_out[vindex1] = ecmFrictionCoeff * v_in[vindex1];

        if ( ! is2D )
            v_out[vindex2] = ecmFrictionCoeff * v_in[vindex2];
        else
            v_out[vindex2] = 0; // in the 2D case the z-component of v must be zero!

        // cell-cell/element-elemnt friction terms:
        for ( int k = 0; k < element->mNumContacts; ++k )
        {

            // get the index of the matrix of this ModelElement-ModelElement interaction:
            unsigned long offset = 9 * element->frictionMatrixEntryNumbers[k];

            // the zeroth element of this interaction's matrix
            double  *fMatrixLine = mpFrictionMatrices + offset;

            // the index of the interacting ModelElement -> index of velocity vector within v_in
            int j = 3 * element->mContacts[k];

            // A.(v_i - v_j)

            double v[3];
            v[0] = v_in[vindex0]  - v_in[j];
            v[1] = v_in[vindex1]  - v_in[j+1];
            v[2] = v_in[vindex2]  - v_in[j+2];

            v_out[vindex0 ] += fMatrixLine[0] * ( v[0] )
                               +  fMatrixLine[1] * ( v[1] )
                               +  fMatrixLine[2] * ( v[2] );
            fMatrixLine += 3;
            v_out[vindex1]  += fMatrixLine[0] * ( v[0] )
                               +  fMatrixLine[1] * ( v[1] )
                               +  fMatrixLine[2] * ( v[2] );
            if( ! is2D )
            {
                fMatrixLine += 3;
                v_out[vindex2]  += fMatrixLine[0] * ( v[0] )
                                   +  fMatrixLine[1] * ( v[1] )
                                   +  fMatrixLine[2] * ( v[2] );
            }
        }
    }
}

void
ModelCellsSpherical::InitialResidualAndPreconditioner( double * v_in, double *residual_out, double *preconditioner_out  )
{
    int problemVecSize = cells2->size();

    for ( int i = 0; i < problemVecSize; ++i )
    {
        int vindex0 = 3*i;
        int vindex1 = vindex0+1;
        int vindex2 = vindex0+2;

        ModelElement * element = cells2->element(i);

        if ( !element )
        {
            residual_out[vindex0] = 0;
            residual_out[vindex1] = 0;
            residual_out[vindex2] = 0;

            preconditioner_out[vindex0] = 0;
            preconditioner_out[vindex1] = 0;
            preconditioner_out[vindex2] = 0;

            continue;
        }

        // ecm friction on the diagonal of the friction matrix:
        if (element->mStatic)
        {
            preconditioner_out[vindex0] = 0;
            preconditioner_out[vindex1] = 0;
            preconditioner_out[vindex2] = 0;

            residual_out[vindex0] = 0;
            residual_out[vindex1] = 0;
            residual_out[vindex2] = 0;

            continue;
        }

		// added by Jieling
		if (element->mType == ModelElement::TypeVesselSphere)
		{
			preconditioner_out[vindex0] = 0;
			preconditioner_out[vindex1] = 0;
			preconditioner_out[vindex2] = 0;

			residual_out[vindex0] = 0;
			residual_out[vindex1] = 0;
			residual_out[vindex2] = 0;

			continue;
		}

		// revised by Jieling
		double gammaECMCoeff = gammaECM;
		if (element->mType == ModelElement::TypeCellLatticeNode)
			gammaECMCoeff = gammaECMLattice;

        double ecmFrictionCoeff = gammaECMCoeff * element->freeSurfaceArea;

        preconditioner_out[vindex0] = ecmFrictionCoeff;
        preconditioner_out[vindex1] = ecmFrictionCoeff;
        preconditioner_out[vindex2] = ecmFrictionCoeff;

        residual_out[vindex0] = element->directedForce.x - ecmFrictionCoeff * v_in[vindex0];
        residual_out[vindex1] = element->directedForce.y - ecmFrictionCoeff * v_in[vindex1];
        if ( ! is2D )
            residual_out[vindex2] = element->directedForce.z - ecmFrictionCoeff * v_in[vindex2];
        else
            residual_out[vindex2] = 0; // in the 2D case the z-component of v must be zero!

		/*if (element->mType == ModelElement::TypeECMSphere)
		{
			cout << "		ECM  sphere " << ((ModelElementECMSphere*)element)->mIndex << "_";
			cout << element->position.y << "_directedForce_x: " << element->directedForce.x;
			cout << ", residual_out_x: " << residual_out[vindex0];
			cout << ", preconditioner_out_x: " << preconditioner_out[vindex0] << endl;
		}*/

        // cell-cell friction terms:
        for ( int k = 0; k < element->mNumContacts; ++k )
        {

            // get the index of the matrix of this ModelElement-ModelElement interaction:
            unsigned long offset = 9 * element->frictionMatrixEntryNumbers[k];

            // the zeroth element of this interaction's matrix
            double  *fMatrixLine = mpFrictionMatrices + offset;

            // the index of the interacting ModelElement -> index of velocity vector within v_in
            int j = 3 * element->mContacts[k];

            preconditioner_out[vindex0] += *(fMatrixLine);
            preconditioner_out[vindex1] += *(fMatrixLine+4);
            preconditioner_out[vindex2] += *(fMatrixLine+8);

            int v[3];
            v[0] = v_in[vindex0]  - v_in[j];
            v[1] = v_in[vindex1]  - v_in[j+1];
            v[2] = v_in[vindex2]  - v_in[j+2];

            // A.(v_i - v_j)
            residual_out[vindex0 ] -= fMatrixLine[0] * ( v[0] )
                                      +  fMatrixLine[1] * ( v[1] )
                                      +  fMatrixLine[2] * ( v[2] );
            fMatrixLine += 3;
            residual_out[vindex1]  -= fMatrixLine[0] * ( v[0] )
                                      +  fMatrixLine[1] * ( v[1] )
                                      +  fMatrixLine[2] * ( v[2] );
            if( ! is2D )
            {
                fMatrixLine += 3;
                residual_out[vindex2]  -= fMatrixLine[0] * ( v[0] )
                                          +  fMatrixLine[1] * ( v[1] )
                                          +  fMatrixLine[2] * ( v[2] );
            }
        }

        preconditioner_out[vindex0] = 1/preconditioner_out[vindex0];
        preconditioner_out[vindex1] = 1/preconditioner_out[vindex1];
        preconditioner_out[vindex2] = 1/preconditioner_out[vindex2];
    }
}

//end Solve equation of motion

#pragma region Simulation control

void ModelCellsSpherical::SimulateInThread()
{
    if (biolink->getTimeInDays(time) <= simulateUntilDays)
    {
        if (enableObservation)
        {
            if (biolink->getTimeInDays(time) > nextObservationTime)
            {
                // Measure
                observe->ObserveModelMeasures();

                writeVesselSphereVTP();

                mpParaview->exec(this);

                // Update next observation time
                nextObservationTime = biolink->getTimeInDays(time) + observeEveryDays;

                ++mOutputCounter;
            }
        }

        Simulate();
    }
    else
    {
        enableSimulation = false;
    }
}

float ModelCellsSpherical::GetSimulationProgress()
{
    double simTime = simulateUntilDays - simulateFromDays;

    if (simTime > 0)
    {
        return (float)CSModelTools::ClampDouble(((biolink->getTimeInDays(time) - simulateFromDays) / simTime), 0, 1);
    }
    else
    {
        return 1.0;
    }
}

#pragma endregion

#pragma region Parameters

void ModelCellsSpherical::SwitchSubmodel(int index)
{
    if (index == 0)
    {
        is2D = true;
    }
    else if (index == 1)
    {
        is2D = false;
    }
    // Reset(); // has to be called afterwards!!
    //
    // comment by t1m:
    // Here, an enum would be handy, especially, if you plan to extend your group of submodels.
    // In this case, we can even equip Reset() with an argument:  void Reset( ModelCellsSpherical::SubModel ).
    //
}


//this method is only used when one click on "apply parameters to all existing cells" in the GUI
void ModelCellsSpherical::UpdateParametersForAllCells()
{
    UpdateParametersFromBioLink();

    for (unsigned int i=0; i<cells.size(); i++)
    {
        UpdateParametersForCell( cells[i] );
    }
}


void ModelCellsSpherical::UpdateParametersForCell(CellSpherical *cell)
{
    // t1m:  divisionRadius and initialRadius have to be initialized after a
    //       data read.  Will override any individual setting.

    #pragma region Size and growth specifics

    cell->divisionRadius = defaultDivisionCellRadius;
    cell->initialRadius  = defaultInitialCellRadius;
    cell->defaultDivisionDistance = defaultDivisionDistance;
    cell->UpdateDeltaRadius(timeStep);

    #pragma endregion

    #pragma region Biophysics

    cell->youngModulus = defaultYoungModulusCells;
    cell->poissonRatio = defaultPoissonRatioCells;
    cell->mQuiescentMin = mQuiescentMin;
    cell->mQuiescentMax = mQuiescentMax;
    cell->mUseDumbbell = mUseDumbbell;

    if ( cell->mType == ModelElement::TypeCellSphericalPolar )
    {
        // Angle of the adhesive pole regions (spherical caps) calculated from a
        // surface area percentage (mParmAdhesivesurface).
        // Two caps with the area 2*Pi*R*h, where R is the cell's radius and h the
        // height of the spherical cap.  With h=R(1-cos(alpha)), the area percentage
        // is
        //   Q = 4 Pi R^2 (1-cos(alpha)) / 4 Pi R^2 = (1-cos(alpha)),
        // and therefore
        //   alpha = acos( 1 - Q ), where Q is mParmAdhesivesurface/100.

        static_cast<CellSphericalPolar *>(cell)->SetPoleRegionAngle( std::acos(1. - .01*mParmAdhesiveSurface) );
        
        //! max rotation angle for a polar cell is given in degree, yet the angles are needed in radians, therefore a conversion is required: *M_PI/180
        static_cast<CellSphericalPolar *> (cell)->setMaxRotationAngleForMetropolis(mMaxRotationAngleForMetropolis*M_PI/180.);
    }

    #pragma endregion

    #pragma region Dumbbell specifics

    if ( mUseDumbbell )
    {
        cell->mDumbbellPeriod = this->mDefaultDumbbellPeriod;
        cell->mDumbbellInitialRadius =
            defaultDivisionCellRadius - (defaultDivisionCellRadius-defaultInitialCellRadius) * cell->mDumbbellPeriod/cell->cycleTime;
        if ( cell->mpDivisionPartner )
            cell->mElongationFactor = defaultDivisionDistance / cell->mDumbbellPeriod;
    }

    #pragma endregion
}


void ModelCellsSpherical::UpdateParametersFromBioLink()
{
    // Cycle time
    defaultCellCycleTime
    = biolink->scaleBiologyToInternal( biolink->cycletime_bio,
                                       BiologyLink::ScaleTime);

    // Cycle time deviation (Gaussian distribution)
    defaultCellCycleTimeStandardDeviation
    = biolink->scaleBiologyToInternal( biolink->cycletime_stddev_bio,
                                       BiologyLink::ScaleTime );

    //min and max pressure for proliferation
    this->mQuiescentMin
    = biolink->scaleBiologyToInternal( biolink->quiescentMinPressure_bio,
                                       BiologyLink::ScalePressure);
    this->mQuiescentMax
    = biolink->scaleBiologyToInternal( biolink->quiescentMaxPressure_bio,
                                       BiologyLink::ScalePressure);

    mCenterPullForce = biolink->scaleBiologyToInternal( mParmCenterPullForce, BiologyLink::ScaleForce );
    std::cout << "# Active Motion Center Pull:  " << mCenterPullForce << "/" << mParmCenterPullForce << " N\n";

	// added by Jieling
	mECMPullForce = mCenterPullForce * 10000.;
   
    // and all the rest...
    UpdateDimensionlessParametersFromBioLink();
}


void ModelCellsSpherical::UpdateDimensionlessParametersFromBioLink()
{
    defaultCellCycleTime = biolink->scaleBiologyToInternal(biolink->cycletime_bio, BiologyLink::ScaleTime);
    defaultCellCycleTimeStandardDeviation = biolink->scaleBiologyToInternal(biolink->cycletime_stddev_bio, BiologyLink::ScaleTime);

    // remark: the defaultcelldiameter don't need to be update, it is always 0.5 by definition (we set the characteristic length scale so that cell diameter is 1)

    // Young modulus
    defaultYoungModulusCells = biolink->getYoungModuleDimensionless(biolink->youngModulus_bioCells);
	defaultYoungModulusECM = biolink->getYoungModuleDimensionless(450); // added by Jieling

    // Poisson ratio
    defaultPoissonRatioCells = biolink->poisson_ratio_bioCells;
	defaultPoissonRatioECM = 0.4; // added by Jieling

    // Young modulus
    defaultYoungModulusSinusoids = biolink->getYoungModuleDimensionless(biolink->youngModulus_bioSinusoids);

    // Poisson ratio
    defaultPoissonRatioSinusoids = biolink->poisson_ratio_bioSinusoids;

    // Diffusion constant
    defaultDiffusionConstantCells
    = biolink->scaleBiologyToInternal( biolink->diffusion_constant_cells_bio,
                                       BiologyLink::ScaleDiffusivity );

    // Get single bond energy
    singleBondEnergy
    = biolink->scaleBiologyToInternal( biolink->single_bond_energy_bio,
                                       BiologyLink::ScaleEnergy );

    // Set cell-cell adhesion density
    adhesionDensity
    = biolink->scaleBiologyToInternal( biolink->adhesion_to_cells_bio,
                                       BiologyLink::ScaleAdhesionDensity );

    // Cell-cell friction coefficient
    gammaCellsPerpendicular
    = biolink->scaleBiologyToInternal( biolink->cell_gamma_perpendicular_bio,
                                       BiologyLink::ScaleFrictionCoefficient );

    gammaCellsParallel
    = biolink->scaleBiologyToInternal( biolink->cell_gamma_parallel_bio,
                                       BiologyLink::ScaleFrictionCoefficient );


    // Cell-ECM friction coefficient
    gammaECM
    = biolink->scaleBiologyToInternal( biolink->ecm_gamma_bio,
                                       BiologyLink::ScaleFrictionCoefficient );

	// added by Jieling, temporarily for lattice node
	gammaECMLattice
	= biolink->scaleBiologyToInternal( 1e+7,
									   BiologyLink::ScaleFrictionCoefficient );

    //set values depend on other variables
	
    #pragma endregion
}


void
ModelCellsSpherical::RegisterParameters()
{
    if ( mpParameters )
        return;

    mpParameters = new CSParameterContext( name );

    CSParameterChoice * choiceContactModel = new CSParameterChoice( mContactModelNames, 0 );

    mpParameters->setParameter("Contact Model",   CSParameter::Choice, choiceContactModel, "" );

    // simulation start time in relation to point 0.
    mpParameters->setParameter("Time Shift",   CSParameter::Double, &mStartTime, "s" );
	cout << "	mStartTime: " << mStartTime << endl; // tested by Jieling

    // length scale
    mpParameters->addParameter("Cell Diameter", CSParameter::Double, &biolink->length_scale, "m");

    // cycle time
    mpParameters->addParameter("Cell Cycle Time", CSParameter::Double, &biolink->cycletime_bio, "s");

    // std deviation of cycle times
    mpParameters->addParameter("Cycle Time SD", CSParameter::Double, &biolink->cycletime_stddev_bio, "s");

    // If to use a dumbbell of two spheres for a period of time before the cell cleaves into two
    mpParameters->addParameter("Use Dumbbell Division", CSParameter::Bool, &mUseDumbbell, "");

    // dumbbell time
    mpParameters->addParameter("Dumbbell time", CSParameter::Double, &this->mDefaultDumbbellPeriod, "s");

    // Young Modulus Cells
    mpParameters->addParameter("Young Modulus Cells", CSParameter::Double, &biolink->youngModulus_bioCells, "Pa");

    // Poisson Ratio Cells
    mpParameters->addParameter("Poisson Ratio Cells", CSParameter::Double, &biolink->poisson_ratio_bioCells, "");

    // Young Modulus Sinusoids
    mpParameters->addParameter("Young Modulus Sinusoids", CSParameter::Double, &biolink->youngModulus_bioSinusoids, "Pa");

    // Poisson Ratio Sinusoids
    mpParameters->addParameter("Poisson Ratio Sinusoids", CSParameter::Double, &biolink->poisson_ratio_bioSinusoids, "");

    // diffusion constant for cells
    mpParameters->addParameter("Diffusion Constant", CSParameter::Double, &biolink->diffusion_constant_cells_bio, "m^2/s" );

    // cell-cell adhesion
    mpParameters->addParameter("Cell-Cell Adhesion", CSParameter::Double, &biolink->adhesion_to_cells_bio, "m^-2");
    
    // If cells should be allocated as CellSphericalPolar or simply CellSpherical objects
    mpParameters->addParameter("Use Cell Polarity", CSParameter::Bool, &mUsePolarCells, "");

    mpParameters->addParameter("Percentage of Adhesive Surface", CSParameter::Double, &mParmAdhesiveSurface, "%");
    
    mpParameters->addParameter("Maximum Angle of Rotation for Metropolis", CSParameter::Double, &mMaxRotationAngleForMetropolis, "degree");

    // cell-cell gamma for shear friction
    mpParameters->addParameter("Cell-Cell gamma for shear friction", CSParameter::Double, &biolink->cell_gamma_perpendicular_bio, "Ns/m^3");

    // cell-cell gamma for collision friction
    mpParameters->addParameter("Cell-Cell gamma for normal friction", CSParameter::Double, &biolink->cell_gamma_parallel_bio, "Ns/m^3");

    // Cell-vs-Extracellular Matrix gamma
    mpParameters->addParameter("Cell-Medium gamma", CSParameter::Double, &biolink->ecm_gamma_bio, "Ns/m^3");

    // for ScenarioEmbeddingmedium:
    mpParameters->addParameter("Cell Population Radius (Embedding Medium)", CSParameter::Double, &mPopulationRadius, "Cell Diameters");

    mpParameters->addParameter("Cell Distance in Initial Population (Embedding Medium)", CSParameter::Double, &mPopulationInitialDistance, "Cell Diameters");

    CSParameterChoice * choiceLobuleShape = new CSParameterChoice( mLobuleShapeNames, 0 );
    mpParameters->setParameter("Lobule Shape",   CSParameter::Choice, choiceLobuleShape, "" );
    mpParameters->addParameter("Lobule shape : radius" , CSParameter::Double , &this->mlobule_radius, "");
    mpParameters->addParameter("Lobule shape : height" , CSParameter::Double , &this->mlobule_height, "");

    mpParameters->addParameter("Blood Vessel Network" , CSParameter::Bool , &this->mBloodVesselNetwork, "");
    mpParameters->addParameter("Path to Blood Vessel Network", CSParameter::FileName, &this->mBloodVesselNetworkPath, "");

    mpParameters->addParameter("Read Cell Position", CSParameter::Bool, &this->mReadCells, "");

    mpParameters->addParameter("Pressure Quiescent->Proliferating" , CSParameter::Double , &biolink->quiescentMinPressure_bio, "Pa");
    mpParameters->addParameter("Pressure Proliferating->Quiescent" , CSParameter::Double , &biolink->quiescentMaxPressure_bio, "Pa");

    mpParameters->addParameter("dynamic timestep",         CSParameter::Bool,   &this->mUseDynamicalTimeSteps,         "" );
    mpParameters->addParameter("time step",                CSParameter::Double, &this->timeStep,                       "s");
    mpParameters->addParameter("time step scaling factor", CSParameter::Int,    &this->mTimeScalingFactor,             "" );
    mpParameters->addParameter("max Displacement",         CSParameter::Double, &this->mMaximumDisplacementCellRadius, "" );
    mpParameters->addParameter("dynamic timestep min",     CSParameter::Double, &this->mDynamicTimeStepMin,            "s");
    mpParameters->addParameter("dynamic timestep max",     CSParameter::Double, &this->mDynamicTimeStepMax,            "s");

    mpParameters->addParameter("output path", CSParameter::String, &this->mOutputPath, "");
    mpParameters->addParameter("output suffix", CSParameter::String, &this->mOutputPrefix, "");
    mpParameters->addParameter("Lesion Radius Due to Intoxication", CSParameter::Double, &mKillZoneRadius,  "&mum");
    mpParameters->addParameter("Active Motion toward Central Vein", CSParameter::Double, &mParmCenterPullForce, "N" );
    mpParameters->addParameter("Scale Factor for Data to Proliferation Rates (experimental)", CSParameter::Double, &mScaleProliferationRates, "" );
}


void
ModelCellsSpherical::InitParameters( CSParameterContext *parms )
{
    if (!parms)
    {
        if ( !mpParameters )
            RegisterParameters();

        DefaultParameters( mpParameters );
        return;
    }

    std::vector<CSParameter *> parameters = parms->getParameters();

    std::vector<CSParameter *>::const_iterator parmsIt;

    for ( parmsIt=parameters.begin(); parmsIt!=parameters.end(); ++parmsIt )
    {
        if ( (*parmsIt)->name() == "Contact Model" )
        {
            CSParameter * foundParm = mpParameters->findParameter("Contact Model");
            ((CSParameterChoice *)foundParm->dataPointer())->setCurrentIndex(((CSParameterChoice *)(*parmsIt)->dataPointer())->currentIndex());
        }
        else if ( (*parmsIt)->name() == "Time Shift" )
            mStartTime = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Cell Diameter" )
            biolink->length_scale = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Cell Cycle Time" )
            biolink->cycletime_bio = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Cycle Time SD" )
            biolink->cycletime_stddev_bio = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Use Dumbbell Division" )
            this->mUseDumbbell = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Dumbbell time" )
            this->mDefaultDumbbellPeriod = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Young Modulus Cells" )
            biolink->youngModulus_bioCells = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Poisson Ratio Cells" )
            biolink->poisson_ratio_bioCells = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Young Modulus Sinusoids" )
            biolink->youngModulus_bioSinusoids = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Poisson Ratio Sinusoids" )
            biolink->poisson_ratio_bioSinusoids = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Diffusion Constant" )
            biolink->diffusion_constant_cells_bio = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Cell-Cell Adhesion" )
            biolink->adhesion_to_cells_bio = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Use Cell Polarity" )
            this->mUsePolarCells = *(bool*)(*parmsIt)->dataPointer();
        else if ( (*parmsIt)->name() == "Percentage of Adhesive Surface" )
            this->mParmAdhesiveSurface = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Maximum Angle of Rotation for Metropolis" )
            this->mMaxRotationAngleForMetropolis = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Cell-Cell gamma for shear friction" )
            biolink->cell_gamma_perpendicular_bio = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Cell-Cell gamma for normal friction" )
            biolink->cell_gamma_parallel_bio = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Cell-Medium gamma" )
            biolink->ecm_gamma_bio = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Cell Population Radius (Embedding Medium)")
            mPopulationRadius = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Cell Distance in Initial Population (Embedding Medium)")
            mPopulationInitialDistance = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Lobule shape : radius" )
            this->mlobule_radius = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Lobule shape : height" )
            this->mlobule_height = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Blood Vessel Network" )
            this->mBloodVesselNetwork = *(bool *)(*parmsIt)->dataPointer();
        else if ( (*parmsIt)->name() == "Path to Blood Vessel Network" )
            this->mBloodVesselNetworkPath = (*parmsIt)->dataString();
        else if ( (*parmsIt)->name() == "Read Cell Position" )
            this->mReadCells= *(bool *)(*parmsIt)->dataPointer();
        else if ( (*parmsIt)->name() == "Pressure Quiescent->Proliferating" )
          biolink->quiescentMinPressure_bio = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Pressure Proliferating->Quiescent" )
          biolink->quiescentMaxPressure_bio = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "time step" )
            this->timeStep= (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "dynamic timestep" )
            this->mUseDynamicalTimeSteps= (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "time step scaling factor" )
            this->mTimeScalingFactor= (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "max Displacement" )
            this->mMaximumDisplacementCellRadius= (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "dynamic timestep min" )
            this->mDynamicTimeStepMin= (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "dynamic timestep max" )
            this->mDynamicTimeStepMax= (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "output suffix" )
            this->mOutputPrefix = (*parmsIt)->dataString();
        else if ( (*parmsIt)->name() == "output path" )
            this->mOutputPath = (*parmsIt)->dataString();
        else if ( (*parmsIt)->name() == "Lobule Shape" )
        {
            CSParameter * foundParm = mpParameters->findParameter("Lobule Shape");
            ((CSParameterChoice *)foundParm->dataPointer())->setCurrentIndex(((CSParameterChoice *)(*parmsIt)->dataPointer())->currentIndex());
        }
        else if ( (*parmsIt)->name() == "Lesion Radius Due to Intoxication" )
            this->mKillZoneRadius = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Active Motion toward Central Vein" )
            this->mParmCenterPullForce = (*parmsIt)->value();
        else if ( (*parmsIt)->name() == "Scale Factor for Data to Proliferation Rates (experimental)" )
            this->mScaleProliferationRates = (*parmsIt)->value();
        else{
          std::cerr << "ModelCellsSpherical::InitParameters:  Unknown parameter given:" << "\t";
          std::cerr << (*parmsIt)->name() << std::endl;
        }
    }

}


void
ModelCellsSpherical::DefaultParameters( CSParameterContext *emptyContext )
{
    if ( !emptyContext )
        return;

    // what interaction model to use:
    CSParameterChoice * choiceContactModel = new CSParameterChoice( mContactModelNames, 0 );
    emptyContext->setParameter("Contact Model",   CSParameter::Choice, choiceContactModel, "" );

    // Simulation start time in relation to point 0
    emptyContext->setParameter( "Time Shift",     CSParameter::Double, new double( 0. ), "s");

    // length scale
    emptyContext->setParameter( "Cell Diameter",      CSParameter::Double, new double( 1.e-5 ), "m");

    // cycle time
    emptyContext->setParameter( "Cell Cycle Time",    CSParameter::Double, new double( 10.e+21 ), "s");

    // std deviation of cycle times
    emptyContext->setParameter( "Cycle Time SD",      CSParameter::Double, new double( 7200. ), "s");

    // use an extending dumbbell for a period of time before the cell cleaves into two
    emptyContext->setParameter( "Use Dumbbell Division",    CSParameter::Bool, new bool(false), "" );

    //dumbbell time period
    emptyContext->addParameter( "Dumbbell time", CSParameter::Double, new double( 1800. ), "s");

    // Young Modulus
    emptyContext->setParameter( "Young Modulus Cells",      CSParameter::Double, new double( 450. ), "Pa");

    // Poisson Ratio
    emptyContext->setParameter( "Poisson Ratio Cells",      CSParameter::Double, new double( .4 ), "");

    // Young Modulus
    emptyContext->setParameter( "Young Modulus Sinusoids",      CSParameter::Double, new double( 600. ), "Pa");

    // Poisson Ratio
    emptyContext->setParameter( "Poisson Ratio Sinusoids",      CSParameter::Double, new double( .4 ), "");

    // diffusion constant for cells
    emptyContext->setParameter( "Diffusion Constant", CSParameter::Double, new double( 1.e-20 ), "m^2/s" );

    // cell-cell adhesion
    emptyContext->setParameter( "Cell-Cell Adhesion", CSParameter::Double, new double( 1.e15 ), "m^-2");
    
    // use simple CellSphericals or CellSphericalPolar
    emptyContext->setParameter( "Use Cell Polarity", CSParameter::Bool, new bool(true), "" );

    emptyContext->setParameter("Percentage of Adhesive Surface", CSParameter::Double, new double(3.4), "%");
    
    emptyContext->setParameter("Maximum Angle of Rotation for Metropolis", CSParameter::Double, new double(10), "degree");

    // cell-cell gamma
    emptyContext->setParameter( "Cell-Cell gamma for shear friction",  CSParameter::Double, new double( 1.e8 ), "Ns/m^3");
    emptyContext->setParameter( "Cell-Cell gamma for normal friction", CSParameter::Double, new double( 1.e8 ), "Ns/m^3");

    // Cell-vs-Extracellular Matrix gamma
    emptyContext->setParameter( "Cell-Medium gamma",     CSParameter::Double, new double( 1.e8 ), "Ns/m^3");

    CSParameterChoice * choiceLobuleShape = new CSParameterChoice( mLobuleShapeNames, 0 );
        // for ScenarioEmbeddingmedium:
    emptyContext->setParameter("Cell Population Radius (Embedding Medium)", CSParameter::Double, new double(20), "Cell Diameters");

    emptyContext->setParameter("Cell Distance in Initial Population (Embedding Medium)", CSParameter::Double, new double(.95), "Cell Diameters");


    emptyContext->setParameter("Lobule Shape",   CSParameter::Choice, choiceLobuleShape, "" );

    emptyContext->setParameter( "Lobule shape : radius",     CSParameter::Double, new double( 10 ), "");
    emptyContext->setParameter( "Lobule shape : height",     CSParameter::Double, new double( 10 ), "");

    emptyContext->setParameter( "Blood Vessel Network",      CSParameter::Bool,   new bool( true ), "");
    emptyContext->setParameter( "Path to Blood Vessel Network",     CSParameter::FileName, new std::string("/Users/nboissie/Documents/_Work/5_Soft/CellSys/CellSys-dev/tests/testDiscvor2.mxf"), "");

    emptyContext->setParameter( "Read Cell Position",     CSParameter::Bool, new bool(true), "");

    emptyContext->setParameter( "Pressure Quiescent->Proliferating",     CSParameter::Double, new double( 80 ), "Pa");
    emptyContext->setParameter( "Pressure Proliferating->Quiescent",     CSParameter::Double, new double( 100 ), "Pa");

    emptyContext->setParameter("time step",                 CSParameter::Double, new double( 1 ),    "s");
    emptyContext->setParameter("dynamic timestep",          CSParameter::Bool,   new bool(true),      "" );
    emptyContext->setParameter("time step scaling factor",  CSParameter::Int,    new unsigned int(2), "" );
    emptyContext->setParameter("max Displacement",          CSParameter::Double, new double( 0.1 ),   "" );
    emptyContext->setParameter("dynamic timestep min",      CSParameter::Double, new double( 1.e-8 ),    "s");
    emptyContext->setParameter("dynamic timestep max",      CSParameter::Double, new double( 256. ),  "s");


    emptyContext->addParameter("output path", CSParameter::DirName, new std::string("/Users/nboissie/Documents/_Work/5_Soft/CellSys/CellSys-dev/tests"), "");
    emptyContext->addParameter("output suffix", CSParameter::String, new std::string("sphericalCells_sim1000"), "");
    emptyContext->setParameter("Lesion Radius Due to Intoxication", CSParameter::Double, new double( 0.00001 ),  "um");
    emptyContext->setParameter("Active Motion toward Central Vein", CSParameter::Double, new double( 1e-14), "N" );
    emptyContext->setParameter("Scale Factor for Data to Proliferation Rates (experimental)", CSParameter::Double, new double( 1. ), "" );
}


CSParameterContext *
ModelCellsSpherical::GetParameters( std::string contextName )
{
    if ( !mpParameters)
    {
        RegisterParameters();
        return mpParameters;
    }

    if (!contextName.size())
        return mpParameters;

    if ( contextName == name )
        return mpParameters;

    return mpParameters->findContext(contextName);
}

#pragma endregion


#pragma region Scenarios

void
ModelCellsSpherical::InitScenarioLobuleRegeneration( bool createInitialConditions )
{
	//std::string outputFileName1 = "debugFile.txt";
	//std::ofstream FO;
	//FO.open(outputFileName1.c_str(), std::ios::out | std::ios::app);
	//FO << "ModelCellsSpherical.cpp / InitScenarioLobuleRegeneration()	starts" << endl;
    // if ( createInitialConditions )
    // {
        double radius = biolink->scaleBiologyToInternal( 1e-6 * mKillZoneRadius, BiologyLink::ScaleLength );

        //mark all cells near central vein to delete
        for( unsigned int i = 0 ; i < this->cells.size() ; i++)
        {
            double dist = std::sqrt(this->cells[i]->position.x*this->cells[i]->position.x+this->cells[i]->position.y*this->cells[i]->position.y);

            if( dist < radius )
            {
                this->mCellsToDelete.push_back(this->cells[i]);
            }

            // on the fly init all cells proliferation control to
            // Cell::ProliferationControlCustom
            cells[i]->mProliferationControl = Cell::ProliferationControlCustom;
        }
    // }
    // else
    // {
    //     // mCellsToDelete's contents are to be read from data file
    //     for( unsigned int i = 0 ; i < this->cells.size() ; i++)
    //     {
    //         cells[i]->mProliferationControl = Cell::ProliferationControlCustom;
    //     }
    // }

    // Cells die with a constant rate within 24h.  The rate is
    // calculated with the chosen number of doomed cells.

    // Cells are removed every 5min.  The death rate therefore is
    // r=n_k(T_0)/k*DeltaT and the probability for any cell to die
    // within a period of DeltaT is
    //   r*DeltaT/n_k(t) = n_k(T_0)/(n_k(t)*k),
    // where n_k(t) is the number of remaining 'doomed' cells at the
    // given time point.
    mKillInterval = 360.; // 8640.
    mKillPeriod = 3600.; // 86400.
    mKillEnd = mKillPeriod; // start at time=0.
    mKillCellsToDiePerInterval = mCellsToDelete.size() * mKillInterval / mKillPeriod;

    mNextKillTime = mKillInterval;

	// added by Jieling
	while (!mCellsToDelete.empty())
	{
		RemoveCell(mCellsToDelete.back());
		mCellsToDelete.pop_back();
	}

    // Proliferation is controlled by the (hard-coded) mvProliferationData.
    // Since these data were acquired as percentage of cells, per day, in
    // S-Phase for each layer in a lobule, they have to be transformed to 1/100
    // per second and multiplied by the update interval.
    mvProliferationProfiles.clear();

    std::deque< std::pair<double, std::vector<double> > >::iterator prolifDataIt;
    for ( prolifDataIt=mvProliferationData.begin(); prolifDataIt!=mvProliferationData.end(); ++prolifDataIt )
    {
        // multiply rates in data by a scaling factor (data is from BrdU stained
        // cells - i.e. cells in S-phase, not those which just passed G1).
        // Also, scale rates to 1/s (from 1/d)
        std::vector<double> scaledRates;
        std::vector<double>::const_iterator rateIterator;
        for ( rateIterator=prolifDataIt->second.begin(); rateIterator!= prolifDataIt->second.end(); ++rateIterator )
        {
            scaledRates.push_back( mScaleProliferationRates * (*rateIterator) );
            std::cout << mScaleProliferationRates * (*rateIterator) << "  ";
        }
                std::cout << std::endl;
        mvProliferationProfiles.push_back( std::pair< double, std::vector<double> >(86400*prolifDataIt->first, scaledRates ) );
    }

    // Do we need to parametrize this?
    mProliferationUpdateInterval = 300.;

    // If we loaded from data and the saved time is not the original start time,
    // we need to search forward through the proliferation data to find the nextProliferationRateUpdate
    if ( time != mStartTime )
    {
        if ( mvProliferationProfiles.size() != 0 )
        {
            std::vector<double> lastProliferationRateUpdate;
            double progressedTime = time - mStartTime;
            while ( mvProliferationProfiles.front().first < progressedTime )
            {
                lastProliferationRateUpdate = mvProliferationProfiles.front().second;
                mvProliferationProfiles.pop_front();
            }

            // get the current rates
            unsigned int vectorSize = (lastProliferationRateUpdate.size()<=14) ? lastProliferationRateUpdate.size() : 14;

            for ( unsigned int i=0; i<vectorSize; ++i )
                mpProliferationRatesPerLayer[i] = lastProliferationRateUpdate[i];

            // ToDo: set nextObservablesUpdate accordingly in GrowAndDivide().
            nextObservablesUpdate = mStartTime + mProliferationUpdateInterval*std::ceil((time-mStartTime)/mProliferationUpdateInterval);
        }
    }
    else
    {
        nextObservablesUpdate = mStartTime;
    }

    nextProliferationRateUpdate = /*mStartTime +*/ mvProliferationProfiles.front().first;

    nextLesionDiscretization = time;

    if ( mpLesionVoxelization )
        delete mpLesionVoxelization;

   // mpLesionVoxelization = new LobuleLesionVoxelization(this, .5);
	//FO << "ModelCellsSpherical.cpp / InitScenarioLobuleRegeneration()	ends" << endl;
	//FO.close();
}

#pragma endregion


#pragma region Save and load model

#  pragma region XML

ModelCellsSpherical *
ModelCellsSpherical::createFromXML( QXmlStreamReader * xml,
                                    std::stringstream & /* errors */,
                                    std::stringstream & warnings )
{
    Q_ASSERT( xml->name() == "Model" && xml->isStartElement() );

    ModelCellsSpherical * model = new ModelCellsSpherical();

    // to be on the safe side query 'dimensions' again (even though done in Model::createFromXML())
    unsigned int dimensions = xml->attributes().value("dimensions").toString().toUInt();
    model->is2D = (dimensions == 2) ? true : false;

    while ( xml->readNextStartElement() )
    {
        CSParameterContextTemporary *parms = NULL;

        QStringRef elementName = xml->name();

        if ( elementName == "parameters" )
        {
            // load defaults
            CSParameterContextTemporary defaultParms("");
            DefaultParameters( &defaultParms );
            model->InitParameters( &defaultParms );

            parms = (CSParameterContextTemporary *) CSParameterContext::createFromXML( xml );
            if (parms)
            {
                model->InitParameters( (CSParameterContext *)parms );
                delete parms;
            }
        }
        else if ( elementName == "simulation")
        {
            parms = (CSParameterContextTemporary *) CSParameterContext::createFromXML( xml );
            if ( parms )
            {
                CSParameter * currentTime = parms->findParameter( "Current time" );
                if (currentTime)
                    model->time = currentTime->value();
                CSParameter * untilDays = parms->findParameter( "Simulation time" );
				if (untilDays)
				{
					model->simulateUntilDays = untilDays->value();
					// added by Jieling
					model->mECMNetworkTime += model->simulateUntilDays * 86400;
				}
                CSParameter * observe = parms->findParameter( "Observation interval" );
                model->enableObservation = false;
                if ( observe )
                {
                    model->enableObservation = true;
                    model->observeEveryDays = observe->value();
                }
                delete parms;
            }
        }
        else if ( elementName == "scenario" )
        {
            std::string scenario = xml->attributes().value("type").toString().toStdString();
            std::string scenarioSearch;
            int i;
            for ( i=0, scenarioSearch = mScenarioNames[i];
                  scenarioSearch != "";
                  ++i, scenarioSearch=mScenarioNames[i] )
                if ( scenarioSearch == scenario )
                { model->mScenario = (ModelCellsSpherical::Scenario) i; break; }

            if ( model->mScenario == 0)
                warnings << "Warning: "
                         << "Scenario section could not be parsed:" << std::endl
                         << "Scenario of Name \"" << scenario << "\" not known" << std::endl;

            parms = (CSParameterContextTemporary *) CSParameterContext::createFromXML( xml );
        }
        else if ( elementName == "cells" )
        {
            // legacy xml format element.
            xml->skipCurrentElement();
        }
        else
        {
            warnings << "Warning: "
                     << "Unknown element in XML description: \""
                     << xml->name().toString().toStdString()
                     << "\"" << std::endl;
            xml->skipCurrentElement();
        }
    }

    return model;
}


// a shortcut for writing Parameter entries directly w/o a CSParameterContext
#define WRITE_XML_DOUBLE_PARM( name, value, unit )              \
    {                                                           \
    xml->writeStartElement( "Parameter" );                      \
    xml->writeAttribute( "name", name );                        \
    xml->writeAttribute( "type", "Double" );                    \
    xml->writeAttribute( "value", QString("%1").arg(value) );   \
    xml->writeAttribute( "unit", unit );                        \
    xml->writeEndElement(); }{}

void
ModelCellsSpherical::writeXML(QXmlStreamWriter *xml) const
{
#undef QT_NO_CAST_FROM_ASCII

    xml->writeStartElement( "Model" );
    xml->writeAttribute( "type", xmlType.c_str() );
    xml->writeAttribute( "name", name.c_str() );
    xml->writeAttribute( "dimensions", is2D?"2":"3" );

    xml->writeStartElement( "simulation" );
    WRITE_XML_DOUBLE_PARM( "Current time", time, "s" );
    WRITE_XML_DOUBLE_PARM( "Simulation time", simulateUntilDays, "d" );
    if ( enableObservation )
        WRITE_XML_DOUBLE_PARM( "Observation interval", observeEveryDays, "d" );
    xml->writeEndElement();

    if ( mScenario != ScenarioReadFromData )
    {
        xml->writeStartElement( "scenario" );
        xml->writeAttribute( "type", mScenarioNames[mScenario].c_str() );
        xml->writeEndElement();
    }

    // parameters:
    xml->writeStartElement( "parameters" );

    if ( mpParameters )
        mpParameters->writeXML(xml, true);

    xml->writeEndElement();

    // cells:  left here for debugging purposes only - OBSOLETE
    xml->writeStartElement("cells");

    std::vector<CellSpherical *>::const_iterator cell = cells.begin();

    while ( cell != cells.end() )
    {
        (*cell)->writeXML( xml );
        ++cell;
    }

    xml->writeEndElement(); // cells

    xml->writeEndElement();
}

#  pragma endregion


#  pragma region HDF5

void
ModelCellsSpherical::readModelData( H5::H5File * inputFile,
                                    std::stringstream & errors,
                                    std::stringstream & warnings )
{
    cells.clear();

    if ( cells2 )
        delete cells2;
    cells2 = new BoundingBoxList( is2D ? 2 : 3 );


    hsize_t numElements;  // used for reading data extents and allocating space


    std::string groupString = "/" + name;


    try
    {
        // read in the state of the Marsaglia random number generator:
        std::string dataSetRNGString = groupString + "/MarsagliaRNG";

        H5::DataSet dataSetRNG =
            inputFile->openDataSet( dataSetRNGString );

        H5::CompType rngDataType;
        Random::HDF5DataFormat( rngDataType );

        dataSetRNG.read( &mRandom, rngDataType );
    }
    catch (...)
    {
        std::cout << "No Random Number Generator data in Data file.  Using default seed\n";
        mRandom.Init();
    }



    // reading in barriers:
    std::string dataSetBarriersString =
        groupString + "/ModelElementBarrierTriangles";

    try
    {
        H5::DataSet dataSetBarriers =
            inputFile->openDataSet(dataSetBarriersString);

        H5::CompType dataTypeBarriers =
            dataSetBarriers.getCompType();

        H5::CompType readTypeBarrier =
            ModelElementBarrierTriangle::ParseHDF5DataFormat( dataTypeBarriers,
                    errors,
                    warnings );

        H5::DataSpace outputSpaceBarriers = dataSetBarriers.getSpace();

        outputSpaceBarriers.getSimpleExtentDims( &numElements, NULL );

        ModelElementBarrierTriangle * readBarriers =
            new ModelElementBarrierTriangle[numElements];

        dataSetBarriers.read( readBarriers, readTypeBarrier );

        for ( unsigned int i=0; i<numElements; ++i )
        {
            ModelElementBarrierTriangle * barrier = readBarriers +i;
            barrier->setNormalVector();
            barrier->setBoundingBox( .1 );
            mBarriers.push_back(barrier);
            cells2->add(barrier);
            mpArena->addObject(barrier->GLObject());
        }
    }
    catch (...)
    {
        std::cout << "No data set found named \"" << dataSetBarriersString << "\".\n";
    }

    // reading in vessel graph / vessel spheres:
    std::string dataSetVesselSpheresString =
        groupString + "/VesselGraph/ModelElementVesselSpheres";

    try
    {
        H5::DataSet dataSetVesselSpheres =
            inputFile->openDataSet(dataSetVesselSpheresString);

        H5::CompType dataTypeVesselSphere =
            dataSetVesselSpheres.getCompType();

        H5::CompType readTypeVesselSphere =
            ModelElementVesselSphere::ParseHDF5DataFormat( dataTypeVesselSphere,
                    errors,
                    warnings );

        H5::DataSpace outputSpaceVesselSpheres = dataSetVesselSpheres.getSpace();

        outputSpaceVesselSpheres.getSimpleExtentDims( &numElements, NULL );

        ModelElementVesselSphere * readVesselSpheres =
            new ModelElementVesselSphere[numElements];

        dataSetVesselSpheres.read( readVesselSpheres, readTypeVesselSphere );

        if ( mpGraphBloodVesselNetwork )
            delete mpGraphBloodVesselNetwork;

        GraphSphere * network = new GraphSphere();

        for ( unsigned int i=0; i<numElements; ++i )
        {
            ModelElementVesselSphere * vesselSphere = readVesselSpheres +i;
			// added by Jieling
			vesselSphere->mIndex = i;
			vesselSphere->color.red = 0.;
			vesselSphere->color.green = 0.;
			vesselSphere->color.blue = 1.;
            network->mvNode.push_back(vesselSphere);
        }

        struct vesselConnection { hsize_t indexStartNode; hsize_t indexEndNode; };

        // reading in vessel sphere connections
        std::string dataSetVesselConnectionString =
            groupString + "/VesselGraph/VesselSphereConnections";

        H5::DataSet dataSetVesselConnections =
            inputFile->openDataSet(dataSetVesselConnectionString);

        H5::CompType dataTypeVesselConnection( sizeof(vesselConnection) );
        dataTypeVesselConnection.insertMember( "startIndex", HOFFSET(vesselConnection, indexStartNode), H5::PredType::STD_U64LE);
        dataTypeVesselConnection.insertMember( "endIndex", HOFFSET(vesselConnection, indexEndNode), H5::PredType::STD_U64LE);

        H5::DataSpace outputSpaceVesselConnections = dataSetVesselConnections.getSpace();
        outputSpaceVesselConnections.getSimpleExtentDims( &numElements, NULL );

        vesselConnection * readVesselConnections =
            new vesselConnection[numElements];

        dataSetVesselConnections.read( readVesselConnections, dataTypeVesselConnection );

        for (unsigned long i=0; i<numElements; ++i)
        {
            CSGraphEdge * edge = new CSGraphEdge();
            edge->mpStart = network->mvNode[readVesselConnections[i].indexStartNode];
            edge->mpEnd   = network->mvNode[readVesselConnections[i].indexEndNode];
            edge->mIndex  = i;
            edge->mRadius_start = edge->mpStart->mRadius;
            edge->mRadius_end   = edge->mpEnd->mRadius;
            network->mvEdge.push_back(edge);
        }

        AddBloodVesselNetwork( network );
    }
    catch (...)
    {
        std::cout << "No data set found named \"" << dataSetVesselSpheresString << "\".\n";
    }

	// added by Jieling
	// reading in ECM graph
	std::string dataSetECMSpheresString =
		groupString + "/ECMGraph/ModelElementLatticeNodes";
	try
	{
		//std::cout << "	-> test read ECMSpheres " << std::endl;
		H5::DataSet dataSetECMSpheres =
			inputFile->openDataSet(dataSetECMSpheresString);

		H5::CompType dataTypeECMSphere =
			dataSetECMSpheres.getCompType();

		H5::CompType readTypeECMSphere =
			ModelElementLatticeNode::ParseHDF5DataFormat(dataTypeECMSphere,
				errors,
				warnings);

		H5::DataSpace outputSpaceECMSpheres = dataSetECMSpheres.getSpace();

		outputSpaceECMSpheres.getSimpleExtentDims(&numElements, NULL);

		ModelElementLatticeNode * readECMSpheres =
			new ModelElementLatticeNode[numElements];

		dataSetECMSpheres.read(readECMSpheres, readTypeECMSphere);

		if (mpCollagenNetwork)
			delete mpCollagenNetwork;

		mpCollagenNetwork = new Lattice();
		mpCollagenNetwork->setArena(mpArena);

		for (unsigned int i = 0; i < numElements;++i)
		{
			ModelElementLatticeNode * ECMSphere = readECMSpheres + i;
			ECMSphere->mIndex = i;
			ECMSphere->mYoung = defaultYoungModulusECM;
			ECMSphere->mPoisson = defaultPoissonRatioECM; // temporarily the same
			mpCollagenNetwork->addNode(ECMSphere);
		}
		mpCollagenNetwork->setBoundingBoxList(cells2);

		struct ECMConnection {
			hsize_t indexStartNode; hsize_t indexEndNode;
			hsize_t StartNodeFlag; hsize_t EndNodeFlag;
			hsize_t StartNodeAlongFlag; hsize_t EndNodeAlongFlag;
		};

		// reading in ECM sphere connections
		std::string dataSetECMConnectionString =
			groupString + "/ECMGraph/ECMSphereConnections";

		H5::DataSet dataSetECMConnections =
			inputFile->openDataSet(dataSetECMConnectionString);

		H5::CompType dataTypeECMConnection(sizeof(ECMConnection));
		dataTypeECMConnection.insertMember("startIndex", HOFFSET(ECMConnection, indexStartNode), H5::PredType::STD_U64LE);
		dataTypeECMConnection.insertMember("endIndex", HOFFSET(ECMConnection, indexEndNode), H5::PredType::STD_U64LE);
		dataTypeECMConnection.insertMember("startFlag", HOFFSET(ECMConnection, StartNodeFlag), H5::PredType::STD_U64LE);
		dataTypeECMConnection.insertMember("endFlag", HOFFSET(ECMConnection, EndNodeFlag), H5::PredType::STD_U64LE);
		dataTypeECMConnection.insertMember("startAlongFlag", HOFFSET(ECMConnection, StartNodeAlongFlag), H5::PredType::STD_U64LE);
		dataTypeECMConnection.insertMember("endAlongFlag", HOFFSET(ECMConnection, EndNodeAlongFlag), H5::PredType::STD_U64LE);

		H5::DataSpace outputSpaceECMConnections = dataSetECMConnections.getSpace();
		outputSpaceECMConnections.getSimpleExtentDims(&numElements, NULL);

		ECMConnection * readECMConnections =
			new ECMConnection[numElements];

		dataSetECMConnections.read(readECMConnections, dataTypeECMConnection);

		for (unsigned long i = 0; i < numElements; ++i)
		{
			ModelElementLatticeNode *sStart = mpCollagenNetwork->getNodes().at(readECMConnections[i].indexStartNode);
			ModelElementLatticeNode *sEnd = mpCollagenNetwork->getNodes().at(readECMConnections[i].indexEndNode);
			LinearSpring * edge = new LinearSpring(sStart, sEnd);
			mpCollagenNetwork->getSprings().push_back(edge);
			sStart->addNeighbor(sEnd);
			sEnd->addNeighbor(sStart);
			sStart->addSpring(edge);
			sEnd->addSpring(edge);
			sStart->vesselNeighbor = mpGraphBloodVesselNetwork->mvNode[readECMConnections[i].StartNodeFlag];
			sStart->vesselNeighbor->highlight = 1; // ECM
			sEnd->vesselNeighbor = mpGraphBloodVesselNetwork->mvNode[readECMConnections[i].EndNodeFlag];
			sEnd->vesselNeighbor->highlight = 1; // ECM
												 // add fiber object into the arena
			mpArena->addObject(edge->GLObject());
		}
	}
	catch (...)
	{
		std::cout << "No data set found named \"" << dataSetECMSpheresString << "\".\n";
	}

	// added by Jieling
	try
	{
		// read in HSCs' data:
		std::string dataSetHSCString =
			groupString + "/HSC/HSCSpheres";

		H5::DataSet dataSetHSC =
			inputFile->openDataSet(dataSetHSCString);

		H5::CompType dataTypeHSC =
			dataSetHSC.getCompType();

		H5::CompType readType
			= ModelElementECMSphere::ParseHDF5DataFormat(dataTypeHSC,
				errors,
				warnings);

		H5::DataSpace outputSpace = dataSetHSC.getSpace();
		outputSpace.getSimpleExtentDims(&numElements, NULL);

		ModelElementECMSphere * readHSCs = (ModelElementECMSphere *)malloc(numElements * sizeof(ModelElementECMSphere));
		for (unsigned int i = 0; i < numElements; ++i)
			new (readHSCs + i) ModelElementECMSphere();

		// the actual reading of the data into the readout array
		dataSetHSC.read(readHSCs, readType);

		for (unsigned int i = 0; i < numElements; ++i)
		{
			ModelElementECMSphere *HSCi = readHSCs + i;
			HSCi->setElementType(ModelElementECMSphere::HSC);
			HSCs.push_back(HSCi);
		}

		struct HSCConnection { hsize_t indexHSC; hsize_t flagHSC; };
		// reading in HSC sphere connections
		std::string dataSetHSCConnectionString =
			groupString + "/HSC/HSCConnections";

		H5::DataSet dataSetHSCConnections =
			inputFile->openDataSet(dataSetHSCConnectionString);

		H5::CompType dataTypeHSCConnection(sizeof(HSCConnection));
		dataTypeHSCConnection.insertMember("HSCIndex", HOFFSET(HSCConnection, indexHSC), H5::PredType::STD_U64LE);
		dataTypeHSCConnection.insertMember("HSCFlag", HOFFSET(HSCConnection, flagHSC), H5::PredType::STD_U64LE);

		H5::DataSpace outputSpaceHSCConnections = dataSetHSCConnections.getSpace();
		outputSpaceHSCConnections.getSimpleExtentDims(&numElements, NULL);

		HSCConnection * readHSCConnections =
			new HSCConnection[numElements];

		dataSetHSCConnections.read(readHSCConnections, dataTypeHSCConnection);

		for (unsigned long i = 0; i < numElements; ++i)
		{
			ModelElementECMSphere * sHSC = HSCs[readHSCConnections[i].indexHSC];
			sHSC->vesselNeighbor = mpGraphBloodVesselNetwork->mvNode[readHSCConnections[i].flagHSC];
			if (sHSC->vesselNeighbor->highlight != 1)
				sHSC->vesselNeighbor->highlight = 2; // HSC
		}
	}
	catch (...)
	{
	}


    try
    {
        // read in cells' data:
        std::string dataSetCellsSphericalString =
            groupString + "/CellsSpherical";

        H5::DataSet dataSetCellsSpherical =
            inputFile->openDataSet(dataSetCellsSphericalString);

        // ToDo:  set defaults for non-essential data fields

        // build up the data type in a flexible way

        // The data layout of the saved data
        H5::CompType dataTypeCellSpherical =
            dataSetCellsSpherical.getCompType();

        // the 'read type', i.e. the data layout for our read-in buffer:
        H5::CompType readType
        = CellSpherical::ParseHDF5DataFormat( dataTypeCellSpherical,
                                              errors,
                                              warnings );

        // how many data elements of type dataTypeCellSpherical are in the dataSet
        H5::DataSpace outputSpace = dataSetCellsSpherical.getSpace();
        outputSpace.getSimpleExtentDims(&numElements, NULL);

        // allocate the read-out buffer
        //CellSpherical * readCells = new CellSpherical[numElements];
        // Hack to initialize readCells to be able to delete single cells from
        // it afterwards.  lookup 'placement new in c++' if you do not
        // understand the syntax of the 'new()' call.
        CellSpherical * readCells = (CellSpherical *) malloc (numElements*sizeof(CellSpherical));
        for (unsigned int i=0; i<numElements; ++i)
            new (readCells+i) CellSpherical();

        // the actual reading of the data into the readout array
        dataSetCellsSpherical.read( readCells, readType );

        // add the cells to this Model:
        for ( unsigned int i=0; i<numElements; ++i )
        {
            CellSpherical *cell = readCells +i;
            AddCell( cell );

            if ( (long)cell->mDivisionPartnerIndex != -1 )
            {
                cell->mpDivisionPartner = readCells + cell->mDivisionPartnerIndex;
                cell->mDaughterCell = ( cell->mDivisionPartnerIndex < i );
                cell->setState(Cell::StateDividing);
            }
            else
                cell->mpDivisionPartner = NULL;

            cell->mSPhaseStart = cell->cycleTime * 5/12;
            cell->mSPhaseEnd   = cell->cycleTime * 3/4;
        }
    }
    catch (...)
    {}

    try
    {
        // read in cells' data:
        std::string dataSetCellsSphericalString =
            groupString + "/CellsSphericalPolar";

        H5::DataSet dataSetCellsSpherical =
            inputFile->openDataSet(dataSetCellsSphericalString);

        // ToDo:  set defaults for non-essential data fields

        // build up the data type in a flexible way

        // The data layout of the saved data
        H5::CompType dataTypeCellSpherical =
            dataSetCellsSpherical.getCompType();

        // the 'read type', i.e. the data layout for our read-in buffer:
        H5::CompType readType
        = CellSphericalPolar::ParseHDF5DataFormat( dataTypeCellSpherical,
                errors,
                warnings );

        // how many data elements of type dataTypeCellSpherical are in the dataSet
        H5::DataSpace outputSpace = dataSetCellsSpherical.getSpace();
        outputSpace.getSimpleExtentDims(&numElements, NULL);

        //CellSphericalPolar * readCells = new CellSphericalPolar[numElements];
        // Hack to initialize readCells to be able to delete single cells from
        // it afterwards.  lookup 'placement new in c++' if you do not
        // understand the syntax of the 'new()' call.
        CellSphericalPolar * readCells = (CellSphericalPolar *) malloc ( numElements*sizeof(CellSphericalPolar) );
        for (unsigned int i=0; i<numElements; ++i)
            new (readCells+i) CellSphericalPolar();

        // the actual reading of the data into the readout array
        dataSetCellsSpherical.read( readCells, readType );

        // add the cells to this Model:
        for ( unsigned int i=0; i<numElements; ++i )
        {
            CellSphericalPolar *cell = readCells +i;
            AddCell( cell );

            if ( (long)cell->mDivisionPartnerIndex != -1 )
            {
                cell->mpDivisionPartner = readCells + cell->mDivisionPartnerIndex;
                cell->mDaughterCell = ( cell->mDivisionPartnerIndex < i );
                cell->setState(Cell::StateDividing);
            }
            else
                cell->mpDivisionPartner = NULL;

            cell->mSPhaseStart = cell->cycleTime * 5/12;
            cell->mSPhaseEnd   = cell->cycleTime * 3/4;
        }

    }
    catch (...)
    {}


    // Now ALWAYS determined in InitScenarioLobuleRegeneration.
    //
    // // Cells to delete.

    // // For the intoxication in the lobule model, we read in all remaining cells
    // // that were picked to be killed in the original run.  These cells will be
    // // removed in the same order if and only if the run procedes at the same
    // // point, where the data was saved and all data relevant for the order of
    // // calling the RNG is complete or initialized as well.
    // try
    // {
    //     std::string dataSetCellsToDeleteString =
    //         groupString + "/CellsToDelete";

    //     H5::DataSet dataSetCellsToDelete =
    //         inputFile->openDataSet(dataSetCellsToDeleteString);

    //     H5::PredType dataTypeCellsToDelete = H5::PredType::STD_U64LE;

    //     // how many data elements are in the dataSet
    //     H5::DataSpace outputSpace = dataSetCellsToDelete.getSpace();
    //     outputSpace.getSimpleExtentDims(&numElements, NULL);

    //     unsigned long *readIndices = (unsigned long *)malloc(numElements*sizeof(unsigned long));

    //     //CellSphericalPolar * readCells = new CellSphericalPolar[numElements];
    //     // the actual reading of the data into the readout array
    //     dataSetCellsToDelete.read( readIndices, dataTypeCellsToDelete );

    //     for ( unsigned int i=0; i<numElements; ++i )
    //     {
    //         mCellsToDelete.push_back( cells.at(readIndices[i]) );
    //     }

    //     free( readIndices );
    // }
    // catch (...)
    // {}

    UpdateParametersForAllCells();

	// added by Jieling
	std::cout << "	-> scenario: " << mScenario << std::endl;
	if (mScenario == ScenarioLobuleRegeneration)
		InitScenarioLobuleRegeneration(true);
}


void
ModelCellsSpherical::writeHDF5( H5::H5File * outputFile ) const
{
    // Todo: error:
    if ( ! outputFile )
        return;

    // dimensions - have to be handed over as pointers:
    hsize_t dims[] = { 1 };
    hsize_t cellsSphericalMaxDims[] = { cells.size() };
    // the memory space of a single data element
    H5::DataSpace memspace( 1, dims, NULL );


    std::string groupString = "/" + name;
    H5::Group modelGroup( outputFile->createGroup(groupString) );


    // preparing the Marsaglia random generator to be saved:
    H5::CompType rngDataType( sizeof(Random) );
    Random::HDF5DataFormat( rngDataType );

    H5::DataSpace rngDataSpace( 1, dims, NULL );

    H5::DataSet rngDataSet( outputFile->createDataSet( groupString + "/MarsagliaRNG",
                            rngDataType,
                            rngDataSpace ) );

    rngDataSet.write( &mRandom, rngDataType, memspace, rngDataSpace );


    // creating the compound type and defining its composition
    // (by calling CellSpherical::HDF5Dataformat())

    H5::DSetCreatPropList cparms;

    hsize_t offset = 0;
    hsize_t count  = 1;

    std::string DataSetName;
    H5::CompType cellSphericalType;
    // hack for now until std::map< ModelElement::Type,
    // std::vector<ModelElement*> > is finished and merged:
    if ( cells[0]->mType == ModelElement::TypeCellSpherical )
    {
        cellSphericalType = H5::CompType( sizeof( CellSpherical ) );
        CellSpherical::HDF5DataFormat( cellSphericalType );
        DataSetName = "/CellsSpherical";
    }
    else if ( cells[0]->mType == ModelElement::TypeCellSphericalPolar )
    {
        cellSphericalType = H5::CompType( sizeof( CellSphericalPolar ) );
        CellSphericalPolar::HDF5DataFormat( cellSphericalType );
        DataSetName = "/CellsSphericalPolar";
    }

    H5::CompType cellSphericalType_packed;
    cellSphericalType_packed.copy(cellSphericalType);
    cellSphericalType_packed.pack();

    H5::DataSpace cellSphericalDataSpace( 1, cellsSphericalMaxDims, NULL );

    H5::DataSet cellSphericalDataSet( outputFile->createDataSet( groupString + DataSetName,
                                                                 cellSphericalType_packed,
                                                                 cellSphericalDataSpace ) );

    std::vector<CellSpherical *>::const_iterator cIt;

    for ( cIt = cells.begin(); cIt != cells.end(); ++cIt )
    {
        // preliminary:  determine dumbbell divisionPartnerIndex
        if ( (*cIt)->mpDivisionPartner )
        {
            std::vector<CellSpherical *>::const_iterator found = find( cells.begin(), cells.end(), (*cIt)->mpDivisionPartner );
            if ( found != cells.end() )
                (*cIt)->mDivisionPartnerIndex = found - cells.begin();
            else
                (*cIt)->mDivisionPartnerIndex = -1;
        }
        else
            (*cIt)->mDivisionPartnerIndex = -1;

        cellSphericalDataSpace.selectHyperslab( H5S_SELECT_SET, &count, &offset );
        cellSphericalDataSet.write( *cIt, cellSphericalType, memspace, cellSphericalDataSpace );
        outputFile->flush(H5F_SCOPE_LOCAL);
        offset++;
    }


    // write data set of cells to delete: as indices within
    // CellsSpherical/CellsSphericalPolar data set, i.e. their index in member
    // variable
    //   std::vector<CellSpherical *> cells.
    H5::PredType doomedCellsType = H5::PredType::STD_U64LE;
    hsize_t doomedCellsDims[] = { mCellsToDelete.size() };
    H5::DataSpace doomedCellsDataSpace(1, doomedCellsDims, NULL);
    H5::DataSet   doomedCellsDataSet( outputFile->createDataSet( groupString + "/CellsToDelete",
                                                                 doomedCellsType,
                                                                 doomedCellsDataSpace,
                                                                 cparms ) );
    if ( mCellsToDelete.size() != 0 )
    {
        offset = 0;
        std::vector<CellSpherical *>::const_iterator doomedCellsIt;
        for( doomedCellsIt = cells.begin();
             doomedCellsIt != cells.end();
             ++doomedCellsIt )
        {
            if ( find( mCellsToDelete.begin(), mCellsToDelete.end(), *doomedCellsIt )
                 == mCellsToDelete.end() )
                continue;
            hsize_t index = doomedCellsIt - cells.begin();
            doomedCellsDataSpace.selectHyperslab( H5S_SELECT_SET, &count, &offset);
            doomedCellsDataSet.write( &index,
                                      doomedCellsType,
                                      memspace,
                                      doomedCellsDataSpace );
            ++offset;
        }
        outputFile->flush(H5F_SCOPE_LOCAL);
    }


    std::string vesselGraphGroupString = groupString + "/VesselGraph";
    H5::Group vesselGraphGroup( outputFile->createGroup(vesselGraphGroupString) );

    H5::CompType vesselSphereType( sizeof(ModelElementVesselSphere) );
    ModelElementVesselSphere::HDF5DataFormat( vesselSphereType );

    hsize_t vesselSpheresMaxDims[] = { (mpGraphBloodVesselNetwork)
                                       ? mpGraphBloodVesselNetwork->mvNode.size()
                                       : 0
                                     };
    H5::DataSpace vesselSpheresDataSpace(1, vesselSpheresMaxDims, NULL );

    H5::CompType vesselSphereType_packed;
    vesselSphereType_packed.copy( vesselSphereType );
    vesselSphereType_packed.pack();

    H5::DataSet vesselSpheresDataSet( outputFile->createDataSet( vesselGraphGroupString + "/ModelElementVesselSpheres",
                                      vesselSphereType_packed,
                                      vesselSpheresDataSpace,
                                      cparms ) );

    struct vesselConnection { hsize_t indexStartNode; hsize_t indexEndNode; } connector;
    H5::CompType vesselConnectionType( sizeof(vesselConnection) );
    vesselConnectionType.insertMember( "startIndex", HOFFSET(vesselConnection, indexStartNode), H5::PredType::STD_U64LE);
    vesselConnectionType.insertMember( "endIndex", HOFFSET(vesselConnection, indexEndNode), H5::PredType::STD_U64LE);

    hsize_t vesselConnectionsMaxDims[] = { (mpGraphBloodVesselNetwork)
                                           ? mpGraphBloodVesselNetwork->mvEdge.size()
                                           : 0 };
    H5::DataSpace vesselConnectionsDataSpace(1, vesselConnectionsMaxDims, NULL );

    H5::DataSet vesselConnectionsDataSet( outputFile->createDataSet( vesselGraphGroupString + "/VesselSphereConnections",
                                      vesselConnectionType,
                                      vesselConnectionsDataSpace,
                                      cparms ) );

    if (mpGraphBloodVesselNetwork)
    {
        offset = 0;
        for ( auto vesselSphere: mpGraphBloodVesselNetwork->mvNode )
        {
            vesselSpheresDataSpace.selectHyperslab( H5S_SELECT_SET, &count, &offset );
            vesselSpheresDataSet.write( vesselSphere,
                                        vesselSphereType,
                                        memspace,
                                        vesselSpheresDataSpace );
            ++offset;
        }
        outputFile->flush(H5F_SCOPE_LOCAL);

        offset = 0;
        for ( auto connection:  mpGraphBloodVesselNetwork->mvEdge )
        {
            // find mvEdge[i].mpStart/mpEnd in mvNode and get their offset
            connector.indexStartNode = find( mpGraphBloodVesselNetwork->mvNode.begin(),
                                             mpGraphBloodVesselNetwork->mvNode.end(),
                                             connection->mpStart )
                - mpGraphBloodVesselNetwork->mvNode.begin();
            if ( connector.indexStartNode >= mpGraphBloodVesselNetwork->mvNode.size() )
            {
                std::cerr << "Error writing connection information for VesselSpheres"
                          << std::endl
                          << "Start node of GraphEdge not found in node list: " << connection->mpStart
                          << std::endl << "type: " << connection->mpStart->mType
                          << "  Global index: " << connection->mpStart->mGlobalIndex
                          << std::endl << "connected to " << connection->mpEnd
                          << std::endl << std::endl;
                continue;
            }
            connector.indexEndNode = find( mpGraphBloodVesselNetwork->mvNode.begin(),
                                           mpGraphBloodVesselNetwork->mvNode.end(),
                                           connection->mpEnd )
                - mpGraphBloodVesselNetwork->mvNode.begin();
            if ( connector.indexEndNode >= mpGraphBloodVesselNetwork->mvNode.size() )
            {
                std::cerr << "Error writing connection information for VesselSpheres"
                          << std::endl
                          << "Start node of GraphEdge not found in node list: " << connection->mpEnd
                          << std::endl << "type: " << connection->mpEnd->mType
                          << "  Global index: " << connection->mpEnd->mGlobalIndex
                          << std::endl << "connected to " << connection->mpStart
                          << std::endl << std::endl;
                continue;
            }

            vesselConnectionsDataSpace.selectHyperslab( H5S_SELECT_SET, &count, &offset );
            vesselConnectionsDataSet.write( &connector,
                                        vesselConnectionType,
                                        memspace,
                                        vesselConnectionsDataSpace );
            ++offset;
        }
        outputFile->flush(H5F_SCOPE_LOCAL);
    }

	// added by Jieling
	if (HSCs.size() > 0)
	{
		offset = 0;

		std::string HSCGroupString = groupString + "/HSC";
		H5::Group HSCGroup(outputFile->createGroup(HSCGroupString));

		H5::CompType HSCType(sizeof(ModelElementECMSphere));
		ModelElementECMSphere::HDF5DataFormat(HSCType);
		hsize_t HSCmaxDims[] = { HSCs.size() };
		H5::CompType HSCType_packed;
		HSCType_packed.copy(HSCType);
		HSCType_packed.pack();
		H5::DataSpace HSCDataSpace(1, HSCmaxDims, NULL);
		H5::DataSet HSCDataSet(outputFile->createDataSet(HSCGroupString + "/HSCSpheres",
			HSCType_packed,
			HSCDataSpace));

		std::vector<ModelElementECMSphere *>::const_iterator hIt;

		for (hIt = HSCs.begin(); hIt != HSCs.end(); ++hIt)
		{
			HSCDataSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
			HSCDataSet.write(*hIt, HSCType, memspace, HSCDataSpace);
			outputFile->flush(H5F_SCOPE_LOCAL);
			offset++;
		}

		struct HSCConnection { hsize_t indexHSC; hsize_t flagHSC; } HSCConnector;
		H5::CompType HSCConnectionType(sizeof(HSCConnection));
		HSCConnectionType.insertMember("HSCIndex", HOFFSET(HSCConnection, indexHSC), H5::PredType::STD_U64LE);
		HSCConnectionType.insertMember("HSCFlag", HOFFSET(HSCConnection, flagHSC), H5::PredType::STD_U64LE);

		hsize_t HSCConnectionsMaxDims[] = { HSCs.size() };
		H5::DataSpace HSCConnectionsDataSpace(1, HSCConnectionsMaxDims, NULL);
		H5::DataSet HSCConnectionsDataSet(outputFile->createDataSet(HSCGroupString + "/HSCConnections",
			HSCConnectionType,
			HSCConnectionsDataSpace,
			cparms));

		offset = 0;
		for (int i = 0; i < HSCs.size(); i++)
		{
			HSCConnector.indexHSC = find(HSCs.begin(),
				HSCs.end(),
				HSCs[i]) - HSCs.begin();
			HSCConnector.flagHSC = find(mpGraphBloodVesselNetwork->mvNode.begin(),
				mpGraphBloodVesselNetwork->mvNode.end(),
				HSCs[i]->vesselNeighbor) - mpGraphBloodVesselNetwork->mvNode.begin();

			HSCConnectionsDataSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
			HSCConnectionsDataSet.write(&HSCConnector,
				HSCConnectionType,
				memspace,
				HSCConnectionsDataSpace);
			++offset;
		}
		outputFile->flush(H5F_SCOPE_LOCAL);
	}

    H5::CompType barrierType( sizeof(ModelElementBarrierTriangle) );
    ModelElementBarrierTriangle::HDF5DataFormat( barrierType );

    hsize_t barrierMaxDims[] = { mBarriers.size() };
    H5::DataSpace barrierDataSpace(1, barrierMaxDims, NULL );

    H5::DataSet barrierDataSet( outputFile->createDataSet( groupString + "/ModelElementBarrierTriangles",
                                barrierType,
                                barrierDataSpace,
                                cparms ) );

    offset = 0;
    ModelElementBarrierTriangle ** barrIt;
    for ( barrIt = mBarriers.begin(); barrIt != mBarriers.end(); ++barrIt )
    {
        barrierDataSpace.selectHyperslab( H5S_SELECT_SET, &count, &offset );
        barrierDataSet.write( *barrIt, barrierType, memspace, barrierDataSpace );
        outputFile->flush(H5F_SCOPE_LOCAL);
        ++offset;
    }

	// added by Jieling
	std::string ECMGraphGroupString = groupString + "/ECMGraph";
	H5::Group ECMGraphGroup(outputFile->createGroup(ECMGraphGroupString));

	H5::CompType ECMSphereType(sizeof(ModelElementLatticeNode));
	ModelElementLatticeNode::HDF5DataFormat(ECMSphereType);

	hsize_t ECMSpheresMaxDims[] = { (mpCollagenNetwork)
		? mpCollagenNetwork->getNodes().size()
		: 0 };
	H5::DataSpace ECMSpheresDataSpace(1, ECMSpheresMaxDims, NULL);

	H5::CompType ECMSphereType_packed;
	ECMSphereType_packed.copy(ECMSphereType);
	ECMSphereType_packed.pack();

	H5::DataSet ECMSpheresDataSet(outputFile->createDataSet(ECMGraphGroupString + "/ModelElementLatticeNodes",
		ECMSphereType_packed,
		ECMSpheresDataSpace,
		cparms));

	struct ECMConnection {
		hsize_t indexStartNode; hsize_t indexEndNode;
		hsize_t StartNodeFlag; hsize_t EndNodeFlag;
		hsize_t StartNodeAlongFlag; hsize_t EndNodeAlongFlag;
	} ECMConnector;
	H5::CompType ECMConnectionType(sizeof(ECMConnection));
	ECMConnectionType.insertMember("startIndex", HOFFSET(ECMConnection, indexStartNode), H5::PredType::STD_U64LE);
	ECMConnectionType.insertMember("endIndex", HOFFSET(ECMConnection, indexEndNode), H5::PredType::STD_U64LE);
	ECMConnectionType.insertMember("startFlag", HOFFSET(ECMConnection, StartNodeFlag), H5::PredType::STD_U64LE);
	ECMConnectionType.insertMember("endFlag", HOFFSET(ECMConnection, EndNodeFlag), H5::PredType::STD_U64LE);
	ECMConnectionType.insertMember("startAlongFlag", HOFFSET(ECMConnection, StartNodeAlongFlag), H5::PredType::STD_U64LE);
	ECMConnectionType.insertMember("endAlongFlag", HOFFSET(ECMConnection, EndNodeAlongFlag), H5::PredType::STD_U64LE);

	hsize_t ECMConnectionsMaxDims[] = { (mpCollagenNetwork)
		? mpCollagenNetwork->getSprings().size()
		: 0 };
	H5::DataSpace ECMConnectionsDataSpace(1, ECMConnectionsMaxDims, NULL);

	H5::DataSet ECMConnectionsDataSet(outputFile->createDataSet(ECMGraphGroupString + "/ECMSphereConnections",
		ECMConnectionType,
		ECMConnectionsDataSpace,
		cparms));
	if (mpCollagenNetwork)
	{
		offset = 0;
		for (auto ECMSphere : mpCollagenNetwork->getNodes())
		{
			ECMSpheresDataSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
			ECMSpheresDataSet.write(ECMSphere,
				ECMSphereType,
				memspace,
				ECMSpheresDataSpace);
			++offset;
		}
		outputFile->flush(H5F_SCOPE_LOCAL);

		offset = 0;
		for (auto connection : mpCollagenNetwork->getSprings())
		{
			ECMConnector.indexStartNode = find(mpCollagenNetwork->getNodes().begin(),
				mpCollagenNetwork->getNodes().end(),
				connection->nodes().at(0)) - mpCollagenNetwork->getNodes().begin();
			if (ECMConnector.indexStartNode >= mpCollagenNetwork->getNodes().size())
			{
				std::cerr << "Error writing connection information for ECMSpheres" << std::endl;
				continue;
			}
			ECMConnector.StartNodeFlag = find(mpGraphBloodVesselNetwork->mvNode.begin(),
				mpGraphBloodVesselNetwork->mvNode.end(),
				connection->nodes().at(0)->vesselNeighbor) - mpGraphBloodVesselNetwork->mvNode.begin();

			if (connection->nodes().at(0)->alongVesselNeighbor == NULL)
				ECMConnector.StartNodeAlongFlag = -1;
			else
			{
				ECMConnector.StartNodeAlongFlag = find(mpGraphBloodVesselNetwork->mvNode.begin(),
					mpGraphBloodVesselNetwork->mvNode.end(),
					connection->nodes().at(0)->alongVesselNeighbor) - mpGraphBloodVesselNetwork->mvNode.begin();
				if (ECMConnector.StartNodeAlongFlag > mpGraphBloodVesselNetwork->mvNode.size())
				{
					std::cerr << "	-> Error: startNodeAlongFlag is too large!" << std::endl;
				}
			}
			//
			ECMConnector.indexEndNode = find(mpCollagenNetwork->getNodes().begin(),
				mpCollagenNetwork->getNodes().end(),
				connection->nodes().at(1)) - mpCollagenNetwork->getNodes().begin();
			if (ECMConnector.indexStartNode >= mpCollagenNetwork->getNodes().size())
			{
				std::cerr << "Error writing connection information for ECMSpheres" << std::endl;
				continue;
			}
			ECMConnector.EndNodeFlag = find(mpGraphBloodVesselNetwork->mvNode.begin(),
				mpGraphBloodVesselNetwork->mvNode.end(),
				connection->nodes().at(1)->vesselNeighbor) - mpGraphBloodVesselNetwork->mvNode.begin();
			if (connection->nodes().at(1)->alongVesselNeighbor == NULL)
				ECMConnector.EndNodeAlongFlag = -1;
			else
			{
				ECMConnector.EndNodeAlongFlag = find(mpGraphBloodVesselNetwork->mvNode.begin(),
					mpGraphBloodVesselNetwork->mvNode.end(),
					connection->nodes().at(1)->alongVesselNeighbor) - mpGraphBloodVesselNetwork->mvNode.begin();
				if (ECMConnector.EndNodeAlongFlag > mpGraphBloodVesselNetwork->mvNode.size())
				{
					std::cerr << "	-> Error: endNodeAlongFlag is too large!" << std::endl;
				}
			}
			ECMConnectionsDataSpace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
			ECMConnectionsDataSet.write(&ECMConnector,
				ECMConnectionType,
				memspace,
				ECMConnectionsDataSpace);
			++offset;
		}
		outputFile->flush(H5F_SCOPE_LOCAL);
	}
}

#  pragma endregion

#pragma endregion


#pragma region Visualization

void ModelCellsSpherical::UpdateCellsStaining(int mode)
{
    #pragma region 0: All white
    if (mode == 0)
    {
        for (unsigned int i=0; i<cells.size(); i++)
        {
			if (cells.at(i)->mSubType==1) cells.at(i)->SetColor(227./255., 179./255., 94./255.);
            else 
			{
            cells.at(i)->SetColor(1, 1, 1);
			// added by Jieling, for testing
			cells.at(i)->SetAlpha(0.);
			}
        }
    }
    #pragma endregion

    #pragma region 1: Cell cycle state

    else if (mode == 1)
    {
        for (unsigned int i=0; i<cells.size(); i++)
        {
			if (cells.at(i)->mSubType==1) cells.at(i)->SetColor(227./255., 179./255., 94./255.);
            else 
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
    }
    #pragma endregion

    #pragma region 2: Last absolute force

    else if (mode == 2)
    {
        for (unsigned int i=0; i<cells.size(); i++)
        {
			if (cells.at(i)->mSubType==1) cells.at(i)->SetColor(227./255., 179./255., 94./255.);
            else 
			{
            (*core->tools).color->CreateSuperTrafficHue(cells.at(i)->lastPressure, 0., 20000.);
            cells.at(i)->SetColor( (*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b );
			}
        }
    }

    #pragma endregion

    #pragma region 3: Cell volume

    else if (mode == 3)
    {
        for (unsigned int i=0; i<cells.size(); i++)
        {
			if (cells.at(i)->mSubType==1) cells.at(i)->SetColor(227./255., 179./255., 94./255.);
            else 
			{
            // Prelim: Later: Actual volume calculation
            (*core->tools).color->CreateTrafficHue(cells.at(i)->mRadius, defaultInitialCellRadius, defaultDivisionCellRadius);
            cells.at(i)->SetColor( (*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b );
			}
        }
    }


    #pragma endregion

    else if (mode == 4) // relative position within a lobule
    {
        for (unsigned int i=0; i<cells.size(); ++i)
        {
            (*core->tools).color->CreateTrafficHue(cells[i]->mLayerIndicator, 0, 13);
            cells[i]->SetColor( (*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b );
        }
    }
    else if (mode == 5)
    {
        for (unsigned int i=0; i<cells.size(); ++i)
        {
            if (cells[i]->mLesionEdge)
                cells[i]->SetColor( 1, 0, 0 );
            else
                cells[i]->SetColor( 1, 1, 1 );
        }
    }
	// Remember colormode
	lastColormodeCells= mode;
	// added by Jieling
	//mpCollagenNetwork->getNodes().at(3)->SimulationObject::SetColor(1., 1., 1., 1.);
	for (int i = 0; i < mpCollagenNetwork->getSprings().size(); i++)
	{
		if (mpCollagenNetwork->getSprings().at(i)->mSpringType == LatticeSpring::TypeLinearSpring)
		{
			LinearSpring *lsI = (LinearSpring*)(mpCollagenNetwork->getSprings().at(i));
			double strain = lsI->curStrainStaining;
			double Bth = 0.1;
			double BScale = 1. / Bth;
			if (strain < 0) // to blue
			{
				double col = std::sqrt(-BScale * strain);
				if (col > 1.) col = 1.;
				lsI->SetColor(0., 1. - col, col, 1.);
			}
			else // to red
			{
				double col = std::sqrt(BScale * strain);
				if (col > 1.) col = 1.;
				lsI->SetColor(col, 1. - col, 0., 1.);
			}
		}
	}
}


void ModelCellsSpherical::UpdateCellsStaining()
{
    UpdateCellsStaining(lastColormodeCells);
}


void ModelCellsSpherical::UpdateCellsStaining(CellSpherical *cell)
{
    #pragma region 0: All white

  switch (lastColormodeCells)
  {
  case 0:
    if (cell->mSubType==1) cell->SetColor(227./255., 179./255., 94./255.);
    else
      cell->SetColor(1, 1, 1);
    break;
#pragma endregion

#pragma region 1: Cell cycle state

  case 1:
    if (cell->mSubType==1) cell->SetColor(227./255., 179./255., 94./255.);
    else
	{
	  if (cell->getState(Cell::StateQuiescent))
      {
        cell->SetColor(0.5, 0.5, 0.5);
      }
      else
      {
        cell->SetColor(1, 1, 1);
      }
	}
    break;

#pragma endregion

#pragma region 2: Last absolute force

  case 2:
    if (cell->mSubType==1) cell->SetColor(227./255., 179./255., 94./255.);
    else
    {
      (*core->tools).color->CreateSuperTrafficHue(cell->lastPressure, 0., 20000.);
      cell->SetColor( (*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b );
    }
    break;

#pragma endregion

#pragma region 3: Cell volume

  case 3:
    if (cell->mSubType==1) cell->SetColor(227./255., 179./255., 94./255.);
    else
    {
      // Prelim: Later: Actual volume calculation
      (*core->tools).color->CreateTrafficHue(cell->mRadius, defaultInitialCellRadius, defaultDivisionCellRadius);
      cell->SetColor( (*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b );
    }
    break;

  case 4:
    // relative position within a Lobule
    (*core->tools).color->CreateTrafficHue(cell->mLayerIndicator, 0, 13);
    cell->SetColor( (*core->tools).color->r, (*core->tools).color->g, (*core->tools).color->b );
    break;

  case 5:
    if (cell->mLesionEdge)
      cell->SetColor( 1, 0, 0 );
    else
      cell->SetColor( 1, 1, 1 );
    break;
  }
    #pragma endregion
}


void ModelCellsSpherical::setVisible( bool visible )
{

    for( unsigned int i = 0 ; i < this->cells.size() ; i++ )
    {
        if( this->cells[i]->mVisible == true)
        {
            mpArena->removeObject(cells[i]->GLObject() );
        }
        else
        {
            mpArena->addObject(this->cells[i]->GLObject());
        }

        this->cells[i]->mVisible = visible;
    }

    mpArena->draw();

}

void ModelCellsSpherical::writeXML2()
{

    std::ofstream f;
    f.open("../../output/default_changed.xml");

    f << "<?xml version=\"1.0\"?>" << std::endl;
    f << "<Model>" << std::endl;
    f << "<VesselGraph dimensions=\"3\" x=\""
      << ceil( (this->mpGraphBloodVesselNetwork->mXmax-this->mpGraphBloodVesselNetwork->mXmin)*23.3)
      << "\" y=\""
      << ceil( (this->mpGraphBloodVesselNetwork->mYmax-this->mpGraphBloodVesselNetwork->mYmin)*23.3)
      << "\" z=\""
      << ceil( (this->mpGraphBloodVesselNetwork->mZmax-this->mpGraphBloodVesselNetwork->mZmin)*23.3)
      << "\" latticeConstant=\""
      << 5
      << "\">"
      << std::endl;

    for( unsigned int i = 0 ; i < this->mpGraphBloodVesselNetwork->mvNode.size() ; i++)
    {
        f << "<Node id=\""
          << i
          << "\" x=\""
          << (this->mpGraphBloodVesselNetwork->mvNode[i]->position.x - this->mpGraphBloodVesselNetwork->mXmin)*23.3
          << "\" y=\""
          << (this->mpGraphBloodVesselNetwork->mvNode[i]->position.y - this->mpGraphBloodVesselNetwork->mYmin)*23.3
          << "\" z=\""
          << (this->mpGraphBloodVesselNetwork->mvNode[i]->position.z - this->mpGraphBloodVesselNetwork->mZmin)*23.3
          << "\" type=\"";


        switch(this->mpGraphBloodVesselNetwork->mvNode[i]->mType)
        {
        case 1:
            f << 1;
            f << "\" pressure=\"0\"/>";
            break;
        case 2:
            f << 2;
            f << "\" pressure=\"\"/>";
            break;
        case 3:
            f << 1;
            f << "\" pressure=\"10\"/>";
            break;
        default:
            fprintf(stderr,"error:%i" , this->mpGraphBloodVesselNetwork->mvNode[i]->mType);
            break;
        }
        f << std::endl;

    }

    for( unsigned int i = 0 ; i < this->mpGraphBloodVesselNetwork->mvEdge.size() ; i++)
    {
        f << "<Segment id=\""
          << i
          << "\" node1=\""
          << this->mpGraphBloodVesselNetwork->mvEdge[i]->mpStart->mIndex
          << "\" node2=\""
          << this->mpGraphBloodVesselNetwork->mvEdge[i]->mpEnd->mIndex
          << "\" radiusStatic=\""
          << 2.5
          << "\"/>"
          << std::endl;
    }
    f << "</VesselGraph>" << std::endl;
    f << "</Model>" << std::endl;


    f.close();



}


void ModelCellsSpherical::spherePacking(double radiusSphere){

  int counter=0;

  double radius = 0.5;

  int max = (radiusSphere + 3) * 2;
  double max_h = radiusSphere;

  int I = max;
  int J = max;
  int K = max;

  Vector3f pos;

  for( int i = 0 ; i < I ; i++ )
    for( int j = 0 ; j < J ; j++ )
      for( int k = 0 ; k < K ; k++ ){

         pos.x = (2.*i+((j+k)%2))         * radius - max_h;
         pos.y = sqrt(3.)*(j+1./3.*(k%2)) * radius - max_h;
         pos.z = 2.*sqrt(6.)/3.*k         * radius - max_h;

        if( pos.Norm() < radiusSphere ){
          AddCell(pos.x,pos.y,pos.z);
          std::cerr << "add cell : " << this->cells.size() << std::endl;
          counter++;
        }
        else
          std::cerr << "does not add cell (cell outside)" << std::endl;
      }


}


void ModelCellsSpherical::CalcRadiusGyration(){

  Vector3f mean;

  for( unsigned int c = 0 ; c < this->cells.size() ; c++ )
    mean += this->cells[c]->position;

  mean /= this->cells.size();

  this->mRadiusGyration = 0;

  for( unsigned int c = 0 ; c < this->cells.size() ; c++ ){
    double n = (this->cells[c]->position - mean ).Norm() * this->biolink->length_scale;
    this->mRadiusGyration += n*n;
  }

  this->mRadiusGyration /= this->cells.size();

  this->mRadiusGyration = sqrt(this->mRadiusGyration);

}


void
ModelCellsSpherical::voxelize()
{
    if ( !mpLesionVoxelization )
        mpLesionVoxelization = new LobuleLesionVoxelization(this, .5);

    mpLesionVoxelization->exec();
}

#pragma endregion

// preliminary vtp writing mechanism for ModelElementVesselSpheres
void ModelCellsSpherical::writeVesselSphereVTP()
{
    if ( !mpGraphBloodVesselNetwork )
        return;

    // construct output filename
    std::stringstream outfilenamestream;
    outfilenamestream << mOutputPrefix << "sinusoids_"
                      << setfill('0') << setw(7) << mOutputCounter << ".vtp";

    std::string outputFilePath = mOutputPath + outfilenamestream.str();

    size_t numElements = mpGraphBloodVesselNetwork->mvNode.size();

    ofstream vtp;
    vtp.open(outputFilePath, std::ios::out|std::ios::trunc );

    vtp << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
    vtp << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"BigEndian\">" << std::endl;
    vtp << "<PolyData GhostLevel=\"0\">" << std::endl;

    vtp << "<Piece NumberOfPoints=\""
        << numElements
        << "\" Simtime=\""
        << time
        << "\" TimeStep=\""
        << biolink->getTimeInSeconds(timeStep)
        << "\" NumberOfVerts=\""
        << 0
        << "\" NumberOfLines=\""
        << 0
        << "\" NumberOfStrips=\"0\" NumberOfPolys=\""
        << 0
        << "\">" << std::endl;

    vtp << "<Points>" << std::endl;
    vtp << "<DataArray Name=\"xyz\" NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">" << std::endl;

    for( auto vsphere: mpGraphBloodVesselNetwork->mvNode )
    {
        vtp << vsphere->position.x << " "
            << vsphere->position.y << " "
            << vsphere->position.z << std::endl;
    }

    vtp << "</DataArray>" << std::endl
        << "</Points>"    << std::endl;


    vtp << "<PointData>" << std::endl
        << "<DataArray Name=\"Radius\" NumberOfComponents=\"1\" type=\"Float64\" format=\"ascii\">" << std::endl;

    for ( auto vsphere: mpGraphBloodVesselNetwork->mvNode )
    {
        vtp << vsphere->mRadius << std::endl;
    }

    vtp << "</DataArray>" << std::endl;

    vtp << "</PointData>" << std::endl;

    vtp << "</Piece>"     << std::endl;
    vtp << "</PolyData>"  << std::endl;
    vtp << "</VTKFile>"   << std::endl;

    vtp.close();
}
