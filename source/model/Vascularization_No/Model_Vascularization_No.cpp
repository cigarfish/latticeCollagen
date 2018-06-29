#include "Model_Vascularization_No.h"

#include "../Model/ModelVascularization/CSVesselGraph.h"

#include "../gui/QCSSimulationThread.h"
#include "../gui/CSGLArena.h"

#include "../Model/ModelVascularization/Perfusion.cpp"
#include "../gui/GLTools/CSGLVascularization.h"

Model_Vascularization_No::Model_Vascularization_No(): CSModel (3)
{
  this->labelCells = 0;
  this->vornoi = NULL;
  this->vg = NULL;
  this->concentration = 0;

}

void Model_Vascularization_No::Reset(){
	
}

void Model_Vascularization_No::SimulateInThread(){

	//  float *du=0, *u=0;
	  CONCENTRATION_T *du = 0, *u = 0;
	  const float dt = 0.001;

	  float borderCondition = 10.;

	  for( float time = 0; time <= 10; time+=dt){

		  fprintf( stderr, "\r%1.6f\b", time);
if( time > 0.1)
	borderCondition = 0.;
	   // float borderCondition = 10.;//(1+sin(time))*3;
	 //   vascularization->vg->updateMarkerVesselsExplicit( dt, borderCondition, du, u, 0);
	   this->vg->updateMarkerVesselsAndInterstitialExplicit(
	        dt ,
	        borderCondition ,
	        du ,
	        u ,
	        0 ,
	        this->concentration ,//not used
	        this->labelCells ,
	        this->vornoi);
	  }
}

void Model_Vascularization_No::SetupSimulation(){

  Tumor *tumor = NULL;
  this->vg =new CSVesselGraph("../../input/sin_graph0.xml",tumor );

  //this->vg =new CSVesselGraph("../../input/default_changed.xml",tumor);
  this->vg->simplifyNetwork(3.);

  this->vg->updatePressureNEW2();
  this->vg->updateFlow();

  for( int i = 0 ; i < this->vg->countVesselSegments ; i++ )
    this->vg->vesselSegments[i]->permeability = 0.1;

  this->mpArena->addObject( this->vg->GLObject() );

  this->vornoi = createVoronoiTesselation( this->vg, 5. );

  CSGLVascularization* tmp = (CSGLVascularization*) this->vg->GLObject() ;
  tmp->setVoronoi(this->vornoi);
  this->vornoi->setVolume();

  this->vg->diffusionCoefficient = 1000.;

}

