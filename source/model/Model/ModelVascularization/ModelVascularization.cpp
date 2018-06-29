#include "ModelVascularization.h"

#include "CSVesselGraph.h"

#include "../gui/QCSSimulationThread.h"
#include "../gui/CSGLArena.h"

#include "Perfusion.cpp"
#include "../gui/GLTools/CSGLVascularization.h"

Model_Vascularization::Model_Vascularization(): CSModel (3)
{
  this->labelCells = 0;
  this->vornoi = NULL;
  this->vg = NULL;
  this->concentration = 0;

}

void Model_Vascularization::Reset(){
	
}

void Model_Vascularization::SimulateInThread(){

	//  float *du=0, *u=0;
	  CONCENTRATION_T *du = 0, *u = 0;
	  const float dt = 0.0001;




		  float borderCondition = 10.01;

	  for( float time = 0; time <= 100000; time+=dt){

		  fprintf( stderr, "\r%1.6f\b", time);
if( time > 20.0)
	borderCondition = 0.;
	   // float borderCondition = 10.;//(1+sin(time))*3;
	 //   vascularization->vg->updateMarkerVesselsExplicit( dt, borderCondition, du, u, 0);
	   this->vg->updateMarkerVesselsAndInterstitialExplicit(
	        dt ,
	        borderCondition ,
	        du ,
	        u ,
	        2 ,
	        this->concentration ,//not used
	        this->labelCells ,
	        this->vornoi);
	  }
	  this->enableSimulation = false;
}

void Model_Vascularization::SetupSimulation(){

  Tumor *tumor = NULL;
//  this->vg =new CSVesselGraph("../../input/sin_graph0.xml",tumor );

 // this->vg =new CSVesselGraph("../../input/default_changed.xml",tumor);

  this->vg =new CSVesselGraph("../../input/sin_graph_minimal4.xml",tumor);


  //this->vg->simplifyNetwork(3.);

  this->vg->updatePressureNEW2();
  this->vg->updateFlow();

  for( int i = 0 ; i < this->vg->countVesselSegments ; i++ )
    this->vg->vesselSegments[i]->permeability = 1.1;

  this->mpArena->addObject( this->vg->GLObject() );

  this->vornoi = createVoronoiTesselationHexagonal( this->vg, this->vg->LATTICE_CONSTANT );
//  this->vornoi = createVoronoiTesselation( this->vg, this->vg->LATTICE_CONSTANT );
  char outputFileName[20] = "../../input/vtk.vtk";
  char outputFileName2[23] = "../../input/vtk_vg.vtk";

  //this->vornoi->printVoronoiDiagramToVTK( outputFileName, NULL);
  //this->vg->printVesselGraphToVTK(outputFileName2);

  CSGLVascularization* tmp = (CSGLVascularization*) this->vg->GLObject() ;
  tmp->setVoronoi(this->vornoi);
  this->vornoi->setVolume();

  this->vg->diffusionCoefficient = 1000.;

}

