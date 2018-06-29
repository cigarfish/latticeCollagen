#ifndef VESSELGRAPH_HPP_
#define VESSELGRAPH_HPP_

#include <float.h>
#include "../../../tools/new_triangulation/Geometry.hpp"
#include "../../../tools/new_triangulation/LinearAlgebra.hpp"

#include "Tumor.hpp"
//#include "VoronoiDiagram.hpp"

#define RAND01 ((double)rand()/(double)(RAND_MAX+1.))
//#define MAX(a,b) (a>b?a:b)
//#define MIN(a,b) (a<b?a:b)

//#define UNIDIRECTIONAL_DIFFUSION
//#define UNIDIRECTIONAL_TRANSPORT

#define DIMENSIONS_		3

/*#if DIMENSIONS_ == 1
	// 1D SETTINGS
	#define DOMAIN_SIZE_X 60//100//601//103//301 //103
	#define DOMAIN_SIZE_Y 1//100//601//103//301 //103
	#define DOMAIN_SIZE_Z 1//50//601//103//301 //103
	#define ROOT_DISTANCE	(BRANCH_LENGTH*3)//50)

#elif DIMENSIONS_ == 2
	// 2D SETTINGS
	#define DOMAIN_SIZE_X 100//601//103//301 //103
	#define DOMAIN_SIZE_Y 100//601//103//301 //103
	#define DOMAIN_SIZE_Z 1//50//601//103//301 //103
	#define ROOT_DISTANCE	(BRANCH_LENGTH*50)
#elif DIMENSIONS_ == 3
	// 3D SETTINGS
	#define DOMAIN_SIZE_X 100//601//103//301 //103
	#define DOMAIN_SIZE_Y 100//601//103//301 //103
	#define DOMAIN_SIZE_Z 50//601//103//301 //103
	#define ROOT_DISTANCE	(BRANCH_LENGTH*50)
#endif*/

#define ROOT_DISTANCE	(BRANCH_LENGTH*50)

//#define LATTICE_CONSTANT	7.//60. //mym
#define BRANCH_LENGTH		1  //lattice sites

#define MIN_PRESSURE	0.
#define MAX_PRESSURE	10.//12.//12.

#define RADIUS_EXPONENT_ARTERIES 	3
#define RADIUS_EXPONENT_VENES 		2.7//.4

#define MIN_VESSEL_RADIUS 4.//2.5//4.
//#define MAX_VESSEL_RADIUS 200.
#define MAX_VESSEL_RADIUS (LATTICE_CONSTANT*0.4)

#define ARTERIAL_MARKER_CONC	1.
#define EXCHANGE_RATE_VESSELS_INTERSTITIAL	200.//10000.
#define MARKER_DIFFUSION_COEFFICIENT	1000.//0
//#define PERMEABILITY
 
// VesselNode Types
#define ROOT	1
#define VESSEL	2
#define SPROUT	3
#define TIP		4

// VesselSegment Types
#define UNKNOWN	0
#define ARTERIE	1
#define VENE	2
#define CAPILLARY	3

#define IMPLICIT		0
#define UPWIND			1
#define CENTERED_EULER	2
#define LAX_WENDROFF	3
#define MACCORMACK		4
//#define BEAM_WARMING	4
//#define FROMM 		5
#define MINMOD		6
#define SUPERBEE		7


/*#if DIMENSIONS_==3
	#define DISTANCE( a, b) sqrt( pow(a->position[0]-b->position[0],2) + pow(a->position[1]-b->position[1],2) + pow(a->position[2]-b->position[2],2))
#elif DIMENSIONS_==2
	#define DISTANCE( a, b) sqrt( pow(a->position[0]-b->position[0],2) + pow(a->position[1]-b->position[1],2))
#elif DIMENSIONS_==1
	#define DISTANCE( a, b) fabs( a->position[0]-b->position[0])
#endif*/

typedef double CONCENTRATION_T;

class VesselNode;
class VesselSegment;
class VoronoiDiagram;

class VesselSegment {
	public:
		char type;
		int index;
		VesselNode *vesselNodes[2];

		bool  radiusStatic;
		double radius;
		double countParallelVessels;
		double wallThickness;
		double permeability;

		double flow;
		double shear;
	public:
		VesselSegment(VesselNode *node1, VesselNode *node2);
		~VesselSegment();
		double getViscosity();
		void set( VesselNode *node1, VesselNode *node2){
			vesselNodes[0] = node1;
			vesselNodes[1] = node2;
		};
	};


class VesselNode {
		char type;
	public:
		int index;
		double position[3];
		VesselNode **neighbors;
		VesselSegment **branches;
		char countNeighbors;
		double pressure;
		double time;
		double marker;
	public:
		VesselNode(double x, double y, double z);
		~VesselNode();
		void addNeighbor(VesselNode *neighbor);
		void addNeighbor(VesselNode *neighbor, VesselSegment *vs);
		void removeNeighbor(VesselNode *neighbor);
		VesselNode *getNeighborUpstream( int n);
		VesselNode *getNeighborDownstream( int n);
		char getType();
		double getFlowFromNeighbor( int i);
		double getFlowToNeighbor( int i);
		void setType( char type);
		bool isNeighbor( VesselNode *neighbor){
			for( int n=0; n<countNeighbors; n++)
				if( neighbors[n] == neighbor)
					return true;
			return false;
		};
	};


class VesselGraph {
public:
	enum Types {RegularLattice, IrregularLattice};
private:
	Types _type;
	double _latticeConstant;
public:
	int dimensions;
	double LATTICE_CONSTANT;
	int DOMAIN_SIZE_X;
	int DOMAIN_SIZE_Y;
	int DOMAIN_SIZE_Z;
	double permeability;
	double diffusionCoefficient;

	VesselNode **vesselNodes;
	int countVesselNodes;
	int maxVesselNodes;
	
	VesselSegment **vesselSegments;
	int countVesselSegments;
	int maxVesselSegments;
	
	Octree<VesselNode*> *octree;

public:
	//VesselGraph();
	VesselGraph( int x, int y, int z, int dimensions);
	VesselGraph( const char *pFilename, Tumor *&tumor);
	~VesselGraph();
	Types type() { /*fprintf(stderr, "type();\n");*/ return _type;};
	void  setType( Types newType) { _type=newType;};
	double latticeConstant() { return _latticeConstant;};
	void addVesselNode( VesselNode *vn);
	void removeVesselNode( VesselNode *vn);
	void addVesselSegment(VesselSegment *vs);
	void addVesselSegmentAndNeighbors( VesselSegment *vs);
	void addVesselSegmentAndNeighbors( VesselSegment *vs, double radius);
	void removeVesselSegmentAndNeighbors( VesselSegment *vs);
	void removeVesselSegment( VesselSegment *vs);
	void removeVesselSegmentAndNodes( VesselSegment *vs);
	void printToVRML(const char *filename, char mode);
	void printToPovray(const char *filename, char mode);
	void printToPovray(const char *filename, CONCENTRATION_T ***marker);
	void printToPovray(const char *filename, char mode, double ***marker);
	void printToEPS(const char *filename, char mode, double ***marker);
	void writeToXML(const char *filename, double ***marker, Tumor *tumor);
	void writeStatisticsToFile(const char *filename, char mode);
	void updateFlow();
	void updateShear();
	void updatePressure( VesselNode **borderVesselNodes, double *borderPressure, int borderSize);
	void updatePressureNEW();
	void updatePressureNEW2();
	void updateRadius();
	void updateSegmentTypes();
	void updateMarkerVessels( double dt, double **markerVesselsIn, double **markerVesselsOut);
	void updateMarkerVesselsAndInterstitialSpace( double dt, double **markerVesselsIn, double **markerVesselsOut, double ***markerIntSpaceIn, double ***markerIntSpaceOut);
	void updateMarkerVesselsAndInterstitialSpaceNEW( double dt, double **markerVesselsIn, double **markerVesselsOut, double ***markerIntSpaceIn, double ***markerIntSpaceOut);
	void updateMarkerVesselsAndInterstitialSpaceNEW2( double dt, double border, double ***markerIntSpaceIn, double ***markerIntSpaceOut, SparseMatrix<double> *&sA, double *&b, double *&x);
	void updateMarkerVesselsAndInterstitialSpaceNEW3( double dt, double border, double ***markerIntSpaceIn, double ***markerIntSpaceOut);
	void updateMarkerVesselsAndInterstitialSpaceNEW4( double dt, double border, double ***markerIntSpaceIn, double ***markerIntSpaceOut);
	void updateMarkerVesselsAndInterstitialSpaceNEW5( double dt, double border, CONCENTRATION_T ***markerIntSpaceIn, CONCENTRATION_T ***markerIntSpaceOut);

	void updateMarkerVesselsExplicit( double dt, double border, double *&x);
	void updateMarkerVesselsExplicitRefined( double dt, double border, double *&dx, double *&ds, double *&s, char schema, int refinement);
	void updateMarkerVesselsExplicitAntiDiffusion( double dt, double border, double *&x);
	void updateMarkerVesselsExplicit( double dt, double border, double *&dx, double *&x, int);
	void updateMarkerVesselsAndInterstitialExplicit( double dt, double border, CONCENTRATION_T *&du, CONCENTRATION_T *&u, int countVesselSegmentSubNodes, CONCENTRATION_T ***markerIntSpace, int*** cells, bool UpdateInterstitialCompartment = true, bool UpdatePlasmaCompartment = true);
	void updateMarkerVesselsAndInterstitialExplicit( double dt, double border, CONCENTRATION_T *&du, CONCENTRATION_T *&u, int countVesselSegmentSubNodes, CONCENTRATION_T ***markerIntSpace, int*** cells, VoronoiDiagram *t);


	double getVascularVolume( int &x, int &y, int &z);
	double getExtraVascularVolume( int &x, int &y, int &z);
	double getExtraVascularVolumeNoRef( int x, int y, int z);
	double getVascularSurface( VesselNode *vn);
	double getVascularSurface( VesselSegment *vs);
	double getVascularVolume( VesselNode *vn);
	double getVascularVolume( VesselSegment *vs);
	double getExtraVascularVolume( VesselNode *vn);
	double getVascularPermeabilitySurfaceProduct( VesselNode *vn);
	double getVascularPermeabilitySurfaceProduct( VesselSegment *vs);
	double getVascularPermeabilitySurfaceProduct( int &x, int &y, int &z);

	void updateMarkerVesselsAndInterstitialSpaceBrix( double dt, double **markerVesselsIn, double **markerVesselsOut, double ***markerIntSpaceIn, double ***markerIntSpaceOut);

	void printMarkerIntensityToBinary(char *filename);
	void readMarkerIntensityFromBinary(char *filename);
	void printMarkerIntensityToBinary(char *filename, VesselGraph *vg, CONCENTRATION_T ***markerIntSpace, VoronoiDiagram *t, bool printAIF, CONCENTRATION_T markerAIF);
	void printMarkerIntensityMap(char *filename, VesselGraph *vg, CONCENTRATION_T ***markerIntSpace, double voxelSizeX, double voxelSizeY,	double voxelSizeZ, bool printAIF, double markerAIF);
	void printFlowMap(char *filename, VesselGraph *vg, double voxelSizeX, double voxelSizeY, double voxelSizeZ);
	void printPermeabilityMap(char *filename, VesselGraph *vg, double voxelSizeX, double voxelSizeY, double voxelSizeZ);
	void printBloodVolumeFractionMap(char *filename, VesselGraph *vg, double voxelSizeX, double voxelSizeY, double voxelSizeZ);
	void printInterstitialSpaceVolumeFractionMap(char *filename, VesselGraph *vg, double voxelSizeX, double voxelSizeY, double voxelSizeZ);
	void printVesselGraphToVTK(char *filename);

	VesselNode *getClosestVesselNode(double *pos, double &sqrDist);

	void setSingleVessel(int length, int height, int depth, double vessel_radius = MIN_VESSEL_RADIUS);
	void setSingleRandomVessel(int length, int height, int depth);
	void setSingleVessel2(int length);
	void setSymmetricVesselGraph(int height);
	void setInitialVesselGraph( int mode);

	void setInterTipConnections(Tumor *tumor, int &countInterTipConnections);
	void setInterTipConnectionsAll( int &countInterTipConnections);
	void removeInterTipConnections( int &countInterTipConnections);

	void simplifyNetwork( double minSegmentLength, double maxSegmentLength = FLT_MAX);

	void remodel(Tumor *tumor);

	double distance( VesselNode *a, VesselNode *b);
	double distanceSqr(VesselNode *a, VesselNode *b);

	void move( VesselNode *node, double x, double y, double z){
		octree->erase(node->position[0], node->position[1], node->position[2]);
		node->position[0]=x;
		node->position[1]=y;
		node->position[2]=z;
		octree->set(node->position[0], node->position[1], node->position[2], node);
	};


	class AIF {
	public:
		// types
		static const char CONSTANT = 1;
		static const char PARKER   = 2;
		static double ArterialInputFunction( double time, double time0, const char type);
	};
};

//#include "Vascularization.tcc"

#endif /*VESSELGRAPH_HPP*/
