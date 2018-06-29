#ifndef __VORONOI_DIAGRAM_EXTENDED_H
#define __VORONOI_DIAGRAM_EXTENDED_H

// DEFINES
#define FILENAMESIZE 512
#define READBUFFERSIZE 20000

#define MAX_TRIANGLES 6727210
#define MAX_NNS 30

#define INIT_FACES_ARRAY_SIZE 10

#define FACES_ARRAY_EXTENSION_SIZE 10
#define TETRAHEDRA_ARRAY_EXTENSION_SIZE 10
#define POINTS_ARRAY_EXTENSION_SIZE 100

#define DIMENSIONS 3
#define NR_FACE_POINTS        DIMENSIONS
#define NR_TETRAHEDRON_POINTS (DIMENSIONS + 1)

#define TRACK_DELAUNAY_REGION_NEIGHBORHOOD 1

#define LIMITED_NEIGHBOR_DISTANCE	3

#define POINT_POSITION_PRECISION	(1e-5)

//#define USE_AGENT

// TYPES
typedef double POSITION_T;
typedef double DISTANCE_T;

// DATA STRUCTURES
//typedef struct _Action Action;
//typedef struct _GridPoint GridPoint;
//typedef struct _VoronoiCell Vertex;
//typedef struct _Tetrahedron Tetrahedron;
//typedef struct _VoronoiDiagram Triangulation;
typedef struct _SphericalCoordinates SphericalCoordinates;
class Triangulation;

struct _SphericalCoordinates{
	double phi;		// 0 <= phi   <  2*PI
	double theta;	// 0 <= theta <= PI
};

// voronoi cell type
class GridPoint {
public:
	// methods
	//static newGridPoint();
	
	// index
	int index;
	
	// spatial position
	double position[DIMENSIONS];

	// neighbor points
	GridPoint **neighborPoints;
	int countNeighborPoints;
};

// voronoi cell type

class Vertex {
public:
	// methodes
	Vertex(double x, double y, double z);
	~Vertex();

	void actualizeExtendedNeighborhood( Triangulation *voronoiDiagram, int base_index, int generations, int* nr_kr_neighbors_gen, int radius);
	void actualizeFreeNeighbors();
	void checkFreeNeighbors();
	int isFree();
	int getState();
	void validate();
	double getDistanceTo( double[DIMENSIONS]);
	double getDistanceTo( Vertex* cell);
	bool isDomainBorder( Triangulation *vd);

	int index;

	// spatial position
	POSITION_T position[DIMENSIONS];

	// neighborhood
	char neighborCellsInitialized;
	int countNeighborCells;
	int countFreeNeighborCells;
	Vertex **neighborCells;

	// neighborhood
	char extendedNeighborCellsInitialized;
	int countExtendedNeighborCells;
	int countFreeExtendedNeighborCells;
	Vertex **extendedNeighborhood;

	// meta data
	void*    agent;
	bool	refined;
	int refinement;
	Vertex *coarseParent;

	//Action** actions;		// each cell knows its possible actions
	//int 	 actionsInitialized;
	double	glucose;
	double	oxygen;
	double	growthfactors;
	double	lactate;

	double  dglucose;
	double  doxygen;

	double	CPT11in;
	double	CPT11out;
	SphericalCoordinates sphericalCoordinates;
};

// tetrahedron type
class Tetrahedron {
public:
	// methodes
	int addTetrahedronNeighbor( Tetrahedron* neighborTet);
	void getCircumCenter( double center[DIMENSIONS]);

	int index;

	// vertices
	Vertex *vertices[DIMENSIONS+1];

	// neighbor tetrahedra
	int countNeighborTetrahedra;
	Tetrahedron *neighborTetrahedra[DIMENSIONS+1];

	// circumsphere
	char circumsphereInitialized;
	double circumsphere[DIMENSIONS+1];
	double radius;
//	double numericError;
};


// voronoi diagram type
//template <int DIMENSIONS>
class Triangulation {
public:
	// methodes
	Triangulation();
	~Triangulation();

	void setVoronoiGrid();
	void setDomain();
	void addVertex( Vertex* newVoronoiCell);
	void triangulate();

	void printToPovray( const char *filename, bool printPoints,  bool printFramePoints, bool printNeighborship, bool printTetrahedra, bool printConvexHull);

	void setPointsOnSphericalSurface( double radius, int N, double variance, double center[DIMENSIONS]);

	void getConvexHull( double thresholdDistance);
	void setFramePoints();

	int getCountConvexFaces();
	Vertex * getCentralCell();
//private:
	void sortTetrahedra();
	void setConvexHullNeighborhood();
	void NEWsetExtendedNeighborhood( int radius);
	void setExtendedNeighborhoodAroundVoronoiCell( int radius, int explorationRadius, Vertex *explorationCenter);
	void setExtendedNeighborhoodWithinSphere( int radius, double sphereRadius, double *sphereCenter);
	void setExtendedNeighborhood( double radius);
	void setExtendedNeighborhood( int cells);
	
	static Vertex *searchClosestVoronoiCell( Vertex *explorationCenter, double targetPosition[DIMENSIONS]);
	Vertex *searchClosestFreeVoronoiCell( Vertex *explorationCenter);
	//static Vertex *searchForVoronoiCell( int countIgnoredPoints, Vertex** ignoredPoints, int countSourcePoints, Vertex** sourcePoints, int &countExploredPoints, Vertex** exploredPoints);
	static Vertex *searchForVoronoiCell( int countIgnoredPoints, Vertex** ignoredPoints, int countSourcePoints, Vertex** sourcePoints, int &countExploredPoints, Vertex** exploredPoints);
	
	static Triangulation* newVoronoiDiagram();
	static Triangulation* newVoronoiDiagram( int x, int y, int z);
	static Triangulation* newVoronoiDiagramFromFile( char* filename);

	int readExtendedNeighborhoodToFile( char* filename, int radius);
	int writeExtendedNeighborhoodToFile( char* filename, int radius);

	bool tetrahedronContainsFramePoints( Tetrahedron* tet);
	Tetrahedron *getTetrahedronContainingPointInCircumSphere( Vertex* point);

	// voronoi cell
	int countVoronoiCells;
	int maxVoronoiCells;
	Vertex **vertices;
	
	// frame
	int countFramePoints;
	Vertex **framePoints;

	// tetrahedra
	int countTetrahedra;
	int maxTetrahedra;
	Tetrahedron **tetrahedra;

	// voronoi grid point
	char voronoiGridSet;
	GridPoint *voronoiGridPoints;

	// link between voronoi grid and voronoi cells
	GridPoint ***voronoiCellToVoronoiGridPoints;
	int        *voronoiCellCountVoronoiGridPoints;

	// faces
	int countFaces;
	int maxFaces;
	Vertex ***faces;
	
	// convex faces
	int countConvexFaces;
	int maxConvexFaces;
	Vertex ***convexFaces;

	// domain
	char domainSet;
	float boundaryThickness;
	float xMin[DIMENSIONS];
	float xMax[DIMENSIONS];
	int xN[DIMENSIONS];
	
	
	bool PointInsideConvexHull( double *p);
	void CheckAndChangeOrientationOfConvexHullTriangles();
};


bool PointInTriangle2D( double *p, double *a, double *b, double *c);




// FUNCTION PROTOTYPES

void getCircumCenter( Tetrahedron *tet, double center[DIMENSIONS]);

// functions for voronoi diagram
Triangulation* newVoronoiDiagram();
Triangulation* newVoronoiDiagramFromFile( char* filename);
void deleteVoronoiDiagram(  Triangulation* voronoiDiagram);
Vertex* findNearestCell3D( Triangulation* voronoiDiagram, Vertex* newVoronoiCell);
void insertVoronoiCell( Triangulation* voronoiDiagram, Vertex* newVoronoiCell);
void removeVoronoiCell( Triangulation* voronoiDiagram, Vertex* removedVoronoiCell);
void printPoints( const char * title, Vertex ** allPoints, int api );
void setTetrahedronNeighbors( Triangulation* voronoiDiagram);
void printTetrahedronNeighbors( Tetrahedron* tetrahedron);
void printTetrahedraNeighbors( Triangulation* voronoiDiagram);
Vertex * getCentralCell( Triangulation *voronoiDiagram);
void setDomain( Triangulation *voronoiDiagram);
double getDistanceOfPointToCircumsphere( Vertex* point, Tetrahedron* tet);
Tetrahedron *getTetrahedronContainingPoint( Triangulation *voronoiDiagram, Vertex* point);
double getCircumsphereRadius( Tetrahedron* tet);
int isElementOf( Vertex **pointList, int nrOfPoints, Vertex *point);

// functions for voronoi cell
//Vertex* newVoronoiCell( double x, double y, double z);
void deleteVoronoiCell( Vertex* voronoiCell);

// functions for tetrahedron
Tetrahedron* newTetrahedron();

// testing
void checkDelaunayCondition( Triangulation *voronoiDiagram, Vertex **newCell, int countNewCells);

void printToPovray( const char *filename, Triangulation *voronoiDiagram, Vertex **newCell, int countNewCells, int printPoints, int printNeighborship, int printTetrahedra, int printConvexHull);




#endif
