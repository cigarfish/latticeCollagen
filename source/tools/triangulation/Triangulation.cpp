#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#ifdef _OPENMP
    #include <omp.h>
#endif


// This should be moved to a central location ////////////////
#include <float.h>
#ifdef _WIN32
#define CS_ISNAN _isnan
#else
#define CS_ISNAN isnan
#endif
//////////////////////////////////////////////////////////////


#include "Triangulation.hpp"
#include "Mathematix.h"

// TEST FROM NICK

double **_A;
double **_B;
double *_x;

#define signof( a) (a==0 ? 0 : (a<0 ? -1 : 1))
#define dot2d(a,b) (a[0]*b[0] + a[1]*b[1])

void Tetrahedron::getCircumCenter( double center[DIMENSIONS])
{
	Tetrahedron *tet = this;
#if DIMENSIONS == 3
//        |d-a|^2 [(b-a)x(c-a)] + |c-a|^2 [(d-a)x(b-a)] + |b-a|^2 [(c-a)x(d-a)]
//m = a + ---------------------------------------------------------------------.
//                                | bx-ax  by-ay  bz-az |
//                              2 | cx-ax  cy-ay  cz-az |
//                                | dx-ax  dy-ay  dz-az |
	POSITION_T *a = tet->vertices[0]->position,
	       *b = tet->vertices[1]->position, 
	       *c = tet->vertices[2]->position, 
	       *d = tet->vertices[3]->position;
	
	double ba[DIMENSIONS], ca[DIMENSIONS], da[DIMENSIONS];
	vectorDifference3D( b, a, ba);
	vectorDifference3D( c, a, ca);
	vectorDifference3D( d, a, da);
	
	double bc[DIMENSIONS], cd[DIMENSIONS], db[DIMENSIONS];
	crossProduct3D( ba, ca, bc);
	crossProduct3D( ca, da, cd);
	crossProduct3D( da, ba, db);
	
	double bcd = 2. * dotProduct3D( ba, cd);

	vectorScale3D( bc, dotProduct3D( da, da), bc);
	vectorScale3D( cd, dotProduct3D( ba, ba), cd);
	vectorScale3D( db, dotProduct3D( ca, ca), db);
		
	center[0] = a[0] + (bc[0] + cd[0] + db[0])/bcd;
	center[1] = a[1] + (bc[1] + cd[1] + db[1])/bcd;
	center[2] = a[2] + (bc[2] + cd[2] + db[2])/bcd;

#elif DIMENSIONS == 2

	   double *a = tet->vertices[0]->position,
	           *b = tet->vertices[1]->position,
	           *c = tet->vertices[2]->position;

	    double D = 2.*( a[0] * (b[1] - c[1]) + b[0] * ( c[1] - a[1] ) + c[0] * ( a[1] - b[1] ) );

	    center[0]  = ( a[1]*a[1] + a[0]*a[0] )*(b[1] - c[1]);
	    center[0] += ( b[1]*b[1] + b[0]*b[0] )*(c[1] - a[1]);
	    center[0] += ( c[1]*c[1] + c[0]*c[0] )*(a[1] - b[1]);
	    center[0] /= D;

	    center[1]  = ( a[1]*a[1] + a[0]*a[0] )*(c[0] - b[0]);
	    center[1] += ( b[1]*b[1] + b[0]*b[0] )*(a[0] - c[0]);
	    center[1] += ( c[1]*c[1] + c[0]*c[0] )*(b[0] - a[0]);
	    center[1] /= D;
#endif
}


GridPoint *newGridPoint()
{
	return (GridPoint *) malloc( sizeof(GridPoint));
}


void Triangulation::setVoronoiGrid()
{
	fprintf( stderr, "Triangulation::setVoronoiGrid()\n");

	// init points
	this->voronoiGridPoints = (GridPoint *) malloc( this->countTetrahedra * sizeof(GridPoint));
	
	// init links
	this->voronoiCellToVoronoiGridPoints = (GridPoint ***) malloc( this->countVoronoiCells * sizeof(GridPoint**));
	this->voronoiCellCountVoronoiGridPoints = (int *) malloc( this->countVoronoiCells * sizeof(int));
	for( int i=0; i<this->countVoronoiCells; i++){
		this->voronoiCellToVoronoiGridPoints[i] = NULL;
		this->voronoiCellCountVoronoiGridPoints[i] = 0;
	}

	// set points
	for( int i=0; i<this->countTetrahedra; i++){
		// set index
		this->voronoiGridPoints[i].index = i;
		
		// set point position
		this->tetrahedra[i]->getCircumCenter( this->voronoiGridPoints[i].position);

		// set neighborhoods of points
		// ...
	}
	
	// set links
	for( int t=0; t<this->countTetrahedra; t++){
		for( int v=0; v<DIMENSIONS+1; v++){
			// for each tet vertice
			int i = this->tetrahedra[t]->vertices[v]->index;
			int l;
			for( l=0; l<this->voronoiCellCountVoronoiGridPoints[i] && this->voronoiCellToVoronoiGridPoints[i][l]!=&(this->voronoiGridPoints[t]); l++);
			if( l==this->voronoiCellCountVoronoiGridPoints[i]){
				this->voronoiCellCountVoronoiGridPoints[i]++;
				this->voronoiCellToVoronoiGridPoints[i] = (GridPoint **)realloc( this->voronoiCellToVoronoiGridPoints[i], this->voronoiCellCountVoronoiGridPoints[i]*sizeof(GridPoint*));
				this->voronoiCellToVoronoiGridPoints[i][l] = &(this->voronoiGridPoints[t]);
			}
		}
	}

	// test output of links
	/*for( int i=0; i<this->countVoronoiCells; i++){
		fprintf(stderr, "%i: ", i);
		for( int l=0; l<this->voronoiCellCountVoronoiGridPoints[i]; l++)
			fprintf(stderr, "%i ", this->voronoiCellToVoronoiGridPoints[i][l]->index);
		fprintf(stderr, "\n");
	}*/

	
}

void Vertex::validate()
{
/*	int i;
	
	int countFreeNeighbors = 0;
	for( i=0; i<this->countNeighborCells; i++){
		if( this->neighborCells[i]->isFree())
			countFreeNeighbors++;
	}
	if( countFreeNeighbors!=this->countFreeNeighborCells){
		fprintf( stderr, "ERROR: Vertex %i has wrong number of free neighbors: %i (correct: %i)\n", this->index, this->countFreeNeighborCells, countFreeNeighbors);
		exit( 0);
	}

	countFreeNeighbors = 0;
	for( i=0; i<this->countExtendedNeighborCells; i++){
		if( this->extendedNeighborhood[i]->isFree())
			countFreeNeighbors++;
	}
	if( countFreeNeighbors!=this->countFreeExtendedNeighborCells){
		if( !this->isFree())
			fprintf( stderr, "ERROR: Vertex %i (agent %i) has wrong number of free extended neighbors: %i/%i (correct: %i)\n",
			         this->index, GetAgent( this)->index, this->countFreeExtendedNeighborCells, this->countExtendedNeighborCells, countFreeNeighbors);
		else{
			fprintf( stderr, "ERROR: Vertex %i has wrong number of free extended neighbors: %i/%i (correct: %i)\n",
			         this->index, this->countFreeExtendedNeighborCells, this->countExtendedNeighborCells, countFreeNeighbors);			
		}

			for( i=0; i<this->countExtendedNeighborCells; i++){
				int found = FALSE;
				int ii;
				for( ii=0; ii<this->extendedNeighborhood[i]->countExtendedNeighborCells; ii++){
					if( this->extendedNeighborhood[i]->extendedNeighborhood[ii] == this){
						found = TRUE;
					//fprintf( stderr, "ERROR: Vertex %i -> extended neighbor %i -> extended neighbor %i :-)\n",
					//         this->index, this->extendedNeighborhood[i]->index, this->extendedNeighborhood[i]->extendedNeighborhood[ii]->index);
					}
				}			
				if( !found){
					fprintf( stderr, "ERROR: Vertex %i has a extended neighbor %i which doesn't contain %i as extended neighbor!\n", this->index, this->extendedNeighborhood[i]->index, this->index);
				}
			}

		exit( 0);
	}

*/
}
/*****************************************************************************/

#define FREE	0


#ifdef USE_AGENT
int Vertex::getState()
{
	return (this->agent==NULL ? FREE : GetAgent( this)->state);	
}
/*****************************************************************************/


int Vertex::isFree()
{
	return (GetAgent( this)==NULL ? TRUE : GetAgent( this)->isFree());	
}
#endif
/*****************************************************************************/



int isElementOf( Vertex **pointList, int nrOfPoints, Vertex *point){
	int i;

	//printf("{to search point:%i,list(%i): ", point->nr, nrOfPoints);

	for( i=0; i<nrOfPoints; i++){
		//printf("%i, ", pointNrList[i]);
		if( pointList[i] == point) {
			//printf("[found]}");
			return 1;
		}
	}
	//printf("}");

	return 0;
}
/*****************************************************************************/



Vertex *Triangulation::searchForVoronoiCell( int countIgnoredPoints, Vertex** ignoredPoints, int countSourcePoints, Vertex** sourcePoints, int &countExploredPoints, Vertex** exploredPoints){
	int i, j;

	/*fprintf( stderr, "Ignored Points: [");
	for( j=0; j<countIgnoredPoints; j++)
		fprintf( stderr, "%i ", ignoredPoints[j]->index);
	fprintf( stderr, "\n");

	fprintf( stderr, "Source Points: [");
	for( j=0; j<countSourcePoints; j++)
		fprintf( stderr, "%i ", sourcePoints[j]->index);
	fprintf( stderr, "\n");
	*/

	for( i = 0; i < countSourcePoints; i++){ // for every element in the previous list
		for( j = 0; j < sourcePoints[i]->countNeighborCells; j++){ // for every neighbor of a element in the previous list
        //printf("  candidates:%i\n", (voronoiDiagram[ point->reachableNeighbors[ base_index + i ] ].neighbors[j])->nr);
			Vertex *candidate = sourcePoints[i]->neighborCells[j];
			if( !isElementOf( ignoredPoints,  // pointer to list of all generations
			                  countIgnoredPoints,     // number of points in the list
			                  candidate) &&
			    !isElementOf( sourcePoints,  // pointer to list of all generations
			                  countSourcePoints,     // number of points in the list
			                  candidate) &&
			    !isElementOf( exploredPoints,  // pointer to list of all generations
			                  countExploredPoints,     // number of points in the list
			                  candidate)
			    ){                                                                       // actual neighbor of neighbor 
			    exploredPoints[countExploredPoints] = candidate;
				countExploredPoints++;
			}
		}
	}
  
	return NULL;
}
/*****************************************************************************/



void Vertex::actualizeExtendedNeighborhood( Triangulation *voronoiDiagram, int base_index, int generations, int* nr_kr_neighbors_gen, int radius){
	int i, j;
	int count = 0;
	int base_gen = generations - 1; // number of last gotten generation

	//printf("get_grid_points_below_kr( Grid *point(%i), base_index(%i), generations(%i))\n", point->nr, base_index, generations);

	//fprintf( stderr, "%ith generation:{", base_gen+1);
	for( i = 0; i < nr_kr_neighbors_gen[ base_gen ]; i++){ // for every element in the previous list
		/*double prevCandidateRadius = 
			sqrt( pow(this->position[0]-this->extendedNeighborhood[ base_index + i ]->position[0], 2) + 
			      pow(this->position[1]-this->extendedNeighborhood[ base_index + i ]->position[1], 2) + 
			      pow(this->position[2]-this->extendedNeighborhood[ base_index + i ]->position[2], 2));*/
		//printf("  point:%i\n", voronoiDiagram[ point->reachableNeighbors[ base_index + i ] ].nr);
		for( j = 0; j < this->extendedNeighborhood[ base_index + i ]->countNeighborCells; j++){ // for every neighbor of a element in the previous list
        //printf("  candidates:%i\n", (voronoiDiagram[ point->reachableNeighbors[ base_index + i ] ].neighbors[j])->nr);
			Vertex *candidate = this->extendedNeighborhood[ base_index + i ]->neighborCells[j];
			/*double candidateRadius = sqrt( pow( this->position[0] - candidate->position[0], 2) + 
            	                           pow( this->position[1] - candidate->position[1], 2) + 
                	                       pow( this->position[2] - candidate->position[2], 2));
			if( candidateRadius < radius &&*/
			if( generations-1 <  radius && 
			    !isElementOf( this->extendedNeighborhood,  // pointer to list of all generations
			                  base_index + nr_kr_neighbors_gen[ base_gen ] + count,     // number of points in the list
			                  candidate)){                                                                       // actual neighbor of neighbor 
			    /*!isElementOf( &(this->extendedNeighborhood[ base_index - nr_kr_neighbors_gen[ base_gen - 1 ]]),// pointer to list of last, actual and new generation
			                  nr_kr_neighbors_gen[ base_gen - 1 ] + nr_kr_neighbors_gen[ base_gen ] + count,     // number of points in the descendant list
			                  candidate)){    */                                                                   // actual neighbor of neighbor 

				//fprintf( stderr, " %i,", candidate->index);
				this->extendedNeighborhood[base_index + nr_kr_neighbors_gen[ base_gen ] + count] = candidate;
				/*if( candidate->countExtendedNeighborCells==0)
					candidate->extendedNeighborhood = (Vertex**) realloc( candidate->extendedNeighborhood, (voronoiDiagram->countVoronoiCells) * sizeof(Vertex*));
				//candidate->extendedNeighborhood = (Vertex**) realloc( candidate->extendedNeighborhood, (candidate->countExtendedNeighborCells + 1) * sizeof(Vertex*));
				candidate->extendedNeighborhood[ (candidate->countExtendedNeighborCells)++] = this;*/
			
				count++;
			}
		}
	}
	//fprintf( stderr, "\b}\n");

	if( count != 0){

		/* Update depending parameters */
		generations++;                             // number of set generations
		this->countExtendedNeighborCells           += count; // actualize number of neighbors
		//point->countReachableFreeNeighbors       += count; // actualize number of free neighbors
		nr_kr_neighbors_gen[base_gen + 1] = count; // set number of neighbors

		/* get next generation */
		this->actualizeExtendedNeighborhood( voronoiDiagram, base_index + nr_kr_neighbors_gen[ base_gen ], generations, nr_kr_neighbors_gen, radius);
	} 
  
	return;
}
/*****************************************************************************/

char *long2str( char *str, int number)
{
	int base = 128 - 3;
	int i = 0;
	// 4 = [100]
	int mod;
	
	do{
		mod = number%base;
		number = number/base;
		str[i] = mod+3;
		//fprintf( stderr, "%i|", mod+3);
		i++;
	}while(number>0);
	str[i] = '\0';
	
	// resort
	int ii;
	int temp;
	for( ii=0; ii<i/2; ii++){
		temp = str[ii];
		str[ii] = str[i-ii-1];
		str[i-ii-1] = temp;
	}
	
	return str;
}

int str2long( char *str)
{
	int base = 128 - 3;
	int i=0;
	int number = 0;
	
	while( str[i]!='\0'){
		number *= base;
		number += (str[i]-3);
		i++;
	}
	
	return number;

}

int strlen( char *str)
{
	int i=0;
	
	while( str[i]!='\0'){
		i++;
	}
	
	return i;
}

int Triangulation::readExtendedNeighborhoodToFile( char* filename, int radius){
	FILE *fp;
	char completeFilename[ FILENAMESIZE ];
	char buffer[ READBUFFERSIZE ];
    char* ptr;

	int i, ii;
	//int hash = 0;
	
	sprintf( completeFilename,"%s.k%i", filename, radius);
	
    fp = fopen( completeFilename, "r" );
    if( fp == NULL ){
      fprintf( stderr, "Error opening file %s\n", completeFilename);
      return FALSE;
    }
    else fprintf( stderr, ">%s eingelesen\n", completeFilename);

	char c = fgetc( fp );
	char word[ 10 ];
	int  word_length = 0;
	int countNeighbors = -1;
	i = 0;
	while( c != EOF){

		word[word_length++] = c;
		switch(c){
			case 0: // TAB
				//fprintf( stderr, "_____");
				word[word_length-1] = '\0';
				ptr = word;
				countNeighbors = str2long( word );
				this->vertices[i]->extendedNeighborhood = (Vertex**) calloc( countNeighbors, sizeof( Vertex*));
				assert (this->vertices[i]->extendedNeighborhood);
				this->vertices[i]->countExtendedNeighborCells = countNeighbors;
				this->vertices[i]->countFreeExtendedNeighborCells = countNeighbors;
				this->vertices[i]->extendedNeighborCellsInitialized = TRUE;
				word_length = 0;
				//fprintf( stderr, "%i\t", countNeighbors);

				break;	
				
			case 1: // ENTER
				//fprintf( stderr, "\n");	
				word_length = 0;
				i++;
#if _COMMENTS_ > 2
				fprintf( stderr, "\rRead Extended Neighborhood %.3lf%%                              \b", 100.*(i+1.)/this->countVoronoiCells);
#endif
				break;	
				
			case 2: // SPACE
				//fprintf( stderr, "_");		
				word[word_length-1] = '\0';
				ptr = word;
				this->vertices[i]->extendedNeighborhood[--countNeighbors] = this->vertices[str2long( word )];
				word_length = 0;
				//fprintf( stderr, "%i ", this->vertices[i]->extendedNeighborhood[countNeighbors]->index);
				break;	
				
			default:
			//fprintf( stderr, "%i", c-48);
			break;
		}
		//word_length++;
		c = fgetc( fp );
		//if( c!='\n');
		//fprintf( stderr, "|%i", c);
	}
	fclose( fp);
	return TRUE;

	// OLD NEW
	while( c != EOF){
		word[word_length++] = c;
		switch(c){
			case 9: // TAB
				//fprintf( stderr, "_____");
				word[word_length-1] = '\0';
				ptr = word;
				countNeighbors = strtol( ptr, &ptr, 0 );
				this->vertices[i]->extendedNeighborhood = (Vertex**) calloc( countNeighbors, sizeof( Vertex*));
				assert(this->vertices[i]->extendedNeighborhood);
				this->vertices[i]->countExtendedNeighborCells = countNeighbors;
				this->vertices[i]->countFreeExtendedNeighborCells = countNeighbors;
				word_length = 0;
				//fprintf( stderr, "%i\t", countNeighbors);

				break;	
				
			case 10: // ENTER
				//fprintf( stderr, "\n");	
				word_length = 0;
				i++;
				break;	
				
			case 32: // SPACE
				//fprintf( stderr, "_");		
				word[word_length-1] = '\0';
				ptr = word;
				this->vertices[i]->extendedNeighborhood[--countNeighbors] = this->vertices[strtol( ptr, &ptr, 0 )];
				word_length = 0;
				//fprintf( stderr, "%i ", this->vertices[i]->extendedNeighborhood[countNeighbors]->index);
				break;	
				
			default:
			//fprintf( stderr, "%i", c-48);
			break;
		}
		//word_length++;
		c = fgetc( fp );
		//if( c!='\n');
		//fprintf( stderr, "|%i", c);
	}
	
	return TRUE;
	
	// OLD
	for( i=0; i<this->countVoronoiCells; i++){
		// read line
		if( fgets( buffer, READBUFFERSIZE, fp ) == NULL ){
			perror( "Error reading line" );
			return FALSE;
		}

		ptr = buffer;

		// read number of extended neighbors
		int countNeighbors = strtol( ptr, &ptr, 0 );
		this->vertices[i]->extendedNeighborhood = (Vertex**) calloc( countNeighbors, sizeof( Vertex*));
		assert(this->vertices[i]->extendedNeighborhood);
		this->vertices[i]->countExtendedNeighborCells = countNeighbors;
		this->vertices[i]->countFreeExtendedNeighborCells = countNeighbors;

		// read all extended neighbors
		for( ii=0; ii<countNeighbors; ii++){
			int index = strtol( ptr, &ptr, 0 );
			//hash += index;
			this->vertices[i]->extendedNeighborhood[ii] = this->vertices[index];
		}
	}

	fclose( fp);

	//fprintf( stderr, "hash: %i\n", hash);

	return TRUE;

}
/*****************************************************************************/



int Triangulation::writeExtendedNeighborhoodToFile( char* filename, int radius){
	FILE *fp;
	char completeFilename[ FILENAMESIZE ];
	int i, ii;
	//int hash = 0;
	char buffer[ 10];
	int byte = 0;
	
	sprintf( completeFilename,"%s.k%i", filename, radius);
	
	fp = fopen( completeFilename, "w");
	if(fp==NULL){
		fprintf(stderr, "Error opening file %s for writing!\n",completeFilename);
		return FALSE;
	}//else
		//fprintf(stderr, "Write file %s\n",completeFilename);

	for( i=0; i<this->countVoronoiCells; i++){
		byte += 2 + strlen( long2str( buffer, this->vertices[i]->countExtendedNeighborCells));
		fprintf( fp, "%s%c", buffer, 0);
		for( ii=0; ii<this->vertices[i]->countExtendedNeighborCells; ii++){
			byte += 1 + strlen( long2str( buffer, this->vertices[i]->extendedNeighborhood[ii]->index));
			fprintf( fp, "%s%c", buffer, 2);
			//hash += this->vertices[i]->extendedNeighborhood[ii]->index;
		}
		fprintf( fp, "%c", 1);
#if _COMMENTS_ > 2
		fprintf( stderr, "\rWrite Extended Neighborhood to File %.3lf%% (%.1lf MBytes)\t\t\b", 100.*(i+1.)/this->countVoronoiCells, (double)this->countVoronoiCells*(double)byte/1024./1024./(i+1.));
#endif

	}
	fprintf( stderr, "\n");
	
/*	for( i=0; i<this->countVoronoiCells; i++){
		fprintf( fp, "%i\t", this->vertices[i]->countExtendedNeighborCells);
		for( ii=0; ii<this->vertices[i]->countExtendedNeighborCells; ii++){
			fprintf( fp, "%i ", this->vertices[i]->extendedNeighborhood[ii]->index);
			hash += this->vertices[i]->extendedNeighborhood[ii]->index;
		}
		fprintf( fp, "\n");
	}
*/
	fclose( fp);
	
	//fprintf( stderr, "hash: %i\n", hash);

	return TRUE;

}
/*****************************************************************************/



void Triangulation::setExtendedNeighborhoodWithinSphere( int radius, double sphereRadius, double *sphereCenter){

	Triangulation *voronoiDiagram = this;
	int i, j, count;
	Vertex **tempNeighborhoodList;
	int generations;
	int countNeighborsGen[30];
	//int TEST_countExtendedNeighbors = 0;

	//printf("radius=%lf\n", radius);

	Vertex** tempExtendedNeighborhoodList = ( Vertex**) calloc ( voronoiDiagram->countVoronoiCells, sizeof( Vertex*));
	assert(tempExtendedNeighborhoodList);

	for( i = 0; i < voronoiDiagram->countVoronoiCells; i++){
		
		double distanceToSphereCenter = 0.;
		for( j=0; j<DIMENSIONS; j++)
			distanceToSphereCenter += pow( voronoiDiagram->vertices[i]->position[j] - sphereCenter[j], 2);
		distanceToSphereCenter = sqrt( distanceToSphereCenter);
		if( distanceToSphereCenter<sphereRadius){

#if _COMMENTS_ > 2
			fprintf( stderr, "\rCalculate Extended Neighborhood %.3lf%%\t\t\b", 100.*(i+1.)/voronoiDiagram->countVoronoiCells);
#endif
			tempNeighborhoodList = voronoiDiagram->vertices[i]->extendedNeighborhood;
			voronoiDiagram->vertices[i]->extendedNeighborhood = tempExtendedNeighborhoodList;


			// 0st generation (Set point itself)
			voronoiDiagram->vertices[i]->extendedNeighborhood[0]  = voronoiDiagram->vertices[i];        // neighbors in 0th generation
			countNeighborsGen[0] = 1;                   // number of neighbors in 0th generation
			//fprintf( stderr, "0st generation:{ %i}\n", voronoiDiagram->vertices[i]->extendedNeighborhood[0]->index);


			// 1st generation (Set direct neighbors)
			count = 0;
			//fprintf( stderr, "1st generation:{");
			for(j = 1; j <= voronoiDiagram->vertices[i]->countNeighborCells; j++){
				//printf("neighbor:%i, ", voronoiDiagram->vertices[i].neighbors[j-1]->nr);
				if(1 <= radius){  
					voronoiDiagram->vertices[i]->extendedNeighborhood[j] = voronoiDiagram->vertices[i]->neighborCells[j-1]; // neighbors in 1st generation
					//fprintf( stderr, " %i,", voronoiDiagram->vertices[i]->extendedNeighborhood[j]->index);
					count++;
				}
			}
			countNeighborsGen[1] = count;       // number of neighbors in 1st generation
			//printf("\ncountNeighborsGen[0]:%i, countNeighborsGen[1]:%i\n", countNeighborsGen[0], countNeighborsGen[1]);

			/* Initialize depending parameters*/
			generations                 = 2;  // number of set generations
			voronoiDiagram->vertices[i]->countExtendedNeighborCells   = 1+count;  // point itself and direct

			/* get all point below a certain radius */
			voronoiDiagram->vertices[i]->actualizeExtendedNeighborhood( this, 1, generations, countNeighborsGen, radius); // beginning with the 2nd generation of neighbors



			voronoiDiagram->vertices[i]->countExtendedNeighborCells   -= 1+voronoiDiagram->vertices[i]->countNeighborCells;  // point itself and direct neighbors aren't important
			voronoiDiagram->vertices[i]->countFreeExtendedNeighborCells = voronoiDiagram->vertices[i]->countExtendedNeighborCells;  // point itself and direct neighbors aren't important

			/* Delete point itself and direct neighbors from list */
			tempNeighborhoodList = ( Vertex**) calloc ( voronoiDiagram->vertices[i]->countExtendedNeighborCells + 1, sizeof( Vertex*)); // memoryallocation
			assert(tempNeighborhoodList);
			for(j = 0; j < voronoiDiagram->vertices[i]->countExtendedNeighborCells; j++)
				tempNeighborhoodList[j] = voronoiDiagram->vertices[i]->extendedNeighborhood[1 + voronoiDiagram->vertices[i]->countNeighborCells + j ];                //

			//free( voronoiDiagram->vertices[i]->extendedNeighborhood);		// free unused allocated memory
			voronoiDiagram->vertices[i]->extendedNeighborhood = tempNeighborhoodList; // set pointer
		}
		else{
			voronoiDiagram->vertices[i]->countExtendedNeighborCells = 0;
			voronoiDiagram->vertices[i]->countFreeExtendedNeighborCells = 0;
			voronoiDiagram->vertices[i]->extendedNeighborhood = NULL;
		}
	}
	free( tempExtendedNeighborhoodList);
}
/*****************************************************************************/



void Triangulation::setExtendedNeighborhoodAroundVoronoiCell( int radius, int explorationRadius, Vertex *explorationCenter){

	Triangulation *voronoiDiagram = this;
	int j, count;
	Vertex **tempNeighborhoodList;
	int generations;
	int countNeighborsGen[30];
	
	double test_countCells = 0.;
	double test_countNeighbors = 0.;
	//int TEST_countExtendedNeighbors = 0;

	//printf("radius=%lf\n", radius);

	// STACKS
	Vertex **voronoiCellQueue;
	int          *depthQueue;
	//int maxQueueLength = (int)pow( explorationRadius*2, DIMENSIONS);
	//maxQueueLength = (maxQueueLength>this->countVoronoiCells ? this->countVoronoiCells : maxQueueLength);
	int maxQueueLength = this->countVoronoiCells;
	int queueLength    = 0;
	int firstQueueElement = 0;
	
	// allocation
	voronoiCellQueue = (Vertex **) calloc( maxQueueLength, sizeof( Vertex *));
	assert(voronoiCellQueue);
	depthQueue       = (int *)          calloc( maxQueueLength, sizeof( int));
	assert(depthQueue);
	printf("maxQueueLength=%i\n", maxQueueLength);

	// initialization
	voronoiCellQueue[queueLength] = explorationCenter;
	depthQueue[queueLength]       = 0;
	queueLength++;
	depthQueue[(firstQueueElement+maxQueueLength-1)%maxQueueLength] = 0;
	Vertex** tempExtendedNeighborhoodList = ( Vertex**) calloc ( voronoiDiagram->countVoronoiCells, sizeof( Vertex*));
	assert(tempExtendedNeighborhoodList);

	do{
		//fprintf( stderr, "first: %i, queueLength: %i, maxQueueLength: %i\n", firstQueueElement, queueLength, maxQueueLength);

		
		// actualize extended neighborhood of actual queue element
		//if( depthQueue[firstQueueElement]<=explorationRadius){

			//fprintf( stderr, "Process Vertex %i, Depth: %i => first: %i, queueLength: %i, maxQueueLength: %i\n", voronoiCellQueue[firstQueueElement]->index, depthQueue[firstQueueElement], firstQueueElement, queueLength, maxQueueLength);
#if _COMMENTS_ > 2
			if( firstQueueElement>0 && depthQueue[firstQueueElement-1]<depthQueue[firstQueueElement])
			fprintf( stderr, "\rCalculate Extended Neighborhood %.3lf%%\t\t\b", 100.*(depthQueue[firstQueueElement])/explorationRadius);
#endif
			tempNeighborhoodList = voronoiCellQueue[firstQueueElement]->extendedNeighborhood;
			voronoiCellQueue[firstQueueElement]->extendedNeighborhood = tempExtendedNeighborhoodList; 


			// 0st generation (Set point itself)
			voronoiCellQueue[firstQueueElement]->extendedNeighborhood[0]  = voronoiCellQueue[firstQueueElement];        // neighbors in 0th generation
			countNeighborsGen[0] = 1;                   // number of neighbors in 0th generation
			//fprintf( stderr, "0st generation:{ %i}\n", voronoiCellQueue[firstQueueElement]->extendedNeighborhood[0]->index);


			// 1st generation (Set direct neighbors)
			count = 0;
			//fprintf( stderr, "1st generation:{");
			for(j = 1; j <= voronoiCellQueue[firstQueueElement]->countNeighborCells; j++){
				//printf("neighbor:%i, ", voronoiCellQueue[firstQueueElement].neighbors[j-1]->nr);
				if(1 <= radius){  
					voronoiCellQueue[firstQueueElement]->extendedNeighborhood[j] = voronoiCellQueue[firstQueueElement]->neighborCells[j-1]; // neighbors in 1st generation
					//fprintf( stderr, " %i,", voronoiCellQueue[firstQueueElement]->extendedNeighborhood[j]->index);
					count++;
					
					// add neighbors to queue
					//printf( "\r%i \b", voronoiCellQueue[firstQueueElement]->extendedNeighborhood[j]->countExtendedNeighborCells);
					//printf( "\r%i \b", depthQueue[firstQueueElement]);
					if(  voronoiCellQueue[firstQueueElement]->extendedNeighborhood[j]->countExtendedNeighborCells==0
					     && depthQueue[firstQueueElement]<explorationRadius){
						// mark as add to queue
						voronoiCellQueue[firstQueueElement]->extendedNeighborhood[j]->countExtendedNeighborCells = -1;
						
						// add to queue
						voronoiCellQueue[(firstQueueElement+queueLength)%maxQueueLength] = voronoiCellQueue[firstQueueElement]->extendedNeighborhood[j];
						depthQueue[(firstQueueElement+queueLength)%maxQueueLength]       = depthQueue[firstQueueElement]+1;
						//fprintf( stderr, "Add To Queue: %i, Depth: %i\n", voronoiCellQueue[(firstQueueElement+queueLength)%maxQueueLength]->index, depthQueue[(firstQueueElement+queueLength)%maxQueueLength]);
						queueLength++;
						
						// memory allocation check
						if( queueLength==maxQueueLength){
							fprintf( stderr, "INFO: Queue reached maximal length!\n");
							exit( 0);
						}
					}
				}
			}
			countNeighborsGen[1] = count;       // number of neighbors in 1st generation
			//printf("\ncountNeighborsGen[0]:%i, countNeighborsGen[1]:%i\n", countNeighborsGen[0], countNeighborsGen[1]);

			/* Initialize depending parameters*/
			generations                 = 2;  // number of set generations
			voronoiCellQueue[firstQueueElement]->countExtendedNeighborCells   = 1+count;  // point itself and direct

			/* get all point below a certain radius */
			voronoiCellQueue[firstQueueElement]->actualizeExtendedNeighborhood( this, 1, generations, countNeighborsGen, radius); // beginning with the 2nd generation of neighbors



			voronoiCellQueue[firstQueueElement]->countExtendedNeighborCells   -= 1+voronoiCellQueue[firstQueueElement]->countNeighborCells;  // point itself and direct neighbors aren't important
			voronoiCellQueue[firstQueueElement]->countFreeExtendedNeighborCells = voronoiCellQueue[firstQueueElement]->countExtendedNeighborCells;  // point itself and direct neighbors aren't important

			/* Delete point itself and direct neighbors from list */
			tempNeighborhoodList = ( Vertex**) calloc ( voronoiCellQueue[firstQueueElement]->countExtendedNeighborCells, sizeof( Vertex*)); // memoryallocation
			if( tempNeighborhoodList==NULL){
				fprintf( stderr, "ERROR: Couldn't allocate the memory for the extended neighborhood list of Voronoi cell %i!\n", voronoiCellQueue[firstQueueElement]->index);	
				fprintf( stderr, "Voronoi cell %i: depth: %i, desired neighborhood list length = %i elements!\n", voronoiCellQueue[firstQueueElement]->index, depthQueue[firstQueueElement], voronoiCellQueue[firstQueueElement]->countExtendedNeighborCells);	
				fprintf( stderr, "average number of extended neighbors = %lf\n", test_countNeighbors/test_countCells);	
			}
			for(j = 0; j < voronoiCellQueue[firstQueueElement]->countExtendedNeighborCells; j++)
				tempNeighborhoodList[j] = voronoiCellQueue[firstQueueElement]->extendedNeighborhood[1 + voronoiCellQueue[firstQueueElement]->countNeighborCells + j ];                // 

			//free( voronoiCellQueue[firstQueueElement]->extendedNeighborhood);		// free unused allocated memory
			voronoiCellQueue[firstQueueElement]->extendedNeighborhood = tempNeighborhoodList; // set pointer
			voronoiCellQueue[firstQueueElement]->extendedNeighborCellsInitialized = TRUE;
	
			test_countCells++;
			test_countNeighbors += voronoiCellQueue[firstQueueElement]->countExtendedNeighborCells;
			
			firstQueueElement = (firstQueueElement+1)%maxQueueLength;
			queueLength--;
		//}
		/*else{
			voronoiCellQueue[firstQueueElement]->countExtendedNeighborCells = 0;
			voronoiCellQueue[firstQueueElement]->countFreeExtendedNeighborCells = 0;
			voronoiCellQueue[firstQueueElement]->extendedNeighborhood = NULL;
		}*/
 		
	}while( queueLength!=0);
	
	free( tempExtendedNeighborhoodList);
	free( voronoiCellQueue);
	free( depthQueue);
}
/*****************************************************************************/



Vertex *Triangulation::searchClosestFreeVoronoiCell( Vertex *explorationCenter){

/*	int j;
	Vertex* tempExplorationList[this->countVoronoiCells];// = ( Vertex**) calloc ( this->countVoronoiCells, sizeof( Vertex*));
	//Vertex** tempExplorationList = ( Vertex**) calloc ( this->countVoronoiCells, sizeof( Vertex*));
	int countExploredCells=0, countIgnoredPoints=0, countSourcePoints=0, countExploredPoints=0;
	int indexIgnoredPoints, indexSourcePoints, indexExploredPoints;

	// 0st generation (Set point itself)
	tempExplorationList[0]  = explorationCenter;        // neighbors in 0th generation
	countExploredCells = 1;                   // number of neighbors in 0th generation
	indexIgnoredPoints = 0;
	countIgnoredPoints = 1;
	//printf("Add: %i to Ignored\n", tempExplorationList[0]->index);


	// 1st generation (Set direct neighbors)
	for(j = 0; j < explorationCenter->countNeighborCells; j++){
		tempExplorationList[countExploredCells++]  = explorationCenter->neighborCells[j];        // neighbors in 0th generation
		if(explorationCenter->countExtendedNeighborCells>0){
			//printf("Add: %i to Ignored\n", tempExplorationList[countExploredCells-1]->index);
			countIgnoredPoints++;
		}else{
			//printf("Add: %i to Source\n", tempExplorationList[countExploredCells-1]->index);
			countSourcePoints++;
		}
	}

	// following generations (Set extended neighbors)
	for(j = 0; j < explorationCenter->countExtendedNeighborCells; j++){
		tempExplorationList[countExploredCells++]  = explorationCenter->extendedNeighborhood[j];        // neighbors in 0th generation
		//printf("Add: %i to Source\n", tempExplorationList[countExploredCells-1]->index);
		countSourcePoints++;
	}
	indexSourcePoints = countIgnoredPoints;
	
	// get all point below a certain radius
	countExploredPoints = 0;
	indexExploredPoints = countIgnoredPoints + countSourcePoints;	

	do{
		//printf("Ignored: count %i, index %i; Source: count %i, index %i;\n", countIgnoredPoints, indexIgnoredPoints, countSourcePoints, indexSourcePoints);
		//fprintf( stderr, ". \b");
		
		// explore neighborhood
		Triangulation::searchForVoronoiCell( countIgnoredPoints,  &tempExplorationList[indexIgnoredPoints],
		                            countSourcePoints,   &tempExplorationList[indexSourcePoints], 
		                            countExploredPoints, &tempExplorationList[indexExploredPoints]);
		
		// search closest free point
		Vertex* closestPoint = NULL;
		double distanceToClosestPoint = 0.;
		double tempDistance;
		for( j=0; j<countExploredPoints; j++){
			if( tempExplorationList[indexExploredPoints+j]->getState()==FREE){
				tempDistance = explorationCenter->getDistanceTo( tempExplorationList[indexExploredPoints+j]);
				if( tempDistance < distanceToClosestPoint || closestPoint==NULL){
					distanceToClosestPoint = tempDistance;
					closestPoint = tempExplorationList[indexExploredPoints+j];
				}
			}
		}
		if( closestPoint!=NULL){
			//fprintf( stderr, "clostest cell to cell (%i) is %i\n", explorationCenter->index, closestPoint->index);
			return closestPoint;
		}
		
		// initialize next
		indexIgnoredPoints+= countIgnoredPoints;
		countIgnoredPoints = countSourcePoints;
		
		indexSourcePoints+= countSourcePoints;
		countSourcePoints = countExploredPoints;
		
		indexExploredPoints+= countExploredPoints;
		countExploredPoints = 0;
	}while( countSourcePoints > 0);

	
	//free( tempExplorationList);
	*/
	fprintf( stderr, "WARNING: No free cell was found which is close to cell %i\n", explorationCenter->index);
	return NULL;
}
/*****************************************************************************/



Vertex *Triangulation::searchClosestVoronoiCell( Vertex *explorationCenter, double targetPosition[DIMENSIONS]){

	double min_dist = explorationCenter->getDistanceTo( targetPosition);
	Vertex *candidate = explorationCenter;

	do{
		explorationCenter = candidate;
		
		for( int i=0; i<explorationCenter->countNeighborCells; i++){
			double dist = explorationCenter->neighborCells[i]->getDistanceTo( targetPosition);
			if( dist < min_dist){
				min_dist = dist;
				candidate = explorationCenter->neighborCells[i];
			}
		}	
	}while( candidate!=explorationCenter);
	
	return candidate;
	
	//free( tempExplorationList);
	fprintf( stderr, "WARNING: No free cell was found which is close to cell %i\n", explorationCenter->index);
	return NULL;
}
/*****************************************************************************/



double Vertex::getDistanceTo( Vertex* cell)
{
	double distance = 0.;
	
	for( int i=0; i<DIMENSIONS; i++)
		distance += pow( this->position[i] - cell->position[i], 2);
	
	return sqrt( distance);
}
/*****************************************************************************/



double Vertex::getDistanceTo( double targetPosition[DIMENSIONS])
{
	double distance = 0.;
	
	for( int i=0; i<DIMENSIONS; i++)
		distance += pow( this->position[i] - targetPosition[i], 2);
	
	return sqrt( distance);
}
/*****************************************************************************/



void Triangulation::NEWsetExtendedNeighborhood( int radius){

	Triangulation *voronoiDiagram = this;
	int i, j, count;
	Vertex **tempNeighborhoodList;
	int generations;
	int countNeighborsGen[30];
	//int TEST_countExtendedNeighbors = 0;

	//printf("radius=%lf\n", radius);

	Vertex** tempExtendedNeighborhoodList = ( Vertex**) calloc ( voronoiDiagram->countVoronoiCells, sizeof( Vertex*));
	assert(tempExtendedNeighborhoodList);

	for( i = 0; i < voronoiDiagram->countVoronoiCells; i++){

		fprintf( stderr, "\rCalculate Extended Neighborhood %.3lf%%\t\t\b", 100.*(i+1.)/voronoiDiagram->countVoronoiCells);
		//fprintf( stderr, "ACTUALIZE EXTENDED NEIGHBORHOOD OF POINT %i, actual count=%i\n", i, this->vertices[i]->countExtendedNeighborCells);
		//if( this->vertices[i]->countExtendedNeighborCells == 0){
			/* Memoryallocation for list of neighbors*/
		tempNeighborhoodList = voronoiDiagram->vertices[i]->extendedNeighborhood;
		voronoiDiagram->vertices[i]->extendedNeighborhood = tempExtendedNeighborhoodList;


		// 0st generation (Set point itself)
		voronoiDiagram->vertices[i]->extendedNeighborhood[0]  = voronoiDiagram->vertices[i];        // neighbors in 0th generation
		countNeighborsGen[0] = 1;                   // number of neighbors in 0th generation
		//fprintf( stderr, "0st generation:{ %i}\n", voronoiDiagram->vertices[i]->extendedNeighborhood[0]->index);


		// 1st generation (Set direct neighbors)
		count = 0;
		//fprintf( stderr, "1st generation:{");
		for(j = 1; j <= voronoiDiagram->vertices[i]->countNeighborCells; j++){
			//printf("neighbor:%i, ", voronoiDiagram->vertices[i].neighbors[j-1]->nr);
			if(1 <= radius){  
				voronoiDiagram->vertices[i]->extendedNeighborhood[j] = voronoiDiagram->vertices[i]->neighborCells[j-1]; // neighbors in 1st generation
				//fprintf( stderr, " %i,", voronoiDiagram->vertices[i]->extendedNeighborhood[j]->index);
				count++;
			}
		}
		/*if( this->vertices[i]->countExtendedNeighborCells != 0){
			// set already added extended neighbors
			int k=j;
			for(j=0; j < this->vertices[i]->countExtendedNeighborCells; j++){
				//printf("neighbor:%i, ", voronoiDiagram->vertices[i].neighbors[j-1]->nr);
				//if(1 <= radius){  
					voronoiDiagram->vertices[i]->extendedNeighborhood[k+j] = tempNeighborhoodList[j]; // neighbors in 1st generation
					//fprintf( stderr, "1st generation=%i [ADDITION]\n", voronoiDiagram->vertices[i]->extendedNeighborhood[k+j]->index);
					//fprintf( stderr, " %i*,", voronoiDiagram->vertices[i]->extendedNeighborhood[k+j]->index);
					count++;
				//}
			}
			free( tempNeighborhoodList);
		}*/
		//fprintf( stderr, "\b}\n");
		countNeighborsGen[1] = count;       // number of neighbors in 1st generation
		//printf("\ncountNeighborsGen[0]:%i, countNeighborsGen[1]:%i\n", countNeighborsGen[0], countNeighborsGen[1]);

		/* Initialize depending parameters*/
		generations                 = 2;  // number of set generations
		voronoiDiagram->vertices[i]->countExtendedNeighborCells   = 1+count;  // point itself and direct
		//GetAgent( voronoiDiagram->vertices[i])->countReachableFreeNeighbors = 0;  // point itself and direct neighbors aren't important
		//}
		/*else{
			// backup old list
			tempNeighborhoodList = voronoiDiagram->vertices[i]->extendedNeighborhood;
			voronoiDiagram->vertices[i]->extendedNeighborhood = tempExtendedNeighborhoodList;
			
			// 1st generation (Set direct neighbors)
			count = 0;
			for(j = 1; j <= voronoiDiagram->vertices[i]->countNeighborCells; j++){
				//printf("neighbor:%i, ", voronoiDiagram->vertices[i].neighbors[j-1]->nr);
				if(1 <= radius){  
					voronoiDiagram->vertices[i]->extendedNeighborhood[j] = voronoiDiagram->vertices[i]->neighborCells[j-1]; // neighbors in 1st generation
					count++;
				}
			}
			// set already added extended neighbors
			int k=j;
			for(j=0; j < this->vertices[i]->countExtendedNeighborCells; j++){
				//printf("neighbor:%i, ", voronoiDiagram->vertices[i].neighbors[j-1]->nr);
				//if(1 <= radius){  
					voronoiDiagram->vertices[i]->extendedNeighborhood[k+j] = tempNeighborhoodList[j]; // neighbors in 1st generation
					count++;
				//}
			}
			countNeighborsGen[1] = count;       // number of neighbors in 1st generation

			generations = 2;
			voronoiDiagram->vertices[i]->countExtendedNeighborCells   = 1+count;  // point itself and direct
			
			free( tempNeighborhoodList);

		}*/

		/* get all point below a certain radius */
		voronoiDiagram->vertices[i]->actualizeExtendedNeighborhood( this, 1, generations, countNeighborsGen, radius); // beginning with the 2nd generation of neighbors



		voronoiDiagram->vertices[i]->countExtendedNeighborCells   -= 1+voronoiDiagram->vertices[i]->countNeighborCells;  // point itself and direct neighbors aren't important
		voronoiDiagram->vertices[i]->countFreeExtendedNeighborCells = voronoiDiagram->vertices[i]->countExtendedNeighborCells;  // point itself and direct neighbors aren't important

		/* Delete point itself and direct neighbors from list */
		tempNeighborhoodList = ( Vertex**) calloc ( voronoiDiagram->vertices[i]->countExtendedNeighborCells + 1, sizeof( Vertex*)); // memoryallocation
		assert(tempNeighborhoodList);
		for(j = 0; j < voronoiDiagram->vertices[i]->countExtendedNeighborCells; j++)
			tempNeighborhoodList[j] = voronoiDiagram->vertices[i]->extendedNeighborhood[1 + voronoiDiagram->vertices[i]->countNeighborCells + j ];                //

		//free( voronoiDiagram->vertices[i]->extendedNeighborhood);		// free unused allocated memory
		voronoiDiagram->vertices[i]->extendedNeighborhood = tempNeighborhoodList; // set pointer
		voronoiDiagram->vertices[i]->extendedNeighborCellsInitialized = TRUE;


		//printf ("%li \n", sizeof(voronoiDiagram->vertices[i]->extendedNeighborhood));

		//test_count += GetAgent( voronoiDiagram->vertices[i])->countReachableNeighbors;
		//TEST_countExtendedNeighbors += voronoiDiagram->vertices[i]->countExtendedNeighborCells;
	}
	free( tempExtendedNeighborhoodList);

  //printf("###### k_r: %lf \n", k_r);
  //printf("###### Xtended neighbors: total %i, average %lf \n", TEST_countExtendedNeighbors, TEST_countExtendedNeighbors/(double)voronoiDiagram->countVoronoiCells);
  //exit(0);
}
/*****************************************************************************/



void getShiftingNeighborhood( Triangulation *voronoiDiagram, Vertex *point, int base_index, int generations, int* nr_kr_neighbors_gen, double radius){
	int i, j;
	int count = 0;
	int base_gen = generations - 1; // number of last gotten generation

	//printf("get_grid_points_below_kr( Grid *point(%i), base_index(%i), generations(%i))\n", point->nr, base_index, generations);

	for( i = 0; i < nr_kr_neighbors_gen[ base_gen ]; i++){ // for every element in the previous list
		//printf("  point:%i\n", voronoiDiagram[ point->reachableNeighbors[ base_index + i ] ].nr);
		for( j = 0; j < point->extendedNeighborhood[ base_index + i ]->countNeighborCells; j++){ // for every neighbor of a element in the previous list
        //printf("  candidates:%i\n", (voronoiDiagram[ point->reachableNeighbors[ base_index + i ] ].neighbors[j])->nr);
        if( generations-1 < (int) radius
/*            get_surface_distance(voronoiDiagram[ point->reachableNeighbors[ base_index + i ] ].neighbors[j], point) <= Divide_cell_depth_radius  // distance between given point and choosen neighbor
*/            && 
            /*!isElementOf( &(point->reachableNeighbors[ base_index - nr_kr_neighbors_gen[ base_gen - 1 ]]), // pointer to list of last, actual and new generation
                          nr_kr_neighbors_gen[ base_gen - 1 ] + nr_kr_neighbors_gen[ base_gen ] + count, // number of points in the descendant list
                          voronoiDiagram[ point->reachableNeighbors[ base_index + i ] ].neighbors[j])){                    // actual neighbor of neighbor 
            *//*!isElementOf( point->reachableNeighbors,                // pointer to list of last, actual and new generation
                         point->countReachableNeighbors + count, // number of points in the descendant list
                         GetAgent( GetVoronoiCell(point->reachableNeighbors[ base_index + i ])->neighborCells[j]))){                   // actual neighbor of neighbor 
            */!isElementOf( &(point->extendedNeighborhood[ base_index - nr_kr_neighbors_gen[ base_gen - 1 ]]),                // pointer to list of last, actual and new generation
                         nr_kr_neighbors_gen[ base_gen - 1 ] + nr_kr_neighbors_gen[ base_gen ] + count, // number of points in the descendant list
                         point->extendedNeighborhood[ base_index + i ]->neighborCells[j])){                   // actual neighbor of neighbor 
	  point->extendedNeighborhood[base_index + nr_kr_neighbors_gen[ base_gen ] + count] = point->extendedNeighborhood[ base_index + i ]->neighborCells[j];

          count++;
        }
     }
  }

  if( count != 0){

    /* Update depending parameters */
    generations++;                             // number of set generations
    point->countExtendedNeighborCells           += count; // actualize number of neighbors
    //point->countReachableFreeNeighbors       += count; // actualize number of free neighbors
    nr_kr_neighbors_gen[base_gen + 1] = count; // set number of neighbors

    /* get next generation */
    getShiftingNeighborhood( voronoiDiagram, point, base_index + nr_kr_neighbors_gen[ base_gen ], generations, nr_kr_neighbors_gen, radius);
  } 
  
  return;
}
/*****************************************************************************/



void Triangulation::setExtendedNeighborhood( double radius){

	Triangulation *voronoiDiagram = this;
	int i, j, count;
	Vertex **tempNeighborhoodList;
	int generations;
	int countNeighborsGen[30];

	//printf("radius=%lf\n", radius);

	for( i = 0; i < voronoiDiagram->countVoronoiCells; i++){
		/* Memoryallocation for list of neighbors*/
		voronoiDiagram->vertices[i]->extendedNeighborhood  = ( Vertex**) calloc ( voronoiDiagram->countVoronoiCells, sizeof( Vertex*));
		assert(voronoiDiagram->vertices[i]->extendedNeighborhood);

		// 0st generation (Set point itself)
		voronoiDiagram->vertices[i]->extendedNeighborhood[0]  = voronoiDiagram->vertices[i];        // neighbors in 0th generation
		countNeighborsGen[0] = 1;                   // number of neighbors in 0th generation


		// 1st generation (Set direct neighbors)
		count = 0;
		for(j = 1; j <= voronoiDiagram->vertices[i]->countNeighborCells; j++){
			//printf("neighbor:%i, ", voronoiDiagram->vertices[i].neighbors[j-1]->nr);
			if(1 <= radius){  
				voronoiDiagram->vertices[i]->extendedNeighborhood[j] = voronoiDiagram->vertices[i]->neighborCells[j-1]; // neighbors in 1st generation
				count++;
			}
		}
		countNeighborsGen[1] = count;       // number of neighbors in 1st generation
		//printf("\ncountNeighborsGen[0]:%i, countNeighborsGen[1]:%i\n", countNeighborsGen[0], countNeighborsGen[1]);

		/* Initialize depending parameters*/
		generations                 = 2;  // number of set generations
		voronoiDiagram->vertices[i]->countExtendedNeighborCells   = 1+count;  // point itself and direct
		//GetAgent( voronoiDiagram->vertices[i])->countReachableFreeNeighbors = 0;  // point itself and direct neighbors aren't important

		/* get all point below a certain radius */
		getShiftingNeighborhood( this, voronoiDiagram->vertices[i], 1, generations, countNeighborsGen, radius); // beginning with the 2nd generation of neighbors

		voronoiDiagram->vertices[i]->countExtendedNeighborCells   -= 1+count;  // point itself and direct neighbors aren't important
		voronoiDiagram->vertices[i]->countFreeExtendedNeighborCells = voronoiDiagram->vertices[i]->countExtendedNeighborCells;  // point itself and direct neighbors aren't important

		/* Delete point itself and direct neighbors from list */
		tempNeighborhoodList = ( Vertex**) calloc ( voronoiDiagram->vertices[i]->countExtendedNeighborCells + 1, sizeof( Vertex*)); // memoryallocation
		assert(tempNeighborhoodList);
		for(j = 0; j < voronoiDiagram->vertices[i]->countExtendedNeighborCells; j++)
			tempNeighborhoodList[j] = voronoiDiagram->vertices[i]->extendedNeighborhood[1 + count + j ];                //

		free( voronoiDiagram->vertices[i]->extendedNeighborhood);		// free unused allocated memory
		voronoiDiagram->vertices[i]->extendedNeighborhood = tempNeighborhoodList; // set pointer
		voronoiDiagram->vertices[i]->extendedNeighborCellsInitialized = TRUE;

		//test_count += GetAgent( voronoiDiagram->vertices[i])->countReachableNeighbors;
	}

  //printf("###### k_r: %lf \n", k_r);
  //printf("###### average k_r neighbors: %lf \n", test_count/(double)voronoiDiagram->countVoronoiCells);
  //exit(0);
}
/*****************************************************************************/


void Triangulation::setDomain()
{
	if( !domainSet){
		if( countVoronoiCells == 0){
			for( int dim = 0; dim<DIMENSIONS; dim++){
				xMin[dim] = 0;
				xMax[dim] = 0;
			}
			return;
		}

		int i = 0, dim;
	
		// domain borders
		for( dim = 0; dim<DIMENSIONS; dim++){
			xMin[dim] = vertices[i]->position[dim];
			xMax[dim] = vertices[i]->position[dim];
		}

		for( i = 1; i < countVoronoiCells; i++){
			for( dim = 0; dim<DIMENSIONS; dim++){
				if( xMin[dim] > vertices[i]->position[dim])
					xMin[dim] = vertices[i]->position[dim];
				if( xMax[dim] < vertices[i]->position[dim])
					xMax[dim] = vertices[i]->position[dim];
			}
		}
	
		// boundary thickness
		boundaryThickness = 1.;
		for( dim = 0; dim<DIMENSIONS; dim++)
			boundaryThickness *= (xMax[dim]-xMin[dim]); 
		boundaryThickness = pow( boundaryThickness / (double)countVoronoiCells ,1./DIMENSIONS);
					
		// elements in each direction
		for( dim = 0; dim<DIMENSIONS; dim++)
			//xN[dim] = (int)(xMax[dim]-xMin[dim] + 1.);
			xN[dim] = (int)(ceil(xMax[dim])-ceil(xMin[dim])+1);

		domainSet = TRUE;
	}
}
/*****************************************************************************/
  


Vertex * Triangulation::getCentralCell(){
	double min[DIMENSIONS], max[DIMENSIONS], mid[DIMENSIONS];
	
	// minima, maxima and average of x and y
	for( int ii=0; ii<DIMENSIONS; ii++)
		min[ii] = max[ii] = this->vertices[0]->position[ii];
  
 	for( int i = 1; i < this->countVoronoiCells; i++){
		for( int ii=0; ii<DIMENSIONS; ii++){
			if( max[ii] < this->vertices[i]->position[ii])
				max[ii] = this->vertices[i]->position[ii];
			if( min[ii] > this->vertices[i]->position[ii] )
				min[ii] = this->vertices[i]->position[ii];
		}
	}
  
	for( int ii=0; ii<DIMENSIONS; ii++)
		mid[ii] = (min[ii] + max[ii]) / 2;
  
	// central cell
	Vertex* central_cell = this->vertices[0];
	double distance = 0.;
	for( int ii=0; ii<DIMENSIONS; ii++)
		distance += pow( mid[ii] - this->vertices[0]->position[ii], 2.);
	distance = sqrt( distance);
	
	for( int i = 1; i < this->countVoronoiCells; i++){
		double temp_distance = 0.;
		for( int ii=0; ii<DIMENSIONS; ii++)
			temp_distance += pow( mid[ii] - this->vertices[i]->position[ii], 2.);
		temp_distance = sqrt( temp_distance);
		if( distance > temp_distance){
			central_cell = this->vertices[i];
			distance = temp_distance;
		}
	}
  
	//fprintf( stderr, "central cell: %i, %lf\n", central_cell->index, central_cell->position[0]);
  
	return central_cell;
}
/*****************************************************************************/


void initCircumsphereRadius( Tetrahedron* tet)
{
	int k, l;
	
	//double A[NR_TETRAHEDRON_POINTS][NR_TETRAHEDRON_POINTS]/*={{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.}}*/;
	//double **A = newDoubleMatrix( NR_TETRAHEDRON_POINTS, NR_TETRAHEDRON_POINTS);
	//double x[NR_TETRAHEDRON_POINTS];

	// initialize matrix and vector of linear system
	for( k=0; k<NR_TETRAHEDRON_POINTS; k++){
		tet->circumsphere[k]=0.;
		for( l=0; l<DIMENSIONS; l++){
			tet->circumsphere[k] -= pow( tet->vertices[k]->position[l], 2);
			//tet->circumsphere[k] -= myPow( tet->vertices[k]->position[l], 2);
			_A[k][l] = tet->vertices[k]->position[l];
			//A[k][l] = tet->vertices[k]->position[l];
		}
		_A[k][l] = 1.;
		//A[k][l] = 1.;
	}
	solveLinearSystemB( (double**)_A, tet->circumsphere, (double*)_x, NR_TETRAHEDRON_POINTS, (double**)_B);
	//solveLinearSystem( (double**)A, tet->circumsphere, x, NR_TETRAHEDRON_POINTS);
	tet->radius = 0.;
	for( l=0; l<DIMENSIONS; l++)
		tet->radius += tet->circumsphere[l]*tet->circumsphere[l];
	tet->radius = tet->radius/4. - tet->circumsphere[l];
	tet->circumsphereInitialized = 1;

	//deleteDoubleMatrix( A, NR_TETRAHEDRON_POINTS);

	// NUMERIC ERROR

	/*double R  = 0.;
	for( l=0; l<DIMENSIONS; l++)
		//R += pow(tet->vertices[0]->position[l] + tet->circumsphere[l]/2., 2);
		R += myPow(tet->vertices[0]->position[l] + tet->circumsphere[l]/2., 2);
	//double maxNumericError = R - tet->radius;
	double avNumericError = R - tet->radius;
	for( k=1; k<NR_TETRAHEDRON_POINTS; k++){
		R  = 0.;
		for( l=0; l<DIMENSIONS; l++)
			//R += pow(tet->vertices[k]->position[l] + tet->circumsphere[l]/2., 2);
			R += myPow(tet->vertices[k]->position[l] + tet->circumsphere[l]/2., 2);
		tet->numericError = R - tet->radius;
		fprintf( stderr, "INFO: numerical error of tet %i for vertice %i: %e\n", tet->index, tet->vertices[k]->index, tet->numericError );
		//if( maxNumericError < tet->numericError)
		//	maxNumericError = tet->numericError;
		avNumericError += R - tet->radius;
	}
	//tet->numericError = maxNumericError;
	tet->numericError = avNumericError/NR_TETRAHEDRON_POINTS;
*/
	// END NUMERIC ERROR
}
/*****************************************************************************/


void getCircumsphereRadius( Tetrahedron* tet, long double *circumsphere, long double &radius)
{
	int k, l;

	//long double A[NR_TETRAHEDRON_POINTS][NR_TETRAHEDRON_POINTS]/*={{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.}}*/;
	//long double B[NR_TETRAHEDRON_POINTS][NR_TETRAHEDRON_POINTS]/*={{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.}}*/;
	long double **A = newLongDoubleMatrix(NR_TETRAHEDRON_POINTS,NR_TETRAHEDRON_POINTS);
	long double **B = newLongDoubleMatrix(NR_TETRAHEDRON_POINTS,NR_TETRAHEDRON_POINTS);
	long double x[NR_TETRAHEDRON_POINTS];

	// initialize matrix and vector of linear system
	for( k=0; k<NR_TETRAHEDRON_POINTS; k++){
		circumsphere[k]=0.;
		for( l=0; l<DIMENSIONS; l++){
			circumsphere[k] -= pow( tet->vertices[k]->position[l], 2);
			A[k][l] = tet->vertices[k]->position[l];
		}
		A[k][l] = 1.;
	}
	solveLinearSystemB( (long double**)A, circumsphere, (long double*)x, NR_TETRAHEDRON_POINTS, (long double**)B);

	for( l=0; l<DIMENSIONS; l++)
		radius += circumsphere[l]*circumsphere[l];
	radius = radius/4. - circumsphere[l];

	deleteLongDoubleMatrix( A, NR_TETRAHEDRON_POINTS);
	deleteLongDoubleMatrix( B, NR_TETRAHEDRON_POINTS);
}
/*****************************************************************************/


double getCircumsphereRadius( Tetrahedron* tet)
{
	
	// initialize circumsphere
	if( tet->circumsphereInitialized == 0)
		initCircumsphereRadius( tet);

	return tet->radius;
}
/*****************************************************************************/


double getDistanceOfPointToCircumsphereCenter( Vertex* point, Tetrahedron* tet)
{
	int l;

	// initialize circumsphere
	if( tet->circumsphereInitialized == 0)
		initCircumsphereRadius( tet);

	double R  = 0.;
	for( l=0; l<DIMENSIONS; l++)
		R += pow(point->position[l] + tet->circumsphere[l]/2., 2);
		//R += myPow(point->position[l] + tet->circumsphere[l]/2., 2);

	return R;
}
/*****************************************************************************/


double getDistanceOfPointToCircumsphere( Vertex* point, Tetrahedron* tet)
{
	int l;

	// initialize circumsphere
	if( tet->circumsphereInitialized == 0)
		initCircumsphereRadius( tet);

	double R  = 0.;
	for( l=0; l<DIMENSIONS; l++)
		//R += pow((float)point->position[l] + (float)tet->circumsphere[l]/2., 2);
		R += myPow(point->position[l] + tet->circumsphere[l]/2., 2);

	if( CS_ISNAN( R - tet->radius )){
		fprintf( stderr, "dist = %e\n", (R - tet->radius));
		fprintf( stderr, "R = %e\n", R);
		fprintf( stderr, "tet->radius = %e\n", tet->radius);
		for( l=0; l<DIMENSIONS; l++){
			//R += myPow(point->position[l] + tet->circumsphere[l]/2., 2);
			fprintf( stderr, "point->position[%i] = %e\n", l, point->position[l]);
			fprintf( stderr, "tet->circumsphere[%i] = %e\n", l, tet->circumsphere[l]);
		}

		exit( 0);
	}

	//if( fabs(R - tet->radius)<1e-7)
	if( fabs(R - tet->radius)<1e-9)
	{

		long double circumsphere[DIMENSIONS+1];
		long double radius = 0.;
		getCircumsphereRadius( tet, circumsphere, radius);
		long double ldR  = 0.;
		for( l=0; l<DIMENSIONS; l++)
			ldR += ((long double)point->position[l] + circumsphere[l]/2.)*((long double)point->position[l] + circumsphere[l]/2.);
		//fprintf( stderr, "WARNING: dist = %.10e ", (R - tet->radius));
		//fprintf( stderr, "----> ld dist = %.10e ", (double)(ldR - radius));
		//fprintf( stderr, "------> error = %.10e\n", fabs((double)(ldR - radius)-(R - tet->radius)));
		//if( signof(ldR - radius) != signof(R - tet->radius))
		if( (ldR - radius)*(R - tet->radius) < 0.)
		{
			fprintf( stderr, "ERROR: wrong sign!\n");
			fprintf( stderr, "WARNING: dist = %.10e\n", (R - tet->radius));
			fprintf( stderr, "----> ld dist = %.10e ", (double)(ldR - radius));
			fprintf( stderr, "------> error = %.10e\n", fabs((double)(ldR - radius)-(R - tet->radius)));
			return (double)(ldR - radius);
			//exit(0);
		}

		/*long double Rl = 0.;
		for( l=0; l<DIMENSIONS; l++)
			//R += pow((float)point->position[l] + (float)tet->circumsphere[l]/2., 2);
			R += ((long double)point->position[l] + (long double)tet->circumsphere[l]/2.) * ((long double)point->position[l] + (long double)tet->circumsphere[l]/2.);
		fprintf( stderr, "-------> dist = %.10e\n", (Rl - (long double)tet->radius));*/
	}
	return (R - tet->radius)/* - tet->numericError*/;
}
/*****************************************************************************/


Tetrahedron *Triangulation::getTetrahedronContainingPointInCircumSphere( Vertex* point)
{
	int i;

	// look for first tet containing new point
	//fprintf( stderr, "look for first tet containing new point\n");
	Tetrahedron *actualTet = NULL,
	            *candidateTet = this->tetrahedra[(int)(myRandE(1.)*(double)this->countTetrahedra)];
	double distanceToCircumsphere;
	double minDistanceToCircumsphere = getDistanceOfPointToCircumsphere( point, candidateTet);

	// chose closest random point as start point
	for( i=0; i<5; i++){
		actualTet = this->tetrahedra[(int)(myRandE(1.)*(double)this->countTetrahedra)];
		distanceToCircumsphere = getDistanceOfPointToCircumsphere( point, actualTet);
		if( minDistanceToCircumsphere > distanceToCircumsphere){
			minDistanceToCircumsphere = distanceToCircumsphere;
			candidateTet = actualTet;
		}
	}
	actualTet = NULL;

	// follow shortest way to first tetrahedron containing new point
	do{
	//	fprintf( stderr, "INFO: actualTet!=candidateTet ti=%i\n", ti);
		actualTet=candidateTet;
		//fprintf( stderr, "INFO: actualTet(%i)! ti=%i\n", actualTet->index, ti);
		//fprintf( stderr, "INFO: actualTet(%i)==candidateTet(%i)! ti=%i\n", actualTet->index, candidateTet->index, ti);
		for( i=0; i<actualTet->countNeighborTetrahedra; i++){
			// circumsphere contains point i
			distanceToCircumsphere = getDistanceOfPointToCircumsphere( point, actualTet->neighborTetrahedra[i]);
			if( minDistanceToCircumsphere > distanceToCircumsphere){
				minDistanceToCircumsphere = distanceToCircumsphere;
				candidateTet = actualTet->neighborTetrahedra[i];
			}
		}
	}while( minDistanceToCircumsphere > 0.);	

	return candidateTet;
}
/*****************************************************************************/


double getDistanceOfPointToTetrahedron( Vertex* point, Tetrahedron* tet)
{
	int i=0, ii;

	double distance;
	double minDistance = 0;
	if( tet==NULL)
		fprintf(stderr, "Tet == NULL\n");
	if( point==NULL)
		fprintf(stderr, "Point == NULL\n");
	//fprintf(stderr, "Point %i (%.3lf, %.3lf, %.3lf)\n", point->index, point->position[0], point->position[1], point->position[2]);
	for( ii=0; ii<DIMENSIONS; ii++)
		minDistance += pow( tet->vertices[i]->position[ii] - point->position[ii], 2);
		//minDistance += myPow( tet->vertices[i]->position[ii] - point->position[ii], 2);

	for( i=1; i<NR_TETRAHEDRON_POINTS; i++){
		distance = 0;
		for( ii=0; ii<DIMENSIONS; ii++)
			distance += pow( tet->vertices[i]->position[ii] - point->position[ii], 2);
			//distance += myPow( tet->vertices[i]->position[ii] - point->position[ii], 2);

		if( minDistance > distance)
			minDistance = distance;
	}

	return minDistance;
}
/*****************************************************************************/


Tetrahedron *getTetrahedronContainingPoint( Triangulation *voronoiDiagram, Vertex* point)
{
	int i;

	// look for first tet containing new point
	//fprintf( stderr, "look for first tet containing new point\n");
	i=(int)(myRandE((double)voronoiDiagram->countTetrahedra));
	if(i>=voronoiDiagram->countTetrahedra){
		fprintf( stderr, "WRONG1! ;)\n");
		exit(0);
	}
	Tetrahedron *actualTet = NULL,
            *candidateTet = voronoiDiagram->tetrahedra[i];
			//*candidateTet = voronoiDiagram->tetrahedra[(int)(myRand()*(double)voronoiDiagram->countTetrahedra)];

	double distance;
	double candidateDistance = getDistanceOfPointToTetrahedron( point, candidateTet);

	// chose closest random point as start point
	for( i=0; i<5; i++){
		int temp = (int)(myRandE((double)voronoiDiagram->countTetrahedra));
		if(temp>=voronoiDiagram->countTetrahedra){
			fprintf( stderr, "WRONG2! ;)\n");
			exit(0);
		}
		//actualTet = voronoiDiagram->tetrahedra[(int)(myRand()*(double)voronoiDiagram->countTetrahedra)];
		actualTet = voronoiDiagram->tetrahedra[temp];
		distance = getDistanceOfPointToTetrahedron( point, actualTet);
		if( candidateDistance > distance){
			candidateDistance = distance;
			candidateTet = actualTet;
		}
	}


	// follow shortest way to first tetrahedron containing new point
	while( candidateTet->vertices[0]!=point && candidateTet->vertices[1]!=point && candidateTet->vertices[2]!=point
#if DIMENSIONS == 3
			&& candidateTet->vertices[3]!=point
#endif
			){
		//fprintf( stderr, "INFO: candidateTet=%i\n", candidateTet->index);
		actualTet=candidateTet;
		//fprintf( stderr, "INFO: actualTet(%i), (%i %i %i %i), dist=%lf (%lf) to point %i ( %lf, %lf, %lf)\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index, /*minDistance*/0, getDistanceOfPointToTetrahedron( point, actualTet), point->index, point->position[0], point->position[1], point->position[2]);


		for( i=0; i<actualTet->countNeighborTetrahedra; i++){
			distance = getDistanceOfPointToTetrahedron( point, actualTet->neighborTetrahedra[i]);
			//fprintf( stderr, "INFO: %i. neighbor tet(%i), (%i %i %i %i), dist=%lf\n", i+1, actualTet->neighborTetrahedra[i]->index, actualTet->neighborTetrahedra[i]->vertices[0]->index, actualTet->neighborTetrahedra[i]->vertices[1]->index, actualTet->neighborTetrahedra[i]->vertices[2]->index, actualTet->neighborTetrahedra[i]->vertices[3]->index, distance);
			//fprintf( stderr, "INFO: dist=%lf", distance);
			if( candidateDistance > distance){
			//if( minDistance > distance = getDistanceOfPointToTetrahedron( point, actualTet->neighborTetrahedra[i])){
				candidateTet      = actualTet->neighborTetrahedra[i];
				candidateDistance = distance;
			}
		}

		if( actualTet==candidateTet){
			int temp = (int)(myRandE((double)actualTet->countNeighborTetrahedra));
			if(temp>=actualTet->countNeighborTetrahedra){
				fprintf( stderr, "WRONG3: %i / %i ;)\n", temp, actualTet->countNeighborTetrahedra);
				exit(0);
			}
			candidateTet      = actualTet->neighborTetrahedra[temp];

			//candidateTet      = actualTet->neighborTetrahedra[(int)(myRand()*(double)actualTet->countNeighborTetrahedra)];
			candidateDistance = getDistanceOfPointToTetrahedron( point, candidateTet);
		}
		
	}	

	return candidateTet;
}
/*****************************************************************************/


int countCommonPointOfTetrahedra( Tetrahedron* tetA, Tetrahedron* tetB)
{
	int l, m;
	int nrCommonPoints = 0;
		
	for( l=0; l<DIMENSIONS+1; l++) // for all points of j
		for( m=0; m<DIMENSIONS+1; m++) // for all points of k
			if( tetA->vertices[l] == tetB->vertices[m])
				nrCommonPoints++;
				
	return nrCommonPoints;
}
/*****************************************************************************/


int removeTetrahedronNeighbor( Tetrahedron* tet, Tetrahedron* neighborTet)
{
	int i;

	//fprintf( stderr, "INFO: removeNeighbor(): remove %i as neighbor of %i!\n", neighborCell->index, voronoiCell->index);
	
	for( i=0; i<tet->countNeighborTetrahedra; i++){
		if( tet->neighborTetrahedra[i] == neighborTet){
			tet->countNeighborTetrahedra--;
			tet->neighborTetrahedra[i] = tet->neighborTetrahedra[tet->countNeighborTetrahedra];
			//tet->neighborCells = (Vertex**) realloc( tet->neighborCells, tet->countNeighborCells*sizeof(Vertex*));
			return 1;
		}	
	}

	fprintf( stderr, "ERROR: removeTetrahedronNeighbor() can not find the neighbor to delete!\n");
	exit( 0);
	return 0;
}
/*****************************************************************************/


int Tetrahedron::addTetrahedronNeighbor( Tetrahedron* neighborTet)
{
	this->neighborTetrahedra[this->countNeighborTetrahedra] = neighborTet;
	this->countNeighborTetrahedra++;	

	return 1;
}
/*****************************************************************************/


int addOrReplaceTetrahedronNeighbor( Tetrahedron* tet, Tetrahedron* neighborTet)
{
	int i;
	for( i=0; i<tet->countNeighborTetrahedra; i++){
		if( tet->neighborTetrahedra[i] == neighborTet){
			//fprintf( stderr, "ERROR: addTetrahedronNeighbor() found neighbor being already included!\n");
			//exit( 0);
			return 0;
		}	
	}
	

	//printf("realloc: %i -> %i\n", voronoiCell->countNeighborCells, voronoiCell->countNeighborCells+1);
	//if( voronoiCell->countNeighborCells != 0)
		//tet->neighborCells = (Vertex**) realloc( voronoiCell->neighborCells, (voronoiCell->countNeighborCells+1)*sizeof(Vertex*));
	//else
	//	voronoiCell->neighborCells = (Vertex**) malloc( sizeof());
	//printf("set new neighbor at the end of neighbor list\n");
	tet->neighborTetrahedra[tet->countNeighborTetrahedra] = neighborTet;
	//printf("increase neighbor counter\n");
	tet->countNeighborTetrahedra++;	
	//printf("done!\n");

	return 1;
}
/*****************************************************************************/


int tetrahedronContainsFace( Vertex** tet, Vertex** face)
{
	int i, j = 0;
	
	for( i=0; i<DIMENSIONS && j!=DIMENSIONS+1; i++)
		for( j=0; j<DIMENSIONS+1 && face[i]!=tet[j]; j++);	

	if( i==DIMENSIONS && j!=DIMENSIONS+1){
		//fprintf( stderr, "( %i, %i, %i, %i) contains ( %i, %i, %i)\n", tet[0]->index, tet[1]->index, tet[2]->index, tet[3]->index, face[0]->index, face[1]->index, face[2]->index);
		return 1;
	}
	else
		return 0;
}
/*****************************************************************************/


void setTetrahedronNeighbors( Triangulation* voronoiDiagram)
{
	int i, j;

	for( i=0; i<voronoiDiagram->countTetrahedra-1; i++){
		// find neighbors of tetrahedron i
		
		for( j=i+1; j<voronoiDiagram->countTetrahedra; j++){
			// check if tetrahedron j is a neighbor of tetrahedron i

			if( countCommonPointOfTetrahedra( voronoiDiagram->tetrahedra[i], voronoiDiagram->tetrahedra[j]) == DIMENSIONS){
				voronoiDiagram->tetrahedra[i]->addTetrahedronNeighbor( voronoiDiagram->tetrahedra[j]);
				voronoiDiagram->tetrahedra[j]->addTetrahedronNeighbor( voronoiDiagram->tetrahedra[i]);
			}
		}
	}

	// TEST: check if all neighborhood relations are reversible
	/*int k;
	for( i=0; i<voronoiDiagram->countTetrahedra; i++){
		for( j=0; j<voronoiDiagram->tetrahedra[i]->countNeighborTetrahedra; j++){
			for( k=0; k<voronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->countNeighborTetrahedra &&
				voronoiDiagram->tetrahedra[i] != voronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->neighborTetrahedra[k]; k++);
			
			if( k==voronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->countNeighborTetrahedra){
				fprintf( stderr, "ERROR: neighborship relation between tets %p -> %p!!\n", voronoiDiagram->tetrahedra[i], voronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->neighborTetrahedra[k]);
				exit( 0);
			}
		}
	}*/
}
/*****************************************************************************/


void printTetrahedronNeighbors( Tetrahedron* tetrahedron)
{
	int j;
	
	/*for( i=0; i<voronoiDiagram->countTetrahedra; i++){
		fprintf( stderr, "Tetrahedron %p has %i neighbors: ", voronoiDiagram->tetrahedra[i], voronoiDiagram->tetrahedra[i]->countNeighborTetrahedra);
		for( j=0; j<voronoiDiagram->tetrahedra[i]->countNeighborTetrahedra; j++){
			fprintf( stderr, "%p ", voronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]);
		}
		fprintf( stderr, "\n");
	}*/
	fprintf( stderr, "Tetrahedron %p has %i neighbors: ", tetrahedron, tetrahedron->countNeighborTetrahedra);
	for( j=0; j<tetrahedron->countNeighborTetrahedra; j++){
		fprintf( stderr, "%p ", tetrahedron->neighborTetrahedra[j]);
	}
	fprintf( stderr, "\n");		
}
/*****************************************************************************/


void printTetrahedraNeighbors( Triangulation* voronoiDiagram)
{
	int i;
	
	for( i=0; i<voronoiDiagram->countTetrahedra; i++){
		printTetrahedronNeighbors( voronoiDiagram->tetrahedra[i]);
	}		
}
/*****************************************************************************/


int removeNeighbor( Vertex* voronoiCell, Vertex* neighborCell)
{
	int i;

	//fprintf( stderr, "INFO: removeNeighbor(): remove %i as neighbor of %i!\n", neighborCell->index, voronoiCell->index);
	
	for( i=0; i<voronoiCell->countNeighborCells; i++){
		if( voronoiCell->neighborCells[i] == neighborCell){
			voronoiCell->countNeighborCells--;
			voronoiCell->neighborCells[i] = voronoiCell->neighborCells[voronoiCell->countNeighborCells];
			voronoiCell->neighborCells = (Vertex**) realloc( voronoiCell->neighborCells, voronoiCell->countNeighborCells*sizeof(Vertex*));
			assert(voronoiCell->countNeighborCells==0 || voronoiCell->neighborCells);
			return 1;
		}	
	}

	//fprintf( stderr, "ERROR: removeNeighbor() can not find the neighbor to delete!\n");
	//exit( 0);
	return 0;
}
/*****************************************************************************/


int addNeighbor( Vertex* voronoiCell, Vertex* neighborCell)
{
	int i;
	for( i=0; i<voronoiCell->countNeighborCells; i++){
		if( voronoiCell->neighborCells[i] == neighborCell){
	//		fprintf( stderr, "ERROR: addNeighbor() found neighbor being already included!\n");
			//exit( 0);
			return 0;
		}	
	}
	

	//printf("realloc: %i -> %i\n", voronoiCell->countNeighborCells, voronoiCell->countNeighborCells+1);
	//if( voronoiCell->countNeighborCells != 0)
		voronoiCell->neighborCells = (Vertex**) realloc( voronoiCell->neighborCells, (voronoiCell->countNeighborCells+1)*sizeof(Vertex*));
		assert(voronoiCell->neighborCells);
	//else
	//	voronoiCell->neighborCells = (Vertex**) malloc( sizeof());
	//printf("set new neighbor at the end of neighbor list\n");
	voronoiCell->neighborCells[voronoiCell->countNeighborCells] = neighborCell;
	//printf("increase neighbor counter\n");
	voronoiCell->countNeighborCells++;	
	//printf("done!\n");

	return 1;
}
/*****************************************************************************/



void removeVoronoiCell( Triangulation* voronoiDiagram, Vertex* removedVoronoiCell)
{
	int i, j, k, l, m;

	/*if(removedVoronoiCell->index==4674){
		double volume = 0.;
		for( int i=0; i<voronoiDiagram->countTetrahedra; i++){
			volume += getVolume(
					voronoiDiagram->tetrahedra[i]->vertices[0]->position,
					voronoiDiagram->tetrahedra[i]->vertices[1]->position,
					voronoiDiagram->tetrahedra[i]->vertices[2]->position,
					voronoiDiagram->tetrahedra[i]->vertices[3]->position);
		}
		fprintf(stderr, "volume of all tets: %.10lf\n", volume);
	/ *FILE* fp = fopen("removedCell.dat","w+");
	fprintf(fp, "\n%lf %lf %lf %i\n", removedVoronoiCell->position[0], removedVoronoiCell->position[1], removedVoronoiCell->position[2], removedVoronoiCell->index);
	for(int i=0; i<removedVoronoiCell->countNeighborCells; i++){
		fprintf(fp, "%lf %lf %lf %i\n", removedVoronoiCell->neighborCells[i]->position[0], removedVoronoiCell->neighborCells[i]->position[1], removedVoronoiCell->neighborCells[i]->position[2], removedVoronoiCell->neighborCells[i]->index);
	}
	fclose(fp);
	checkDelaunayCondition(voronoiDiagram, 0, 0);* /
	}*/
	// tetrahedra
	int ti = voronoiDiagram->countTetrahedra;
	int timax = voronoiDiagram->maxTetrahedra;
	Tetrahedron **tetrahedra = voronoiDiagram->tetrahedra;

	// faces
	int fi = 0;			// face counter
	int fimax = voronoiDiagram->maxFaces;			// maximal number of faces in list
	Vertex ***faces = voronoiDiagram->faces;		// list of faces

	// points
	int pi = 0;			// face counter
	int api = 0;			// face counter
	int pimax = INIT_FACES_ARRAY_SIZE;	// maximal number of faces in list
	Vertex **points;	// list of faces
	points = (Vertex**) calloc( sizeof(Vertex*), INIT_FACES_ARRAY_SIZE);
	assert(points);
	int *pointsInvolved;
	pointsInvolved = (int*) calloc( sizeof(int), INIT_FACES_ARRAY_SIZE);
	assert(pointsInvolved);
	Vertex **allPoints;	// list of faces



	// DELETE POINT
	//voronoiDiagram->countVoronoiCells--;



	// SEARCH INVOLVED TETRAHEDRA

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	Tetrahedron **correspondingTetrahedron = (Tetrahedron**) calloc( sizeof(Tetrahedron*), fimax);

	//fprintf( stderr, "look for first tetrahedra containing deleted point\n");

	Tetrahedron *actualTet = getTetrahedronContainingPoint( voronoiDiagram, removedVoronoiCell);

	// look for all tetrahedra deleted by the point insertion	
	int maxTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;		
	Tetrahedron **tetStack = (Tetrahedron**) calloc( sizeof(Tetrahedron*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(tetStack);
	int maxUsedTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;		
	Tetrahedron **usedTetStack = (Tetrahedron**) calloc( sizeof(Tetrahedron*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(usedTetStack);

	tetStack[0] = actualTet; //tetrahedra[0];
	int tetStackLength = 1;	
	int usedTetStackLength = 0;	


	// TEST:
	/*int test_count_tetrahedra_deleted = 0;
	for( j=0; j<ti; j++){
		if( getDistanceOfPointToCircumsphereCenter( newVoronoiCell, tetrahedra[j]) < getCircumsphereRadius( tetrahedra[j]))
			test_count_tetrahedra_deleted++;
	}*/
	int countIts = 0;
	// END TEST

	//fprintf( stderr, "look for all tetrahedra deleted by the point deletion\n");


	do{
		// TAKE NEXT TETRAHEDRON FROM STACK
		tetStackLength--;
		actualTet = tetStack[tetStackLength];
		if( usedTetStackLength == maxUsedTetStackLength){
			maxUsedTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
			usedTetStack = (Tetrahedron**) realloc( usedTetStack, sizeof(Tetrahedron*) * maxUsedTetStackLength);
			assert(usedTetStack);
		}
		usedTetStack[usedTetStackLength] = actualTet;
		usedTetStackLength++;

		//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

		// ADD NEIGHBOR TETRAHEDRA TO STACK
		//int temp_tetStackLength = tetStackLength;
		for( j=0; j<actualTet->countNeighborTetrahedra /*&& temp_tetStackLength == tetStackLength*/; j++){
			if( tetStackLength + actualTet->countNeighborTetrahedra < maxTetStackLength){
				maxTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
				tetStack = (Tetrahedron**) realloc( tetStack, sizeof(Tetrahedron*) * maxTetStackLength);
				assert(tetStack);
			}
			//distanceToCircumsphereCenter = getDistanceOfPointToCircumsphereCenter( newVoronoiCell, actualTet->neighborTetrahedra[j]);
			for( k=0; k<usedTetStackLength && usedTetStack[k]!=actualTet->neighborTetrahedra[j]; k++);
			for( l=0; l<tetStackLength     && tetStack[l]!=actualTet->neighborTetrahedra[j];     l++);
			// circumsphere contains new point ?
			if( k==usedTetStackLength && l==tetStackLength && 
			    ( actualTet->neighborTetrahedra[j]->vertices[0]==removedVoronoiCell || actualTet->neighborTetrahedra[j]->vertices[1]==removedVoronoiCell || actualTet->neighborTetrahedra[j]->vertices[2]==removedVoronoiCell
#if DIMENSIONS == 3
			    		|| actualTet->neighborTetrahedra[j]->vertices[3]==removedVoronoiCell
#endif
			    		)){
				tetStack[tetStackLength]=actualTet->neighborTetrahedra[j];
				tetStackLength++;
			}
		}

		// SET j TO POSITION OF VERTEX EQUAL TO REMOVED POINT
		for( j=0; j<NR_TETRAHEDRON_POINTS && actualTet->vertices[j] != removedVoronoiCell; j++);
		j++;

#else //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	// for all tetrahedra
	for( i=0; i<ti; i++)
	{
		Tetrahedron *actualTet = tetrahedra[i];

		// for all vertices
		int containsRemovedPoint = 0;
		for( j=0; j<NR_TETRAHEDRON_POINTS && !containsRemovedPoint; j++){
			// tetrehedron contains removed point 
			if( actualTet->vertices[j] == removedVoronoiCell){
				//fprintf( stderr, "Tetrahedra %i contains removed point %i\n", i, removedVoronoiCell->index);
				containsRemovedPoint++;
			}
		}

		if( containsRemovedPoint){
			//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

#endif //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

			// ADD FACES TO LIST
			if( fi+NR_TETRAHEDRON_POINTS >= fimax){
				//fprintf( stderr, "INFO: realloc: %i -> %i!\n", fimax, fimax + FACES_ARRAY_EXTENSION_SIZE);
				faces = (Vertex***) realloc( faces, sizeof(Vertex**) * (fimax + FACES_ARRAY_EXTENSION_SIZE));
				assert(faces);
				for( l=fimax; l<fimax + FACES_ARRAY_EXTENSION_SIZE; l++){
					faces[l] = (Vertex**) calloc( sizeof(Vertex*), NR_FACE_POINTS);
					assert(faces[l]);
				}
				fimax += FACES_ARRAY_EXTENSION_SIZE;
			}
			//k=j;
			//fprintf( stderr, "j=%i, Add faces %i (", j, fi);
			for( l=0; l<NR_FACE_POINTS; l++){
				k=0;
				for( k=0; k<pi && points[k]!=actualTet->vertices[(j+l)%(NR_TETRAHEDRON_POINTS)]; k++);
				//if( pi==0 || points[k]!=actualTet->vertices[(j+l)%(NR_TETRAHEDRON_POINTS)]){
				if( k==pi){
					// realloc
					if( pi==pimax){
						pimax += FACES_ARRAY_EXTENSION_SIZE;
						points = (Vertex**) realloc( points, sizeof(Vertex*) * pimax);
						assert(points);
						pointsInvolved = (int*)    realloc( pointsInvolved,    sizeof(int) * pimax);
						assert(pointsInvolved);
					}
					//j=
					points[pi] = actualTet->vertices[(j+l)%(NR_TETRAHEDRON_POINTS)];
					pointsInvolved[pi] = 1;
					pi++;
				}else
					pointsInvolved[k]++;
				faces[fi][l%NR_FACE_POINTS]=actualTet->vertices[(j+l)%(NR_TETRAHEDRON_POINTS)];
				//fprintf( stderr, " %i", faces[fi][l]->index);

				//pointsInvolved[k]++;


				// REMOVE NEIGHBORSHIP RELATIONS TO REMOVED POINT
				removeNeighbor( removedVoronoiCell, faces[fi][l]);
				removeNeighbor( faces[fi][l], removedVoronoiCell);
			}
			//fprintf( stderr, ")\n");
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

			// DELETE TETRAHEDRA NEIGHBORSHIP RELATIONS

			for( l=actualTet->countNeighborTetrahedra-1; l>=0; l--){
				if( tetrahedronContainsFace( actualTet->neighborTetrahedra[l]->vertices, faces[fi]))
					correspondingTetrahedron[fi] = actualTet->neighborTetrahedra[l];
				removeTetrahedronNeighbor( actualTet->neighborTetrahedra[l], actualTet);
				removeTetrahedronNeighbor( actualTet, actualTet->neighborTetrahedra[l]);
			}

#endif	
			fi++;


			// DELETE TETRAHEDRON
				
			//fprintf( stderr, "INFO: delete tetrahedron!\n");
			//actualTet->countNeighborTetrahedra = 0;

			ti--;
			int tempIndex = actualTet->index;
			Tetrahedron * temp = tetrahedra[tempIndex];
	
			tetrahedra[tempIndex] = tetrahedra[ti];
			tetrahedra[tempIndex]->index = tempIndex;
	
			tetrahedra[ti] = temp;
			tetrahedra[ti]->index = ti;
				
			//printPoints( "Points: ", points, pi );


#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	}while( tetStackLength != 0);

	free( tetStack);
	free( usedTetStack);

#else // TRACK_DELAUNAY_REGION_NEIGHBORHOOD
			
			if( i>=tempIndex)
				i--;
		}
	}

#endif // TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	api = pi;
	allPoints = (Vertex**) calloc( sizeof(Vertex*), api);
	assert(allPoints);
	for( i=0; i<pi; i++)
		allPoints[i] = points[i];
		//pointsInvolved[i] = 2;

	//printPoints( "Points: ", allPoints, api );




	// CREATE NEW TETRAHEDRA FROM FACES
		
	int countAddedTetrahedra = 0;
	//int last_minCountValidTetrahedra = 1;
	//int minCountValidTetrahedra = fi;
	//int countNotAddedFacesInAChain = 0;

	for( i=0; pi!=0; i = (fi==0 ? 0 : (i+1)%fi)){ // for all faces

		// STOP ENDLESS LOOP
		if( countIts==2000){
			fprintf( stderr, "WARNING: Apparently there is a Problem in removeVoronoiCell()\n");
			//fprintf( stderr, "Algorithm imprisoned in endless loop while trying to remove point %i!\n", removedVoronoiCell->index);
			fprintf( stderr, "Print Point coordinates to file\n");
			/*FILE *fp;
			fp = fopen( "errorPointSet.nodes", "w+");
			fprintf( fp, "%i %i\n", voronoiDiagram->countVoronoiCells, DIMENSIONS);
			int index = 0;
			for( int v=0; v<voronoiDiagram->countVoronoiCells; v++)
			if(voronoiDiagram->vertices[v]->countNeighborCells>0 || voronoiDiagram->vertices[v] == removedVoronoiCell){
				if( voronoiDiagram->vertices[v] == removedVoronoiCell)
					fprintf( stderr, "Algorithm imprisoned in endless loop while trying to remove point %i!\n", index);

				for( int p=0; p<DIMENSIONS; p++)
					fprintf( fp, "%.20lf ", voronoiDiagram->vertices[v]->position[p]);
				fprintf( fp, "\n");
				index++;
			}
			fclose( fp);*/

			FILE *pFile = fopen( "errorPointSet.nodes", "w+");
			int index = 0;
			for( int v=0; v<voronoiDiagram->countVoronoiCells; v++)
			if(voronoiDiagram->vertices[v]->countNeighborCells>0 || voronoiDiagram->vertices[v] == removedVoronoiCell)
			{
				if( voronoiDiagram->vertices[v] == removedVoronoiCell){
					fprintf( stderr, "Algorithm imprisoned in endless loop while trying to remove point %i (in output file: %i)!\n", removedVoronoiCell->index, index);
					for( int p=0; p<DIMENSIONS; p++)
						fprintf( stderr, "%lf ", voronoiDiagram->vertices[v]->position[p]);
					fprintf( stderr, "\n");
				}
				fwrite( voronoiDiagram->vertices[v]->position , sizeof(voronoiDiagram->vertices[v]->position), 1, pFile );
				index++;
			}
			fclose(pFile);

			exit( 0);
		}
		countIts++;
		// STOP ENDLESS LOOP


		if( ti == timax){
			//fprintf( stderr, "INFO: realloc: %i -> %i!\n", timax, timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE);
			tetrahedra = (Tetrahedron **) realloc( tetrahedra, sizeof(Tetrahedron*)*(timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE));
			assert(tetrahedra);
			for( k=timax; k<timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE;k++)
				tetrahedra[k] = newTetrahedron();
			timax += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
		}
		//fprintf( stderr, "Complete new tetrahedron %i: (", ti);	
		// set face as points of new tetrahedron
		for( k=0; k<NR_FACE_POINTS; k++){
			//fprintf( stderr, " %i", faces[i][k]->index);	
			tetrahedra[ti]->vertices[k] = faces[i][k];
		}
		//fprintf( stderr, " ?)\n");	

		// look for point to complete tetrahedron
		int tetrahedronValid = 0;
		int countValidTetrahedra = 0;
		Vertex * validPoint = NULL;


		for( j=0; j<pi /*TEST && !tetrahedronValid*/; j++){ // for all points
			// check if point is not contained yet by face
			for( k=0; k<NR_FACE_POINTS && points[j]!=faces[i][k]; k++);
			if( k==NR_FACE_POINTS){
				//fprintf( stderr, "point %i is candidate to form a tetrahedron with face %i\n", points[j]->index, i);

				// set point as point of new tetrahedron
				tetrahedra[ti]->vertices[k] = points[j];
				tetrahedra[ti]->circumsphereInitialized = 0;
				
				// check validy of new tetrahedron
				tetrahedronValid = 1;

				// removed point lies in the circumsphere?
				if( getDistanceOfPointToCircumsphere( removedVoronoiCell, tetrahedra[ti]) < 0.){
					// other points lie in the circumsphere?
					/*for( k=0; k<pi; k++){ // for all points
						// check if point is contained tetrahedron
						for( l=0; l<NR_TETRAHEDRON_POINTS && points[k]!=tetrahedra[ti]->vertices[l]; l++);
						if( l==NR_TETRAHEDRON_POINTS){
							// is point in the new tetrahedrons circumsphere ???
							if( getDistanceOfPointToCircumsphere( points[k], tetrahedra[ti]) < 0.){
								tetrahedronValid = 0;
							}					
						}
					}*/
					for( k=0; k<api; k++){ // for all points
						// check if point is contained tetrahedron
						for( l=0; l<NR_TETRAHEDRON_POINTS && allPoints[k]!=tetrahedra[ti]->vertices[l]; l++);
						if( l==NR_TETRAHEDRON_POINTS){
							// is point in the new tetrahedrons circumsphere ???
							if( getDistanceOfPointToCircumsphere( allPoints[k], tetrahedra[ti]) < 0.){
								tetrahedronValid = 0;
							}					
						}
					}

//#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
					// tet exists already?
					if( correspondingTetrahedron[i]!=NULL){
						for( k=0; k<NR_TETRAHEDRON_POINTS && correspondingTetrahedron[i]->vertices[k]!=points[j]; k++);
						if( k<NR_TETRAHEDRON_POINTS)
							tetrahedronValid = 0;
					}
//#endif
					
				}else{
					//fprintf( stderr, "WARNING: new tetrahedron %i is not valid: removed point is not within the circumsphere\n", ti);
					//TEST 
					tetrahedronValid = 0;
				}

				if( tetrahedronValid){
					countValidTetrahedra++;
					validPoint = points[j];
				}
				//printPoints( "test1 points:      ", allPoints, api );
			}
		}

		//fprintf( stderr, "INFO: fi=%i, i=%i, ti=%i:  countValidTetrahedra=%i   < last_minCountValidTetrahedra=%i ?\n", fi, i, ti, countValidTetrahedra, last_minCountValidTetrahedra);

		/*if( countValidTetrahedra > 1){
			fprintf( stderr, "INFO: countValidTetrahedra = %i\n", countValidTetrahedra);
			for( k=0; k<NR_FACE_POINTS; k++){
				fprintf( stderr, " %i", faces[i][k]->index);
				//tetrahedra[ti]->vertices[k] = faces[i][k];
			}
			fprintf( stderr, "\n");
			if( countIts>1000)
				countValidTetrahedra=1;
		}
		if( countValidTetrahedra < 1){
			fprintf( stderr, "INFO: countValidTetrahedra = %i\n", countValidTetrahedra);
			for( k=0; k<NR_FACE_POINTS; k++){
				fprintf( stderr, " %i", faces[i][k]->index);
				//tetrahedra[ti]->vertices[k] = faces[i][k];
			}
			fprintf( stderr, "\n");
			exit( 0);
		}*/

		if( countValidTetrahedra == 1){
		//if( countValidTetrahedra && countValidTetrahedra == last_minCountValidTetrahedra){
			tetrahedra[ti]->vertices[DIMENSIONS] = validPoint;
			tetrahedra[ti]->circumsphereInitialized = 0;

			countAddedTetrahedra++;

			// print new tetrahedron
			//fprintf( stderr, "new tetrahedron %i is (", ti);	
			for( k=0; k<NR_FACE_POINTS; k++){
				addNeighbor( tetrahedra[ti]->vertices[k], validPoint);
				addNeighbor( validPoint, tetrahedra[ti]->vertices[k]);
				//fprintf( stderr, " %i", tetrahedra[ti]->vertices[k]->index);	
			}
			//fprintf( stderr, ")\n");	
			

			// DELETE USED FACE FROM LIST
			fi--;
			//fprintf( stderr, "Remove Face %i: (", i);	
			for( k=0; k<NR_FACE_POINTS; k++){
				//fprintf( stderr, " %i", faces[i][k]->index);	
				for( l=0; faces[i][k]!=points[l]; l++);
				pointsInvolved[l]--;
				faces[i][k] = faces[fi][k];
			}

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

			if(correspondingTetrahedron[i]!=NULL){
				addOrReplaceTetrahedronNeighbor( tetrahedra[ti], correspondingTetrahedron[i]);
				addOrReplaceTetrahedronNeighbor( correspondingTetrahedron[i], tetrahedra[ti]);
			}
			correspondingTetrahedron[i] = correspondingTetrahedron[fi];

#endif
			
			//fprintf( stderr, ")\n");	
			i--;
			// ADD OTHER FACES TO LIST
			for( j=1; j<NR_TETRAHEDRON_POINTS; j++){ // for all new faces
				int commonVertices = 0; //TODO: optimize: replace faceAlreadyInList by commonVertices!!!
				for( k=0; commonVertices<NR_FACE_POINTS && k<fi; k++){ // check all faces for a duplicate
					commonVertices = 0;
					for( l=0; l<NR_FACE_POINTS; l++){ // for all vertices of new face
						//commonVertices = 0;
						for( m=0; m<NR_FACE_POINTS; m++){ // for all vertices of face in list
							if( tetrahedra[ti]->vertices[(j+l)%NR_TETRAHEDRON_POINTS] == faces[k][m])
								commonVertices++;
						}
					}
				}
				// face already in the list ?
				if( commonVertices>=NR_FACE_POINTS){
					// already in the list => delete face
					k--;
					fi--;
					//fprintf( stderr, "Remove Face %i: (", k);	
					for( l=0; l<NR_FACE_POINTS; l++){
						//fprintf( stderr, " %i", faces[k][l]->index);
						// deleted face has common points with previous deleted faces?	
						for( m=0; faces[k][l]!=points[m]; m++);
						/*if(m>=pi){
							fprintf( stderr, "\n ERROR: face point not in point list!!!\n");	
							exit( 0);
						}*/
						pointsInvolved[m]--;
						if( pointsInvolved[m] == 0){
							// delete from point list
							pi--;
							points[m] = points[pi];
							pointsInvolved[m] = pointsInvolved[pi];
								
						}

						faces[k][l]=faces[fi][l];
					}
					//fprintf( stderr, ")\n");

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

					if(correspondingTetrahedron[k]!=NULL){
						addOrReplaceTetrahedronNeighbor( tetrahedra[ti], correspondingTetrahedron[k]);
						addOrReplaceTetrahedronNeighbor( correspondingTetrahedron[k], tetrahedra[ti]);
					}
					correspondingTetrahedron[k] = correspondingTetrahedron[fi];

#endif

				}else{
					// not yet in the list => add face	
					//fprintf( stderr, "Add Face %i: (", fi);	
					for( l=0; l<NR_FACE_POINTS; l++){
						faces[fi][l]=tetrahedra[ti]->vertices[(j+l)%NR_TETRAHEDRON_POINTS];
						//fprintf( stderr, " %i", faces[fi][l]->index);	
						for( m=0; faces[fi][l]!=points[m]; m++);
						/*if(m>=pi){
							fprintf( stderr, "\n ERROR: face point not in point list!!!\n");	
							exit( 0);
						}*/
						pointsInvolved[m]++;
					}
					//fprintf( stderr, ")\n");

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

					correspondingTetrahedron[fi] = tetrahedra[ti];

#endif

					fi++;
				}
	
			}
			//fprintf( stderr, "\n");	
			tetrahedra[ti]->index = ti;
			ti++;

			//countNotAddedFacesInAChain = 0;				
			//last_minCountValidTetrahedra = 1;
		}
		/*else{
			countNotAddedFacesInAChain++;
			if( countValidTetrahedra == 0){
				fprintf( stderr, "ERROR: countValidTetrahedra == 0!!!\n");
				exit( 0);
			}
		}
		if( fi>0 && countNotAddedFacesInAChain == fi){
			last_minCountValidTetrahedra = minCountValidTetrahedra;
			minCountValidTetrahedra = fi;
			countNotAddedFacesInAChain = 0;
			fprintf( stderr, "\n WARNING: More than one posibility to include new tetrahedra: %i!\n", last_minCountValidTetrahedra);
			fprintf( stderr, "INFO: ti=%i, fi=%i\n", ti, fi);
			printPoints( "all points:  ", allPoints, api );
			printPoints( "points left: ", points, pi );
			voronoiDiagram->countTetrahedra = ti;
			//checkDelaunayCondition( voronoiDiagram, NULL, 0);
			//exit( 0);
			//occurenceOfAmbiguity++;
		}*/

	}


	// DELETE POINT
	//free( removedVoronoiCell);
	

	// FREE ALLOCATED MEMORY

	//deleteMatrix( A, NR_TETRAHEDRON_POINTS);
	free( points);
	free( pointsInvolved);
	free( allPoints);

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	free( correspondingTetrahedron);

#endif	

	// RESET VALUES OF VORONOI DIAGRAM

	voronoiDiagram->countFaces = fi;			// face counter
	voronoiDiagram->maxFaces = fimax;			// maximal number of faces in list
	voronoiDiagram->faces = faces;		// list of faces

	voronoiDiagram->countTetrahedra = ti; 
	voronoiDiagram->maxTetrahedra = timax;
	voronoiDiagram->tetrahedra = tetrahedra;



	// FINAL TESTING
	/*if( occurenceOfAmbiguity)
		checkDelaunayCondition( voronoiDiagram, NULL, 0);
	*/
	
}
/*****************************************************************************/


void Triangulation::addVertex( Vertex* newVoronoiCell)
{
	newVoronoiCell->index = this->countVoronoiCells;

	// ADD NEW POINT TO POINT ARRAY
	if( this->countVoronoiCells==this->maxVoronoiCells){
		this->maxVoronoiCells += POINTS_ARRAY_EXTENSION_SIZE;
		this->vertices = (Vertex **) realloc( this->vertices, sizeof(Vertex*) * this->maxVoronoiCells);
		assert(this->vertices);
	}
	this->vertices[this->countVoronoiCells] = newVoronoiCell;
	this->countVoronoiCells++;
}
/*****************************************************************************/


void triangulateVoronoiCell( Triangulation* voronoiDiagram, Vertex* newVoronoiCell)
{
	int j, k, l;

	//double **A = newMatrix( DIMENSIONS + 1, DIMENSIONS + 1);

	int count_faces_added = 0;
	int count_common_faces_deleted = 0;
	int count_tetrahedra_deleted = 0;
	int count_tetrahedra_added = 0;
	//int npr = 0;

	// faces
	int fi = 0;			// face counter
	int fimax = voronoiDiagram->maxFaces;			// maximal number of faces in list
	Vertex ***faces = voronoiDiagram->faces;		// list of faces


	// tetrahedra
	int ti = voronoiDiagram->countTetrahedra;
	int timax = voronoiDiagram->maxTetrahedra;
	Tetrahedron **tetrahedra = voronoiDiagram->tetrahedra;

	//newVoronoiCell->index = voronoiDiagram->countVoronoiCells;


	// DELAUNAY TRIANGULATION (WATSON)

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	Tetrahedron **correspondingTetrahedron = (Tetrahedron**) calloc( sizeof(Tetrahedron*), fimax);
	assert(correspondingTetrahedron);

	Tetrahedron *actualTet = voronoiDiagram->getTetrahedronContainingPointInCircumSphere( newVoronoiCell);

	// look for all tetrahedra deleted by the point insertion
	int maxTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;
	Tetrahedron **tetStack = (Tetrahedron**) calloc( sizeof(Tetrahedron*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(tetStack);
	int maxUsedTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;
	Tetrahedron **usedTetStack = (Tetrahedron**) calloc( sizeof(Tetrahedron*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(usedTetStack);

	tetStack[0] = actualTet; //tetrahedra[0];
	int tetStackLength = 1;
	int usedTetStackLength = 0;


	// TEST:
	/*int test_count_tetrahedra_deleted = 0;
	for( j=0; j<ti; j++){
		if( getDistanceOfPointToCircumsphereCenter( newVoronoiCell, tetrahedra[j]) < getCircumsphereRadius( tetrahedra[j]))
			test_count_tetrahedra_deleted++;
	}*/
	// END TEST

	//fprintf( stderr, "look for all tetrahedra deleted by the point insertion\n");


	do{

		// END TEST

		// TAKE NEXT TETRAHEDRON FROM STACK
		tetStackLength--;
		actualTet = tetStack[tetStackLength];
		if( usedTetStackLength == maxUsedTetStackLength){
			maxUsedTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
			usedTetStack = (Tetrahedron**) realloc( usedTetStack, sizeof(Tetrahedron*) * maxUsedTetStackLength);
			assert(usedTetStack);
		}
		usedTetStack[usedTetStackLength] = actualTet;
		usedTetStackLength++;


		// TEST
		/*int i;
		for( i=0; i<tetStackLength; i++)
		{
			for( j=0; j<actualTet->countNeighborTetrahedra; j++){
				if( actualTet->neighborTetrahedra[j] == NULL)
					exit( 0);

				if( actualTet->neighborTetrahedra[j]->circumsphereInitialized == 1){
					if( actualTet->neighborTetrahedra[j]->circumsphere[0] == 0.123456789)
						exit( 0);
				}
			}

			if( tetStack[i]->circumsphereInitialized == 1){
				if( tetStack[i]->circumsphere[0] == 0.123456789)
					exit( 0);
			}
		}*/


		//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

		// ADD NEIGHBOR TETRAHEDRA TO STACK
		//int temp_tetStackLength = tetStackLength;
		for( j=0; j<actualTet->countNeighborTetrahedra /*&& temp_tetStackLength == tetStackLength*/; j++){
			if( tetStackLength + actualTet->countNeighborTetrahedra < maxTetStackLength){
				maxTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
				tetStack = (Tetrahedron**) realloc( tetStack, sizeof(Tetrahedron*) * maxTetStackLength);
				assert(tetStack);
			}
			//distanceToCircumsphereCenter = getDistanceOfPointToCircumsphereCenter( newVoronoiCell, actualTet->neighborTetrahedra[j]);
			for( k=0; k<usedTetStackLength && usedTetStack[k]!=actualTet->neighborTetrahedra[j]; k++);
			for( l=0; l<tetStackLength     && tetStack[l]!=actualTet->neighborTetrahedra[j];     l++);
			// circumsphere contains new point ?
			if( k==usedTetStackLength && l==tetStackLength &&
			    getDistanceOfPointToCircumsphere( newVoronoiCell, actualTet->neighborTetrahedra[j]) < 0./*getCircumsphereRadius( actualTet->neighborTetrahedra[j])*/){
				tetStack[tetStackLength] = actualTet->neighborTetrahedra[j];
				tetStackLength++;
				//fprintf( stderr, "INFO: Add neighbor tet %i (%i, %i, %i, %i) as candidate\n", actualTet->neighborTetrahedra[j]->index, actualTet->neighborTetrahedra[j]->vertices[0]->index, actualTet->neighborTetrahedra[j]->vertices[1]->index, actualTet->neighborTetrahedra[j]->vertices[2]->index, actualTet->neighborTetrahedra[j]->vertices[3]->index);

			}//else
				//fprintf( stderr, "INFO: Neighbor tet %i (%i, %i, %i, %i) doesn't contain new point\n", actualTet->neighborTetrahedra[j]->index, actualTet->neighborTetrahedra[j]->vertices[0]->index, actualTet->neighborTetrahedra[j]->vertices[1]->index, actualTet->neighborTetrahedra[j]->vertices[2]->index, actualTet->neighborTetrahedra[j]->vertices[3]->index);

		}

#else //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	for( j=0; j<ti; j++){

		Tetrahedron *actualTet = tetrahedra[j];
		if( getDistanceOfPointToCircumsphere( newVoronoiCell, actualTet) < 0. /*getCircumsphereRadius( actualTet)*/){
			//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

#endif //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

		// ADD FACES

		// reallocation if nessessary
		if( fi+NR_TETRAHEDRON_POINTS >= fimax){
			//fprintf( stderr, "INFO: realloc: %i -> %i!\n", fimax, fimax + FACES_ARRAY_EXTENSION_SIZE);
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
			correspondingTetrahedron = (Tetrahedron**) realloc( correspondingTetrahedron, sizeof(Tetrahedron*) * (fimax + FACES_ARRAY_EXTENSION_SIZE));
#endif
			faces = (Vertex***) realloc( faces, sizeof(Vertex**) * (fimax + FACES_ARRAY_EXTENSION_SIZE));
			assert(faces);
			for( l=fimax; l<fimax + FACES_ARRAY_EXTENSION_SIZE; l++){
				faces[l] = (Vertex**) calloc( sizeof(Vertex*), NR_FACE_POINTS);
				assert(faces[l]);
			}
			fimax += FACES_ARRAY_EXTENSION_SIZE;
		}

		for( k=0; k<NR_TETRAHEDRON_POINTS; k++){ // for all faces of the tetrahedron
			// add face
			for( l=0; l<NR_FACE_POINTS; l++){ // for all points of one face
				faces[fi][l]=actualTet->vertices[(l+k)%(NR_TETRAHEDRON_POINTS)];
			}

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

			//fprintf( stderr, "INFO: save neighbor tet!\n");
			for( l=0; l<actualTet->countNeighborTetrahedra &&
				!tetrahedronContainsFace( actualTet->neighborTetrahedra[l]->vertices,
				faces[fi]); l++);
			if( l!=actualTet->countNeighborTetrahedra){
				correspondingTetrahedron[fi] = actualTet->neighborTetrahedra[l];
				removeTetrahedronNeighbor( actualTet->neighborTetrahedra[l], actualTet);
				removeTetrahedronNeighbor( actualTet, actualTet->neighborTetrahedra[l]);
			}else{
				correspondingTetrahedron[fi] = NULL;
			}

#endif

			// sort point labels of face
			Vertex * temp;
			int dim;
			for( dim=NR_FACE_POINTS; dim>1; dim--){
				for( l=1; l<dim; l++){
					// (memory address)
					//if( faces[fi][l-1] > faces[fi][l])
					// (point index)
					if( faces[fi][l-1]->index > faces[fi][l]->index)
					{
						temp = faces[fi][l-1];
						faces[fi][l-1] = faces[fi][l];
						faces[fi][l] = temp;
					}
				}
			}


			fi++;
			count_faces_added++;
		}


		// DELETE TETRAHEDRA FROM LIST

		ti--;
		int tempIndex = actualTet->index;
		Tetrahedron * temp = tetrahedra[tempIndex];

		tetrahedra[tempIndex] = tetrahedra[ti];
		tetrahedra[tempIndex]->index = tempIndex;

		tetrahedra[ti] = temp;
		tetrahedra[ti]->index = ti;

		count_tetrahedra_deleted++;

		/*for( k=0; k<ti; k++){
			if( k != tetrahedra[k]->index){
				fprintf( stderr, "ERROR: index of tet %i wrong: %i ! ti=%i\n", k, tetrahedra[k]->index, ti);
				exit( 0);
			}
		}*/

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	}while( tetStackLength != 0);

	free( tetStack);
	free( usedTetStack);

#else // TRACK_DELAUNAY_REGION_NEIGHBORHOOD

		j--;
	}}

#endif // TRACK_DELAUNAY_REGION_NEIGHBORHOOD



		// DELETE COMMON FACES FROM LIST

		//fprintf( stderr, "Number of Faces: fi=%d", fi);
		for( j=0; j<fi-1; j++){
			//fprintf( stderr, "Duplicate Check for Face %d\n", j);

			// check for copy of face j
			int duplicateFound = 0;
			for( l=j+1; l<fi && !duplicateFound; l++){ // for all following faces
				// check if faces j and l are equal
				for( k=0; k<DIMENSIONS && faces[j][k] == faces[l][k]; k++); // for all points of the faces
				if( k==DIMENSIONS && faces[j][k-1] == faces[l][k-1]){
					duplicateFound++;
					//fprintf( stderr, "%d. Face Duplicate found!!! Face %d (%d %d %d) and %d (%d %d %d) are equal. \n", duplicateFound, j, faces[j][0], faces[j][1], faces[j][2],  l, faces[l][0], faces[l][1], faces[l][2]);
					//exit(0);

					// remove neighborships
					for( k=0; k<NR_FACE_POINTS; k++){
						removeNeighbor( faces[l][k], faces[l][(k+1)%NR_FACE_POINTS]);
						removeNeighbor( faces[l][(k+1)%NR_FACE_POINTS], faces[l][k]);
					}

					// remove face l
					//int faceContainsInsertedPoint = 0;
					for( k=0; k<NR_FACE_POINTS; k++){
						faces[l][k]=faces[fi-1][k];
					}
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
					correspondingTetrahedron[l]=correspondingTetrahedron[fi-1];
#endif
					//if
					//facesTetrahedronIndex[l]=facesTetrahedronIndex[fi-1];
					fi--;
					count_common_faces_deleted++;

					// remove face j
					for( k=0; k<NR_FACE_POINTS; k++){
						faces[j][k]=faces[fi-1][k];
					}
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
					correspondingTetrahedron[j]=correspondingTetrahedron[fi-1];
#endif
					fi--;
					count_common_faces_deleted++;
					j--;
				}
			}

		}
		//fprintf( stderr, "( after deleting common faces: fi=%d\n", fi);



		// ADD NEW TETRAHEDRA TO LIST

		for( j=0; j<fi; j++){ // for all faces
			if( ti == timax){
				//fprintf( stderr, "INFO: realloc: %i -> %i!\n", timax, timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE);
				tetrahedra = (Tetrahedron **) realloc( tetrahedra, sizeof(Tetrahedron*)*(timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE));
				assert(tetrahedra);
				for( k=timax; k<timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE;k++)
					tetrahedra[k] = newTetrahedron();
				timax += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
			}

			//tetrahedra[ti] = newTetrahedron();
			//fprintf( stderr, "Add tetrahedra %d: (", ti);
			for( k=0; k<NR_FACE_POINTS; k++){ // for all points of face
				tetrahedra[ti]->vertices[k] = faces[j][k]; // set point of face
				//fprintf( stderr, "%d ", tetrahedra[ti]->vertices[k]->index);
				addNeighbor( newVoronoiCell, faces[j][k]);
				addNeighbor( faces[j][k], newVoronoiCell);
				addNeighbor( faces[j][k], faces[j][(k+1)%NR_FACE_POINTS]);
				addNeighbor( faces[j][(k+1)%NR_FACE_POINTS], faces[j][k]);
			}
			tetrahedra[ti]->vertices[k] = newVoronoiCell; // set added point
			//fprintf( stderr, "%d )\n", tetrahedra[ti]->vertices[k]->index);
			//tetrahedra[ti]->included=1;
			tetrahedra[ti]->countNeighborTetrahedra=0;
			tetrahedra[ti]->circumsphereInitialized=0;
			tetrahedra[ti]->index=ti;

			count_tetrahedra_added++;
			ti++;

		}
		//fprintf( stderr, ", add %d Tetrehedra\n", fi);
		//fprintf( stderr, "to %d already existing Tetrehedra\n", ti);


#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

		// DESTINE NEIGHBORSHIPS OF NEW TETRAHEDRA
		for( j=ti-fi/*-1*/; j<ti-1; j++){ // for each new tetrahedra j
		//for( j=ti-fi-1; j<ti-1; j++){ // for each new tetrahedra j


			for( k=j+1; k<ti; k++){ // for all other new tetrahedra k
				//fprintf( stderr, "Compare tetrahedra %i vs %i, max %i\n", j, k, ti);
				if( countCommonPointOfTetrahedra( tetrahedra[j], tetrahedra[k]) == DIMENSIONS){
					// j and k are neighbors
					tetrahedra[j]->addTetrahedronNeighbor( tetrahedra[k]);
					tetrahedra[k]->addTetrahedronNeighbor( tetrahedra[j]);
				}
			}
		}

		for( j=ti-fi/*-1*/; j<ti; j++){ // for each new tetrahedra j
			for( k=0; k<fi; k++){ // for all surrounding tetrahedra k
				//fprintf( stderr, "INFO: (%i)\n", k);
				if( correspondingTetrahedron[k]!=NULL && countCommonPointOfTetrahedra( tetrahedra[j], correspondingTetrahedron[k]) == DIMENSIONS){
					// j and k are neighbors
					addOrReplaceTetrahedronNeighbor( tetrahedra[j], correspondingTetrahedron[k]);
					addOrReplaceTetrahedronNeighbor( correspondingTetrahedron[k], tetrahedra[j]);
				}
			}
		}

#endif

	// ADD NEW POINT TO POINT ARRAY
	/*if( voronoiDiagram->countVoronoiCells==voronoiDiagram->maxVoronoiCells){
		voronoiDiagram->maxVoronoiCells += POINTS_ARRAY_EXTENSION_SIZE;
		voronoiDiagram->vertices = (Vertex **) realloc( voronoiDiagram->vertices, sizeof(Vertex*) * voronoiDiagram->maxVoronoiCells);
	}
	voronoiDiagram->vertices[voronoiDiagram->countVoronoiCells] = newVoronoiCell;
	voronoiDiagram->countVoronoiCells++;
	*/

	// FREE MEMORY

	//deleteMatrix( A, DIMENSIONS+1);

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
	free( correspondingTetrahedron);
#endif

	// RESET VALUES OF VORONOI DIAGRAM

	voronoiDiagram->faces = faces;
	voronoiDiagram->maxFaces = fimax;
	voronoiDiagram->countFaces = fi;

	voronoiDiagram->tetrahedra = tetrahedra;
	voronoiDiagram->maxTetrahedra = timax;
	voronoiDiagram->countTetrahedra = ti;
}


void insertVoronoiCell( Triangulation* voronoiDiagram, Vertex* newVoronoiCell)
{
	int j, k, l;

	//double **A = newMatrix( DIMENSIONS + 1, DIMENSIONS + 1);

	int count_faces_added = 0;
	int count_common_faces_deleted = 0;
	int count_tetrahedra_deleted = 0;
	int count_tetrahedra_added = 0;
	//int npr = 0;

	// faces
	int fi = 0;			// face counter
	int fimax = voronoiDiagram->maxFaces;			// maximal number of faces in list
	Vertex ***faces = voronoiDiagram->faces;		// list of faces


	// tetrahedra
	int ti = voronoiDiagram->countTetrahedra;
	int timax = voronoiDiagram->maxTetrahedra;
	Tetrahedron **tetrahedra = voronoiDiagram->tetrahedra;



	// DELAUNAY TRIANGULATION (WATSON)

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	Tetrahedron **correspondingTetrahedron = (Tetrahedron**) calloc( sizeof(Tetrahedron*), fimax);
	assert(correspondingTetrahedron);

	Tetrahedron *actualTet = voronoiDiagram->getTetrahedronContainingPointInCircumSphere( newVoronoiCell);

	// look for all tetrahedra deleted by the point insertion	
	int maxTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;		
	Tetrahedron **tetStack = (Tetrahedron**) calloc( sizeof(Tetrahedron*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(tetStack);
	int maxUsedTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;		
	Tetrahedron **usedTetStack = (Tetrahedron**) calloc( sizeof(Tetrahedron*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(usedTetStack);

	tetStack[0] = actualTet; //tetrahedra[0];
	int tetStackLength = 1;	
	int usedTetStackLength = 0;	


	// TEST:
	/*int test_count_tetrahedra_deleted = 0;
	for( j=0; j<ti; j++){
		if( getDistanceOfPointToCircumsphereCenter( newVoronoiCell, tetrahedra[j]) < getCircumsphereRadius( tetrahedra[j]))
			test_count_tetrahedra_deleted++;
	}*/
	// END TEST

	//fprintf( stderr, "look for all tetrahedra deleted by the point insertion\n");


	do{
		// TAKE NEXT TETRAHEDRON FROM STACK
		tetStackLength--;
		actualTet = tetStack[tetStackLength];
		if( usedTetStackLength == maxUsedTetStackLength){
			maxUsedTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
			usedTetStack = (Tetrahedron**) realloc( usedTetStack, sizeof(Tetrahedron*) * maxUsedTetStackLength);
			assert(usedTetStack);
		}
		usedTetStack[usedTetStackLength] = actualTet;
		usedTetStackLength++;

		//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

		// ADD NEIGHBOR TETRAHEDRA TO STACK
		//int temp_tetStackLength = tetStackLength;
		for( j=0; j<actualTet->countNeighborTetrahedra /*&& temp_tetStackLength == tetStackLength*/; j++){
			if( tetStackLength + actualTet->countNeighborTetrahedra < maxTetStackLength){
				maxTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
				tetStack = (Tetrahedron**) realloc( tetStack, sizeof(Tetrahedron*) * maxTetStackLength);
				assert(tetStack);
			}
			//distanceToCircumsphereCenter = getDistanceOfPointToCircumsphereCenter( newVoronoiCell, actualTet->neighborTetrahedra[j]);
			for( k=0; k<usedTetStackLength && usedTetStack[k]!=actualTet->neighborTetrahedra[j]; k++);
			for( l=0; l<tetStackLength     && tetStack[l]!=actualTet->neighborTetrahedra[j];     l++);
			// circumsphere contains new point ?
			if( k==usedTetStackLength && l==tetStackLength &&
			    getDistanceOfPointToCircumsphere( newVoronoiCell, actualTet->neighborTetrahedra[j]) < 0./*getCircumsphereRadius( actualTet->neighborTetrahedra[j])*/){
				tetStack[tetStackLength] = actualTet->neighborTetrahedra[j];
				tetStackLength++;
				//fprintf( stderr, "INFO: Add neighbor tet %i (%i, %i, %i, %i) as candidate\n", actualTet->neighborTetrahedra[j]->index, actualTet->neighborTetrahedra[j]->vertices[0]->index, actualTet->neighborTetrahedra[j]->vertices[1]->index, actualTet->neighborTetrahedra[j]->vertices[2]->index, actualTet->neighborTetrahedra[j]->vertices[3]->index);

			}//else
				//fprintf( stderr, "INFO: Neighbor tet %i (%i, %i, %i, %i) doesn't contain new point\n", actualTet->neighborTetrahedra[j]->index, actualTet->neighborTetrahedra[j]->vertices[0]->index, actualTet->neighborTetrahedra[j]->vertices[1]->index, actualTet->neighborTetrahedra[j]->vertices[2]->index, actualTet->neighborTetrahedra[j]->vertices[3]->index);

		}

#else //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	for( j=0; j<ti; j++){

		Tetrahedron *actualTet = tetrahedra[j];
		if( getDistanceOfPointToCircumsphere( newVoronoiCell, actualTet) < 0. /*getCircumsphereRadius( actualTet)*/){
			//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

#endif //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

		// ADD FACES

		// reallocation if nessessary 
		if( fi+NR_TETRAHEDRON_POINTS >= fimax){
			//fprintf( stderr, "INFO: realloc: %i -> %i!\n", fimax, fimax + FACES_ARRAY_EXTENSION_SIZE);
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
			correspondingTetrahedron = (Tetrahedron**) realloc( correspondingTetrahedron, sizeof(Tetrahedron*) * (fimax + FACES_ARRAY_EXTENSION_SIZE));
			assert(correspondingTetrahedron);
#endif
			faces = (Vertex***) realloc( faces, sizeof(Vertex**) * (fimax + FACES_ARRAY_EXTENSION_SIZE));
			assert(faces);
			for( l=fimax; l<fimax + FACES_ARRAY_EXTENSION_SIZE; l++){
				faces[l] = (Vertex**) calloc( sizeof(Vertex*), NR_FACE_POINTS);
				assert(faces[l]);
			}
			fimax += FACES_ARRAY_EXTENSION_SIZE;
		}

		for( k=0; k<NR_TETRAHEDRON_POINTS; k++){ // for all faces of the tetrahedron
			// add face
			for( l=0; l<NR_FACE_POINTS; l++){ // for all points of one face
				faces[fi][l]=actualTet->vertices[(l+k)%(NR_TETRAHEDRON_POINTS)];
			}

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

			//fprintf( stderr, "INFO: save neighbor tet!\n");
			for( l=0; l<actualTet->countNeighborTetrahedra && 
				!tetrahedronContainsFace( actualTet->neighborTetrahedra[l]->vertices, 
				faces[fi]); l++);
			if( l!=actualTet->countNeighborTetrahedra){
				correspondingTetrahedron[fi] = actualTet->neighborTetrahedra[l];
				removeTetrahedronNeighbor( actualTet->neighborTetrahedra[l], actualTet);
				removeTetrahedronNeighbor( actualTet, actualTet->neighborTetrahedra[l]);
			}else{
				correspondingTetrahedron[fi] = NULL;
			}

#endif	

			// sort point labels of face (memory address)
			Vertex * temp;
			int dim;
			for( dim=NR_FACE_POINTS; dim>1; dim--){
				for( l=1; l<dim; l++){
					// (memory address)
					//if( faces[fi][l-1] > faces[fi][l])
					// (point index)
					if( faces[fi][l-1]->index > faces[fi][l]->index)
					{
						temp = faces[fi][l-1];
						faces[fi][l-1] = faces[fi][l];
						faces[fi][l] = temp;
					}
				}
			}

			fi++;			
			count_faces_added++;
		}


		// DELETE TETRAHEDRA FROM LIST

		ti--;
		int tempIndex = actualTet->index;
		Tetrahedron * temp = tetrahedra[tempIndex];

		tetrahedra[tempIndex] = tetrahedra[ti];
		tetrahedra[tempIndex]->index = tempIndex;

		tetrahedra[ti] = temp;
		tetrahedra[ti]->index = ti;
				
		count_tetrahedra_deleted++;

		/*for( k=0; k<ti; k++){
			if( k != tetrahedra[k]->index){
				fprintf( stderr, "ERROR: index of tet %i wrong: %i ! ti=%i\n", k, tetrahedra[k]->index, ti);
				exit( 0);
			}
		}*/

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	}while( tetStackLength != 0);

	free( tetStack);
	free( usedTetStack);

#else // TRACK_DELAUNAY_REGION_NEIGHBORHOOD

		j--;
	}}

#endif // TRACK_DELAUNAY_REGION_NEIGHBORHOOD



		// DELETE COMMON FACES FROM LIST

		//fprintf( stderr, "Number of Faces: fi=%d", fi);
		for( j=0; j<fi-1; j++){
			//fprintf( stderr, "Duplicate Check for Face %d\n", j);

			// check for copy of face j
			int duplicateFound = 0;
			for( l=j+1; l<fi && !duplicateFound; l++){ // for all following faces
				// check if faces j and l are equal
				for( k=0; k<DIMENSIONS && faces[j][k] == faces[l][k]; k++); // for all points of the faces
				if( k==DIMENSIONS && faces[j][k-1] == faces[l][k-1]){
					duplicateFound++;
					//fprintf( stderr, "%d. Face Duplicate found!!! Face %d (%d %d %d) and %d (%d %d %d) are equal. \n", duplicateFound, j, faces[j][0], faces[j][1], faces[j][2],  l, faces[l][0], faces[l][1], faces[l][2]);
					//exit(0);

					// remove neighborships
					for( k=0; k<NR_FACE_POINTS; k++){
						removeNeighbor( faces[l][k], faces[l][(k+1)%NR_FACE_POINTS]);
						removeNeighbor( faces[l][(k+1)%NR_FACE_POINTS], faces[l][k]);
					}
					
					// remove face l
					//int faceContainsInsertedPoint = 0;
					for( k=0; k<NR_FACE_POINTS; k++){
						faces[l][k]=faces[fi-1][k];
					}
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
					correspondingTetrahedron[l]=correspondingTetrahedron[fi-1];
#endif
					//if
					//facesTetrahedronIndex[l]=facesTetrahedronIndex[fi-1];
					fi--;
					count_common_faces_deleted++;
					
					// remove face j
					for( k=0; k<NR_FACE_POINTS; k++){
						faces[j][k]=faces[fi-1][k];
					}
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
					correspondingTetrahedron[j]=correspondingTetrahedron[fi-1];
#endif
					fi--;
					count_common_faces_deleted++;
					j--;
				}
			}		
			
		}
		//fprintf( stderr, "( after deleting common faces: fi=%d\n", fi);
			


		// ADD NEW TETRAHEDRA TO LIST

		for( j=0; j<fi; j++){ // for all faces
			if( ti == timax){
				//fprintf( stderr, "INFO: realloc: %i -> %i!\n", timax, timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE);
				tetrahedra = (Tetrahedron **) realloc( tetrahedra, sizeof(Tetrahedron*)*(timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE));
				if( tetrahedra==NULL){
					fprintf( stderr, "ERROR reallocating memory for tetrahedron array!\n");
					exit(0);
				}
				for( k=timax; k<timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE;k++)
					tetrahedra[k] = newTetrahedron();
				timax += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
			}

			//tetrahedra[ti] = newTetrahedron();
			//fprintf( stderr, "Add tetrahedra %d: (", ti);
			for( k=0; k<NR_FACE_POINTS; k++){ // for all points of face
				tetrahedra[ti]->vertices[k] = faces[j][k]; // set point of face
				//fprintf( stderr, "%d ", tetrahedra[ti]->vertices[k]->index);
				addNeighbor( newVoronoiCell, faces[j][k]);
				addNeighbor( faces[j][k], newVoronoiCell);
				addNeighbor( faces[j][k], faces[j][(k+1)%NR_FACE_POINTS]);
				addNeighbor( faces[j][(k+1)%NR_FACE_POINTS], faces[j][k]);
			}
			tetrahedra[ti]->vertices[k] = newVoronoiCell; // set added point
			//fprintf( stderr, "%d )\n", tetrahedra[ti]->vertices[k]->index);
			//tetrahedra[ti]->included=1;
			tetrahedra[ti]->countNeighborTetrahedra=0;
			tetrahedra[ti]->circumsphereInitialized=0;
			tetrahedra[ti]->index=ti;
			
			count_tetrahedra_added++;
			ti++;

		}
		//fprintf( stderr, ", add %d Tetrehedra\n", fi);


#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
           		
		// DESTINE NEIGHBORSHIPS OF NEW TETRAHEDRA
		for( j=ti-fi; j<ti; j++){ // for each new tetrahedra j
		//for( j=ti-fi-1; j<ti-1; j++){ // for each new tetrahedra j

           		
			for( k=j+1; k<ti; k++){ // for all other new tetrahedra k
				if( countCommonPointOfTetrahedra( tetrahedra[j], tetrahedra[k]) == DIMENSIONS){
					// j and k are neighbors
					tetrahedra[j]->addTetrahedronNeighbor( tetrahedra[k]);
					tetrahedra[k]->addTetrahedronNeighbor( tetrahedra[j]);	
				}
			}
		}
           		
		for( j=ti-fi-1; j<ti; j++){ // for each new tetrahedra j
			for( k=0; k<fi; k++){ // for all surrounding tetrahedra k
				//fprintf( stderr, "INFO: (%i)\n", k);
				if( correspondingTetrahedron[k]!=NULL && countCommonPointOfTetrahedra( tetrahedra[j], correspondingTetrahedron[k]) == DIMENSIONS){
					// j and k are neighbors
					addOrReplaceTetrahedronNeighbor( tetrahedra[j], correspondingTetrahedron[k]);
					addOrReplaceTetrahedronNeighbor( correspondingTetrahedron[k], tetrahedra[j]);
				}
			}
		}

#endif

	// ADD NEW POINT TO POINT ARRAY
	if( voronoiDiagram->countVoronoiCells==voronoiDiagram->maxVoronoiCells){
		voronoiDiagram->maxVoronoiCells += POINTS_ARRAY_EXTENSION_SIZE;
		voronoiDiagram->vertices = (Vertex **) realloc( voronoiDiagram->vertices, sizeof(Vertex*) * voronoiDiagram->maxVoronoiCells);
		assert(voronoiDiagram->vertices);
	}
	voronoiDiagram->vertices[voronoiDiagram->countVoronoiCells] = newVoronoiCell;
	voronoiDiagram->countVoronoiCells++;

	// FREE MEMORY

	//deleteMatrix( A, DIMENSIONS+1);

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
	free( correspondingTetrahedron);
#endif

	// RESET VALUES OF VORONOI DIAGRAM

	voronoiDiagram->faces = faces;
	voronoiDiagram->maxFaces = fimax;
	voronoiDiagram->countFaces = fi; 

	voronoiDiagram->tetrahedra = tetrahedra;
	voronoiDiagram->maxTetrahedra = timax;
	voronoiDiagram->countTetrahedra = ti; 
}
/*****************************************************************************/


Tetrahedron* newTetrahedron()
{
	// allocate memory
	Tetrahedron* newTet = (Tetrahedron*) malloc( sizeof(Tetrahedron));
	assert(newTet);
	
	// initialize values
	newTet->circumsphereInitialized = 0;
	
	newTet->countNeighborTetrahedra = 0;

	// return pointer
	return newTet;
}
/*****************************************************************************/


Vertex::~Vertex()
{
	free(this->neighborCells);
}
/*****************************************************************************/


Vertex::Vertex( double x, double y, double z)
{
	// initialize values
	this->position[0] = x;
	this->position[1] = y;
	this->position[2] = z;

	// neighborhood
	this->neighborCells = 0;
	this->countNeighborCells = 0;
	this->neighborCellsInitialized = FALSE;
	this->extendedNeighborCellsInitialized = FALSE;

	// agent
	this->agent = NULL;

	this->refined = false;
}
/*****************************************************************************/


void deleteVoronoiCell( Vertex* voronoiCell)
{
	free( voronoiCell->neighborCells);
	free( voronoiCell);
}
/*****************************************************************************/


Triangulation::~Triangulation()
{
	fprintf( stderr, "Triangulation::deleteVoronoiDiagram()\n");
	int i;
	
	// free tetrahedra
	for( i=0; i<this->maxTetrahedra; i++)
		free( (this->tetrahedra[i]));
	free( this->tetrahedra);

	// free faces
	for( i=0; i<this->maxFaces; i++)
		free( (this->faces[i]));
	free( this->faces);

	// free points
	for( i=0; i<this->countVoronoiCells; i++){
		delete this->vertices[i];
		//free( this->vertices[i]->neighborCells);
		//free( this->vertices[i]);
	}
	free( this->vertices);
	
	// free frame points
	for( i=0; i<this->countFramePoints; i++){
		delete this->framePoints[i];
		//free( this->framePoints[i]->neighborCells);
		//free( this->framePoints[i]);
	}
	free( this->framePoints);

	// free voronoi grid
	if( voronoiGridSet == TRUE){
		free( this->voronoiGridPoints);
		free( this->voronoiCellCountVoronoiGridPoints);
		for( int i=0; i<this->countVoronoiCells; i++)
			free( this->voronoiCellToVoronoiGridPoints[i]);
		free( this->voronoiCellToVoronoiGridPoints);
	}
	
	deleteDoubleMatrix(_A, NR_TETRAHEDRON_POINTS);
	deleteDoubleMatrix(_B, NR_TETRAHEDRON_POINTS);
	free(_x);
}
/*****************************************************************************/



Triangulation::Triangulation()
{
	// init vertices
	vertices = 0;
	this->countVoronoiCells = 0;
	this->maxVoronoiCells = 0;

	// init tets
	this->tetrahedra = 0;
	this->countTetrahedra = 0;
	this->maxTetrahedra = 0;

	// voronoiDiagram->faces;
	this->faces = 0;
	this->maxFaces = 0;
	this->countFaces = 0;
	
	// init Voronoi diagram
	voronoiGridSet = false;

	// init domain information
	domainSet = false;

	_A = newDoubleMatrix( NR_TETRAHEDRON_POINTS, NR_TETRAHEDRON_POINTS);
	_B = newDoubleMatrix( NR_TETRAHEDRON_POINTS, NR_TETRAHEDRON_POINTS);
	_x = (double*) malloc( sizeof(double) * NR_TETRAHEDRON_POINTS);
	assert( _x);
	for( int i=0; i<NR_TETRAHEDRON_POINTS; i++){
		for( int ii=0; ii<NR_TETRAHEDRON_POINTS; ii++){
			_A[i][ii]=0;
			_B[i][ii]=0;
		}
		_x[i]=0;
	}
}
/*****************************************************************************/


Triangulation* Triangulation::newVoronoiDiagramFromFile( char* filename){
/****************************************************************************
 * Definieren der Variablen                                                 *
 ****************************************************************************/

	//int DIMENSIONS = 3;

	Triangulation	*newVoronoiDiagram = new Triangulation();

    char fin[ FILENAMESIZE ];
    char f_node[ FILENAMESIZE ];
	sprintf( fin,	"%s.qele", filename );
	sprintf( f_node,"%s.qnode", filename );
    FILE*	pin;
    FILE*	fp_node;

    int		i, j, l, new_point,
			t, ta, tx,
			tetrahedra	= 0,
			last		= 0,
			lastpoint	= 0;
    long	triangle[ DIMENSIONS+2 ];

    char	buffer[ READBUFFERSIZE ], 
            cTri[ READBUFFERSIZE ],
		    buffer_node[ READBUFFERSIZE ];
    char*	ptr = buffer;

/**************************************************************************** 
 * Einlesen der Dreiecke			                                        *
 ****************************************************************************/

    

	/* open input file (*.qele)*/
    pin = fopen( fin, "r" );
    if( pin == NULL ){
      fprintf( stderr, "Error opening file %s\n", fin);
      exit( 1);
    }
    else fprintf( stderr, ">%s eingelesen\n", fin);

	/* open input file (*.qnode)*/
    fp_node = fopen( f_node, "r" );
    if( fp_node == NULL ){
      fprintf( stderr, "Error opening file %s\n", f_node);
      exit( 1 );
    }
    else fprintf( stderr, ">%s eingelesen\n", f_node);

	/* erste Zeile einlesen und Anzahl der Dreiecke auslesen */
    if( fgets( cTri, READBUFFERSIZE, pin ) == NULL )
      perror( "Error reading first line" );

    /* maximale Punkte = 3*Dreiecke;
     * jeder Punkt kann maximal doppelt so viele Nachbarn haben,
     * wie es Dreiecke gibt */
    tetrahedra = atoi( cTri );
    newVoronoiDiagram->countTetrahedra = tetrahedra;
    if (tetrahedra > MAX_TRIANGLES) {
       printf("\n Too many triangle (%d)", tetrahedra);
       exit(0);
       }
    //static int point[MAX_TRIANGLES][MAX_NNS]; /* * 2; * 3 */
    	fgets( buffer_node, READBUFFERSIZE, fp_node ); 

	//printf("Allocate Memory!\n");
	lastpoint = atoi( buffer_node);
	//printf("%i Points to initialize!\n", lastpoint);
	int **point;
	int *pointAllocSize;
	point = (int**) calloc( sizeof(int*), lastpoint);
	assert(point);
	pointAllocSize = (int*) calloc( sizeof(int), lastpoint);
	assert(pointAllocSize);
	for( i=0; i<lastpoint; i++){
		point[i] = (int*) calloc( sizeof(int), MAX_NNS);
		assert(point[i]);
		pointAllocSize[i] = MAX_NNS;
	}
	//printf("Finished: Allocate Memory!\n");

/****************************************************************************
 * Dreieck- in Punktbeziehung wandeln                                       *
 ****************************************************************************/
	newVoronoiDiagram->tetrahedra = ( Tetrahedron **) calloc ( tetrahedra, sizeof( Tetrahedron *));
	assert(newVoronoiDiagram->tetrahedra);
	for( t = 0; t < tetrahedra; t++ ){
		newVoronoiDiagram->tetrahedra[t] = newTetrahedron(); //( Tetrahedron *) calloc ( DIMENSIONS+1, sizeof( Tetrahedron));
		newVoronoiDiagram->tetrahedra[t]->index = t;
	}
	newVoronoiDiagram->maxTetrahedra =      tetrahedra;
	newVoronoiDiagram->countTetrahedra =    tetrahedra;

	/* Anzahl der Punkte speichern */
	newVoronoiDiagram->countVoronoiCells = lastpoint;            // Anzahl der Punkte speichern
	newVoronoiDiagram->maxVoronoiCells = lastpoint;            // Anzahl der Punkte speichern
	newVoronoiDiagram->vertices = ( Vertex**) malloc ( newVoronoiDiagram->countVoronoiCells * sizeof( Vertex*));
	assert(newVoronoiDiagram->vertices);
	for( i=0; i<newVoronoiDiagram->countVoronoiCells; i++){
		newVoronoiDiagram->vertices[i] = ( Vertex*) malloc ( sizeof( Vertex));
		assert(newVoronoiDiagram->vertices[i]);
	}

	/* alle Dreiecke verarbeiten */
	for( t = 0; t < tetrahedra; t++ ){
#if _COMMENTS_ > 2
		fprintf( stderr, "\rRead and Process Tetrahedra Data %.3lf%%                         \b", 100.*(t+1.)/tetrahedra);
#endif    	
		if( fgets( buffer, READBUFFERSIZE, pin ) == NULL )
			perror( "Error reading line1" );

		ptr = buffer;
		triangle[ 0 ] = t;

		/* alle Punkte eines Dreiecks einlesen */
		for( tx = 1; tx < DIMENSIONS+2; tx++){
			triangle[ tx ]  = strtol( ptr, &ptr, 0 );
			newVoronoiDiagram->tetrahedra[t]->vertices[tx-1] = newVoronoiDiagram->vertices[triangle[ tx ]];
		}

		/* zu jedem Punkt die Nachbarn dieses Dreiecks sichern */
		for( ta = 1; ta < DIMENSIONS+2; ta++ ){
			if( point[ triangle[ ta ] ][ 0 ] == 0 ){	// Punkt neu?
				point[ triangle[ ta ] ][ 0 ] = 1;
				last = 1;
			if( lastpoint < triangle[ ta ] )
				lastpoint = triangle[ ta ];
			}else
				last = point[ triangle[ ta ] ][ 0 ];	// Speicherposition
		for( tx = 1; tx < DIMENSIONS+2; tx++ ){
			new_point = 1;
			if( triangle[ tx ] == triangle[ ta ] )	// identischer Punkt
				new_point = 0;
			else 
				for( l = 1; l < last; l++ )		// alle Nachbarn pruefen
					if( point[ triangle[ ta ] ][ l ] == triangle[ tx ] )
						new_point = 0;	// Punkt ist schon Nachbar
				if( new_point != 0 ){
					if( pointAllocSize[triangle[ ta ]] == last){
						pointAllocSize[triangle[ ta ]] += MAX_NNS;
						point[triangle[ ta ]] = (int*) realloc( point[triangle[ ta ]], sizeof(int) * pointAllocSize[triangle[ ta ]]);
						assert(point[triangle[ ta ]]);
					}
					point[ triangle[ ta ] ][ last++ ] = triangle[ tx ];
				}
			}
			point[ triangle[ ta ] ][ 0 ] = last;
		}
	}

/****************************************************************************
 * Dreieckbeziehung einlesen                                                *
 ****************************************************************************/
	if( fgets( buffer, READBUFFERSIZE, pin ) != NULL )
	for( t = 0; t < tetrahedra; t++ ){
#if _COMMENTS_ > 2
		fprintf( stderr, "\rDetermine Neighborhood Relations between VoronoiCells %.3lf%%                    \b", 100.*(t+1.)/tetrahedra);
#endif    	
		if( fgets( buffer, READBUFFERSIZE, pin ) == NULL )
			perror( "Error reading line1" );

		ptr = buffer;

		/* alle Punkte eines Dreiecks einlesen */
		strtol( ptr, &ptr, 0 );
		for( tx = 0; tx < DIMENSIONS+1; tx++){
			int index  = strtol( ptr, &ptr, 0 );
			if( index >= 0)
				newVoronoiDiagram->tetrahedra[t]->addTetrahedronNeighbor( newVoronoiDiagram->tetrahedra[index]);
			//newVoronoiDiagram->tetrahedra[t]->vertices[tx-1] = &newVoronoiDiagram->vertices[triangle[ tx ]];
		}
	}else{
		perror( "Error reading line1" );
	}

/*	for( i=0; i<newVoronoiDiagram->countTetrahedra-1; i++){
		// find neighbors of tetrahedron i
		
		for( j=i+1; j<newVoronoiDiagram->countTetrahedra; j++){
			// check if tetrahedron j is a neighbor of tetrahedron i

			if( countCommonPointOfTetrahedra( newVoronoiDiagram->tetrahedra[i], newVoronoiDiagram->tetrahedra[j]) == DIMENSIONS){
				addTetrahedronNeighbor( newVoronoiDiagram->tetrahedra[i], newVoronoiDiagram->tetrahedra[j]);
				addTetrahedronNeighbor( newVoronoiDiagram->tetrahedra[j], newVoronoiDiagram->tetrahedra[i]);
			}
		}
	}
	int k;
	for( i=0; i<newVoronoiDiagram->countTetrahedra; i++){
		for( j=0; j<newVoronoiDiagram->tetrahedra[i]->countNeighborTetrahedra; j++){
			for( k=0; k<newVoronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->countNeighborTetrahedra &&
				newVoronoiDiagram->tetrahedra[i] != newVoronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->neighborTetrahedra[k]; k++);
			
			if( k==newVoronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->countNeighborTetrahedra){
				fprintf( stderr, "ERROR: neighborship relation between tets %p -> %p!!\n", newVoronoiDiagram->tetrahedra[i], newVoronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->neighborTetrahedra[k]);
				exit( 0);
			}
		}
	}*/


/****************************************************************************
 *  alle Punkte in Datei schreiben                                          *
 ****************************************************************************/

	/* Gitter initialisieren*/

	
	/* Punkte mit Nachbarn in Datei sichern */
	for( i = 0; i<lastpoint; i++ ){
#if _COMMENTS_ > 2
		fprintf( stderr, "\rInitialize Triangulation %.3lf%%                              \b", 100.*(i+1.)/newVoronoiDiagram->countVoronoiCells);
#endif		
		if( fgets( buffer_node, READBUFFERSIZE, fp_node ) == NULL )
			perror( "Error reading line4" );
		ptr = buffer_node;
	  
	  //printf("initialize point %i\n", i);
		/* set number of voronoi cell */
		newVoronoiDiagram->vertices[i]->index = i;

		newVoronoiDiagram->vertices[i]->agent = NULL;

		newVoronoiDiagram->vertices[i]->refined = false;

		/* set coordinates */
		int dim;
		for( dim=0; dim<DIMENSIONS; dim++)
			newVoronoiDiagram->vertices[i]->position[dim] = strtod( ptr, &ptr );
		
		/* number of neighbours */
		newVoronoiDiagram->vertices[i]->countNeighborCells = point[ i ][ 0 ] - 1;
		newVoronoiDiagram->vertices[i]->countFreeNeighborCells = point[ i ][ 0 ] - 1;
		//newVoronoiDiagram->vertices[i].free_neighbors = point[ i ][ 0 ] - 1;
		
		/* set neighbors of voronoi cell */
		newVoronoiDiagram->vertices[i]->neighborCells = ( Vertex**) calloc ( point[ i ][ 0 ] - 1, sizeof( Vertex*)); // allocate memory for list of neighbors
		assert(newVoronoiDiagram->vertices[i]->neighborCells);
		for( j = 0; j < point[i][0] - 1; j++) {
			newVoronoiDiagram->vertices[i]->neighborCells[j] = newVoronoiDiagram->vertices[point[ i ][ j + 1]]; // Nachbarn in Liste eifimaxuegen
		}
		newVoronoiDiagram->vertices[i]->neighborCellsInitialized = TRUE;

		// extended neighbors
		newVoronoiDiagram->vertices[i]->countExtendedNeighborCells = 0;
		newVoronoiDiagram->vertices[i]->countFreeExtendedNeighborCells = 0;
		newVoronoiDiagram->vertices[i]->extendedNeighborhood = NULL;
		newVoronoiDiagram->vertices[i]->extendedNeighborCellsInitialized = FALSE;

	}
	fprintf( stderr, "\n");
	
#ifdef LIMITED_NEIGHBOR_DISTANCE
	for( int i=0; i<newVoronoiDiagram->countVoronoiCells; i++){
		int countReachable = 0;
		for( int ii=0; ii<newVoronoiDiagram->vertices[i]->countNeighborCells; ii++){
			if( abs( (int)(newVoronoiDiagram->vertices[i]->position[0]) - (int)(newVoronoiDiagram->vertices[i]->neighborCells[ii]->position[0])) < LIMITED_NEIGHBOR_DISTANCE
#if DIMENSIONS >= 2
			    && abs( (int)(newVoronoiDiagram->vertices[i]->position[1]) - (int)(newVoronoiDiagram->vertices[i]->neighborCells[ii]->position[1])) < LIMITED_NEIGHBOR_DISTANCE
#endif
#if DIMENSIONS == 3
			    && abs( (int)(newVoronoiDiagram->vertices[i]->position[2]) - (int)(newVoronoiDiagram->vertices[i]->neighborCells[ii]->position[2])) < LIMITED_NEIGHBOR_DISTANCE
#endif
			){
			
				//countReachable++;
				newVoronoiDiagram->vertices[i]->neighborCells[countReachable++] = newVoronoiDiagram->vertices[i]->neighborCells[ii];
			}
		}
		
		newVoronoiDiagram->vertices[i]->countNeighborCells     = countReachable;
		newVoronoiDiagram->vertices[i]->countFreeNeighborCells = countReachable;
	}
#endif		

	/*for( i = 0; i<lastpoint; i++ ){
		printf( "point %i: ", newVoronoiDiagram->vertices[i].index);
		for( j = 0; j < newVoronoiDiagram->vertices[i].countNeighborCells; j++)
			printf( " %i", newVoronoiDiagram->vertices[i].neighborCells[j]->index);
		printf( "\n");
	}*/
	
	/* Ein- und Ausgabedatei schließen */
    if( fclose( pin ) )
      printf( "can't close %s\n", fin );
    if( fclose( fp_node ) )
      printf( "can't close %s\n", f_node );

	for( i=0; i<newVoronoiDiagram->countVoronoiCells; i++)
		free( point[i]);
	free( point);
	free( pointAllocSize);
	return newVoronoiDiagram;
}
/*****************************************************************************/


Triangulation* Triangulation::newVoronoiDiagram( int x, int y, int z){

	Triangulation	*newVoronoiDiagram = new Triangulation();

	newVoronoiDiagram->countVoronoiCells = x*y*z;            // Anzahl der Punkte speichern
	newVoronoiDiagram->maxVoronoiCells = x*y*z;            // Anzahl der Punkte speichern
	newVoronoiDiagram->vertices = ( Vertex**) malloc ( newVoronoiDiagram->countVoronoiCells * sizeof( Vertex*));
	assert(	newVoronoiDiagram->vertices);
	//for( i=0; i<newVoronoiDiagram->countVoronoiCells; i++)
	//	newVoronoiDiagram->vertices[i] = ( Vertex*) malloc ( sizeof( Vertex));



	/* Punkte mit Nachbarn in Datei sichern */
	int i=0;
	// Regular Mesh
	//double alpha = 0.0001;
	// Irregular Mesh
	double alpha = 1.;
	for( int ix = 0; ix<x; ix++ )
#if DIMENSIONS >= 2
		for( int iy = 0; iy<y; iy++ )
#endif
#if DIMENSIONS == 3
			for( int iz = 0; iz<z; iz++ )
#endif
		{
#if _COMMENTS_ > 2
		fprintf( stderr, "\rInitialize Triangulation %.3lf%%                              \b", 100.*(i+1.)/newVoronoiDiagram->countVoronoiCells);
#endif

#if DIMENSIONS == 3
		newVoronoiDiagram->vertices[i] = new Vertex( (float)(ix+alpha*myRandE(1.)), (float)(iy+alpha*myRandE(1.)), (float)(iz+alpha*myRandE(1.)));
#elif DIMENSIONS == 2
		if( ix<=1 || ix>=x-2 || iy<=1 || iy>=y-2 )
			alpha = 0.0001;
		else
			alpha = 0.5;
		newVoronoiDiagram->vertices[i] = newVoronoiCell( (float)(ix+alpha*myRandE(1.)), (float)(iy+alpha*myRandE(1.)), 0);
#elif DIMENSIONS == 1
		newVoronoiDiagram->vertices[i] = newVoronoiCell( (float)(ix+alpha*myRandE(1.)), 0, 0);
#endif

		/* set number of voronoi cell */
		newVoronoiDiagram->vertices[i]->index = i;

		newVoronoiDiagram->vertices[i]->agent = NULL;

		newVoronoiDiagram->vertices[i]->refined = false;

		i++;
	}
	fprintf( stderr, "\n");

	newVoronoiDiagram->setFramePoints();
	newVoronoiDiagram->triangulate();

	// Post-Triangulation Initialization
	for( int i=0; i<newVoronoiDiagram->countVoronoiCells; i++ )
		newVoronoiDiagram->vertices[i]->countFreeNeighborCells = newVoronoiDiagram->vertices[i]->countNeighborCells;

	return newVoronoiDiagram;
}
/*****************************************************************************/


void printPoints( const char * title, Vertex ** allPoints, int api )
{
	int j;
	int negativeValues = 0;
		// print points
		fprintf( stderr, "%s ", title);	
		for( j=0; j<api; j++){
			fprintf( stderr, "%i ", allPoints[j]->index);
			if( allPoints[j]->index<0)
				negativeValues++;	
		}
		fprintf( stderr, "\n");
		if( negativeValues)
			exit( 0);	
}
/*****************************************************************************/


void checkDelaunayCondition( Triangulation *voronoiDiagram, Vertex **newCell, int countNewCells)
{
	int ii, l;
	fprintf( stderr, "Check Delaunay Condition of Voronoi Diagram...\n");

	for( int i=0; i<voronoiDiagram->countTetrahedra; i++)
		if( voronoiDiagram->tetrahedra[i]->circumsphereInitialized == 0)
			initCircumsphereRadius( voronoiDiagram->tetrahedra[i]);

#pragma omp parallel for private(l,ii)
	//for( int i=voronoiDiagram->countTetrahedra-1; i>=0; i--){
	for( int i=0; i<voronoiDiagram->countTetrahedra; i++){

		for( ii=0; ii<voronoiDiagram->countVoronoiCells; ii++){
			if( voronoiDiagram->vertices[ii]->countNeighborCells > 0){

				for( l=0; l<NR_TETRAHEDRON_POINTS && voronoiDiagram->vertices[ii]!=voronoiDiagram->tetrahedra[i]->vertices[l]; l++);
				if( 	l==NR_TETRAHEDRON_POINTS &&
					getDistanceOfPointToCircumsphere( voronoiDiagram->vertices[ii], voronoiDiagram->tetrahedra[i]) < 0.){
					fprintf( stderr, "... Delaunay Condition is not valid!!! :-(\n");
					/*fprintf( stderr, "%i\n", voronoiDiagram->tetrahedra[i]->index);
					fprintf( stderr, "%i\n", voronoiDiagram->tetrahedra[i]->vertices[0]->index);
					fprintf( stderr, "%i\n", voronoiDiagram->tetrahedra[i]->vertices[1]->index);
					fprintf( stderr, "%i\n", voronoiDiagram->tetrahedra[i]->vertices[2]->index);
					//fprintf( stderr, "%i\n", voronoiDiagram->tetrahedra[i]->vertices[3]->index);
					fprintf( stderr, "%i\n", voronoiDiagram->vertices[ii]->index);
					fprintf( stderr, "%f\n", voronoiDiagram->vertices[ii]->position[0]);
					fprintf( stderr, "%f\n", voronoiDiagram->vertices[ii]->position[1]);*/
					//fprintf( stderr, "%f\n", voronoiDiagram->vertices[ii]->position[2]);

					/*fprintf( stderr, "tet %i (%i [%p] %i [%p] %i [%p]) should not contain point %i [%p] (%lf %lf)\n",
										voronoiDiagram->tetrahedra[i]->index,
										voronoiDiagram->tetrahedra[i]->vertices[0]->index,
										voronoiDiagram->tetrahedra[i]->vertices[0],
										voronoiDiagram->tetrahedra[i]->vertices[1]->index,
										voronoiDiagram->tetrahedra[i]->vertices[1],
										voronoiDiagram->tetrahedra[i]->vertices[2]->index,
										voronoiDiagram->tetrahedra[i]->vertices[2],
										voronoiDiagram->vertices[ii]->index,
										voronoiDiagram->vertices[ii],
										voronoiDiagram->vertices[ii]->position[0],
										voronoiDiagram->vertices[ii]->position[1]);
*/
					fprintf( stderr, "tet %i (%i %i %i %i) should not contain point %i (%lf %lf %lf)\n",
						voronoiDiagram->tetrahedra[i]->index,
						voronoiDiagram->tetrahedra[i]->vertices[0]->index,
						voronoiDiagram->tetrahedra[i]->vertices[1]->index,
						voronoiDiagram->tetrahedra[i]->vertices[2]->index,
						voronoiDiagram->tetrahedra[i]->vertices[3]->index,
						voronoiDiagram->vertices[ii]->index,
						voronoiDiagram->vertices[ii]->position[0],
						voronoiDiagram->vertices[ii]->position[1],
						voronoiDiagram->vertices[ii]->position[2]);
					fprintf( stderr, "vertex %i (%lf %lf %lf)\n",
							voronoiDiagram->tetrahedra[i]->vertices[0]->index,
							voronoiDiagram->tetrahedra[i]->vertices[0]->position[0],
							voronoiDiagram->tetrahedra[i]->vertices[0]->position[1],
							voronoiDiagram->tetrahedra[i]->vertices[0]->position[2]);
					fprintf( stderr, "vertex %i (%lf %lf %lf)\n",
							voronoiDiagram->tetrahedra[i]->vertices[1]->index,
							voronoiDiagram->tetrahedra[i]->vertices[1]->position[0],
							voronoiDiagram->tetrahedra[i]->vertices[1]->position[1],
							voronoiDiagram->tetrahedra[i]->vertices[1]->position[2]);
					fprintf( stderr, "vertex %i (%lf %lf %lf)\n",
							voronoiDiagram->tetrahedra[i]->vertices[2]->index,
							voronoiDiagram->tetrahedra[i]->vertices[2]->position[0],
							voronoiDiagram->tetrahedra[i]->vertices[2]->position[1],
							voronoiDiagram->tetrahedra[i]->vertices[2]->position[2]);
					fprintf( stderr, "vertex %i (%lf %lf %lf)\n",
							voronoiDiagram->tetrahedra[i]->vertices[3]->index,
							voronoiDiagram->tetrahedra[i]->vertices[3]->position[0],
							voronoiDiagram->tetrahedra[i]->vertices[3]->position[1],
							voronoiDiagram->tetrahedra[i]->vertices[3]->position[2]);
					/*fprintf( stderr, "distance of point to circumsphere of tet = %e\n", getDistanceOfPointToCircumsphere( voronoiDiagram->vertices[ii], voronoiDiagram->tetrahedra[i]));
					fprintf( stderr, "distance of point to circumsphere of tet = %e\n", getDistanceOfPointToCircumsphere( voronoiDiagram->tetrahedra[i]->vertices[0], voronoiDiagram->tetrahedra[i]));
					fprintf( stderr, "distance of point to circumsphere of tet = %e\n", getDistanceOfPointToCircumsphere( voronoiDiagram->tetrahedra[i]->vertices[1], voronoiDiagram->tetrahedra[i]));
					fprintf( stderr, "distance of point to circumsphere of tet = %e\n", getDistanceOfPointToCircumsphere( voronoiDiagram->tetrahedra[i]->vertices[2], voronoiDiagram->tetrahedra[i]));
					fprintf( stderr, "distance of point to circumsphere of tet = %e\n", getDistanceOfPointToCircumsphere( voronoiDiagram->tetrahedra[i]->vertices[3], voronoiDiagram->tetrahedra[i]));*/

					exit( 0);
				}
			}else{
				//fprintf( stderr, "Point %i is not part of the triangulation and thus will be ignored!\n", voronoiDiagram->vertices[ii]->index);
			}
			/*for( ii=0; ii<countNewCells; ii++){
				for( l=0; l<NR_TETRAHEDRON_POINTS && newCell[ii]!=voronoiDiagram->tetrahedra[i]->vertices[l]; l++);
				if( l==NR_TETRAHEDRON_POINTS && getDistanceOfPointToCircumsphere( newCell[ii], voronoiDiagram->tetrahedra[i]) < 0){
					fprintf( stderr, "... Delaunay Condition is not valid!!! :-(\n");
					exit( 0);
				}
			}*/
		}
	}

	fprintf( stderr, "... Delaunay Condition is valid :-)\n");
}

/*
Hey salut Luna!
tu as encore reussi hier de prendre le dernier metro? Pour nous c'etait un peu chaud...en plus avec des trajets prolonges de la forme sinusoidaire a cause la biere ;) Je t'ecrit juste pour te proposer a boire quelque chose ce soir/weekend dans le cas que il te tombe le plafond sur la tete dans ta maison model ikea :) Il aura meme Chiara a Paris ce weekend je croit.

 */

/*****************************************************************************/



void Triangulation::getConvexHull( double thresholdDistance)
{
	this->convexFaces = (Vertex ***)malloc( sizeof(Vertex **) * this->countVoronoiCells);
	this->countConvexFaces = 0;
	this->maxConvexFaces = this->countVoronoiCells;
	//int countConvexFaces;
	//int maxConvexFaces;
	//Vertex ***convexFaces;

	Tetrahedron ** innerTets = (Tetrahedron **)malloc( sizeof(Tetrahedron *) * this->countVoronoiCells);
	int countInnerTets = 0;
	int maxInnerTets = this->countVoronoiCells;

	Tetrahedron ** usedTets = (Tetrahedron **)malloc( sizeof(Tetrahedron *) * this->countVoronoiCells);
	int countUsedTets = 0;
	int maxUsedTets = this->countVoronoiCells;

	int i, t, f;

	// ADD ALL FACES CONNECTED TO FRAME

	for( i=0; i<this->countTetrahedra; i++){ // for each existing tet
		// check if tet contains one of the frame points
		f=this->countFramePoints;
		int countFramePointsInTet = 0;
		int tt = 0, ff = 0;
		for( t=0; t<NR_TETRAHEDRON_POINTS; t++){
			//fprintf( stderr, "%i ", this->tetrahedra[i]->vertices[t]->index);
			for( f=0; f<this->countFramePoints && this->framePoints[f]!=this->tetrahedra[i]->vertices[t]; f++);
			if( f<this->countFramePoints){
				countFramePointsInTet++;
				tt = t;
				ff = f;
				//fprintf( stderr, "(f) ");
			}
		}
		//fprintf( stderr, "\n");
		//fprintf( stderr, "INFO: tet %i contains %i frame points\n", i, countFramePoints);


		int j, jj;
		double dist = 0.;
		double maxDist = 0.;
		if( countFramePointsInTet==1){
//		if( countFramePointsInTet>=1){
			/*fprintf( stderr, "INFO: vertice %i of tet %i (%i %i %i %i) contains frame point %i\n",
			         this->tetrahedra[i]->vertices[tt]->index, i,
			         this->tetrahedra[i]->vertices[0]->index,
			         this->tetrahedra[i]->vertices[1]->index,
			         this->tetrahedra[i]->vertices[2]->index,
			         this->tetrahedra[i]->vertices[3]->index,
			         this->framePoints[ff]->index);
			*/
			for( j=0; j<NR_TETRAHEDRON_POINTS-1; j++){
				for( jj=1; jj<NR_TETRAHEDRON_POINTS; jj++){
					if( j!=tt && jj!=tt){
						dist = pow( this->tetrahedra[i]->vertices[j]->position[0] - this->tetrahedra[i]->vertices[jj]->position[0], 2)
						     + pow( this->tetrahedra[i]->vertices[j]->position[1] - this->tetrahedra[i]->vertices[jj]->position[1], 2)
						     + pow( this->tetrahedra[i]->vertices[j]->position[2] - this->tetrahedra[i]->vertices[jj]->position[2], 2);
						dist = sqrt( dist);
						//fprintf( stderr, "INFO: distance of vertice %i and %i is %lf\n", this->tetrahedra[i]->vertices[j]->index, this->tetrahedra[i]->vertices[jj]->index, dist);
						if( maxDist<dist)
							maxDist = dist;
					}
				}
			}
		//}

		//if( countFramePoints==1 && maxDist<thresholdDistance){
			if( maxDist<thresholdDistance){
				if( this->countConvexFaces == this->maxConvexFaces){
					this->maxConvexFaces += FACES_ARRAY_EXTENSION_SIZE;
					this->convexFaces = (Vertex ***)realloc( this->convexFaces, sizeof(Vertex **) * this->maxConvexFaces);
				}

				t = tt+1;
				f = ff;

				//tet contains one of the frame points
				//fprintf( stderr, "INFO: vertice %i of tet %i contains frame point %i\n", this->tetrahedra[i]->vertices[t-1]->index, i, this->framePoints[f]->index);

				// add face
				this->convexFaces[this->countConvexFaces] = (Vertex **)malloc( sizeof(Vertex *) * NR_FACE_POINTS);
				//fprintf( stderr, "Add %ith Convex Face: ", this->countConvexFaces);
				for( f=0; f<NR_FACE_POINTS; f++){
					this->convexFaces[this->countConvexFaces][f] = this->tetrahedra[i]->vertices[t];
					//fprintf( stderr, " %i", this->convexFaces[this->countConvexFaces][f]->index);

					t++;
					t = t%NR_TETRAHEDRON_POINTS;
					//t = (++t)%NR_TETRAHEDRON_POINTS;
				}
				//fprintf( stderr, "\n");
				this->countConvexFaces++;

			}else{
				// find inner neighbors of actual tet
				//fprintf( stderr, "Add inner tets: ");
				for( j=0; j<this->tetrahedra[i]->countNeighborTetrahedra; j++){
					for( jj=0; jj<NR_TETRAHEDRON_POINTS && this->tetrahedra[i]->neighborTetrahedra[j]->vertices[jj]!=this->framePoints[ff]; jj++);
					if(jj==NR_TETRAHEDRON_POINTS){
						// inner tet
						for( jj=0; jj<countInnerTets && innerTets[jj]!=this->tetrahedra[i]->neighborTetrahedra[j]; jj++);
						if(jj==countInnerTets){
							if( countInnerTets == maxInnerTets){
								maxInnerTets += FACES_ARRAY_EXTENSION_SIZE;
								innerTets = (Tetrahedron **)realloc( innerTets, sizeof(Tetrahedron *) * maxInnerTets);
							}
							if( countUsedTets == maxUsedTets){
								maxUsedTets += FACES_ARRAY_EXTENSION_SIZE;
								usedTets = (Tetrahedron **)realloc( usedTets, sizeof(Tetrahedron *) * maxUsedTets);
							}
							innerTets[countInnerTets++] = this->tetrahedra[i]->neighborTetrahedra[j];
							usedTets[countUsedTets++] = this->tetrahedra[i];

							//fprintf( stderr, "%i (used:%i)", this->tetrahedra[i]->neighborTetrahedra[j]->index, this->tetrahedra[i]->index);
						}//else
					}
				}
				//fprintf( stderr, "\n");


			}
		}else{
			//fprintf( stderr, "INFO: tet %i doesn't contain any frame point\n", i);
		}

	}

	// ADD INNER FACES

	for( i=0; i<countInnerTets; i++){ // for each inner tet
		int j, jj, ii;
		Vertex* tempFace[DIMENSIONS];

		for( ii=0; ii<NR_TETRAHEDRON_POINTS; ii++){ // for each face

			// actual face
			for( j=0; j<NR_FACE_POINTS; j++)
				tempFace[j] = innerTets[i]->vertices[(ii+j)%NR_TETRAHEDRON_POINTS];

			// distance of face points
			double dist = 0.;
			double maxDist = 0.;

			for( j=0; j<NR_FACE_POINTS-1; j++){
				for( jj=j+1; jj<NR_FACE_POINTS; jj++){
					dist = pow( tempFace[j]->position[0] - tempFace[jj]->position[0], 2)
					     + pow( tempFace[j]->position[1] - tempFace[jj]->position[1], 2)
					     + pow( tempFace[j]->position[2] - tempFace[jj]->position[2], 2);
					dist = sqrt( dist);
					//fprintf( stderr, "INFO: distance of face points %i and %i is %lf (%i-%i)\n", tempFace[j]->index, tempFace[jj]->index, dist, j, jj);
					if( maxDist<dist)
						maxDist = dist;
				}
			}
			if( maxDist<thresholdDistance){
				if( this->countConvexFaces == this->maxConvexFaces){
					this->maxConvexFaces += FACES_ARRAY_EXTENSION_SIZE;
					this->convexFaces = (Vertex ***)realloc( this->convexFaces, sizeof(Vertex **) * this->maxConvexFaces);
				}

				// add face
				this->convexFaces[this->countConvexFaces] = (Vertex **)malloc( sizeof(Vertex *) * NR_FACE_POINTS);
				//fprintf( stderr, "Add %ith Convex Face: ", this->countConvexFaces);
				for( j=0; j<NR_FACE_POINTS; j++){
					this->convexFaces[this->countConvexFaces][j] = tempFace[j];
					//fprintf( stderr, " %i", this->convexFaces[this->countConvexFaces][j]->index);
				}
				//fprintf( stderr, "\n");
				this->countConvexFaces++;
				//fprintf( stderr, "INFO: distance of vertice %i and %i is %lf\n", this->tetrahedra[i]->vertices[j]->index, this->tetrahedra[i]->vertices[jj]->index, dist);
			}else{
				// add corresponding neighbor tet to inner tet list
				for( j=0; j<innerTets[i]->countNeighborTetrahedra && !tetrahedronContainsFace( innerTets[i]->neighborTetrahedra[j]->vertices, tempFace); j++);
				if( j<innerTets[i]->countNeighborTetrahedra){
					//fprintf( stderr, "Add %ith Neighbor tet: %i\n", j, innerTets[i]->neighborTetrahedra[j]->index);
					// is already used?
					for( jj=0; jj<countUsedTets && usedTets[jj]!=innerTets[i]->neighborTetrahedra[j]; jj++);
					if( jj<countUsedTets){
						//fprintf( stderr, "Neighbor tet %i is already used!!!\n", innerTets[i]->neighborTetrahedra[j]->index);
						//exit( 0);
					}else{
						// insert?
						//fprintf( stderr, "Add inner tets: ");
						for( jj=0; jj<countInnerTets && innerTets[jj]!=innerTets[i]->neighborTetrahedra[j]; jj++);
						if(jj==countInnerTets){
							if( countInnerTets == maxInnerTets){
								maxInnerTets += FACES_ARRAY_EXTENSION_SIZE;
								innerTets = (Tetrahedron **)realloc( innerTets, sizeof(Tetrahedron *) * maxInnerTets);
							}
							if( countUsedTets == maxUsedTets){
								maxUsedTets += FACES_ARRAY_EXTENSION_SIZE;
								usedTets = (Tetrahedron **)realloc( usedTets, sizeof(Tetrahedron *) * maxUsedTets);
							}
							innerTets[countInnerTets++] = innerTets[i]->neighborTetrahedra[j];
							usedTets[countUsedTets++] = innerTets[i];

							//fprintf( stderr, "%i (used:%i)", innerTets[i]->neighborTetrahedra[j]->index, innerTets[i]->index);
						}
						//fprintf( stderr, "\n");
					}
				}
			}
		}
	}
}
/*****************************************************************************/


void Triangulation::setConvexHullNeighborhood()
{
	int i, ii;
	//fprintf( stderr, "countVoronoiCells=%i\n", countVoronoiCells);
	for( i=0; i<this->countVoronoiCells; i++){
		this->vertices[i]->countNeighborCells = 0;
	}

	for( i=0; i<this->countConvexFaces; i++){
		for( ii=0; ii<NR_FACE_POINTS; ii++){
			addNeighbor( this->convexFaces[i][ii], this->convexFaces[i][(ii+1)%NR_FACE_POINTS]);
			addNeighbor( this->convexFaces[i][(ii+1)%NR_FACE_POINTS], this->convexFaces[i][ii]);
		}
	}


}
/*****************************************************************************/


void Triangulation::printToPovray( const char *filename,
		bool printPoints, bool printFramePoints, bool printNeighborship, bool printTetrahedra, bool printConvexHull)
{
	Triangulation *voronoiDiagram=this;
	Vertex **newCell=0;
	int countNewCells=0;

	//int NR_TETRAHEDRON_POINTS = DIMENSIONS+1;
	int i, j, k, l;

	FILE* fp;

	fp = fopen( filename, "w+" );

	// PRINT GRID TO POV-FILE
	fprintf( fp, "#include \"finish.inc\"\n"\
				"#include \"colors.inc\"\n"\
				"#include \"textures.inc\"\n"\
				"background { color White }\n"\
				"camera {location <%lf,%lf,%lf> look_at  <%lf,%lf,%lf>}\n"\
				"light_source {<%lf,%lf,%lf> color rgb <1, 1, 1>}\n",
				(this->xMin[0]+this->xMax[0]) / 2., (this->xMin[1]+this->xMax[1]) / 2., this->xMax[2] + this->xMax[0]-this->xMin[0],
				(this->xMin[0]+this->xMax[0]) / 2., (this->xMin[1]+this->xMax[1]) / 2., this->xMax[2],
				(this->xMin[0]+this->xMax[0]) / 2., (this->xMin[1]+this->xMax[1]) / 2., this->xMax[2] + this->xMax[0]-this->xMin[0]);

	// all points (GRAY)
	if( printPoints)
	for( i=0; i<voronoiDiagram->countVoronoiCells; i++){
		{
		fprintf( fp, "sphere {<");
		for( j=0; j<DIMENSIONS; j++){
			fprintf( fp, " %.3lf", voronoiDiagram->vertices[i]->position[j]);
			if( j+1!=DIMENSIONS)
				fprintf( fp, ",");
		}
		fprintf( fp, ">,0.1 texture {pigment { color rgbf <0.95, 0.95, 0.95, 0.2> }}}\n");
		}

		// all neighborships (YELLOW)
		if( printNeighborship){
		for( j=0; j<voronoiDiagram->vertices[i]->countNeighborCells; j++)
			if(printFramePoints || voronoiDiagram->vertices[i]->neighborCells[j]->index >= 0){
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		voronoiDiagram->vertices[i]->position[0],
				voronoiDiagram->vertices[i]->position[1],
				voronoiDiagram->vertices[i]->position[2],
	       		voronoiDiagram->vertices[i]->position[0]
				+ 0.4 * (voronoiDiagram->vertices[i]->neighborCells[j]->position[0] - voronoiDiagram->vertices[i]->position[0]),
	       		voronoiDiagram->vertices[i]->position[1]
				+ 0.4 * (voronoiDiagram->vertices[i]->neighborCells[j]->position[1] - voronoiDiagram->vertices[i]->position[1]),
				voronoiDiagram->vertices[i]->position[2]
				+ 0.4 * (voronoiDiagram->vertices[i]->neighborCells[j]->position[2] - voronoiDiagram->vertices[i]->position[2]),
				0.05,
	           	"Yellow");
			fprintf( fp, "\n");
		}
		}
	}

	// all frame points (GRAY)
	if( printFramePoints)
	for( i=0; i<voronoiDiagram->countFramePoints; i++){
		{
		fprintf( fp, "sphere {<");
		for( j=0; j<DIMENSIONS; j++){
			fprintf( fp, " %.3lf", voronoiDiagram->framePoints[i]->position[j]);
			if( j+1!=DIMENSIONS)
				fprintf( fp, ",");
		}
		fprintf( fp, ">,0.1 texture {pigment { color rgbf <0.95, 0.95, 0.95, 0.2> }}}\n");
		}

		// all neighborships (YELLOW)
		if( printNeighborship){
		for( j=0; j<voronoiDiagram->framePoints[i]->countNeighborCells; j++)
			if(printPoints || voronoiDiagram->framePoints[i]->neighborCells[j]->index < 0){
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		voronoiDiagram->framePoints[i]->position[0],
				voronoiDiagram->framePoints[i]->position[1],
				voronoiDiagram->framePoints[i]->position[2],
	       		voronoiDiagram->framePoints[i]->position[0]
				+ 0.4 * (voronoiDiagram->framePoints[i]->neighborCells[j]->position[0] - voronoiDiagram->framePoints[i]->position[0]),
	       		voronoiDiagram->framePoints[i]->position[1]
				+ 0.4 * (voronoiDiagram->framePoints[i]->neighborCells[j]->position[1] - voronoiDiagram->framePoints[i]->position[1]),
				voronoiDiagram->framePoints[i]->position[2]
				+ 0.4 * (voronoiDiagram->framePoints[i]->neighborCells[j]->position[2] - voronoiDiagram->framePoints[i]->position[2]),
				0.05,
	           	"Yellow");
			fprintf( fp, "\n");
		}
		}
	}

	// convex hull (BLUE)
	if( printConvexHull)
	for( j=0; j<voronoiDiagram->countConvexFaces; j++){
		for( k=0; k<NR_FACE_POINTS; k++)
		for( l=k+1; l<NR_FACE_POINTS; l++){ // TEST
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		voronoiDiagram->convexFaces[j][k]->position[0],
				voronoiDiagram->convexFaces[j][k]->position[1],
				voronoiDiagram->convexFaces[j][k]->position[2],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[0],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[1],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[2],
				voronoiDiagram->convexFaces[j][l]->position[0],
				voronoiDiagram->convexFaces[j][l]->position[1],
				voronoiDiagram->convexFaces[j][l]->position[2],
				0.04,
           		"rgbf <0, 0, 1, 0.2>");//"Blue");
			// zugehoerigkeitsmarkierung
			/*double centerX, centerY, centerZ;
			centerX = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[0] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[0]);
			centerY = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[1] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[1]);
			centerZ = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[2] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[2]);
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		centerX,
				centerY,
				centerZ,
				0.8 * centerX + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(DIMENSIONS+1)]->position[0],
	       		0.8 * centerY + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(DIMENSIONS+1)]->position[1],
				0.8 * centerZ + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(DIMENSIONS+1)]->position[2],
				0.02,
           		"rgbf <0, 0, 1, 0.2>");//"Blue");*/
		}
		fprintf( fp, "\n");
       		fprintf( fp, "triangle {<%lf, %lf, %lf>, <%lf, %lf, %lf>, <%lf, %lf, %lf> pigment { color %s }}\n",
       		voronoiDiagram->convexFaces[j][0]->position[0],
			voronoiDiagram->convexFaces[j][0]->position[1],
			voronoiDiagram->convexFaces[j][0]->position[2],
			voronoiDiagram->convexFaces[j][1]->position[0],
			voronoiDiagram->convexFaces[j][1]->position[1],
			voronoiDiagram->convexFaces[j][1]->position[2],
			voronoiDiagram->convexFaces[j][2]->position[0],
			voronoiDiagram->convexFaces[j][2]->position[1],
			voronoiDiagram->convexFaces[j][2]->position[2],
			//"rgbf <0, 0, 1, 0.2>"
			"rgbf <0.95, 0.95, 0.95, 0.2> "
			);
	}


	// all tetrahedra (BLUE)
	if( printTetrahedra)
	for( j=0; j<voronoiDiagram->countTetrahedra; j++){
		for( k=0; k<NR_TETRAHEDRON_POINTS; k++)
		for( l=k+1; l<NR_TETRAHEDRON_POINTS; l++){ // TEST
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    			voronoiDiagram->tetrahedra[j]->vertices[k]->position[0],
				voronoiDiagram->tetrahedra[j]->vertices[k]->position[1],
				voronoiDiagram->tetrahedra[j]->vertices[k]->position[2],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[0],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[1],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[2],
				voronoiDiagram->tetrahedra[j]->vertices[l]->position[0],
				voronoiDiagram->tetrahedra[j]->vertices[l]->position[1],
				voronoiDiagram->tetrahedra[j]->vertices[l]->position[2],
				0.04,
           		"rgbf <0, 0, 1, 0.2>");//"Blue");
			// zugehoerigkeitsmarkierung
			/*double centerX, centerY, centerZ;
			centerX = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[0] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[0]);
			centerY = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[1] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[1]);
			centerZ = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[2] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(DIMENSIONS+1)]->position[2]);
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		centerX,
				centerY,
				centerZ,
				0.8 * centerX + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(DIMENSIONS+1)]->position[0],
	       		0.8 * centerY + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(DIMENSIONS+1)]->position[1],
				0.8 * centerZ + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(DIMENSIONS+1)]->position[2],
				0.02,
           		"rgbf <0, 0, 1, 0.2>");//"Blue");*/
		}
		fprintf( fp, "\n");
	}

	for( i=0; i<countNewCells; i++){
	// new point
	fprintf( fp, "sphere {<");
	for( j=0; j<DIMENSIONS; j++){
		fprintf( fp, " %.3lf", newCell[i]->position[j]);
		if( j+1!=DIMENSIONS)
			fprintf( fp, ",");
	}
	fprintf( fp, ">,0.11 texture {pigment { color rgbf <0.95, 0., 0.95, 0.2> }}}\n");

	// all neighborships
	for( j=0; j<newCell[i]->countNeighborCells; j++){
		//for( k=0; k<NR_TETRAHEDRON_POINTS; k++){ // TEST
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		newCell[i]->position[0],
				newCell[i]->position[1],
				newCell[i]->position[2],
				newCell[i]->position[0]
				+ 0.4 * (newCell[i]->neighborCells[j]->position[0] - newCell[i]->position[0]),
	       		newCell[i]->position[1]
				+ 0.4 * (newCell[i]->neighborCells[j]->position[1] - newCell[i]->position[1]),
				newCell[i]->position[2]
				+ 0.4 * (newCell[i]->neighborCells[j]->position[2] - newCell[i]->position[2]),
				0.051,
           		"rgbf <0.95, 0., 0.95, 0.2>");
		//}
		fprintf( fp, "\n");
	}
	}


	fclose( fp);
}
/*****************************************************************************/



void Triangulation::setPointsOnSphericalSurface( double radius, int N, double variance, double center[DIMENSIONS])
{
	int k;
	double height, theta=0., phi=0., alpha, var_angle;
	double x, y, z, x_moved, y_moved, z_moved;

	// extend Voronoi diagram
	this->vertices = (Vertex**) realloc( this->vertices, sizeof(Vertex*)*(this->countVoronoiCells+N));
	assert(this->vertices);

	// polar coordinates of point k of N points
	//alpha = 0.5*2*acos(1-1/(2*PI*radius*radius*density)); // wenn Flächeninhalt der Kaloppe gleich 1/densitiy
	alpha = acos(1.-2./(double)N)*variance;

	/*for( countNewCells=1; countNewCells<=maxNewCells; countNewCells++){
		srand ( countNewCells+randomGeneratorSeed);
		//sprintf( filename, "test%i.pov", countNewCells);
		//testCell = newVoronoiCell( myRand()*8.+1., myRand()*8.+1., myRand()*8.+1.);
		//testCell = newVoronoiCell( myRand()*60., myRand()*60., myRand()*60.);
		testCell = newVoronoiCell( myRand()*2., myRand()*2., myRand()*2.);
		//printf("Insert %i. Voronoi Cell ( %lf, %lf, %lf)\n", countNewCells, testCell->position[0], testCell->position[1], testCell->position[2]);
		testCell->index = testGrid->countVoronoiCells;
		insertVoronoiCell( testGrid, testCell);
		//newCells[countNewCells-1] = testCell;
		//printToPovray( filename, testGrid, newCells, countNewCells);

	}*/
	srand(0);

	for(k=1; k<=N; k++){

		// spherical coordinates
		if( k==1){
			height = -1.;
			theta = acos(height);
			phi = 0.;
		}
		if( 1<k && k<N){
			height = -1. + 2.*((double)k-1.)/((double)N-1.);
			theta = acos(height);
			phi = (phi + 3.6/(sqrt((double)N) * sqrt(1.-height*height)));
		}
		if( k==N ){
			height = 1.;
			theta = acos(height);
			phi = 0.;
		}

		/*x = 1.;
		y = 0.;
		z = 0.;

		//rotateXYZ( x, y, z, x_moved, y_moved, z_moved, phi, 0, theta);
		rotateZ( x, y, z, x, y, z, theta);
		rotateX( x, y, z, x_moved, y_moved, z_moved, phi);*/


		// cartesian coordinates
		x = radius * sin(theta) * cos(phi);
		y = radius * sin(theta) * sin(phi);
		z = radius * cos(theta);

		// random angle of variance
		var_angle = 0.;//asin( sin(alpha) * myRand());

		x_moved = radius * sin(theta - var_angle) * cos(phi);
		y_moved = radius * sin(theta - var_angle) * sin(phi);
		z_moved = radius * cos(theta - var_angle);

		rotate( &x_moved, &y_moved, &z_moved, x, y, z, 2*PI*myRand()-PI);



		// set new voronoi cell
		//this->vertices[this->countVoronoiCells] = newVoronoiCell( x_moved, y_moved, z_moved);
		//this->countVoronoiCells++;
		Vertex* testCell = new Vertex( center[0]+x_moved+myRand()*POINT_POSITION_PRECISION, center[1]+y_moved+myRand()*POINT_POSITION_PRECISION, center[2]+z_moved+myRand()*POINT_POSITION_PRECISION);
		testCell->index = this->countVoronoiCells;
		this->addVertex( testCell);

		// set spherical coordinates
		testCell->sphericalCoordinates.phi = phi;
		testCell->sphericalCoordinates.theta = theta;

		//printf("%lf %lf %lf\n", x_moved, y_moved, z_moved);
		//printf("%lf %lf %lf\n", x, y, z);
	}
}
/*****************************************************************************/



void Triangulation::triangulate()
{
	int i;

	for( i=0; i<this->countVoronoiCells; i++){
		//printf("\rAdd Cell %i \b", i);
		//fprintf( stderr, "Add Cell %i \n", i);
		fprintf( stderr, "\r%.3lf%%\t\t\b", 100.*(i+1.)/this->countVoronoiCells);
		triangulateVoronoiCell( this, this->vertices[i]);
		//checkDelaunayCondition( this, 0, 0);
	}
}
/*****************************************************************************/



void Triangulation::sortTetrahedra()
{
	int i,ii;

	// sort vertex order of tets
	for( i=0; i<this->countTetrahedra; i++){
		for( int v=0; v<DIMENSIONS; v++)
			for( int vv=0; vv<DIMENSIONS-v; vv++)
				if( abs(this->tetrahedra[i]->vertices[vv]->index) > abs(this->tetrahedra[i]->vertices[vv+1]->index))
				{
					Vertex *temp = this->tetrahedra[i]->vertices[vv];
					this->tetrahedra[i]->vertices[vv] = this->tetrahedra[i]->vertices[vv+1];
					this->tetrahedra[i]->vertices[vv+1] = temp;
				}

		/*for( int v=0; v<DIMENSIONS+1; v++)
			fprintf( stderr, "%i ", this->tetrahedra[i]->vertices[v]->index);
		fprintf( stderr, "\n");*/
	}

	// sort tets
	for( i=0; i<this->countTetrahedra-1; i++)
		for( ii=0; ii<this->countTetrahedra-i-1; ii++)
			if( this->tetrahedra[ii]->vertices[0]->index <= this->tetrahedra[ii+1]->vertices[0]->index &&
				this->tetrahedra[ii]->vertices[1]->index <= this->tetrahedra[ii+1]->vertices[1]->index &&
				this->tetrahedra[ii]->vertices[2]->index <= this->tetrahedra[ii+1]->vertices[2]->index &&
				this->tetrahedra[ii]->vertices[3]->index <= this->tetrahedra[ii+1]->vertices[3]->index ){

				Tetrahedron *temp      = this->tetrahedra[ii];
				this->tetrahedra[ii]   = this->tetrahedra[ii+1];
				this->tetrahedra[ii+1] = temp;
			}

	for( i=0; i<this->countTetrahedra; i++){
		for( int v=0; v<DIMENSIONS+1; v++)
			fprintf( stderr, "%i ", this->tetrahedra[i]->vertices[v]->index);
		fprintf( stderr, "\n");
	}

}
/*****************************************************************************/



bool Triangulation::tetrahedronContainsFramePoints( Tetrahedron* tet)
{
	for( int v=0; v<NR_TETRAHEDRON_POINTS; v++){
		if( isElementOf( this->framePoints, this->countFramePoints, tet->vertices[v]))
			return true;
	}
	return false;
}
/*****************************************************************************/



void Triangulation::setFramePoints()
{

#if DIMENSIONS == 3

	//double dimMin[DIMENSIONS], dimMax[DIMENSIONS];
	int i, dim;

	if( this->countVoronoiCells==0)
		return;

	/*for( dim=0; dim<DIMENSIONS; dim++){
		dimMin[dim] = dimMax[dim] = this->vertices[0]->position[dim] ;
	}

	// find border
	for( i=1; i<this->countVoronoiCells; i++){
		for( dim=0; dim<DIMENSIONS; dim++){
			if( dimMin[dim] > this->vertices[i]->position[dim])
				dimMin[dim] = this->vertices[i]->position[dim] ;
			if( dimMax[dim] < this->vertices[i]->position[dim])
				dimMax[dim] = this->vertices[i]->position[dim] ;
		}
	}*/
	this->setDomain();

	// BACKUP CELL LIST
	//Vertex ** tempList = this->vertices;
	//int tempCountCells = this->countVoronoiCells;

	// NEW EMPTY & EXTENDED CELL LIST
	this->framePoints = (Vertex**) malloc( sizeof(Vertex*) * ((int)pow( 2., (int)DIMENSIONS)));
	assert(this->framePoints);
	this->countFramePoints = 0;
	//this->maxFramePoints   = tempCountCells + (int)pow( 2, DIMENSIONS);

	// SET INITIAL FRAME POINTS

	double strech = MAX( this->xMax[0]-this->xMin[0], MAX( this->xMax[1]-this->xMin[1], this->xMax[2]-this->xMin[2]));

	for( i=0; i<(int)pow( 2., DIMENSIONS); i++){
		int temp = i;
		int countOnes = 0;
		this->framePoints[this->countFramePoints] = new Vertex( 0., 0., 0.);
		for( dim=0; dim<DIMENSIONS; dim++){
			countOnes += temp%2;
			//printf( "%lf ", (1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim] + ((temp%2 - 0.5)));
			//printf( "%lf ", (1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim] + (2*(temp%2)-1)*myRand()*1.);

			this->framePoints[this->countFramePoints]->position[dim] = (1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim] + (2*(temp%2)-1)*myRand()*strech;
			temp = temp/2;
		}
		//fprintf( stderr, " (%i)", countOnes);
		this->framePoints[this->countFramePoints]->index = -1-i;

		temp = i;
		if( countOnes%2 == 1)
		for( dim=0; dim<DIMENSIONS; dim++){
			this->framePoints[this->countFramePoints]->position[dim] += (myRand()*strech*0.1+strech)*(2*(temp%2) - 1);
			temp = temp/2;
		}

		//for( dim=0; dim<DIMENSIONS; dim++)
			//printf( "%lf ", this->framePoints[this->countFramePoints]->position[dim]);
		//printf( "\n");

		this->countFramePoints++;
	}

	// SET INITIAL TETS
	this->countTetrahedra = 0;
	this->maxTetrahedra = 1 + 4 + 6;
	this->tetrahedra = (Tetrahedron**) realloc( this->tetrahedra, sizeof(Tetrahedron*)*this->maxTetrahedra);
	assert(this->tetrahedra);

	// central tet
	this->tetrahedra[this->countTetrahedra] = newTetrahedron();
	for( i=0; i<NR_TETRAHEDRON_POINTS; i++){
		int temp = 0;
		for( dim=0; dim<DIMENSIONS; dim++){
			temp = temp*2;
			if(dim==i || DIMENSIONS==i)
				temp++;
		}
		this->tetrahedra[this->countTetrahedra]->vertices[i] = this->framePoints[this->countFramePoints-temp-1];
		//fprintf( stderr, "%i ", 7-temp);
	}
	this->tetrahedra[this->countTetrahedra]->index = this->countTetrahedra;
	//fprintf( stderr, "\n");

	// 4 surrounding tets
	for( i=0; i<NR_TETRAHEDRON_POINTS; i++){

		// first tet vertex
		int temp = 0;
		int temp_others[DIMENSIONS];
		for( dim=0; dim<DIMENSIONS; dim++){
			temp = temp*2;
			if(dim!=i && DIMENSIONS!=i){
				temp++;
				temp_others[dim]=1;
			}else{
				temp_others[dim]=0;
			}
		}

		this->tetrahedra[this->countTetrahedra+1+i] = newTetrahedron();

		//this->tetrahedra[this->countTetrahedra+1+i]->addTetrahedronNeighbor( this->tetrahedra[this->countTetrahedra]    );
		//this->tetrahedra[this->countTetrahedra    ]->addTetrahedronNeighbor( this->tetrahedra[this->countTetrahedra+1+i]);

		this->tetrahedra[this->countTetrahedra+1+i]->vertices[0] = this->framePoints[this->countFramePoints-temp-1];
		//fprintf( stderr, "%i ", 7-temp);

		// following 3 vertices
		for( dim=0; dim<DIMENSIONS; dim++){
			int dim2;
			int temp = 0;
			for( dim2=0; dim2<DIMENSIONS; dim2++){
				temp = temp*2;
				if(dim!=dim2)
					temp += temp_others[dim2];
				else
					temp += 1 - temp_others[dim2];
			}
			this->tetrahedra[this->countTetrahedra+1+i]->vertices[1+dim] = this->framePoints[this->countFramePoints-temp-1];
			//fprintf( stderr, "%i ", this->tetrahedra[this->countTetrahedra+1+i]->vertices[1+dim]->index); //7-temp);
		}
		//fprintf( stderr, "\n");
		this->tetrahedra[this->countTetrahedra+1+i]->index = this->countTetrahedra+1+i;
	}
	this->countTetrahedra += 5;

	// 6 peripheric tets
	for( dim=0; dim<DIMENSIONS; dim++){ // for all dimensions
		for( i=0; i<2; i++){        // for each direction
			//fprintf( stderr, "tet %i:", this->countTetrahedra+dim*2+i);
			int ii;
			this->tetrahedra[this->countTetrahedra+dim*2+i] = newTetrahedron();
			this->tetrahedra[this->countTetrahedra+dim*2+i]->index = this->countTetrahedra+dim*2+i;
			for( ii=0; ii<4; ii++){ // for all tet points
				int temp2 = 0;
				int temp = ii;
				int dim2;
				for( dim2=0; dim2<DIMENSIONS; dim2++){
					temp2 = temp2*2;
					if(dim==dim2)
						temp2 += i;
					else{
						temp2 += temp%2;
						temp = temp/2;
					}
				}
				this->tetrahedra[this->countTetrahedra+dim*2+i]->vertices[ii] = this->framePoints[temp2];
				//fprintf( stderr, "%i ", this->tetrahedra[this->countTetrahedra+dim*2+i]->vertices[ii]->index);//temp2);

			}
			//fprintf( stderr, "\n");

		}

	}
	this->countTetrahedra += 6;
	setTetrahedronNeighbors( this);

	// SET NEIGHBORHOOD RELATIONS
	int ii, iii;
	for( i=0; i<this->countTetrahedra; i++){
		//fprintf( stderr, "%i \n", this->tetrahedra[i]->index);

		for( ii=0; ii<NR_TETRAHEDRON_POINTS; ii++){
			for( iii=ii+1; iii<NR_TETRAHEDRON_POINTS; iii++){
				//fprintf( stderr, "add %i & %i \n", this->tetrahedra[i]->vertices[ii]->index, this->tetrahedra[i]->vertices[iii]->index);
				addNeighbor( this->tetrahedra[i]->vertices[ii], this->tetrahedra[i]->vertices[iii]);
				addNeighbor( this->tetrahedra[i]->vertices[iii], this->tetrahedra[i]->vertices[ii]);
			}
		}
	}

	//addNeighbor(

	// TEST
	//for( i=0; i<tempCountCells; i++)

	// PRINT ALL

	/*fprintf( stderr, "countVoronoiCells:%i, maxVoronoiCells:%i\n", countVoronoiCells, maxVoronoiCells);
	for( i=0; i<this->countVoronoiCells; i++)
		fprintf( stderr, "%i: index=%i, (%lf %lf %lf)\n",
		         i, this->vertices[i]->index,
		         this->vertices[i]->position[0], this->vertices[i]->position[1], this->vertices[i]->position[2]);

	fprintf( stderr, "countFramePoints:%i, maxFramePoints:%i\n", countFramePoints, countFramePoints);
	for( i=0; i<this->countFramePoints; i++)
		fprintf( stderr, "%i: index=%i, (%lf %lf %lf)\n",
		         i, this->framePoints[i]->index,
		         this->framePoints[i]->position[0], this->framePoints[i]->position[1], this->framePoints[i]->position[2]);

	fprintf( stderr, "countTetrahedra:%i, maxTetrahedra:%i\n", countTetrahedra, maxTetrahedra);
	for( i=0; i<this->countTetrahedra; i++){
		fprintf( stderr, "%i: index=%i, ", i, this->tetrahedra[i]->index);
		fprintf( stderr, "(%i %i %i %i)",
		         this->tetrahedra[i]->vertices[0]->index, this->tetrahedra[i]->vertices[1]->index, this->tetrahedra[i]->vertices[2]->index, this->tetrahedra[i]->vertices[3]->index);
		fprintf( stderr, " => %i neighbors: (", this->tetrahedra[i]->countNeighborTetrahedra);
		for( int ii=0; ii<this->tetrahedra[i]->countNeighborTetrahedra; ii++)
			fprintf( stderr, " %i", this->tetrahedra[i]->neighborTetrahedra[ii]->index);
		fprintf( stderr, " )\n");
	}*/
	// TRIANGULATION
	//Vertex ** tempList = this->vertices;
	//int tempCountCells = this->countVoronoiCells;
	/*printToPovray( "test0.pov", this, NULL, 0);
	checkDelaunayCondition( this, NULL, 0);

	Triangulation *testDiagram = Triangulation::newVoronoiDiagramFromFile( "data");
	fprintf( stderr, "countVoronoiCells:%i, maxVoronoiCells:%i\n", countVoronoiCells, maxVoronoiCells);
	for( i=0; i<testDiagram->countVoronoiCells; i++)
		fprintf( stderr, "%i: index=%i, (%lf %lf %lf)\n",
		         i, testDiagram->vertices[i]->index,
		         testDiagram->vertices[i]->position[0], testDiagram->vertices[i]->position[1], testDiagram->vertices[i]->position[2]);

	fprintf( stderr, "countTetrahedra:%i, maxTetrahedra:%i\n", countTetrahedra, maxTetrahedra);
	for( i=0; i<testDiagram->countTetrahedra; i++){
		fprintf( stderr, "%i: index=%i, ", i, testDiagram->tetrahedra[i]->index);
		fprintf( stderr, "(%i %i %i %i)",
		         testDiagram->tetrahedra[i]->vertices[0]->index, testDiagram->tetrahedra[i]->vertices[1]->index, testDiagram->tetrahedra[i]->vertices[2]->index, testDiagram->tetrahedra[i]->vertices[3]->index);
		fprintf( stderr, " => %i neighbors: (", testDiagram->tetrahedra[i]->countNeighborTetrahedra);
		for( int ii=0; ii<testDiagram->tetrahedra[i]->countNeighborTetrahedra; ii++)
			fprintf( stderr, " %i", testDiagram->tetrahedra[i]->neighborTetrahedra[ii]->index);
		fprintf( stderr, " )\n");
	}

	char filename[100];
	for( i=0; i<tempCountCells; i++){
		fprintf( stderr, "Add point:%i  => (%lf %lf %lf)\n", tempList[i]->index = testDiagram->countVoronoiCells, tempList[i]->position[0], tempList[i]->position[1], tempList[i]->position[2]);
		insertVoronoiCell( testDiagram, tempList[i]);
		sprintf( filename, "anim_%i.pov", i+1);
		printToPovray( filename, testDiagram, tempList, i+1);
		//checkDelaunayCondition( testDiagram, NULL, 0);
	}	*/

#elif	DIMENSIONS == 2
	this->setDomain();

	// BACKUP CELL LIST
	//Vertex ** tempList = this->vertices;
	//int tempCountCells = this->countVoronoiCells;

	// NEW EMPTY & EXTENDED CELL LIST
	this->framePoints = (Vertex**) malloc( sizeof(Vertex*) * ((int)pow( 2, DIMENSIONS)));
	assert(this->framePoints);
	this->countFramePoints = 0;
	//this->maxFramePoints   = tempCountCells + (int)pow( 2, DIMENSIONS);

	//fprintf( stderr, " min frame:(%f, %f) - (%f, %f)\n", xMin[0], xMin[1], xMax[0], xMax[1]);

	// SET INITIAL FRAME POINTS
	for( int i=0; i<(int)pow( 2, DIMENSIONS); i++){
		int temp = i;
		int countOnes = 0;
		this->framePoints[this->countFramePoints] = newVoronoiCell( 0., 0., 0.);
		//fprintf( stderr, " (%i, %i) [%p]", i%2, i/2, this->framePoints[this->countFramePoints]);
		for( int dim=0; dim<DIMENSIONS; dim++){
			countOnes += temp%2;
			//fprintf( stderr, " %i,", temp%2);
			//fprintf( stderr, " %f+%f,",
			//		(1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim],
			//		( i%2 == i/2 ? 2.*(temp%2)-1. : (2.*(temp%2)-1.)*(1+1) ));
			//printf( "%lf ", (1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim] + ((temp%2 - 0.5)));
			//printf( "%lf ", (1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim] + (2*(temp%2)-1)*myRand()*1.);
			this->framePoints[this->countFramePoints]->position[dim] = (1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim] + ( i%2 != i/2 ? (2.*(temp%2)-1.)*myRand() : (2.*(temp%2)-1.)*(myRand()+1));
			temp = temp/2;
		}
		//fprintf( stderr, " (%i)", countOnes);
		//fprintf( stderr, " (%i, %i)", countOnes%2, temp);
		this->framePoints[this->countFramePoints]->index = -1-i;

		/*printf( " -> point %i: ", i);
		for( int dim=0; dim<DIMENSIONS; dim++)
			printf( "%lf ", this->framePoints[this->countFramePoints]->position[dim]);
		printf( "\n");*/
		//printf( "%i %i\n", i%2, i%3);

		this->countFramePoints++;
	}



	// SET INITIAL TRIANGLES
	this->countTetrahedra = 0;
	this->maxTetrahedra = 2;
	this->tetrahedra = (Tetrahedron**) realloc( this->tetrahedra, sizeof(Tetrahedron*)*this->maxTetrahedra);
	assert(this->tetrahedra);

	for( int t=0; t<2; t++){
		// alloc triangle
		this->tetrahedra[this->countTetrahedra] = newTetrahedron();

		// set vertices
		for( int i=0; i<NR_TETRAHEDRON_POINTS; i++){
			this->tetrahedra[this->countTetrahedra]->vertices[i] = this->framePoints[t+i];
			//fprintf( stderr, "%i ", 7-temp);
		}

		// set index
		this->tetrahedra[this->countTetrahedra]->index = this->countTetrahedra;
		//fprintf( stderr, "\n");

		this->countTetrahedra++;
	}

	setTetrahedronNeighbors( this);

	// SET NEIGHBORHOOD RELATIONS
	int ii, iii;
	for( int i=0; i<this->countTetrahedra; i++){
		//fprintf( stderr, "%i \n", this->tetrahedra[i]->index);

		for( ii=0; ii<NR_TETRAHEDRON_POINTS; ii++){
			for( iii=ii+1; iii<NR_TETRAHEDRON_POINTS; iii++){
				//fprintf( stderr, "add %i & %i \n", this->tetrahedra[i]->vertices[ii]->index, this->tetrahedra[i]->vertices[iii]->index);
				addNeighbor( this->tetrahedra[i]->vertices[ii], this->tetrahedra[i]->vertices[iii]);
				addNeighbor( this->tetrahedra[i]->vertices[iii], this->tetrahedra[i]->vertices[ii]);
			}
		}
	}

	/*this->countTetrahedra = 0;
	this->maxTetrahedra = 2;
	this->tetrahedra = (Tetrahedron**) realloc( this->tetrahedra, sizeof(Tetrahedron*)*this->maxTetrahedra);

	// central tet
	this->tetrahedra[this->countTetrahedra] = newTetrahedron();
	for( i=0; i<NR_TETRAHEDRON_POINTS; i++){
		int temp = 0;
		for( dim=0; dim<DIMENSIONS; dim++){
			temp = temp*2;
			if(dim==i || DIMENSIONS==i)
				temp++;
		}
		this->tetrahedra[this->countTetrahedra]->vertices[i] = this->framePoints[this->countFramePoints-temp-1];
		//fprintf( stderr, "%i ", 7-temp);
	}
	this->tetrahedra[this->countTetrahedra]->index = this->countTetrahedra;
	//fprintf( stderr, "\n");

	// 4 surrounding tets
	for( i=0; i<NR_TETRAHEDRON_POINTS; i++){

		// first tet vertex
		int temp = 0;
		int temp_others[DIMENSIONS];
		for( dim=0; dim<DIMENSIONS; dim++){
			temp = temp*2;
			if(dim!=i && DIMENSIONS!=i){
				temp++;
				temp_others[dim]=1;
			}else{
				temp_others[dim]=0;
			}
		}

		this->tetrahedra[this->countTetrahedra+1+i] = newTetrahedron();

		//this->tetrahedra[this->countTetrahedra+1+i]->addTetrahedronNeighbor( this->tetrahedra[this->countTetrahedra]    );
		//this->tetrahedra[this->countTetrahedra    ]->addTetrahedronNeighbor( this->tetrahedra[this->countTetrahedra+1+i]);

		this->tetrahedra[this->countTetrahedra+1+i]->vertices[0] = this->framePoints[this->countFramePoints-temp-1];
		//fprintf( stderr, "%i ", 7-temp);

		// following 3 vertices
		for( dim=0; dim<DIMENSIONS; dim++){
			int dim2;
			int temp = 0;
			for( dim2=0; dim2<DIMENSIONS; dim2++){
				temp = temp*2;
				if(dim!=dim2)
					temp += temp_others[dim2];
				else
					temp += 1 - temp_others[dim2];
			}
			this->tetrahedra[this->countTetrahedra+1+i]->vertices[1+dim] = this->framePoints[this->countFramePoints-temp-1];
			//fprintf( stderr, "%i ", this->tetrahedra[this->countTetrahedra+1+i]->vertices[1+dim]->index); //7-temp);
		}
		//fprintf( stderr, "\n");
		this->tetrahedra[this->countTetrahedra+1+i]->index = this->countTetrahedra+1+i;
	}
	this->countTetrahedra += 5;

	// 6 peripheric tets
	for( dim=0; dim<DIMENSIONS; dim++){ // for all dimensions
		for( i=0; i<2; i++){        // for each direction
			//fprintf( stderr, "tet %i:", this->countTetrahedra+dim*2+i);
			int ii;
			this->tetrahedra[this->countTetrahedra+dim*2+i] = newTetrahedron();
			this->tetrahedra[this->countTetrahedra+dim*2+i]->index = this->countTetrahedra+dim*2+i;
			for( ii=0; ii<4; ii++){ // for all tet points
				int temp2 = 0;
				int temp = ii;
				int dim2;
				for( dim2=0; dim2<DIMENSIONS; dim2++){
					temp2 = temp2*2;
					if(dim==dim2)
						temp2 += i;
					else{
						temp2 += temp%2;
						temp = temp/2;
					}
				}
				this->tetrahedra[this->countTetrahedra+dim*2+i]->vertices[ii] = this->framePoints[temp2];
				//fprintf( stderr, "%i ", this->tetrahedra[this->countTetrahedra+dim*2+i]->vertices[ii]->index);//temp2);

			}
			//fprintf( stderr, "\n");

		}

	}
	this->countTetrahedra += 6;
	setTetrahedronNeighbors( this);

	// SET NEIGHBORHOOD RELATIONS
	int ii, iii;
	for( i=0; i<this->countTetrahedra; i++){
		//fprintf( stderr, "%i \n", this->tetrahedra[i]->index);

		for( ii=0; ii<NR_TETRAHEDRON_POINTS; ii++){
			for( iii=ii+1; iii<NR_TETRAHEDRON_POINTS; iii++){
				//fprintf( stderr, "add %i & %i \n", this->tetrahedra[i]->vertices[ii]->index, this->tetrahedra[i]->vertices[iii]->index);
				addNeighbor( this->tetrahedra[i]->vertices[ii], this->tetrahedra[i]->vertices[iii]);
				addNeighbor( this->tetrahedra[i]->vertices[iii], this->tetrahedra[i]->vertices[ii]);
			}
		}
	}*/


#elif	DIMENSIONS == 1
	this->setDomain();

	// NEW EMPTY & EXTENDED CELL LIST
	this->framePoints = (Vertex**) malloc( sizeof(Vertex*) * 2);
	this->countFramePoints = 2;

	// SET INITIAL FRAME POINTS
	for( int i=0; i<DIMENSIONS+1; i++){
		this->framePoints[i] = newVoronoiCell( 0., 0., 0.);
		this->framePoints[i]->index = -i-1;
		this->framePoints[i]->position[0] = (i? this->xMin[0]-myRand() : this->xMax[0]+myRand());
	}

	// SET FRAME POINTS NEIGHBORHOOD RELATIONS
	addNeighbor( this->framePoints[0], this->framePoints[1]);
	addNeighbor( this->framePoints[1], this->framePoints[0]);

	// SET INITIAL TRIANGLES
	this->countTetrahedra = 1;
	this->maxTetrahedra = 1;
	this->tetrahedra = (Tetrahedron**) realloc( this->tetrahedra, sizeof(Tetrahedron*)*this->maxTetrahedra);
	this->tetrahedra[0] = newTetrahedron();
	this->tetrahedra[0]->vertices[0] = this->framePoints[0];
	this->tetrahedra[0]->vertices[1] = this->framePoints[1];
	this->tetrahedra[0]->index = 0;

	// SET TRIANGLES NEIGHBORHOOD RELATIONS
	// no neighbors in 1D

#endif

}


#ifdef USE_AGENT
void Vertex::actualizeFreeNeighbors()
{
	this->countFreeNeighborCells = 0;
	for( int n=0; n<this->countNeighborCells; n++)
		if( this->neighborCells[n]->isFree())
			this->countFreeNeighborCells++;
}


void Vertex::checkFreeNeighbors()
{
	int countFree = 0;
	for( int n=0; n<this->countNeighborCells; n++)
		if( this->neighborCells[n]->isFree())
			countFree++;

	if( this->countFreeNeighborCells != countFree)
	{
		fprintf( stderr, "Vertex %i has a wrong number of free neighbors: %i (!= %i found)\n", this->index, this->countFreeNeighborCells, countFree);
		exit( 0);
	}
}


void Triangulation::refine( Vertex *vc, int scale, ActionTree *actionTree)
{
	int base_x = vc->position[0];
	int base_y = vc->position[1];
#if DIMENSIONS >= 3
	int base_z = vc->position[2];
#endif

	//fprintf( stderr, "REFINE( cell: %i, scale: %i)\n", vc->index, scale);

	// LIST OF NEW VERTICES
	//int countNewVertices = scale*scale*scale;
	//fprintf( stderr, "countNewVertices = %i -> %i\n", countNewVertices, (this->countVoronoiCells+countNewVertices));

	// LIST OF NEIGHBORS
	Vertex *neighbors[vc->countNeighborCells];
	//int countNeighbors = vc->countNeighborCells;
	for( int n=0; n<vc->countNeighborCells; n++)
		neighbors[n] = vc->neighborCells[n];

	// INIT STACK
	int actualStackElement = 0;
	int stackSize = 0;
	Vertex *stack[2000];
	for( stackSize=0; stackSize<vc->countNeighborCells; stackSize++){
		// add neighbors of refined cells
		stack[stackSize] = vc->neighborCells[stackSize];
		//stack[stackSize]->countFreeNeighborCells = -1;
		//stackSize++;
		//fprintf( stderr, "-> add point %ith point to stack (neighbor of original cells)\n", stackSize+1);
	}

	// ADD NEW VERTICES & POINTS
	int i=this->countVoronoiCells;//-countNewVertices;
	int firstIndexNew = this->countVoronoiCells;
	Vertex *newVC;
	for( int x=0; x<scale; x++)
		for( int y=0; y<scale; y++)
#if DIMENSIONS >= 3
			for( int z=0; z<scale; z++)
#endif
			{

			// allocate Voronoi Cell
			double rx = (float)(base_x + (x+myRandE(1.))/(double)scale);
			double ry = (float)(base_y + (y+myRandE(1.))/(double)scale);
#if DIMENSIONS >= 3
			double rz = (float)(base_z + (z+myRandE(1.))/(double)scale);
#endif
			//this->vertices[i] = newVoronoiCell( rx, ry, rz);
#if DIMENSIONS >= 3
			newVC = newVoronoiCell( rx, ry, rz);
#elif DIMENSIONS == 2
			newVC = newVoronoiCell( rx, ry, 0);
#endif

			newVC->index = i;
			newVC->position[0] = rx;
			newVC->position[1] = ry;
#if DIMENSIONS >= 3
			newVC->position[2] = rz;
#endif
			newVC->glucose = vc->glucose;
			newVC->oxygen  = vc->oxygen;

			newVC->refined = true;

			newVC->coarseParent = vc;

			insertVoronoiCell( this, newVC);
			//fprintf( stderr, "-> add point %ith point to stack (refined cell)\n", stackSize+1);
			//checkDelaunayCondition( this, 0, 0);

			//fprintf( stderr, "-> add point %i (%lf, %lf, %lf)\n", i, rx, ry, (DIMENSIONS>=3?rz:0));
			//fprintf( stderr, "-> add point %i (%lf, %lf)\n", i, rx, ry);

			// add refined cells
			stack[stackSize] = newVC;
			//stack[stackSize]->countFreeNeighborCells = -1;
			stackSize++;

			i++;
		}

	// DELETE OLD POINT FROM TRIANGULATION
	//fprintf( stderr, "remove point %i\n", vc->index);
	removeVoronoiCell( this, vc);
	vc->refined = true;
	vc->refinement = firstIndexNew;

	//vc = newVC;

	// ADD NEIGHBORS OF REFINED CELLS
	for( int c=firstIndexNew; c<this->countVoronoiCells; c++){
		// add neighbors
		for( int n=0; n<this->vertices[c]->countNeighborCells; n++){
			if( !isElementOf( stack, stackSize, this->vertices[c]->neighborCells[n])){
				stack[stackSize] = this->vertices[c]->neighborCells[n];
				//stack[stackSize]->countFreeNeighborCells = -1;
				stackSize++;
				//fprintf( stderr, "-> add point %ith point to stack (neighbor of refined cell)\n", stackSize);
			}

			//fprintf( stderr, "STACK to actualize free neighbor cells in actualized lattice is not implemented yet!\n");
			//exit( 0);
		}
	}

	// ACTUALIZE ALL CONSIDERED VERTICES
	for( actualStackElement=0; actualStackElement<stackSize; actualStackElement++ ){
		// actualize number of free neighbor Voronoi cells
		stack[actualStackElement]->actualizeFreeNeighbors();

		// actualize attached agent
		if(	stack[actualStackElement]->agent!=NULL){
			//fprintf( stderr, "Actualize Agent on Voronoi cell %i\n", stack[actualStackElement]->index);
			GetAgent(stack[actualStackElement])->actualize( actionTree);
		}
	}
	//checkDelaunayCondition( this, 0,0);
	//fprintf( stderr, "...finished\n");
}

Vertex *Triangulation::coarsen( Vertex *vc, int scale, ActionTree *actionTree, AgentList *agentList)
{
	double base_x = floor(vc->position[0]);
	double base_y = floor(vc->position[1]);
#if DIMENSIONS == 3
	double base_z = floor(vc->position[2]);
#endif

#if DIMENSIONS == 3
	fprintf( stderr, "COARSEN( cell: %i (%i), scale: %i)\n", vc->index, (int)base_x*this->xN[1]*this->xN[2] + (int)base_y*this->xN[2] + (int)base_z, scale);
#else
	//fprintf( stderr, "COARSEN( cell: %i (%i), scale: %i)\n", vc->index, (int)base_x*this->xN[1] + (int)base_y, scale);
#endif
	assert( vc->refined);
	assert( vc->countNeighborCells>0);


	// SEARCH & DELETE CONSIDERED POINTS
	int actualStackElement = 0;
	int stackSize = 0;
	Vertex *stack[2000];
	stack[stackSize] = vc; stackSize++;

	int stack2Size = 0;
	Vertex *stack2[2000];

	/* do{ // OLD
		for( int n=0; n<stack[actualStackElement]->countNeighborCells; n++){
			if( stack[actualStackElement]->neighborCells[n]->position[0] > base_x &&
				stack[actualStackElement]->neighborCells[n]->position[0] < base_x+1 &&
				stack[actualStackElement]->neighborCells[n]->position[1] > base_y &&
				stack[actualStackElement]->neighborCells[n]->position[1] < base_y+1 &&
#if DIMENSIONS == 3
				stack[actualStackElement]->neighborCells[n]->position[2] > base_z &&
				stack[actualStackElement]->neighborCells[n]->position[2] < base_z+1 &&
#endif
				!isElementOf( stack, stackSize, stack[actualStackElement]->neighborCells[n])){

				stack[stackSize] = stack[actualStackElement]->neighborCells[n];
				stackSize++;
			}else{
				// mark neighbors to be actualized
				//stack[actualStackElement]->neighborCells[n]->countFreeNeighborCells = -1;
				stack2[stack2Size] = stack[actualStackElement]->neighborCells[n];
				stack2Size++;
			}
		}

//		removeVoronoiCell( this, stack[actualStackElement]);

		actualStackElement++;
	}while( actualStackElement<stackSize);*/
	//fprintf( stderr, "found %i cells to fuse\n", stackSize);

	stackSize=0;
	int i_base = vc->coarseParent->refinement;
	for( int i=0; i<scale*scale
#if DIMENSIONS == 3
								*scale
#endif
								; i++){
		//fprintf( stderr, "->coarsen %i\n", i_base+i);
		stack[stackSize] = this->vertices[i_base+i];
		stackSize++;
	}
	for( int i=0; i<stackSize; i++){
		for( int n=0; n<stack[i]->countNeighborCells; n++)
			if( !isElementOf( stack, stackSize, stack[i]->neighborCells[n]) &&
				!isElementOf( stack2, stack2Size, stack[i]->neighborCells[n])	){

				stack2[stack2Size] = stack[i]->neighborCells[n];
				stack2Size++;

			}
	}


	// REMOVE FINE POINTS FROM TRIANGULATION
	//fprintf( stderr, "removeVoronoiCell\n");
	for( int i=0; i<stackSize; i++){
		//fprintf( stderr, "removeVoronoiCell: %i \b", stack[i]->index);
		removeVoronoiCell( this, stack[i]);
		//fprintf( stderr, "...finished\n");

		// Deactivate Agent
		agentList->deactivateAgent( GetAgent(stack[i]));

		// Delete Actions from List
		actionTree->deleteAllActions( GetAgent(stack[i]));

		// Detach from Grid
		GetAgent(stack[i])->detach();

		// Free Memory of coarsened Lattice Point
		this->vertices[stack[i]->index]=NULL;
		free( stack[i]);
		stack[i]=NULL;
	}


	// REINSERT COARSE POINT INTO TRIANGULATION
	//fprintf( stderr, "triangulateVoronoiCell\n");
#if DIMENSIONS == 3
	int index = base_x*this->xN[1]*this->xN[2] + base_y*this->xN[2] + base_z;
#elif DIMENSIONS == 2
	int index = base_x*this->xN[1] + base_y;
#elif DIMENSIONS == 1
	int index = base_x;
#endif
	triangulateVoronoiCell( this, this->vertices[index]);
	this->vertices[index]->refined = false;


	// ACTUALIZE ALL CONSIDERED VERTICES

	// mark direct
	//fprintf( stderr, "mark direct\n");
	for( int n=0; n<this->vertices[index]->countNeighborCells; n++ )
		if( this->vertices[index]->neighborCells[n]->countFreeNeighborCells != -1){
			this->vertices[index]->neighborCells[n]->countFreeNeighborCells = -1;
			stack2[stack2Size] = this->vertices[index]->neighborCells[n];
			stack2Size++;
		}

	// mark indirect
	//fprintf( stderr, "mark indirect\n");
	for( int i=0; i<stack2Size; i++ )
		for( int n=0; n<stack2[i]->countNeighborCells; n++ ){
		//if( stack2[n]->countFreeNeighborCells != -1){
			stack2[i]->neighborCells[n]->countFreeNeighborCells = -1;
		}

	// actualize all
	//fprintf( stderr, "actualize all\n");
	int ii=0;
	for( int i=0; i<stack2Size; i++ )
		for( int n=0; n<stack2[i]->countNeighborCells; n++ ){
			if( stack2[i]->neighborCells[n]->countFreeNeighborCells == -1){

				// actualize number of free neighbor Voronoi cells
				stack2[i]->neighborCells[n]->actualizeFreeNeighbors();
				//stack[ii] = stack2[i]->neighborCells[n];
				//ii++;

				// actualize attached agent
				if(	stack2[i]->neighborCells[n]->agent!=NULL)
					GetAgent(stack2[i]->neighborCells[n])->actualize( actionTree);
			}
		}

//	for( int i=0; i<ii; i++ )
//		if(	stack[i]->agent!=NULL)
//			GetAgent(stack[i])->actualize( actionTree);
	fprintf( stderr, "...finished\n");

	return this->vertices[index];
}


bool Triangulation::isHomogen( Vertex *vc, int scale, int &type)
{
	double base_x = floor(vc->position[0]);
	double base_y = floor(vc->position[1]);
#if DIMENSIONS == 3
	double base_z = floor(vc->position[2]);
#endif

	//fprintf( stderr, "HOMOGEN? cell: %i (%i)\n", vc->index, (int)base_x*this->xN[1] + (int)base_y);
	//assert( vc->refined);
	//assert( vc->countNeighborCells>0);

	if( GetAgent(vc->coarseParent)->countFree == GetAgent(vc->coarseParent)->maxCellCount){
		type = FREE;
		//fprintf( stderr, "FREE\n"); //exit(0);
		return true;
	}else
	if(	GetAgent(vc->coarseParent)->countActive == GetAgent(vc->coarseParent)->maxCellCount){
		type = ACTIVE;
		//fprintf( stderr, "ACTIVE\n"); //exit(0);
		return true;
	}else
	if(	GetAgent(vc->coarseParent)->countNonactive == GetAgent(vc->coarseParent)->maxCellCount ){
		type = NONACTIVE;
		//fprintf( stderr, "NONACTIVE\n"); exit(0);
		return true;
	}else{
		//fprintf( stderr, "NOT HOMOGENOUS!\n");
		return false;
	}
	// TEST HOMOGENITY over CONSIDERED POINTS
	/*int actualStackElement = 0;
	int stackSize = 0;
	Vertex *stack[200];
	stack[stackSize] = vc; stackSize++;

	type = vc->getState();

	do{
		for( int n=0; n<stack[actualStackElement]->countNeighborCells; n++){
			if( stack[actualStackElement]->neighborCells[n]->position[0] > base_x &&
				stack[actualStackElement]->neighborCells[n]->position[0] < base_x+1 &&
				stack[actualStackElement]->neighborCells[n]->position[1] > base_y &&
				stack[actualStackElement]->neighborCells[n]->position[1] < base_y+1 &&
#if DIMENSIONS == 3
				stack[actualStackElement]->neighborCells[n]->position[2] > base_z &&
				stack[actualStackElement]->neighborCells[n]->position[2] < base_z+1 &&
#endif
				!isElementOf( stack, stackSize, stack[actualStackElement]->neighborCells[n])){
				//fprintf( stderr, "> neighbor: %i\n", stack[actualStackElement]->neighborCells[n]->index);
				stack[stackSize] = stack[actualStackElement]->neighborCells[n];
				stackSize++;

				if( stack[actualStackElement]->neighborCells[n]->getState() != type){
					type =-1;
					return false;
				}
			}
		}


//		removeVoronoiCell( this, stack[actualStackElement]);

		actualStackElement++;
	}while( actualStackElement<stackSize);*/

#if DIMENSIONS == 3
	int index = base_x*this->xN[1]*this->xN[2] + base_y*this->xN[2] + base_z;
#elif DIMENSIONS == 2
	int index = base_x*this->xN[1] + base_y;
#elif DIMENSIONS == 1
	int index = base_x;
#endif

	Vertex **refinement = &this->vertices[this->vertices[index]->refinement];
	/*fprintf( stderr, "fine point:   (%lf, %lf %lf)\n", vc->position[0], vc->position[1], vc->position[2]);
	fprintf( stderr, "coarse point: (%lf, %lf %lf)\n", this->vertices[index]->position[0], this->vertices[index]->position[1], this->vertices[index]->position[2]);
	if( this->vertices[index]->refinement==NULL)
		exit( 0);
	if( this->vertices[index]->refinement[0]==NULL)
		exit( 0);
	fprintf( stderr, "1st refined point: (%lf, %lf %lf)\n", this->vertices[index]->refinement[0]->position[0], this->vertices[index]->refinement[0]->position[1], this->vertices[index]->refinement[0]->position[2]);
*/
	/*type = refinement[0]->getState();
	for( int i=1; i<scale*scale
#if DIMENSIONS == 3
								*scale
#endif
								; i++){
		if(	refinement[i]->getState() != type)
			return false;
	}*/
	//fprintf( stderr, "homogeniously %i (%i cells)\n", type, actualStackElement);

	return true;
}
#endif

/*void Triangulation::coarsen( Vertex *vc, int scale)
{
	double base_x = floor(vc->position[0]);
	double base_y = floor(vc->position[1]);
#if DIMENSIONS == 3
	double base_z = floor(vc->position[2]);
#endif

	fprintf( stderr, "COARSEN( cell: %i, scale: %i)\n", vc->index, scale);
	assert( vc->refined);
	assert( vc->countNeighborCells>0);


	// SEARCH & DELETE CONSIDERED POINTS
	int actualStackElement = 0;
	int stackSize = 0;
	Vertex *stack[200];
	stack[stackSize] = vc; stackSize++;
	do{
		for( int n=0; n<stack[actualStackElement]->countNeighborCells; n++){
			if( stack[actualStackElement]->neighborCells[n]->position[0] > base_x &&
				stack[actualStackElement]->neighborCells[n]->position[0] < base_x+1 &&
				stack[actualStackElement]->neighborCells[n]->position[1] > base_y &&
				stack[actualStackElement]->neighborCells[n]->position[1] < base_y+1 &&
				!isElementOf( stack, stackSize, stack[actualStackElement]->neighborCells[n])){

				stack[stackSize] = stack[actualStackElement]->neighborCells[n];
				stackSize++;
			}
		}

		removeVoronoiCell( this, stack[actualStackElement]);

		actualStackElement++;
	}while( actualStackElement<stackSize);
	//fprintf( stderr, "found %i cells to fuse\n", stackSize);


	// REINSERT COARSE POINT INTO TRIANGULATION
#if DIMENSIONS == 3
	int index = base_x*this->xN[1]*this->xN[2] + base_y*this->xN[2] + base_z;
#elif DIMENSIONS == 2
	int index = base_x*this->xN[1] + base_y;
#endif
	triangulateVoronoiCell( this, this->vertices[index]);
	this->vertices[index]->refined = false;


	// ACTUALIZE ALL CONSIDERED VERTICES
	//actualStackElement = 0;
	//stackSize = 0;
	//stack[stackSize] = this->vertices[index]; stackSize++;

/ *
	// mark new cell to be actualized
	this->vertices[index]->countFreeNeighborCells = -1;

	// mark direct neighbors to be actualized
	for( int n=0; n<this->vertices[index]->countNeighborCells; n++ )
		this->vertices[index]->neighborCells[n]->countFreeNeighborCells = -1;

	// actualize indirect neighbors
	for( int n=0; n<this->vertices[index]->countNeighborCells; n++ )
		// direct neighbors
		for( int nn=0; nn<this->vertices[index]->neighborCells[n]->countNeighborCells; nn++ )
			// indirect neighbors
			if( this->vertices[index]->neighborCells[n]->neighborCells[nn]->countFreeNeighborCells != -1)
			{
				int countFree = 0;
				for( int nnn=0; nnn<this->vertices[index]->neighborCells[n]->neighborCells[nn]->countNeighborCells; nnn++ )
					if( this->vertices[index]->neighborCells[n]->neighborCells[nn]->neighborCells[nnn]->isFree())
						countFree++;
				this->vertices[index]->neighborCells[n]->neighborCells[nn]->countFreeNeighborCells = countFree;
				fprintf( stderr, "-> update free neighbors of %ith point\n", this->vertices[index]->neighborCells[n]->neighborCells[nn]->index);
			}

	// actualize direct neighbors
	for( int n=0; n<this->vertices[index]->countNeighborCells; n++ ){
		// direct neighbors
		int countFree = 0;
		for( int nn=0; nn<this->vertices[index]->neighborCells[n]->countNeighborCells; nn++ )
			if( this->vertices[index]->neighborCells[n]->neighborCells[nn]->isFree())
				countFree++;
		this->vertices[index]->neighborCells[n]->countFreeNeighborCells = countFree;
		fprintf( stderr, "-> update free neighbors of %ith point\n", this->vertices[index]->neighborCells[n]->index);
	}

	// actualize new cell
	int countFree = 0;
	for( int n=0; n<this->vertices[index]->countNeighborCells; n++ )
		if( this->vertices[index]->neighborCells[n]->isFree())
			countFree++;
	this->vertices[index]->countFreeNeighborCells = countFree;
	fprintf( stderr, "-> update free neighbors of %ith point\n", this->vertices[index]->index);
* /

	// mark considered cells to be actualized
	for( int n=0; n<this->vertices[index]->countNeighborCells; n++ )
		// direct neighbors
		for( int nn=0; nn<this->vertices[index]->neighborCells[n]->countNeighborCells; nn++ )
			// indirect neighbors
			this->vertices[index]->neighborCells[n]->neighborCells[nn]->countFreeNeighborCells = -1;

	// actualize marked cells
	for( int n=0; n<this->vertices[index]->countNeighborCells; n++ )
		// direct neighbors
		for( int nn=0; nn<this->vertices[index]->neighborCells[n]->countNeighborCells; nn++ )
			// indirect neighbors
			if( this->vertices[index]->neighborCells[n]->neighborCells[nn]->countFreeNeighborCells == -1)
			{
				int countFree = 0;
				for( int nnn=0; nnn<this->vertices[index]->neighborCells[n]->neighborCells[nn]->countNeighborCells; nnn++ )
					if( this->vertices[index]->neighborCells[n]->neighborCells[nn]->neighborCells[nnn]->isFree())
						countFree++;
				this->vertices[index]->neighborCells[n]->neighborCells[nn]->countFreeNeighborCells = countFree;
				fprintf( stderr, "-> update free neighbors of %ith point\n", this->vertices[index]->neighborCells[n]->neighborCells[nn]->index);
			}else{
				fprintf( stderr, "-> keep free neighbors of %ith point\n", this->vertices[index]->neighborCells[n]->neighborCells[nn]->index);

			}



	/ *for( int c=0; c<this->vertices[index]->countNeighborCells; c++ ){
		// mark
		int countFree = 0;
		for( int n=0; n<this->vertices[index]->neighborCells[c]->countNeighborCells; n++)
			if( this->vertices[index]->neighborCells[c]->neighborCells[n]->isFree())
				countFree++;

		//fprintf( stderr, "-> %s free neighbors of %ith point in stack\n",
//				(stack[actualStackElement]->countFreeNeighborCells == countFree ? "keep" : "update"), actualStackElement+1);
		this->vertices[index]->neighborCells[c]->countFreeNeighborCells = countFree;
	}* /

}*/

bool Vertex::isDomainBorder( Triangulation *vd)
{
	return this->position[0] < vd->xMin[0] + vd->boundaryThickness
	    && this->position[0] > vd->xMax[0] - vd->boundaryThickness
#if DIMENSIONS >= 2
		&& this->position[1] < vd->xMin[1] + vd->boundaryThickness
		&& this->position[1] > vd->xMax[1] - vd->boundaryThickness
#endif
#if DIMENSIONS >= 3
		&& this->position[2] < vd->xMin[2] + vd->boundaryThickness
		&& this->position[2] > vd->xMax[2] - vd->boundaryThickness
#endif
		;
}

int Triangulation::getCountConvexFaces()
{return this->countConvexFaces;}



bool PointInTriangle2D( double *p, double *a, double *b, double *c){
	// Compute vectors
	double v0[2] = {c[0]-a[0],c[1]-a[1]};
	double v1[2] = {b[0]-a[0],b[1]-a[1]};
	double v2[2] = {p[0]-a[0],p[1]-a[1]};

	// Compute dot products
	double dot00 = dot2d(v0, v0);
	double dot01 = dot2d(v0, v1);
	double dot02 = dot2d(v0, v2);
	double dot11 = dot2d(v1, v1);
	double dot12 = dot2d(v1, v2);

	// Compute barycentric coordinates
	double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return (u >= 0) && (v >= 0) && (u + v < 1);
}

bool Triangulation::PointInsideConvexHull(double *p) {
	if (countConvexFaces) {

		int countCollisions = 0;

		for (int f = 0; f < countConvexFaces; f++) {
			// Faces collide with Ray originating from p and parallel to z-Axis ?
			double *a = convexFaces[f][0]->position;
			double *b = convexFaces[f][1]->position;
			double *c = convexFaces[f][2]->position;

			if (PointInTriangle2D(p, a, b, c)) {
				//fprintf(stderr,"found collision!\n");
				//fprintf(stderr,"face: (%.3e, %.3e, %.3e) (%.3e, %.3e, %.3e) (%.3e, %.3e, %.3e) \n", a[0],a[1],a[2], b[0],b[1],b[2], c[0],c[1],c[2]);
				// z-Coordinate of Ray-Collision
				double ab[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
				double ac[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };

				//double j = (ac[1]*(a[0]+p[0]) - ac[0]*(a[1]+p[1])) / (ac[0]*ab[1] - ac[1]*ab[0]);
				double j = (ac[1] * (p[0] - a[0]) - ac[0] * (p[1] - a[1]))
						/ (ac[1] * ab[0] - ac[0] * ab[1]);
				double i;
				if (ac[0] != 0)
					i = (p[0] - a[0] - j * ab[0]) / ac[0];
				else
					i = (p[1] - a[1] - j * ab[1]) / ac[1];
				double z = a[2] + i * ac[2] + j * ab[2];

				if (z <= p[2])
					countCollisions++;
				//fprintf(stderr,"%lf <= %lf, %lf, %lf\n", z, p[2], (ac[1]*ab[0] - ac[0]*ab[1]), ac[0]);
			}
		}

		//fprintf(stderr,"%i collisions accepted!\n", countCollisions);

		if (countCollisions % 2 == 1)
			return true;
	}
	return false;
}


void Triangulation::CheckAndChangeOrientationOfConvexHullTriangles(){
	for( int f=0; f<countConvexFaces; f++){

		//calculate point inside the triangle
		double v[3];
		double v01[3];
		double v02[3];

		v01[0] = this->convexFaces[f][1]->position[0] - this->convexFaces[f][0]->position[0];
		v01[1] = this->convexFaces[f][1]->position[1] - this->convexFaces[f][0]->position[1];
		v01[2] = this->convexFaces[f][1]->position[2] - this->convexFaces[f][0]->position[2];

		v02[0] = this->convexFaces[f][2]->position[0] - this->convexFaces[f][0]->position[0];
		v02[1] = this->convexFaces[f][2]->position[1] - this->convexFaces[f][0]->position[1];
		v02[2] = this->convexFaces[f][2]->position[2] - this->convexFaces[f][0]->position[2];


		v[0] = this->convexFaces[f][0]->position[0] + 0.5 * v01[0];
		v[1] = this->convexFaces[f][0]->position[1] + 0.5 * v01[1];
		v[2] = this->convexFaces[f][0]->position[2] + 0.5 * v01[2];

		double vv2[3];

		vv2[0] = 0.5 * (this->convexFaces[f][1]->position[0] - v[0]);
		vv2[1] = 0.5 * (this->convexFaces[f][1]->position[2] - v[1]);
		vv2[2] = 0.5 * (this->convexFaces[f][1]->position[2] - v[2]);

		v[0] = v[0] + vv2[0];
		v[1] = v[1] + vv2[1];
		v[2] = v[2] + vv2[2];

		double normalvector[3];

		normalvector[0] = v01[1]*v02[2] - v01[2]*v02[1];
		normalvector[1] = v01[2]*v02[0] - v01[0]*v02[2];
		normalvector[2] = v01[0]*v02[1] - v01[1]*v02[0];

		double n = sqrt(normalvector[0]*normalvector[0]+normalvector[1]*normalvector[1]+normalvector[2]*normalvector[2]);

		normalvector[0] /= n;
		normalvector[1] /= n;
		normalvector[2] /= n;


		bool tmp = 1;
		double factor = 0.01;

		while( tmp ){

			bool p=0,m=0;

			double p0[3];
			double p1[3];

			p0[0] = factor*normalvector[0];
			p0[1] = factor*normalvector[1];
			p0[2] = factor*normalvector[2];

			p1[0] = factor*normalvector[0];
			p1[1] = factor*normalvector[1];
			p1[2] = factor*normalvector[2];

			p = this->PointInsideConvexHull(p0);
			m = this->PointInsideConvexHull(p1);

			if( p && !m ){

				Vertex *tmp =  this->convexFaces[f][2];
				this->convexFaces[f][2] = this->convexFaces[f][1];
				this->convexFaces[f][1] = tmp;

			}else{
				if( !p && m){
					tmp = 0;
				}else{
					if( (p && m) || (!p&& !m)){

						factor *= 0.1;

					}
				}

			}

		}

	}
}


