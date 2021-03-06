/*
 * Mathematix.tcc
 *
 *  Created on: Mar 1, 2011
 *      Author: jagiella
 */

#include <stdio.h>

#include "math.h"

template <class T> void printVector( T *v, int N, const char* name, const char* format)
{
	fprintf( stderr, "%s = [\n", name);
	for( int i=0; i<N; i++)
	{
		fprintf( stderr, format, v[i]);
	}
	fprintf( stderr, "\n]\n");
}

template <class T> T vectorNorm( T *v, int N, int exponent)
{
	T norm = 0;
	for( int i=0; i<N; i++)
		norm += pow( v[i], exponent);
	return pow( norm, 1./exponent);
}

template <class T> T vectorNormInfinity( T *v, int N)
{
	T norm = 0;
	for( int i=0; i<N; i++)
		norm = MAX( fabs(v[i]), norm);
	return norm;
}
