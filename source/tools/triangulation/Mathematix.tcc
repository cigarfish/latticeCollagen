/*
 * Mathematix.tcc
 *
 *  Created on: Mar 1, 2011
 *      Author: jagiella
 */

#include <stdio.h>

template <class T> void printVector( T *v, int N, const char* name, const char* format)
{
	fprintf( stderr, "%s = [\n", name);
	for( int i=0; i<N; i++)
	{
		fprintf( stderr, format, v[i]);
	}
	fprintf( stderr, "\n]\n");
}
