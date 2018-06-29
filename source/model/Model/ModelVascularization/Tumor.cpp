/*
 * Tumor.cpp
 *
 *  Created on: 17.01.2012
 *      Author: jagiella
 */

#include <math.h>
#include "Tumor.hpp"
#include "../../../tools/new_triangulation/IO.hpp"


const char Tumor::NONE = 0;
const char Tumor::SPHERICAL = 1;

Tumor::Tumor()
{
	topology = NONE;
	radius = 0;
	necroticCoreRadius = 0;
	countParallelVessels = 1;
	permeability = 0;
	center[0]=center[1]=center[2]=0;
}

Tumor::~Tumor(){

}

void Tumor::setTopology(char topology)
{
	this->topology = topology;
}

void Tumor::setRadius(double radius)
{
	this->radius = radius;
}

void Tumor::setRadiusNecroticCore(double radius)
{
	this->necroticCoreRadius = radius;
}

void Tumor::setCenter(double* center)
{
	for( int i=0; i<3; i++)
		this->center[i] = center[i];
}

void Tumor::setX(double value)
{
	this->center[0] = value;
}

void Tumor::setY(double value)
{
	this->center[1] = value;
}

void Tumor::setZ(double value)
{
	this->center[2] = value;
}

bool Tumor::isTumor( double* pos)
{
	double dist = 0;
	for( int i=0; i<3; i++)
		dist += pow( center[i] - pos[i], 2);
	dist = sqrt(dist);

	return ( dist<=radius );
}

bool Tumor::isNecrotic( double* pos)
{
	// TEMP
	/*if( sqrt( pow( 50 - pos[0], 2) + pow( 35 - pos[1], 2)) <= necroticCoreRadius) return true;
	if( sqrt( pow( 45 - pos[0], 2) + pow( 25 - pos[1], 2)) <= necroticCoreRadius) return true;
	if( sqrt( pow( 55 - pos[0], 2) + pow( 25 - pos[1], 2)) <= necroticCoreRadius) return true;
	return false;*/
	// TEMP

	double dist = 0;
	for( int i=0; i<3; i++)
		dist += pow( center[i] - pos[i], 2);
	dist = sqrt(dist);

	return ( dist<=necroticCoreRadius );
}

bool Tumor::isTumor( double x, double y, double z)
{
	double dist = pow( center[0] - x, 2) + pow( center[1] - y, 2) + pow( center[2] - z, 2);
	dist = sqrt(dist);

	return ( dist<=radius );
}

void Tumor::printToEPS(const char *filename)
{
	std::fstream fs;
	fs.open( filename, std::fstream::out | std::fstream::app);

	char color[512];

	// TUMOR
	EPS::PScolor( color, 0,0,0);
	EPS::PSwriteCircle( &fs, center[0], center[1], radius, 1, color);

	// NECROTIC CORE
	EPS::PScolor( color, 1,0,0);
	EPS::PSwriteCircle( &fs, center[0], center[1], this->necroticCoreRadius, 1, color);

	// TEMP
	/*EpsIO::PSwriteCircle( &fs, 50, 35, this->necroticCoreRadius, 1, color);
	EpsIO::PSwriteCircle( &fs, 45, 25, this->necroticCoreRadius, 1, color);
	EpsIO::PSwriteCircle( &fs, 55, 25, this->necroticCoreRadius, 1, color);*/
	// TEMP

	fs.close();
}

double Tumor::getCenter( int dim){
	return center[dim];
}

double Tumor::getRadius()
{
	return radius;
}

double Tumor::getNecroticCoreRadius()
{
	return necroticCoreRadius;
}

void Tumor::setPermeability(double value)
{
	permeability = value;
}
double Tumor::getPermeability()
{
	return permeability;
}

void Tumor::setParallelVessels( double value)
{
	this->countParallelVessels = value;
}

double Tumor::getParallelVessels()
{
	return this->countParallelVessels;
}
