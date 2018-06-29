/*
 * Tumor.hpp
 *
 *  Created on: 17.01.2012
 *      Author: jagiella
 */

#ifndef TUMOR_HPP_
#define TUMOR_HPP_


class Tumor{
private:
	char topology;
	double center[3];
	double radius;
	double necroticCoreRadius;
	double permeability;
	double countParallelVessels;


public:
	// Constants
	static const char NONE;
	static const char SPHERICAL;

	// Constructor / Destructor
	Tumor();
	~Tumor();

	// Methods
	void setTopology(char);
	void setRadius(double);
	void setRadiusNecroticCore(double);
	void setCenter(double*);
	void setX(double value);
	void setY(double value);
	void setZ(double value);
	void setPermeability(double);
	void setParallelVessels(double);

	bool isTumor( double*);
	bool isTumor( double x, double y, double z);
	bool isNecrotic( double* pos);

	double getCenter( int dim);
	double getRadius();
	double getNecroticCoreRadius();
	double getPermeability();
	double getParallelVessels();

	void printToEPS(const char *filename);
};

#endif /* TUMOR_HPP_ */
