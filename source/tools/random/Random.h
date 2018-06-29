#ifndef RANDOM_H
#define RANDOM_H

namespace H5 {class CompType;}


class Random
{
private:

    #pragma region Marsaglia random numbers (Uniform and Gauss)

    // Data for marsaglia initialisation // Data for marsaglia initialisation
    double MARSarray[ 97 ], MARSc, MARScd, MARScm ;
    int    MARSi, MARSj; // Data for marsaglia initialisation
    double GRan_x1, GRan_x2, GRan_w, GRan_y1, GRan_y2; // persistent variables for GRan
    int GRan_use_last; // persistent variables for GRan

    // Init the random number generator
    void InitRandomMarsaglia( int nA1, int nA2, int nA3, int nB1 );

    #pragma endregion

public:


    // Returns a uniformly distributed random number (0,1]
    double GetRandomUniform01(void);
    // Returns a random number [-1,1]
    double GetRandomUniform11(void);
    // Returns a (long int) random number (1,x]
    long int GetRandomUniform1x(long int x);

    // Returns a (unsigned long) random number from [0,x)
    unsigned long int GetRandomUniform0x(unsigned long int x)
    {return (unsigned long int)GetRandomUniform1x(x) -1;};

    // Generates random number that are guassian distributed
    double GetRandomGauss(double mean, double standard_dev);
    // Generates a log-normally distributed random number
    double GetRandomLogNormal ( double mean , double cv ) ;

	//! Returns a random number between min and max (uniformly distributed)
	int GetRandomUniformInt(int min, int max);
	double GetRandomUniformDouble(double min, double max);

	//! Returns a unit vector in 2D with correctly distributed circle-point picking
	void GetRandomUnitVector(double *ex, double *ey); 

	//! Returns a unit vector in 3D with correctly distributed circle-point picking
    void GetRandomUnitVector(double *ex, double *ey, double *ez); 

    static void HDF5DataFormat( H5::CompType & type );

    Random();

    // Init all generators
    void Init();

};

#endif
