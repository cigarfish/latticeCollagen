
#include <cstdio>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <cmath>

#include "H5Cpp.h"

#include "Random.h"
#include "../output/OutputText.h" // ?

Random::Random()
{
    // Init marsaglia random number generator + others
//    Init();

    // Init persitent variable for GRan
    GRan_use_last = 0;
}

void Random::Init()
{
    InitRandomMarsaglia(22,12,67,54);

    // Init another random number generator
    srand( (unsigned) time(NULL)) ;
}

double Random::GetRandomUniform11(void)
{
    return (2.f * GetRandomUniform01()) - 1.f;
}

// Initialization of the Marsaglia Random number Generator according to
//   Marsaglia, G. and A. Zaman and W.W. Tsang. 1990. Toward a Universal
//   Random Number Generator. /Statistics and Probability Letters./ 8: 35-39.
//   North-Holland.
void Random::InitRandomMarsaglia( int nA1, int nA2, int nA3, int nB1 )
{
    // Initializes the MARSAGLIA pseudo random number generator
    std::cout << "Init: Marsaglia random number generator.\n";

    int mA1, mA2, mA3, mANEW, mB1, mHELP ;
    int i1, i2 ;
    float varS, varT ;

    mA1 = nA1 ;
    mA2 = nA2 ;
    mA3 = nA3 ;
    mB1 = nB1 ;
    MARSi  = 96 ;
    MARSj  = 32 ;

    for ( i1 = 0 ; i1 < 97 ; i1++ )
    {
        varS = 0.0 ;
        varT = 0.5 ;
        for ( i2 = 0 ; i2 < 24 ; i2++ )
        {
            mANEW = ((( mA1 * mA2 ) % 179 )*mA3 ) % 179 ;
            mA1 = mA2 ;
            mA2 = mA3 ;
            mA3 = mANEW ;
            mB1 = ( 53*mB1 + 1 ) % 169 ;
            mHELP = ( mB1 * mANEW ) % 64 ;
            if ( mHELP > 31 )
                varS += varT ;
            varT *= 0.5 ;
        }
        MARSarray[ i1 ] = varS ;
    }

    MARSc  =   362436.0 / 16777216.0 ;
    MARScd =  7654321.0 / 16777216.0 ;
    MARScm = 16777213.0 / 16777216.0 ;

    return;
}

long int Random::GetRandomUniform1x(long int x)
{
    if (x<1)
    {
        printf("CellsysRandom::RAN0x ERROR: Parameter must be greater than zero.\n");
        return 1;
    }
    double r = GetRandomUniform01();
    if (r==0) return 1;
    else
    {
        return (long int)ceil(x * r);
    }
}

// Generating a uniform distribution of pseudo-random numbers with the Marsaglia
// Random number Generator according to
//   Marsaglia, G. and A. Zaman and W.W. Tsang. 1990. Toward a Universal
//   Random Number Generator. /Statistics and Probability Letters./ 8: 35-39.
//   North-Holland.
double Random::GetRandomUniform01(void)
{
    // Generates a pseudo random number 0 .. +1 following the proposal of MARSAGLIA
    double ranMARS ;

    ranMARS = MARSarray[ MARSi ] - MARSarray[ MARSj ] ;
    if ( ranMARS < 0.0 )
        ranMARS += 1.0 ;

    MARSarray[ MARSi ] = ranMARS ;

    MARSi-- ;
    if ( MARSi < 0 )
        MARSi = 96 ;

    MARSj-- ;
    if ( MARSj < 0 )
        MARSj = 96 ;

    MARSc -= MARScd ;
    if ( MARSc < 0.0 )
        MARSc += MARScm ;

    ranMARS -= MARSc ;
    if ( ranMARS < 0.0 )
        ranMARS += 1.0 ;

    return ranMARS ;
}

double Random::GetRandomGauss(double mean, double standard_dev)
{
    // Normal random variate generator
    // Mean m, standard deviation s
    // http://www.taygeta.com/random/gaussian.html

    if (GRan_use_last)              // use value from previous call
    {
        GRan_y1 = GRan_y2;
        GRan_use_last = 0;
    }
    else
    {
        do
        {
            GRan_x1 = 2.0 * GetRandomUniform01() - 1.0;
            GRan_x2 = 2.0 * GetRandomUniform01() - 1.0;
            GRan_w = GRan_x1 * GRan_x1 + GRan_x2 * GRan_x2;
        }
        while ( GRan_w >= 1.0 );

        GRan_w = sqrt( (-2.0 * log( GRan_w ) ) / GRan_w );
        GRan_y1 = GRan_x1 * GRan_w;
        GRan_y2 = GRan_x2 * GRan_w;
        GRan_use_last = 1;
    }

    return( mean + GRan_y1 * standard_dev);
}


double Random::GetRandomLogNormal (double mean, double cv)
{
    double mu = log(mean)-0.5*log(1+cv*cv ) ;
    double sig = sqrt( log(cv*cv+1) );
    return exp ( GetRandomGauss (mu,sig) ) ;
}


int Random::GetRandomUniformInt(int min, int max)
{
    return (int)(GetRandomUniform01() * (double)(max-min)) + min;
}

double Random::GetRandomUniformDouble(double min, double max)
{
    return (GetRandomUniform01() * (max-min)) + min;
}

void Random::GetRandomUnitVector(double *ex, double *ey, double *ez)
{
    // Method after Marsaglia72 (http://mathworld.wolfram.com/SpherePointPicking.html)
    float a,b;
    a = GetRandomUniform11();
    b = GetRandomUniform11();
    while (((a*a) + (b*b)) >= 1.)
    {
        a = GetRandomUniform11();
        b = GetRandomUniform11();
    }
    *ex = 2 * a *sqrt(1-(a*a)-(b*b));
    *ey = 2*b*sqrt(1-(a*a)-(b*b));
    *ez = 1-2*((a*a)+(b*b));
}

void Random::GetRandomUnitVector(double *ex, double *ey)
{
    // Method after Marsaglia72 (http://mathworld.wolfram.com/CirclePointPicking.html)

    float a = GetRandomUniform11();
    float b = GetRandomUniform11();
    while (((a*a) + (b*b))>=1.)
    {
        a = GetRandomUniform11();
        b = GetRandomUniform11();
    }
    *ex = ((a*a) - (b*b)) / ((a*a) + (b*b));
    *ey = (2 * a * b)     / ((a*a) + (b*b));
}


void Random::HDF5DataFormat( H5::CompType & rngType )
{
    rngType = H5::CompType( sizeof( Random ) );

    hsize_t arrayLength = 97;
    H5::ArrayType marsArrayType( H5::PredType::NATIVE_DOUBLE, 1, &arrayLength );

    rngType.insertMember( "MARSarray", HOFFSET( Random, MARSarray ), marsArrayType );

    rngType.insertMember( "MARSi", HOFFSET( Random, MARSi), H5::PredType::NATIVE_INT );
    rngType.insertMember( "MARSj", HOFFSET( Random, MARSj), H5::PredType::NATIVE_INT );

    rngType.insertMember( "GRan_x1", HOFFSET( Random, GRan_x1), H5::PredType::NATIVE_DOUBLE );
    rngType.insertMember( "GRan_x2", HOFFSET( Random, GRan_x2), H5::PredType::NATIVE_DOUBLE );
    rngType.insertMember( "GRan_w",  HOFFSET( Random, GRan_w),  H5::PredType::NATIVE_DOUBLE );
    rngType.insertMember( "GRan_y1", HOFFSET( Random, GRan_y1), H5::PredType::NATIVE_DOUBLE );
    rngType.insertMember( "GRan_y2", HOFFSET( Random, GRan_y2), H5::PredType::NATIVE_DOUBLE );

    rngType.insertMember( "GRan_use_last", HOFFSET( Random, GRan_use_last), H5::PredType::NATIVE_INT );
}
