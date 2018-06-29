
#include "BiologyLink.h"

#include <cstdio>
#include <iostream>

BiologyLink::BiologyLink(int model)
{
    #pragma region 0: Init for ModelCellsSpherical Model

    if (model == 0)
    {
        time_scale   = 1;       // seconds
        length_scale = 1e-5;     // meter (= 10 �m) = cell diameter  (formerly: 2D: 10�m  3D: 14 �m)
        energy_scale = 1e-16;    // J
        youngModulus_bioCells = 450;  // Pa
        poisson_ratio_bioCells = 0.4;
        youngModulus_bioSinusoids = 600;  // Pa
        poisson_ratio_bioSinusoids = 0.4;
        single_bond_energy_bio = 1e-19; // J from DrasdoHoehme 2005
        diffusion_constant_cells_bio = 1e-16; // 1e-16; // m^2/s // 1e-16 m^2/s = 1e-12 cm^2/s
        cell_gamma_perpendicular_bio = 1e8; // Ibrahim value: 1e8 // Old default: 3e7; // Ns/m^3 from Galle et al
        cell_gamma_parallel_bio = 1e8; // gamma_parallel = gamma_perpendicular => isotropic friction with value same value as the variables
        ecm_gamma_bio = 1e8; // New ibrahim value // Old default: 3e7; // Ns/m^3 from Galle et al
        single_bond_energy_bio = 1e-19; //
        adhesion_to_cells_bio = 1e13; // m^-2  // Default: 1e15/m� = 1000/�m� // cell cell
        cycletime_bio = 18 * 3600; // s // 18 h = 64800 s
        cycletime_stddev_bio = 2 * 3600; // s // 4 h
        surface_tension_bio = 1e-4;
    }

    #pragma endregion

    #pragma region 1: Init for Liver3DCellsSpherical Model
	/*
    else if (model == 1)
    {
        time_scale   = 1;    // seconds
        length_scale = 2.33e-5;  // meter (= 23.3 �m) = hepatocyte diameter
        energy_scale = 1e-16; // J
        youngModulus_bio = 450; // Pa
		poisson_ratio_bio = 0.4;
        single_bond_energy_bio = 1e-19; // J from DrasdoHoehme 2005
        //timestep_bio = 0.10; // 1; // s (Just a starting point - is automatically changed by the simulation system) unless fixed
        diffusion_constant_cells_bio[0] = 1e-16; // m^2/s // 1e-16 m^2/s = 1e-12 cm^2/s
        diffusion_constant_cells_bio[1] = 1e-16; // m^2/s // 1e-16 m^2/s = 1e-12 cm^2/s
        cell_gamma_bio = 1e8; // New ibrahim value // Old default: 3e7; // Ns/m^3 from Galle et al
        ecm_gamma_bio = 1e8; // New ibrahim value // Old default: 3e7; // Ns/m^3 from Galle et al
        //constant_friction_option = false;
        single_bond_energy_bio = 1e-19; //
        default_adhesion_to_celltype_bio[0] = 0; //1e15; // m^-2  // 1e15/m� = 1000/�m� // cell cell
        default_adhesion_to_celltype_bio[1] = 0; //1e15; // m^-2  // 1e15/m� = 1000/�m� // tissue tissue
        //    cell_cell_adhesion_density_bio[2] = 0; // m^-2  // 1e15/m� = 1000/�m� // cell tissue
        cell_sinusoid_adhesion_density_bio = 0; // m^-2  // 1e15/m� = 1000/�m�
        cycletime_bio = 24 * 3600; // s // 18 h = 64800 s
        cycletime_stddev_bio = 0; //4 * 3600; // s // 4 h
        surface_tension_bio = 1e-4;
    }
	*/
    #pragma endregion

    #pragma region Other: Error: Unknown model

    else
    {
        // Prelim. Output error
    }

    #pragma endregion
}

#pragma region Generic methods to rescale and convert biology <-> dimensionless internal units

double BiologyLink::getTimeInDays(double internal)
{
    return (internal * time_scale) / (24.f * 60.f * 60.f); // returns days
}

double BiologyLink::getTimeInSeconds(double internal)
{
    return (internal * time_scale); // returns seconds
}

double BiologyLink::getLengthInMicrometers(double length_internal)
{
    return length_internal * (length_scale * 1e6);
}

double BiologyLink::getPressureInPascal(double press_intern)
{
    // N / m�
    return press_intern * ((energy_scale / length_scale) / (length_scale * length_scale));
}

// J = [Nm]
double BiologyLink::getForceInNanoNewton(double force_intern)
{
    return force_intern * (energy_scale / length_scale) * 1e9;
}

double BiologyLink::getYoungModuleInPascal(double young_intern)
{
    return young_intern * (energy_scale / (length_scale * length_scale * length_scale));
}

double BiologyLink::getYoungModuleDimensionless(double young_in_pascal)
{
    return young_in_pascal / (energy_scale / (length_scale * length_scale * length_scale));
}

double BiologyLink::getSurfaceTensionDimensionless(double surftensBiounit)
{
    return surftensBiounit * ((length_scale * length_scale) / energy_scale);
}

double BiologyLink::scaleBiologyToInternal(double x, ScaleMode how)
{
    switch (how)
    {
    case ScaleTime:   // time
        x /= time_scale;
        break;
    case ScaleLength:   // length
        x /= length_scale;
        break;
    case ScaleVelocity:   // velocity = m/s
        x /= (length_scale/time_scale);
        break;
    case ScaleYoungModule:   // youngs module N/m^2 = J/m^3
        x /= (energy_scale/(length_scale*length_scale*length_scale));
        break;
    case ScaleDiffusivity:   // diffusivity m^2/s
        x /= (length_scale*length_scale/time_scale);
        break;
    case ScaleEnergy:   // energy J
        x /= (energy_scale);
        break;
    case ScaleAdhesionDensity:   // adhesion density = 1/m^2
        x *= (length_scale*length_scale);
        break;
    case ScaleViscosity:   // viscosity = Ns/m^2 = Js/m^3
        x /= (energy_scale*time_scale/(length_scale*length_scale*length_scale));
        break;
    case ScaleFrictionCoefficient:   // friction coefficients = Ns/m^3 = Js/m^4
        x /= (energy_scale*time_scale/(length_scale*length_scale*length_scale*length_scale));
        break;
    case ScaleForce:
        x /= (energy_scale/length_scale);
        break;
    case ScalePressure:
        x /= (energy_scale/(length_scale*length_scale*length_scale));
        break;
    default:
        x = 0;
        break;
        //printf("ERROR: Unknown conversion mode (%c) in biology.\n",how);
    }

    return x;
}


#pragma endregion
