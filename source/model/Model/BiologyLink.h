#ifndef BIOLOGY_LINK_H
#define BIOLOGY_LINK_H

class BiologyLink
{
public:
    enum ScaleMode {
        ScaleTime,                // time
        ScaleLength,              // length
        ScaleVelocity,            // velocity = m/s
        ScaleYoungModule,         // youngs module N/m^2 = J/m^3
        ScaleDiffusivity,         // diffusivity m^2/s
        ScaleEnergy,              // energy J
        ScaleAdhesionDensity,     // adhesion density = 1/m^2
        ScaleViscosity,           // viscosity = Ns/m^2 = Js/m^3
        ScaleFrictionCoefficient, // friction coefficients = Ns/m^3 = Js/m^4
        ScaleForce,               // Newton = J/m
        ScalePressure,            // pressure Pa
    };

    #pragma region Generic methods to rescale and convert biology <-> dimensionless internal units

    // Convert biological parameters to internal (dimensionless) parameters
    double scaleBiologyToInternal(double x, ScaleMode how);
    double getYoungModuleDimensionless(double young_in_pascal);
    double getSurfaceTensionDimensionless(double surftensBiounit);

    // Convert internal parameters to biological
    double getLengthInMicrometers(double length_internal);
    double getTimeInDays(double internal_time);
    double getTimeInSeconds(double internal_time);
    double getPressureInPascal(double press_intern);
    double getYoungModuleInPascal(double young_intern);
    double getForceInNanoNewton(double force_intern);

    #pragma endregion


    #pragma region constructor

    BiologyLink(int model);

    #pragma endregionngModuleInPasca


    #pragma region Biological (default) values with biological units

    #pragma region Primary scales

    // Primary scales used to make parameters dimensionless (model)
    double time_scale;
    double energy_scale;
    double length_scale; // model length scale, also cell diameter which is accessible via GUI (different in monolayer and spheroid)

    #pragma endregion

    #pragma region Secondary biological parameters

    double youngModulus_bioCells;
    double poisson_ratio_bioCells;
    double youngModulus_bioSinusoids;
    double poisson_ratio_bioSinusoids;
    double diffusion_constant_cells_bio; // cell diffusion constant
    double cell_gamma_perpendicular_bio; // friction coefficient cells vs. cells for shear friction component
    double cell_gamma_parallel_bio; // friction coefficient cells vs. cells for collisional friction component
    double ecm_gamma_bio;  // friction coefficient cells vs. ecm
    double single_bond_energy_bio; // single bond energy
    double adhesion_to_cells_bio; // density of adhesion molecules to other cells
    double cycletime_bio; // cell cycle time in seconds
    double cycletime_stddev_bio; // max. difference to real cycletime (+-) cell cycle time in seconds
    double surface_tension_bio;
    double quiescentMinPressure_bio;
    double quiescentMaxPressure_bio;

    #pragma endregion
};

#endif
