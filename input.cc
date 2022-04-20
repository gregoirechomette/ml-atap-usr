#include <iostream>
#include <vector>

#include "input.h"


Input::Input() = default;

std::vector<float> Input::prepareInput(double velocity, double incidenceAngle, double azimuth){
    
    // Baseline material parameters 
    float diameter = 70.0;
    float density = 3500;
    float strength = 100000;
    float alpha = 0.2;
    float lumEff = 0.003;
    float ablation = 0.000000001;

    // Conversion of parameters from double to float
    float velocity_f = (float) velocity;
    float incidenceAngle_f = (float) incidenceAngle;
    float azimuth_f = (float) azimuth;

    // Group and return
    std::vector<float> inputNN = {diameter, density, strength, alpha, velocity_f, incidenceAngle_f, azimuth_f, lumEff, ablation};
    return inputNN;
}