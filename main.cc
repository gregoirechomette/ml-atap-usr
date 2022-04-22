#include <iostream>
#include <vector>

#include "model.h"
#include "popGrid.h"


// Input structure
struct Input{
    double latitude = 48.8647;
    double longitude = 2.3490;
    double velocity = 10000;
    double angle = 45;
    double azimuth = 180;
};

// Output structure
struct Output{
    double affectedPopulation;
    double velocityDerivative;
    double angleDerivative;
    double azimuthDerivative;
};


// Main function call
int modelCall(Input input, Model model, PopGrid &popGrid, Output &output){

    // Structure the input vector
    std::vector<float> scenario = model.structInputVector(input.velocity, input.angle, input.azimuth);

    // Find the damage radius
    double damageRadius = std::max(model.evaluateOutput(scenario),0.0);

    // Find the number of people affected
    double affectedPop = 0.1 * popGrid.getAffectedPop(input.latitude, input.longitude, damageRadius);

    // The derivatives of the output w.r.t. the inputs
    std::vector<double> globalDerivatives = model.globalDerivatives(affectedPop, damageRadius);

    // Output derivatives of number of people w.r.t. inputs
    output.affectedPopulation = affectedPop;
    output.velocityDerivative = globalDerivatives[4];
    output.angleDerivative = globalDerivatives[5];
    output.azimuthDerivative = globalDerivatives[6];

    return 0;
}


int main() {

    // ************ Instantiations ************
    
    // Instantiate the neural network model
    const std::string folderName = "../models/BlastRad1/";
    Model model(folderName);

    // Instantiate the population grid vector
    const std::string popGridFile = "../pop-grids/popgrid-2020-2pt5arcmin.bin";
    PopGrid popGrid(popGridFile);


    // ************ Evaluations ************
    Input input = {}; Output output;
    for (int i=0; i<20; i++){

        // Update the input (if necessary)
        input.latitude += 0.05 * i;

        // Call the main function
        int functionReturn = modelCall(input, model, popGrid, output);

        // Print output
        std::cout << "The affected population is: " << output.affectedPopulation << std::endl;

    }

    return 0;
}