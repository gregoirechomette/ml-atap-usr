#include <iostream>
#include <vector>

#include "model.h"
#include "popGrid.h"


// Input structure 
struct Input{
    std::vector<double> positionVector;     // [m]
    std::vector<double> velocityVector;     // [m/s]
};

// Output structure
struct Output{
    double affectedPopulation;              // number of people N [-]
    double velocityDerivative;              // dN/dvelocity [s/m]
    double angleDerivative;                 // dN/dangle [1/degrees]
    double azimuthDerivative;               // dN/dazimuth [1/degrees]
};



// Main function call
int neuralNetPredictions(Input input, Model model, PopGrid &popGrid, Output &output){

    // Structure the input vector
    std::vector<double> scenario = model.structInputVector(input.positionVector, input.velocityVector);

    // Find the damage radius
    double damageRadius = std::max(model.evaluateOutput(scenario),0.0);

    // Find the number of people affected
    double affectedPop = 0.1 * popGrid.getAffectedPop(input.positionVector, damageRadius);

    // Find the derivatives of the output w.r.t. the inputs
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



    // ************ Evaluations during optimization loop ************
    Input input; Output output;
    for (int i=0; i<3; i++){

        // Generate an input (position and velocity vectors)
        input.positionVector = { 3654.32, 4381.75108, 2850.44503 };    
        input.velocityVector = { -13252.019, -24529.8, 1302.29 };  

        // Call the main function
        int functionReturn = neuralNetPredictions(input, model, popGrid, output);

        // Print output for demonstration purposes
        std::cout << "The affected population is: " << output.affectedPopulation << std::endl;
        std::cout << "The derivative w.r.t the velocity is: " << output.velocityDerivative << std::endl;

    }

    return 0;
}