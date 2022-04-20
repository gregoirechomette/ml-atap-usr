#include <iostream>
#include <vector>

#include "input.h"
#include "model.h"
#include "popGrid.h"


int main() {

    /**
     *  Instantiation of model and population file (only once)
     */
    
    // Instantiate the neural network model
    const std::string folderName = "../models/BlastRad1/";
    Model model(folderName);

    // Instantiate the population grid vector
    const std::string popGridFile = "../pop-grids/popgrid-2020-2pt5arcmin.bin";
    PopGrid popGrid(popGridFile);

    // Instantiate the input class for the neural network
    Input inputNN;
    

    /**
     * Optimization loop calling the NN model as many times as necessary
     */

    for (int i=0; i<100; i++){

        // Position and motion generated by the mission design team
        double latitude = 48.8647;
        double longitude = 2.3490;
        double velocity = 10000;
        double angle = 45;
        double azimuth = 180;

        //  Preparation of NN input
        std::vector<float> scenario = inputNN.prepareInput(velocity, angle, azimuth);

        // Find the damage radius
        double damageRadius = std::max(model.evaluateOutput(scenario),0.0);

        // Find the number of people affected
        double affectedPop = 0.1 * popGrid.getAffectedPop(latitude, longitude, damageRadius);

        // The derivatives of the output w.r.t. the inputs
        std::vector<double> nnDerivatives = model.backPropagation();
        std::vector<double> globalDerivatives = model.globalDerivatives(nnDerivatives, affectedPop, damageRadius);

        // Output derivatives of number of people w.r.t. inputs
        double velocityDerivative = globalDerivatives[4];
        double angleDerivative = globalDerivatives[5];
        double azimuthDerivative = globalDerivatives[6];

    }

    return 0;
}