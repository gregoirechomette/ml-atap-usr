#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <ctime>
#include <time.h>
#include <algorithm>
#include <vector>

#include "model.h"
#include "popGrid.h"

// Constants
const double pi = 3.14159;
const double deg2rad = 0.01745329;
const double earthRadius = 6.371e6;
const double worldPop = 7.93 * pow(10,9);


class Input {

    public:

        /**
         * @brief Constructor of a new Input object
         * 
         * @param fileName path and name of the .dat file to be read
         */

        Input(const std::string fileName): _fileName(fileName) {
            readFile();
        }

    private:

        /**
         * @brief Method to read and store the dat file containing scenario parameters
         * 
         */
        
        void readFile() {

            // Prepare the vector   
            std::vector<float> scenarioParameters(9);

            // Try to open filename
            std::ifstream source(_fileName.data());
            if (source.fail()) {
            std::cout << "Error: Cannot open parameter input file: " << _fileName.data() 
                        << std::endl;
            return;
            }
            // actually read file
            source >> scenarioParameters[0]; 
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> scenarioParameters[1];    
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> scenarioParameters[2];
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> scenarioParameters[3];
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> scenarioParameters[4]; 
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> scenarioParameters[5];
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> scenarioParameters[6];
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> scenarioParameters[7]; 
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> scenarioParameters[8]; 
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> _latitude; 
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> _longitude; 
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source.close();

            _scenarioParameters = scenarioParameters;

            // all done
            return;
        }

    public:
        // Declare attributes
        double _latitude, _longitude;
        const std::string _fileName;
        std::vector<float> _scenarioParameters;
};



int main() {

    // Material parameters 
    float diameter = 70.0;
    float density = 3500;
    float strength = 100000;
    float alpha = 0.2;
    float lumEff = 0.003;
    float ablation = 0.000000001;

    // Coordinates
    double latitude = 48.8647;
    double longitude = 2.3490;

    // Trajectory
    float velocity = 10000;
    float angle = 45;
    float azimuth = 180;

    std::vector<float> dataInp = {diameter, density, strength, alpha, velocity, angle, azimuth, lumEff, ablation};
    

    // Instantiate the input objects with the scenario properties
    const std::string propertiesFile = "../input.dat";
    Input data(propertiesFile);

    // Instantiate the population grid vector
    const std::string popGridFile = "../pop-grids/popgrid-2020-2pt5arcmin.bin";
    PopGrid popGrid(popGridFile);

    // Instantiate the neural network model
    const std::string folderName = "../models/BlastRad1/";
    Model model(folderName);

    // Find the damage radius
    double damageRadius = std::max(model.evaluateOutput(dataInp),0.0);
    std::cout << "The radius of damage is: " << damageRadius <<  " m" << std::endl;

    // Find the number of people affected
    double affectedPop = 0.1 * popGrid.getAffectedPop(latitude, longitude, damageRadius);
    std::cout << "The number of people affected is: " << affectedPop << std::endl;

    // The derivatives of the output w.r.t. the inputs
    std::vector<double> nnDerivatives = model.backPropagation();
    std::vector<double> globalDerivatives = model.globalDerivatives(nnDerivatives, affectedPop, damageRadius);

    // Print some derivatives
    std::cout << "The derivative with respect to the velocity is: " << globalDerivatives[4] << std::endl;
    std::cout << "The derivative with respect to the incidence angle is: " << globalDerivatives[5] << std::endl;
    std::cout << "The derivative with respect to the azimuth is: " << globalDerivatives[6] << std::endl;

    return 0;
}