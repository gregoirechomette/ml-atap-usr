#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <stdio.h>

#include "./lib/cppflow/include/cppflow/ops.h"
#include "./lib/cppflow/include/cppflow/model.h"

// Constants
const double pi = 3.14159;
double const deg2rad = 0.01745329;
const double earthRadius = 6.371e6;
const double worldPop = 7.93 * pow(10,9);


class Input {

    public:

        Input(const std::string fileName): _fileName(fileName) {
            _readFile();
        }

    private:

        // Method to read the dat file containing scenario parameters
        void _readFile() {
            // Try to open filename
            std::ifstream source(_fileName.data());
            if (source.fail()) {
            std::cout << "Error: Cannot open parameter input file: " << _fileName.data() 
                        << std::endl;
            return;
            }
            // actually read file
            source >> _diameter; 
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> _density;    
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> _strength;
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> _alphaCoeff;
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> _velocity; 
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> _angleInc;
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> _azimuth;
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> _luminousEff; 
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source >> _ablation; 
            source.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            source.close();

            // Create a vector to contain all the parameters
            std::vector<float> scenarioParameters;
            scenarioParameters = {_diameter, _density, _strength, _alphaCoeff, _velocity, _angleInc, _azimuth, _luminousEff, _ablation};
            _scenarioParameters = scenarioParameters;

            // all done
            return;
        }

    public:
        // Declare attributes
        const std::string _fileName;
        std::vector<float> _scenarioParameters;
        float _diameter, _density, _strength, _alphaCoeff, _velocity, _angleInc, _azimuth, _luminousEff, _ablation;
};


class PopGrid {

    public:

        PopGrid(const std::string fileName): _fileName(fileName) {
                readFile();
        }

        /**
         * @brief Method to read the binary file and store the information about densities
         */
        void readFile(){

            std::ifstream gridFile(_fileName, std::ios::binary);

            // Lecture of the file - saving of properties and data
            if (gridFile.is_open()){
                gridFile.read((char *)&nRows, sizeof(nRows));
                gridFile.read((char *)&nCol, sizeof(nCol));
                gridFile.read((char *)&xlCorner, sizeof(xlCorner));
                gridFile.read((char *)&ylCorner, sizeof(ylCorner));
                gridFile.read((char *)&arcLengthDeg, sizeof(arcLengthDeg));
                gridFile.read((char *)&noDat, sizeof(noDat));
                yuCorner = ylCorner + nRows*arcLengthDeg;
                vSize = nRows*nCol;
                maxcellsidekm = earthRadius*arcLengthDeg*deg2rad;
                maxcelldiagkm = maxcellsidekm*sqrt(2);

                popCount.reserve(vSize);
                for (int i=0; i< vSize; i++) {
                    
                    gridFile.read((char *)&dum,sizeof(dum));
                    popCount[i] = dum;
                    
                }
                gridFile.close();
            }
            else{
                std::cout << "Population file read failed!" << std::endl;
            }
            return;
        }

        /**
         * @brief Method to find the density given the latitude and longitude
         * 
         * @param[in] latitude in degrees
         * @param[in] longtitude in degrees
         * @param[out] density in #people/km^2
         */
        double getCellPop(double latitude, double longitude){

            // Retrive the indices on the table
            int i = int((yuCorner - latitude)/arcLengthDeg);
            int j = int((longitude - xlCorner)/arcLengthDeg);

            // Verify that the longitude is in the bounds, otherwise correct it
            if (j<0){
                j = nCol + j;
            }else if(j >= nCol){
                j = j - nCol;
            }

            // Find the density
            if (popCount[i*nCol + j] < 0){
                return 0;
            } else {
                return popCount[i * nCol + j];
            }
        }

        // Atributes
        const std::string _fileName;
        int nRows, nCol, xlCorner, ylCorner, yuCorner, noDat, vSize;
        double arcLengthDeg, maxcellsidekm, maxcelldiagkm, dum;
        std::vector<double> popCount;

};


class Model {

    public:

        Model(const std::string folderName): 
        _folderName(folderName) {
            _readingScalingParameters();
        }

        // Method to read csv file containing scaling parameters
        void _readingScalingParameters(){

            // Open the file and declare reading variables
            std::string scalingFileName = _folderName + "/Scaling_parameters.csv";
            std::fstream file (scalingFileName, std::ios::in);
            std::vector<std::string> row;
            std::string line, word;

            // Read the file and store data
            if(file.is_open()){
                while(getline(file, line)){
                    row.clear();
                    std::stringstream str(line);
                    while(getline(str, word, ',')){
                        row.push_back(word);
                    }
                    _scalingParameters.push_back(row);
                }
            }
            else{
                std::cout << "Error: Cannot open parameter scaling parameter file "  << std::endl;
            }
            return;
        }

        std::vector<float> _normalizeInputs(Input data){
            std::vector<float> normalizedInputs (9);
            for (int i=0; i<9; i++){
                float mean = std::stof(_scalingParameters[1][i+1]);
                float std = std::stof(_scalingParameters[2][i+1]);
                normalizedInputs[i] = (data._scenarioParameters[i] - mean)/ (std);
            }
            return normalizedInputs;
        }

        float _rescaleOutput(float damageNormalized){
            float mean = std::stof(_scalingParameters[1][10]);
            float std = std::stof(_scalingParameters[2][10]);
            return (damageNormalized * std) + mean;
        }


    public:

        // Declare attributes
        const std::string _folderName;
        std::vector<float> _scenarioParameters;
        std::vector<float> _normalizedScenarioParameters;
        std::vector<std::vector<std::string>> _scalingParameters;

};



class ClassificationModel : public Model{

    public:

        ClassificationModel(const std::string folderName): Model(folderName),
        _model(cppflow::model(folderName + "Classification_model")) {}

        // Evaluate scenario damage from data
        float _evaluateOutput(Input data){

            // Rescale the input data
            std::vector<float> normalizedInputs = _normalizeInputs(data);

            // Give the input the good format and feed it to the NN
            auto scenarioTensorFlow = cppflow::tensor(normalizedInputs, {1,9});
            auto threatProbability = _model({{"serving_default_input1:0", scenarioTensorFlow}}, 
                                                            {"StatefulPartitionedCall:0"});
            return threatProbability[0].get_data<float>()[0];
        }

    public:
        cppflow::model _model;
        const std::string _folderName;
        
};

class RegressionModel : public Model{

    public:

        RegressionModel(const std::string folderName): Model(folderName),
        _model(cppflow::model(folderName + "Regression_model")) {}

        // Evaluate scenario damage from data
        float _evaluateOutput(Input data){

            // Normalize the input data
            std::vector<float> normalizedInputs = _normalizeInputs(data);

            // Give the input the good format and feed it to the NN
            auto scenarioTensorFlow = cppflow::tensor(normalizedInputs, {1,9});
            auto damageNormalized = _model({{"serving_default_input1:0", scenarioTensorFlow}, 
                                    {"serving_default_input2:0", cppflow::fill({1, 1}, -1.0f)}, 
                                    {"serving_default_input3:0", cppflow::fill({1, 1}, -1.0f)}}, 
                                    {"StatefulPartitionedCall:0", "StatefulPartitionedCall:1"});

            // Rescale and return the output data
            return _rescaleOutput(damageNormalized[0].get_data<float>()[0]);
        }

    public:
        cppflow::model _model;
        const std::string _folderName;       
};



int main() {

    // Instantiate the input objects with the scenario properties
    const std::string propertiesFile = "../input.dat";
    Input data(propertiesFile);

    // Instantiate the classifiction and regression models
    const std::string folderName = "../models/BlastRad1/";
    ClassificationModel classificationModel(folderName);
    RegressionModel regressionModel(folderName);

    // Instantiate the population grid vector
    const std::string popGridFile = "../pop-grids/popgrid-2020-2pt5arcmin.bin";
    PopGrid popGrid(popGridFile);

    // Define some interesting coordinates
    double latitude = 48.8647;
    double longitude = 2.3490;

    // Find a specific density (ony once)
    double population = popGrid.getCellPop(latitude, longitude);
    std::cout << "The population in this cell is: " << population << std::endl;
    
    // Predict the probability of any sort of damage
    float threatProbability = classificationModel._evaluateOutput(data);

    // Prediction of ground damage radius for dangerous scenarios
    float damageRadius;
    int peopleAffected;
    if  (threatProbability > 0.5){
        damageRadius = regressionModel._evaluateOutput(data);
        peopleAffected = 0.1 * pi * damageRadius * damageRadius * (worldPop/(4 * pi * earthRadius * earthRadius));
        std::cout << "The number of people affected is: " << peopleAffected << std::endl;
    } else {
        std::cout << "The number of people affected is 0" << std::endl;
    }

    // Test: other predictions using the same model objects (e.g. increase diameter by 10m)
    for (int i=0; i<5; i++){
        data._scenarioParameters[0] += 10; 
        damageRadius = regressionModel._evaluateOutput(data);
        peopleAffected = 0.1 * pi * damageRadius * damageRadius * (worldPop/(4 * pi * earthRadius * earthRadius));
        std::cout << "The number of people affected is: " << peopleAffected << std::endl;
    }

    return 0;
}
