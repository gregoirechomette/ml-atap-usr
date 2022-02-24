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
