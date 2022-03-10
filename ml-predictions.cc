#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <stdio.h>

#include "./lib/cppflow/include/cppflow/ops.h"
#include "./lib/cppflow/include/cppflow/model.h"

#define WITH_LOCAL_DENSITIES

// Constants
const double pi = 3.14159;
double const deg2rad = 0.01745329;
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
            _readFile();
        }

    private:

        /**
         * @brief Method to read and store the dat file containing scenario parameters
         * 
         */
        
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

        /**
         * @brief Constructor of a new Pop Grid object
         * 
         * @param fileName path and name of the binary file to be read
         */

        PopGrid(const std::string fileName): _fileName(fileName) {
                readFile();
        }

        /**
         * @brief Method to read the binary file and store the information about densities
         * 
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
         * @param[out] popCount number of people in the grid cell
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

        /**
         * @brief Get the distance between two points in the Earth based on their lat and long
         * 
         * @param lat1 latidude of the first point [degrees]
         * @param lon1 longitude of the first point [in degrees]
         * @param lat2 latidude of the second point [in degrees]
         * @param lon2 longitude of the second point [in degrees]
         * @return double distance [meters]
         */

        double getDistance(double lat1, double lon1, double lat2, double lon2){

            double dTheta, dLat, dLon;
            dLat = deg2rad*(lat1 - lat2);
            dLon = deg2rad * (lon1 - lon2);
            dTheta = 2.0*asin(sqrt(sin(0.5*dLat)*sin(0.5*dLat) + cos(lat1*deg2rad)*cos(lat2*deg2rad)*sin(dLon*0.5)*sin(dLon*0.5)));
            return dTheta*earthRadius;
        }

        /**
         * @brief Get the total affected population given the radius of a circle, together with its lat and long
         * 
         * @param latitude latitude of the center of the circle [degrees]
         * @param longitude longitude of the center of the circle [degrees]
         * @param damagedRadius radius of the circular damaged area [meters]
         * @return double number of people affected [-]
         */

        double getAffectedPop(double latitude, double longitude, double damagedRadius){

            // Number of cells to check on the latitude direction
            int cellLatNumber = 2 * asin(0.5 * damagedRadius / earthRadius) / arcLengthDeg + 1;
            int cellLonNumber = 2 * asin(0.5 * damagedRadius / earthRadius) / arcLengthDeg / cos(deg2rad * latitude) + 1;

            double distanceComp;
            double totalAffPop = 0.0;
            for (int i=-cellLatNumber; i <=cellLatNumber; i++){
                for (int j=-cellLonNumber; j<=cellLonNumber; j++){
                    distanceComp = getDistance(latitude, longitude, latitude + i * arcLengthDeg, longitude + j * arcLengthDeg);
                    if (distanceComp < damagedRadius){
                        totalAffPop += getCellPop(latitude + i * arcLengthDeg, longitude + j * arcLengthDeg);
                    }
                }
            }
            return totalAffPop;

        }

        // Atributes
        const std::string _fileName;
        int nRows, nCol, xlCorner, ylCorner, yuCorner, noDat, vSize;
        double arcLengthDeg, maxcellsidekm, maxcelldiagkm, dum;
        std::vector<double> popCount;

};


class Model {

    public:

        /**
         * @brief Constructor of a new Model object
         * 
         * @param folderName path of the folder containing the models and scaling parameters
         */

        Model(const std::string folderName): 
        _folderName(folderName) {
            _readingScalingParameters();
        }

        /**
         * @brief Method to read and store the scaling parameters
         * 
         */

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

        /**
         * @brief Method to normalize the input parameters before feeding the neural network
         * 
         * @param data Input object containing real-scale asteroid properties and trajectories
         * @return std::vector<float> vector containing normalized asteroid properties and trajectories
         */

        std::vector<float> _normalizeInputs(Input data){
            std::vector<float> normalizedInputs (9);
            for (int i=0; i<9; i++){
                float mean = std::stof(_scalingParameters[1][i+1]);
                float std = std::stof(_scalingParameters[2][i+1]);
                normalizedInputs[i] = (data._scenarioParameters[i] - mean)/ (std);
            }
            return normalizedInputs;
        }

        /**
         * @brief Method to rescale the neural network output to real-life casualties
         * 
         * @param damageNormalized output of the neural network
         * @return float real-life casualty radius from the asteroid(s) impacts(s)
         */

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

        /**
         * @brief Constructor of a new Classification Model object
         * 
         * @param folderName name of the folder containing the saved and trained TF model
         */

        ClassificationModel(const std::string folderName): Model(folderName),
        _model(cppflow::model(folderName + "Classification_model")) {}

        /**
         * @brief Method to evaluate the propability of casuatlies on the Earth
         * 
         * @param data Input object containing real-scale asteroid properties and trajectories
         * @return float probability of casualty (0=no damage; 1=damage) on the Earth for a given level (e.g. BlastRad1)
         */
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

        /**
         * @brief Constructor of a new Regression Model object
         * 
         * @param folderName name of the folder containing the saved and trained TF model
         */

        RegressionModel(const std::string folderName): Model(folderName),
        _model(cppflow::model(folderName + "Regression_model")) {}

        /**
         * @brief Method to evaluate the extent of casuatlies on the Earth
         * 
         * @param data Input object containing real-scale asteroid properties and trajectories
         * @return float radius of the damaged area for a given level (e.g. BlastRad1)
         */

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

    // Boolean to state if local densities of populations should be used
    bool localDensities = true;

    // Instantiate the input objects with the scenario properties
    const std::string propertiesFile = "../input.dat";
    Input data(propertiesFile);

    // Instantiate the classifiction and regression models
    const std::string folderName = "../models/BlastRad1/";
    ClassificationModel classificationModel(folderName);
    RegressionModel regressionModel(folderName);

    #if defined(WITH_LOCAL_DENSITIES)
    // Instantiate the population grid vector
    const std::string popGridFile = "../pop-grids/popgrid-2020-2pt5arcmin.bin";
    PopGrid popGrid(popGridFile);

    // Define some dummy coordinates (Paris)
    double latitude = 48.8647;
    double longitude = 2.3490;
    #endif
    
    // Predict the probability of any sort of damage
    double threatProbability = classificationModel._evaluateOutput(data);

    // Prediction of ground damage radius for dangerous scenarios
    double damageRadius;
    int peopleAffected;
    if  (threatProbability > 0.5){
        damageRadius = regressionModel._evaluateOutput(data);
        #if defined(WITH_LOCAL_DENSITIES)
        peopleAffected = popGrid.getAffectedPop(latitude, longitude, damageRadius);
        #else
        peopleAffected = 0.1 * pi * damageRadius * damageRadius * (worldPop/(4 * pi * earthRadius * earthRadius));
        #endif
        std::cout << "The number of people affected is: " << peopleAffected << std::endl;
        
    } else {
        std::cout << "The number of people affected is 0" << std::endl;
    }

    return 0;
}
