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



int main() {

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
    double damageRadius = std::max(model.evaluateOutput(data._scenarioParameters),0.0);
    std::cout << "The radius of damage is: " << damageRadius <<  " m" << std::endl;

    // Find the number of people affected
    double affectedPop = 0.1 * popGrid.getAffectedPop(data._latitude, data._longitude, damageRadius);
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