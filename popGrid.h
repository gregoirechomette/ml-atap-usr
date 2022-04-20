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


class PopGrid {

    public:

        /**
         * @brief Constructor of a new Pop Grid object
         * 
         * @param fileName path and name of the binary file to be read
         */

        PopGrid(const std::string fileName);

        /**
         * @brief Method to read the binary file and store the information about densities
         * 
         */

        void readFile();

        /**
         * @brief Method to find the density given the latitude and longitude
         * 
         * @param[in] latitude in degrees
         * @param[in] longtitude in degrees
         * @param[out] popCount number of people in the grid cell
         */

        double getCellPop(double latitude, double longitude);

        /**
         * @brief Get the distance between two points in the Earth based on their lat and long
         * 
         * @param lat1 latidude of the first point [degrees]
         * @param lon1 longitude of the first point [in degrees]
         * @param lat2 latidude of the second point [in degrees]
         * @param lon2 longitude of the second point [in degrees]
         * @return double distance [meters]
         */

        double getDistance(double lat1, double lon1, double lat2, double lon2);

        /**
         * @brief Get the total affected population given the radius of a circle, together with its lat and long
         * 
         * @param latitude latitude of the center of the circle [degrees]
         * @param longitude longitude of the center of the circle [degrees]
         * @param damagedRadius radius of the circular damaged area [meters]
         * @return double number of people affected [-]
         */

        double getAffectedPop(double latitude, double longitude, double damagedRadius);

        // Atributes
        const std::string _fileName;
        int nRows, nCol, xlCorner, ylCorner, yuCorner, noDat, vSize;
        double arcLengthDeg, maxcellsidekm, maxcelldiagkm, dum;
        std::vector<double> popCount;

};