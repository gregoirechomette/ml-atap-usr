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

#include "popGrid.h"

// Constants
const double pi = 3.14159;
const double deg2rad = 0.01745329;
const double earthRadius = 6.371e6;
const double worldPop = 7.93 * pow(10,9);


PopGrid::PopGrid(const std::string fileName): _fileName(fileName) {
    readFile();
}


void PopGrid::readFile(){

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



double PopGrid::getCellPop(double latitude, double longitude){

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


double PopGrid::getDistance(double lat1, double lon1, double lat2, double lon2){

    double dTheta, dLat, dLon;
    dLat = deg2rad*(lat1 - lat2);
    dLon = deg2rad * (lon1 - lon2);
    dTheta = 2.0*asin(sqrt(sin(0.5*dLat)*sin(0.5*dLat) + cos(lat1*deg2rad)*cos(lat2*deg2rad)*sin(dLon*0.5)*sin(dLon*0.5)));
    return dTheta*earthRadius;
}


double PopGrid::getAffectedPop(double latitude, double longitude, double damagedRadius){

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