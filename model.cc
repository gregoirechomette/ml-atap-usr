#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#include "model.h"

const double pi = 3.14159;

Model::Model(const std::string folderName):
        _folderName(folderName)
        {
            readWeightsFile(), 
            readingScalingParameters();
        }


void Model::readWeightsFile(){

    // Open the file and declare reading variables
    std::fstream file (_folderName + "weights.csv", std::ios::in);
    std::vector<double> row;
    std::string line, word;

    // Read the file and store data
    if(file.is_open()){
        while(getline(file, line)){
            row.clear();
            std::stringstream str(line);
            while(getline(str, word, ',')){
                row.push_back(std::stod(word));
            }
            _weights.push_back(row);
        }
    }
    else{
        std::cout << "Error: Cannot open weights file "  << std::endl;
    }
    return;
}


void Model::readingScalingParameters(){

    // Open the file and declare reading variables
    std::string scalingFileName = _folderName + "scaling_parameters.csv";
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

std::vector<float> Model::structInputVector(double velocity, double incidenceAngle, double azimuth){
    
    // Baseline material parameters 
    float diameter = 70.0;
    float density = 3500;
    float strength = 100000;
    float alpha = 0.2;
    float lumEff = 0.003;
    float ablation = 0.000000001;

    // Conversion of parameters from double to float
    float velocity_f = (float) velocity;
    float incidenceAngle_f = (float) incidenceAngle;
    float azimuth_f = (float) azimuth;

    // Group and return
    std::vector<float> inputNN = {diameter, density, strength, alpha, velocity_f, incidenceAngle_f, azimuth_f, lumEff, ablation};
    return inputNN;
}

float Model::computeVelocity(std::vector<double> v){
    
    // Compute the absolute velocity
    float velocity = (float) sqrt( pow(v[0],2) + pow(v[1],2) + pow(v[2],2) );
    return velocity;
}

float Model::computeIncidenceAngle(std::vector<double> x, std::vector<double> v){

    // Compute the scalar product
    double scalarProduct = x[0] * v[0] + x[1] * v[1] + x[2] * v[2];
    // Compute the norms of both vectors
    double xNorm = sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2));
    double vNorm = sqrt( pow(v[0],2) + pow(v[1],2) + pow(v[2],2));
    // Compute the incidence angle
    float incidenceAngle = (float) (180/pi) * (acos(scalarProduct/(xNorm * vNorm)) - 0.5 * pi);

    return incidenceAngle;
}

float Model::computeAzimuth(std::vector<double> x, std::vector<double> v){

    // Compute the norms of both vectors
    double xNorm = sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2));
    double vNorm = sqrt( pow(v[0],2) + pow(v[1],2) + pow(v[2],2));

    // Retrieve the spherical coordinates and angle of incidence
    std::vector<double> sphericalCoordinates = cartesianToSpherical(x);
    double incidenceAngle = (double) computeIncidenceAngle(x,v);

    // Compute the horizontal velocity vector
    std::vector<double> velocityHorizontal(3);
    for (int i(0); i<3; i++){
        velocityHorizontal[i] = v[i] + (x[i]/xNorm) * vNorm * sin(incidenceAngle * pi /180);
    }

    // Find the North and East direction in the horizontal plan
    std::vector<double> north(3);
    std::vector<double> east(3);
    north[0] = - cos(sphericalCoordinates[1]) * cos(sphericalCoordinates[2]);
    north[1] = - cos(sphericalCoordinates[1]) * sin(sphericalCoordinates[2]);
    north[2] = sin(sphericalCoordinates[1]);
    east[0] = - sin(sphericalCoordinates[2]);
    east[1] = cos(sphericalCoordinates[2]);
    east[2] = 0;

    // Project horizontal vector on North and East
    double velocityHorizontalNorth = velocityHorizontal[0] * north[0] 
                                    + velocityHorizontal[1] * north[1] 
                                    + velocityHorizontal[2] * north[2];
    double velocityHorizontalEast = velocityHorizontal[0] * east[0] 
                                    + velocityHorizontal[1] * east[1] 
                                    + velocityHorizontal[2] * east[2];
    if (abs(velocityHorizontalEast) < 0.001){
        velocityHorizontalEast = 0.001;
    }

    // Find azimuth with trigonometrical direction, starting from East
    double trigonometricAzimuth;
    if (velocityHorizontalEast > 0){
        trigonometricAzimuth = (180/pi) * atan(velocityHorizontalNorth/velocityHorizontalEast);
    }else{
        trigonometricAzimuth = (180/pi) * atan(velocityHorizontalNorth/velocityHorizontalEast) + pi;
    }

    float pairAzimuth = (float) 90 - trigonometricAzimuth;
    pairAzimuth = ((int(pairAzimuth) % 360) + 360) % 360;

    return pairAzimuth;
}

std::vector<double> Model::cartesianToSpherical(std::vector<double> x){

    // Create the output vector
    std::vector<double> sphericalCoordinates(3);

    // Populate the output vector (r, theta, phi)
    sphericalCoordinates[0] = sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2));
    sphericalCoordinates[1] = atan((sqrt(pow(x[0],2) + pow(x[1],2)))/(pow(x[2],2)));
    if (x[0]>0){
        sphericalCoordinates[2] = atan((x[1])/(x[0]));
    } else{
        sphericalCoordinates[2] = atan((x[1])/(x[0])) + pi;
    }

    return sphericalCoordinates;
}
    
std::vector<double> Model::normalizeInputs(std::vector<float> data){
    std::vector<double> normalizedInputs(9);
    for (int i=0; i<9; i++){
        double mean = std::stof(_scalingParameters[1][i+1]);
        double std = std::stof(_scalingParameters[2][i+1]);
        normalizedInputs[i] = (data[i] - mean)/ (std);
    }
    return normalizedInputs;
}


double Model::rescaleOutput(double damageNormalized){
    double mean = std::stof(_scalingParameters[1][10]);
    double std = std::stof(_scalingParameters[2][10]);
    return (damageNormalized * std) + mean;
}

        
std::vector<double> Model::matrixTransposedVectorMultiplication(std::vector<double> A, std::vector<double> x, size_t mA, size_t nA){
    std::vector<double> output(nA, 0.0);
    for (int j=0; j<nA; j++){
        for (int i=0; i<mA; i++){
            output[j] += A[i * nA + j] * x[i];
        }
    }
    return output;
}


std::vector<double> Model::matrixVectorMultiplication(std::vector<double> A, std::vector<double> x, size_t mA, size_t nA){
    std::vector<double> output(mA, 0.0);
    for (int i=0; i<mA; i++){
        for (int j=0; j<nA; j++){
            output[i] += A[i * nA + j] * x[j];
        }
    }
    return output;
}


std::vector<double> Model::vectorsElementwiseMultiplication(std::vector<double> a, std::vector<double> b){
    std::vector<double> output(a.size(), 0.0);
    for (int i=0; i<a.size(); i++){
        output[i] = a[i] * b[i];
    }
    return output;
}


std::vector<double> Model::vectorsAddition(std::vector<double> a, std::vector<double> b){
    for (int i=0; i<a.size(); i++){
        a[i] += b[i];
    }
    return a;
}


std::vector<double> Model::rectifiedLinearUnit(std::vector<double> a){
    for (int i=0; i<a.size(); i++){
        if(a[i] < 0){
            a[i] = 0;
        }
    }
    return a;
}


std::vector<double> Model::derRectifiedLinearUnit(std::vector<double> a){
    for (int i=0; i<a.size(); i++){
        if (a[i] < 0){
            a[i] = 0;
        }else{
            a[i] = 1;
        }
    }
    return a;
}


std::vector<double> Model::forwardPropagation(std::vector<double> input){

    _z1 = vectorsAddition(matrixTransposedVectorMultiplication(_weights[0], input, _layers[0], _layers[1]), _weights[1]);
    _a1 = rectifiedLinearUnit(_z1);
    _z2 = vectorsAddition(matrixTransposedVectorMultiplication(_weights[2], _a1, _layers[1], _layers[2]), _weights[3]);
    _a2 = rectifiedLinearUnit(_z2);
    _z3 = vectorsAddition(matrixTransposedVectorMultiplication(_weights[4], _a2, _layers[2], _layers[3]), _weights[5]);
    _a3 = rectifiedLinearUnit(_z3);
    _z4 = vectorsAddition(matrixTransposedVectorMultiplication(_weights[6], _a3, _layers[3], _layers[4]), _weights[7]);
    return _z4;
}

double Model::evaluateOutput(std::vector<float> data){

    // Normalize the input data
    std::vector<double> normalizedInputs = normalizeInputs(data);
    // Call fhe forward propagation method
    double damageNormalized = forwardPropagation(normalizedInputs)[0];
    // Rescale and return the output data
    return rescaleOutput(damageNormalized);
}


std::vector<double> Model::backPropagation(){
    _dZ4dZ3 = vectorsElementwiseMultiplication(derRectifiedLinearUnit(_z3), _weights[6]);
    _dZ4dA2 = matrixVectorMultiplication(_weights[4], _dZ4dZ3, _layers[2], _layers[3]);
    _dZ4dZ2 = vectorsElementwiseMultiplication(derRectifiedLinearUnit(_z2), _dZ4dA2);
    _dZ4dA1 = matrixVectorMultiplication(_weights[2], _dZ4dZ2, _layers[1], _layers[2]);
    _dZ4dZ1 = vectorsElementwiseMultiplication(derRectifiedLinearUnit(_z1), _dZ4dA1);
    _dZ4dA0 = matrixVectorMultiplication(_weights[0], _dZ4dZ1, _layers[0], _layers[1]);
    return _dZ4dA0;
}


std::vector <double> Model::globalDerivatives(double pop, double rad){
    std::vector<double> nnDerivatives = backPropagation();
    std::vector<double> derivatives(9);
    for (int i=0; i<derivatives.size(); i++){
        derivatives[i] = 2 * pop * std::stof(_scalingParameters[2][i+1]) * nnDerivatives[i] / (std::stof(_scalingParameters[2][10]) * rad);
    }
    return derivatives;
}