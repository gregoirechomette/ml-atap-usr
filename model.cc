#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "model.h"


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


std::vector <double> Model::globalDerivatives(std::vector<double> nnDerivatives, double pop, double rad){
    std::vector<double> derivatives(9);
    for (int i=0; i<derivatives.size(); i++){
        derivatives[i] = 2 * pop * std::stof(_scalingParameters[2][i+1]) * nnDerivatives[i] / (std::stof(_scalingParameters[2][10]) * rad);
    }
    return derivatives;
}