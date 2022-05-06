# ml-atap-usr
Use accurate machine learning models to predict the number of people affected by an asteroid impact scenario given the entry conditions. Compute derivatives of the number of people affected w.r.t. the entry conditions in order to help mission design teams.

## How to install it
### Clone the current repository
```
git clone https://github.com/gregoirechomette/ml-atap-usr.git
```

### Select the trained ML model to use
In order to predict the number of people affected by an asteroid impact scenario, a machine learning model is used to assess the size of the ground damage (more specifically, the radius of the circular damaged area). Two well-trained machine learning models are uploaded with this project:
- BlastRad 1: outputs the radius of the circular area with at least 1psi blast overpressure
- ThermRad 2: outputs the radius of the circular area with thermal radiations leading to at least 2nd degree burns  

If other models need to be used, either contact gregoire.a.chomette@nasa.gov or generate them with ``` https://github.com/gregoirechomette/ml-atap-dv.git``` - they then need to be converted from a tensorflow format to a csv file using the script ``` python-utils/model2csv.py```

-> The path to the ML model selected needs to be specified with the variable ```folderName``` in the main function of ```main.cc```. **By default, the main function uses the BlastRad 1 model and can be used without modifications.**


### Obtain the world population data
The world population data is not pushed here because of its large size (~300MB). You can reach out directly to me (gregoire.a.chomette@nasa.gov) to receive a binary file containing the world grid data. 

-> The path to the population file needs to be updated with the variable ``` popGridFile ``` in the main function of ```main.cc```.


## How to use it

### Inputs and outputs
The ML model needs to be fed with the position and velocity vector of the asteroid at the time of impact, using Earth-centered ICRF as the reference frame. The position [x;y;z] and velocity [vx;vy;vz] are two 3D vectors with units of the international system ([m] and [m/s]). The structure definition and an example with numerical values are provided in ```main.cc```

```
struct Input{
    std::vector<double> positionVector;     // [m]
    std::vector<double> velocityVector;     // [m/s]
};
```

The outputs of the program can be retrieved through an output structure - an example is provided below:
```
struct Output{
    double affectedPopulation;              // number of people N [-]
    double velocityDerivative;              // dN/dvelocity [s/m]
    double angleDerivative;                 // dN/dangle [1/degrees]
    double azimuthDerivative;               // dN/dazimuth [1/degrees]
};
```


### Compile the C++ file and create the executable
The C++ program is compiled through the CMake build system:
```
cd ml-atap-usr
mkdir build
cd build
cmake ..
make
```

### Run the program (after creating the executable)
```
cd ml-atap-usr/build
./main
```