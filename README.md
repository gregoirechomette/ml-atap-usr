# ml-atap-usr
Use accurate machine learning models to predict the number of people affected by an asteroid impact scenario given the entry conditions. Compute derivatives of the number of people affected w.r.t. the entry conditions in order to help mission design teams.

## How to install it
### Clone the current repository
```
git clone https://github.com/gregoirechomette/ml-atap-usr.git
```

### Select the trained ML model to use
In order to predict the number of people affected by an asteroid impact scenario, a machine learning model is used to assess the size of the ground damage (more specifically, the radius of the circular damaged area). Two well-trained machine learning models are uploaded with this project:
- ./models/BlastRad 1: outputs the radius of the circular area with at least 1psi blast overpressure
- ./models/ThermRad 2: outputs the radius of the circular area with thermal radiations leading to at least 2nd degree burns  

If other models need to be used, either contact gregoire.a.chomette@nasa.gov or generate them with ``` https://github.com/gregoirechomette/ml-atap-dv.git``` - they then need to be converted from a tensorflow format to a csv file using the script ``` python-utils/model2csv.py```

The path to the ML model selected needs to be specified with the variable ```folderName``` in the main function of ```main.cc```:

```
// Instantiate the neural network model
const std::string folderName = "../models/BlastRad1/";
```


### Obtain the world population data
The world population data is not pushed here because of its large size (~300MB). You can reach out directly to me (gregoire.a.chomette@nasa.gov) to receive a binary file containing the world grid data. 

The path to the population file needs to be updated with the variable ``` popGridFile ``` in the main function of ```main.cc```:

```
// Instantiate the population grid vector
const std::string popGridFile = "../pop-grids/popgrid-2020-2pt5arcmin.bin";
```


## How to use it

### Inputs and outputs
The model **inputs** are the position and velocity vectors of the asteroid at the time of impact, using Earth-centered ICRF as the reference frame. The position [x;y;z] and velocity [vx;vy;vz] are two 3D vectors with units of the international system ([m] and [m/s]). The structure definition and an example with numerical values are provided in ```main.cc```

```
struct Input{
    std::vector<double> positionVector;     // [m]
    std::vector<double> velocityVector;     // [m/s]
};
```

The model **outputs** are the number of people affected together with the derivatives w.r.t. the velocity (ICRF reference frame), the incidence angle (w.r.t. the horizontal), and the azimuth (clockwise from North). The structure definition and an example with numerical values are provided in ```main.cc```

```
struct Output{
    double affectedPopulation;              // number of people N [-]
    double velocityDerivative;              // dN/dvelocity [s/m]
    double angleDerivative;                 // dN/dangle [1/degrees]
    double azimuthDerivative;               // dN/dazimuth [1/degrees]
};
```


### Compiling the C++ file and creating the executable
The C++ program is compiled through the CMake build system:
```
cd ml-atap-usr
mkdir build
cd build
cmake ..
make
```

### Running the program (after creating the executable)
```
cd ml-atap-usr/build
./main
```