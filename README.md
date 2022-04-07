# ml-atap-usr
Load asteroid damage ML models from a C++ code. Use the models for inference to predict the ground damage of asteroid impacts as well as the sensitivity (derivative) of the ground damages w.r.t. the input parameters.

## How to install it
### Clone the current repository
```
git clone https://github.com/gregoirechomette/ml-atap-usr.git
```

### Request the grid population file
The world population data is available [online](https://sedac.ciesin.columbia.edu/data/collection/gpw-v4) - not pushed here because of its large size (~300MB). Alternatively, you can reach out directly to me (gregoire.a.chomette@nasa.gov) to receive a binary file containing the world population data. Then, the path to the population file needs to be updated with the variable ``` popGridFile ``` in the main function.


## How to use it

### Enter the correct asteroid entry conditions
The asteroid properties and trajectories can be modified in the ```input.dat``` file. An example is provided:
```
50              Diameter [m] 
3500            Density [kg/m^3]
100000          Strength [Pa]
0.2             Alpha (strength scaling coefficient) [-]
10000           Velocity [m/s]
45              Incidence Angle [degrees]
180             Azimuth [degrees]
0.003           Luminous Efficiency [-]
0.000000001     Ablation [kg m^-3]
48.8647         Latitude [degrees]
2.3490          Longitude [degrees]
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

### Run the program
```
cd ml-atap-usr/build
./ml-predictions
```

<!-- ## Key run time results

Instantiate the TensorFlow C model: ~10^-1 seconds  
Evaluate the damage radius for one scenario: ~10^-4 seconds

Instantiate the population grid vector: ~ 10^0 seconds  
Evaluate the population affected for one scenario: ~ 10^-6 seconds -->