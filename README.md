# ml-atap-usr
Load asteroid damage ML models from a C++ code with [cppflow](https://github.com/serizba/cppflow) and the [TensorFlow C API](https://www.tensorflow.org/install/lang_c). The models can then be used to evaluate the ground damages of asteroids based on atmospheric entry conditions.

## How to install it
### 1) Clone the current repository
```
git clone https://github.com/gregoirechomette/ml-atap-usr.git
```

### 2) Install the [cppflow](https://github.com/serizba/cppflow) library
We recommend to install the *cppflow* library in a *lib* folder inside the *ml-atap-usr* repository.  Suggested code:

```
cd ml-atap-usr
mkdir lib
cd lib
git clone https://github.com/serizba/cppflow.git
```
However, if the cppflow is installed in another location, the ```#include``` instructions of the C++ files must be updated with the correct locations, i.e. ```#include "<path-to-cppflow>/include/cppflow/ops.h" ```.


### 3) Install the [TensorFlow C API](https://www.tensorflow.org/install/lang_c) library
We also recommend to install the *TensorFlow C API* library in the *lib* repository already created to install *cppflow*. The *TensorFlow C API* can be downloaded [here](https://www.tensorflow.org/install/lang_c), and the associated documentation is found [here](https://serizba.github.io/cppflow/installation.html). Suggested code:
```
cd ml-atap-usr/lib
mkdir libtensorflow2
tar -C <path-to-ml-atap-usr>/lib/libtensorflow2 -xf ~/Downloads/libtensorflow.tar
```
If the *TensorFlow C API* is installed in another location, the ```CMakeLists.txt``` configuration file must be modified.

### 4) Obtain the trained TensorFlow ML models
These models are not uploaded here, please reach out directly to me (gregoire.a.chomette@nasa.gov) if you are interested, or see the repository [ml-atap-dv](https://github.com/gregoirechomette/ml-atap-dv) to create and train these models. Then, the path to the model of interest needs to be updated with the variable ``` folderName ```.

### 5) Obtain the grid population file
The world population data is available [online](https://sedac.ciesin.columbia.edu/data/collection/gpw-v4). Alternatively, you can reach out directly to me to receive a binary file containing the world population data. Then, the path to the population file needs to be updated with the variable ``` popGridFile ```.


## How to use it

### 1) Enter the correct asteroid entry conditions
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


### 2) Compile the C++ file and create the executable
The C++ program is compiled through the CMake build system:
```
cd ml-atap-usr
mkdir build
cd build
cmake ..
make
```

### 3) Run the program
```
cd ml-atap-usr/build
./ml-predictions
```

## Key run time results

Instantiate the TensorFlow C model: ~10^-1 seconds  
Evaluate the damage radius for one scenario: ~10^-4 seconds

Instantiate the population grid vector: ~ 10^0 seconds  
Evaluate the population affected for one scenario: ~ 10^-6 seconds
