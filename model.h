class Model{

    public:
        /**
         * @brief Construct the Model object
         * 
         */

        Model(const std::string folderName);

        /**
         * @brief Method to read the weights and store the three layers
         * 
         */
        void readWeightsFile();

        /**
         * @brief Method to read and store the scaling parameters
         * 
         */

        void readingScalingParameters();

        /**
         * @brief Method to structure the input vector with all parameters
         * 
         * @param x position vector of the astroid (fixed, Earth-centered coordinate system)
         * @param v velocity vector of the astroid (fixed, Earth-centered coordinate system)
         * @return double vector containing the 9 inputs to feed to the NN
         * 
         */

        std::vector<double> structInputVector(std::vector<double> x, std::vector<double> v);

        /**
         * @brief Method to compute the absolute velocity of the asteroid 
         * 
         * @param v velocity vector of the astroid (fixed, Earth-centered coordinate system)
         * @return double absolute velocity of the asteroid
         */

        double computeVelocity(std::vector<double> v);

        /**
         * @brief Method to compute the incidence angle of the asteroid
         * 
         * @param x position vector of the astroid (fixed, Earth-centered coordinate system)
         * @param v velocity vector of the astroid (fixed, Earth-centered coordinate system)
         * @return double incidence angle of asteroid (w.r.t. horizontal), in range [0-90]
         */

        double computeIncidenceAngle(std::vector<double> x, std::vector<double> v);

        /**
         * @brief Method to compute the azimuth (PAIR format) of the asteroid
         * 
         * @param x position vector of the astroid (fixed, Earth-centered coordinate system)
         * @param v velocity vector of the astroid (fixed, Earth-centered coordinate system)
         * @return double azimuth angle, clockwise from North, in range [0-360]
         */

        double computeAzimuth(std::vector<double> x, std::vector<double> v);

        /**
         * @brief Method to transform position vector from cartesian to spherical coordinates
         * 
         * @param x position vector of the astroid (fixed, Earth-centered coordinate system) (x, y, z)
         * @return std::vector<double> position vector of the asteroid in spherical system (r, theta, phi)
         */

        std::vector<double> cartesianToSpherical(std::vector<double> x);

        /**
         * @brief Method to normalize the input parameters before feeding the neural network
         * 
         * @param data Input object containing real-scale asteroid properties and trajectories
         * @return std::vector<double> vector containing normalized asteroid properties and trajectories
         */

        std::vector<double> normalizeInputs(std::vector<double> data);

        /**
         * @brief Method to rescale the neural network output to real-life casualties
         * 
         * @param damageNormalized output of the neural network
         * @return double real-life casualty radius from the asteroid(s) impacts(s)
         */

        double rescaleOutput(double damageNormalized);


        /**
         * @brief Method to compute A^Tx = b
         * 
         * @param A matrix stored in a vector of size (mA x nA), row-major-ordered
         * @param x vector of size (mB x 1)
         * @param mA number of rows of matrix A
         * @param nA number of columns of matrix A
         * @return A^Tx
         */
        
        std::vector<double> matrixTransposedVectorMultiplication(std::vector<double> A, std::vector<double> x, size_t mA, size_t nA);

        /**
         * @brief Method to compute Ax = b
         * 
         * @param A matrix stored in a vector of size (mA x nA), row-major-ordered
         * @param x vector of size (mB x 1)
         * @param mA number of rows of matrix A
         * @param nA number of columns of matrix A
         * @return Ax
         */

        std::vector<double> matrixVectorMultiplication(std::vector<double> A, std::vector<double> x, size_t mA, size_t nA);

        /**
         * @brief Method for elementwise multiplication between two vectors
         * 
         * @param a First vector
         * @param b Second vector
         * @return vector containing a[i] * b[i]
         */

        std::vector<double> vectorsElementwiseMultiplication(std::vector<double> a, std::vector<double> b);

        /**
         * @brief Method for elementwise addition between two  vectors
         * 
         * @param a First vector
         * @param b Second vector
         * @return vector containing a[i] + b[i]
         */

        std::vector<double> vectorsAddition(std::vector<double> a, std::vector<double> b);

        /**
         * @brief Method to compute the elementwise max(a,0)
         * 
         * @param a Vector
         * @return vector containing max(a[i],0)
         */

        std::vector<double> rectifiedLinearUnit(std::vector<double> a);

        /**
         * @brief Method to compute the elementwise derivative of max(a,0)
         * 
         * @param a Vector
         * @return vector containing the derivative of max(a[i],0)
         */

        std::vector<double> derRectifiedLinearUnit(std::vector<double> a);

        /**
         * @brief Method to infer the output of a neural network containing 4 layers
         * 
         * @param input vector fed to the neural network
         * @return output vector of the neural network
         */

        std::vector<double> forwardPropagation(std::vector<double> input);


        double evaluateOutput(std::vector<double> data);

        /**
         * @brief Method to compute the derivative of the ouput of a neural network w.r.t the inputs
         * 
         * @return derivatives of the outputs with respect to the inputs (all normalized)
         */

        std::vector<double> backPropagation();

        /**
         * @brief Method to compute the derivatives of full-scale output w.r.t. full-scale inputs
         * 
         * @param pop number of people affected by the asteroid scenario
         * @param rad radius of damage
         * @return std::vector <double> 
         */

        std::vector <double> globalDerivatives(double pop, double rad);

        /**
         * @brief Attributes of the model
         * 
         */

        // Folder path for the neural network model
        const std::string _folderName;
        // Number of units of the different layers (including input and output)
        std::vector<int> _layers = {9,64,128,256,1};
        // Intermediate pre-activation results of the neural network
        std::vector<double> _z1, _z2, _z3, _z4;
        // Intermediate post-activation results of the neural network
        std::vector<double> _a1, _a2, _a3;
        // Intermediate derivatives 
        std::vector<double> _dZ4dZ3, _dZ4dA2, _dZ4dZ2, _dZ4dA1, _dZ4dZ1, _dZ4dA0;
        // Vector of vectors containing all the weights of the neural network
        std::vector<std::vector<double> > _weights;
        // Trajectory and material parameters
        std::vector<double> _scenarioParameters;
        // Scaling parameters for normalization of inputs/outputs
        std::vector<std::vector<std::string> > _scalingParameters;
};