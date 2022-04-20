class Input {

    public:

        /**
         * @brief Constructor of a new Input object
         * 
         */

        Input();

        /**
         * @brief Method to read the binary file and store the information about densities
         * 
         * @param velocity velocity of the asteroid [m/s]
         * @param incidenceAngle angle of incidence of asteroid entering the atmosphere, 0: horizontal, 90: vertical [degrees]
         * @param azimuth azimuth of the asteroid according to PAIR conventions [degrees]
         * @return float vector containing the 9 inputs to feed to the NN
         * 
         */

        std::vector<float> prepareInput(double velocity, double incidenceAngle, double azimuth);
};