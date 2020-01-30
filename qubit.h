#ifndef QUBIT_H
#define QUBIT_H
#include </usr/include/eigen3/Eigen/Dense>
#include <random>

class Qubit {
    public:
        Qubit();
        Qubit(std::complex<double> a, std::complex<double> b);
        Eigen::Vector2cd get_state();
        int get_size();
        void set_state(Eigen::VectorXcd new_state);
        double measure();
    
    private:
        int size;
        Eigen::Vector2cd state;
        std::default_random_engine generator;
};

#endif