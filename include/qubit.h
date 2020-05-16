/*
        Title: qubit.h
        Author: Blake Gerard
        Date: 05/16/2020
        Description: Class to define a single qubit as an Eigen vector of length two. 
        The state of a quantum bit can be represented as a vector
        of the form [a b]^T, where a and b are the complex probabilities that
        the qubit will collapse to 0 or 1, respectively, on measurement.
*/

#ifndef QUBIT_H
#define QUBIT_H
#include <Eigen/Dense>

class Qubit {
    public:
        Qubit();
        Qubit(std::complex<double> a, std::complex<double> b);
        Eigen::Vector2cd get_state();
        void set_state(Eigen::VectorXcd new_state);
    private:
        Eigen::Vector2cd state;
};
#endif