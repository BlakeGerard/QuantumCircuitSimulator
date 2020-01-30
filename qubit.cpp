#include <chrono>
#include <complex>
#include <iostream>
#include "qubit.h"

Qubit::Qubit() {
    state << 1.0, 0.0;
    size = 2;

    // Seed the generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
};

Qubit::Qubit(std::complex<double> a, std::complex<double> b) {
    state << a, b;
    size = 2;
    
    if (state.norm() != 1.0) {
        state.normalize();
    };

    // Seed the generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
};

Eigen::Vector2cd Qubit::get_state() {
    return state;
};

int Qubit::get_size() {
    return size;
};

void Qubit::set_state(Eigen::VectorXcd new_state) {
    state.resize(new_state.size());
    size = new_state.size();
    state = new_state;
}

double Qubit::measure() {
    std::vector<double> amplitudes = std::vector<double>();
    double aProb = pow(abs(state(0)), 2.0);
    double bProb = pow(abs(state(1)), 2.0);
    std::discrete_distribution<int> distribution {aProb, bProb};

    return (distribution(generator));
};