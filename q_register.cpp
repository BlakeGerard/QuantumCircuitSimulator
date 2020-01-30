#include "q_register.h"
#include </usr/include/eigen3/unsupported/Eigen/KroneckerProduct>
#include <iostream>

Q_Register::Q_Register(int n_qubits) {
    this->n_qubits = n_qubits;
};

void Q_Register::initialize() {
    Eigen::Vector2cd ket_zero;
    ket_zero << 1.0, 0.0;
    Eigen::Vector2cd ket_one;
    ket_one << 0.0, 1.0;

    Eigen::VectorXcd temp_state_vector = ket_zero;
    for (int i = 1; i < n_qubits; ++i) {
        temp_state_vector = Eigen::kroneckerProduct(temp_state_vector, ket_zero).eval();
    }
    state = temp_state_vector;
};

void Q_Register::set_bit_string(std::vector<int> bit_string) {
    assert(pow(2, bit_string.size()) == pow(2, n_qubits));

    Eigen::Vector2cd ket_zero;
    ket_zero << 1.0, 0.0;
    Eigen::Vector2cd ket_one;
    ket_one << 0.0, 1.0;

    Eigen::VectorXcd temp_state_vector;
    if (bit_string.at(0) == 0) {
        temp_state_vector = ket_zero;
    } else {
        temp_state_vector = ket_one;
    }

    for (unsigned int i = 1; i < bit_string.size(); ++i) {
        Eigen::Vector2cd next_qubit;
        if (bit_string.at(i) == 0) {
            next_qubit = ket_zero;
        } else {
            next_qubit = ket_one;
        }
        temp_state_vector = Eigen::kroneckerProduct(temp_state_vector, next_qubit).eval();
    }
    state = temp_state_vector;
};

void Q_Register::set_state(Eigen::VectorXcd new_state) {
    assert(new_state.size() == state.size());
    state = new_state;

};

Eigen::VectorXcd Q_Register::get_state() {
    return state;
};