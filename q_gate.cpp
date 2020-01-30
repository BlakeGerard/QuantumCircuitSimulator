#include "q_gate.h"
#include <iostream>

void hadamard(Qubit *qubit) {
    Eigen::Matrix2d hadamard_matrix;
    hadamard_matrix << 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), -1/sqrt(2);

    qubit->set_state(hadamard_matrix * qubit->get_state());
};

Eigen::Vector4cd cnot(Qubit *control, Qubit *target) {
    Eigen::VectorXcd tensor_product = Eigen::kroneckerProduct(control->get_state(), target->get_state());
    Eigen::MatrixXd cnot_matrix;

/*
    cnot_matrix.resize(tensor_product.size());

    if (tensor_product.size() == 2) {
        cnot_matrix << 1.0, 0.0, 0.0, 0.0,
                       0.0, 1.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 1.0,
                       0.0, 0.0, 1.0, 0.0;
    } else if (tensor_product.size() == 4) {
        cnot_matrix << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;
    };
    return cnot_matrix * tensor_product;
    */
    
   return tensor_product;

};