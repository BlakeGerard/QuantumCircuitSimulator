#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <iostream>
#include <complex>
#include "q_circuit.h"

TEST_CASE("Quantum Teleportation") {
    Eigen::Vector2cd ket_zero;
    Eigen::Vector2cd ket_one;
    ket_zero << 1.0, 0.0;
    ket_one << 0.0, 1.0;
    Eigen::Vector2cd qubit_1;
    Eigen::Vector2cd qubit_2;
    Eigen::Vector2cd qubit_m;
    qubit_1 << 1.0, 0.0;
    qubit_2 << 1.0, 0.0;
    qubit_m << 1.0/sqrt(3.0), sqrt(2.0/3.0);
    Eigen::VectorXcd state;
    state.resize(8);
    state = kroneckerProduct(kroneckerProduct(qubit_m, qubit_1).eval(), qubit_2).eval();

    Eigen::Matrix2d H;
    H << 1.0/sqrt(2), 1.0/sqrt(2), 1.0/sqrt(2), -1.0/sqrt(2);
    Eigen::Matrix2d I;
    I << 1.0, 0.0, 0.0, 1.0;
    Eigen::Matrix2d X;
    X << 0.0, 1.0, 1.0, 0.0;

    Eigen::Matrix2d Z;
    Z << 1.0, 0.0, 0.0, -1.0;
    Eigen::Matrix4d CNOT;
    CNOT << 1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 1.0, 0.0;

    Eigen::Matrix4d CNOT2;
    CNOT2 << 1.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 1.0,
             0.0, 0.0, 1.0, 0.0,
             0.0, 1.0, 0.0, 0.0;

    Eigen::MatrixXd I_H_I;
    I_H_I = kroneckerProduct(kroneckerProduct(I, H).eval(), I);
    
    Eigen::MatrixXd I_CNOT;
    I_CNOT = kroneckerProduct(I, CNOT).eval();

    Eigen::MatrixXd CNOT_I;
    CNOT_I = kroneckerProduct(CNOT, I).eval();

    Eigen::MatrixXd H_I_I;
    H_I_I = kroneckerProduct(kroneckerProduct(H, I).eval(), I).eval();

    Eigen::MatrixXd I_I_X;
    I_I_X = kroneckerProduct(kroneckerProduct(I, I).eval(), X).eval();

    Eigen::MatrixXd I_I_Z;
    I_I_Z = kroneckerProduct(kroneckerProduct(I, I).eval(), Z).eval();

    Eigen::MatrixXd H_I_H;
    H_I_H = kroneckerProduct(kroneckerProduct(H, I).eval(), H).eval();
    std::cout << H_I_H * state << std::endl;

    state = H_I_I*(CNOT_I*(I_CNOT*(I_H_I*state)));

    std::default_random_engine generator;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
    double prob = (1.0/4.0)*(pow(abs(qubit_m(0)), 2.0) + pow(abs(qubit_m(1)), 2.0));
    std::discrete_distribution<int> distribution {prob, prob, prob, prob};

    int result = distribution(generator);

    if (result == 0) {
        std::cout << "Measured |00>" << std::endl;
        state = (kroneckerProduct(kroneckerProduct(ket_zero, ket_zero).eval(), state.segment(0, 2).eval()).eval())/std::complex<double>(sqrt(prob), 0.0);
    } else if (result == 1) {
        std::cout << "Measured |01>" << std::endl;
        state = (kroneckerProduct(kroneckerProduct(ket_zero, ket_one).eval(), state.segment(2, 2).eval()).eval())/std::complex<double>(sqrt(prob), 0.0);
        state = I_I_X*state;
    } else if (result == 2) {
        std::cout << "Measured |10>" << std::endl;
        state = (kroneckerProduct(kroneckerProduct(ket_one, ket_zero).eval(), state.segment(4, 2).eval()).eval())/sqrt(prob);
        state = I_I_Z * state;
    } else {
        std::cout << "Measured |11>" << std::endl;
        state = (kroneckerProduct(kroneckerProduct(ket_one, ket_one).eval(), state.segment(6, 2).eval()).eval())/sqrt(prob);
        state = I_I_Z * (I_I_X * state);
    };
    std::cout << "-------------------------------" << std::endl;
    std::cout << state << std::endl;
    std::cout << "-------------------------------" << std::endl;
};

/*
Tested cases:
{1}: Passed in Quantum teleportation circuit
{0, 2}: Passed in Quantum teleportation circuit

*/
TEST_CASE("Evaluating Hadamard gate application to one or multiple qubits in QT circuit") {
    Qubit qubit_m = Qubit(1.0/sqrt(3.0), sqrt(2.0/3.0));
    Qubit qubit_1 = Qubit();
    Qubit qubit_2 = Qubit();
    std::vector<Qubit> qubits = {qubit_m, qubit_1, qubit_2};

    Eigen::VectorXcd initialization_state;
    initialization_state.resize(8);
    initialization_state << 1.0/sqrt(3.0), 0.0, 0.0, 0.0, sqrt(2.0/3.0), 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);
    REQUIRE(circuit.get_state() == initialization_state);

    Eigen::VectorXcd I_H_I_state;
    I_H_I_state.resize(8);
    I_H_I_state << 1.0/sqrt(3.0) * 1.0/sqrt(2.0), 0.0, 1.0/sqrt(3.0) * 1.0/sqrt(2.0), 0.0, 1.0/sqrt(3.0), 0.0, 1.0/sqrt(3.0), 0.0;

    std::vector<int> qubit_indices = {1};
    circuit.H(qubit_indices);
    std::cout << circuit.get_state() << std::endl;
    REQUIRE(fabs(circuit.get_state()(0) - I_H_I_state(0)) < 0.00001);
};

TEST_CASE("Evaluating CNOT gate on adjacent qubits") {
    Qubit qubit_0 = Qubit(0.0, 1.0);
    Qubit qubit_1 = Qubit(1.0, 0.0);
    Qubit qubit_2 = Qubit(0.0, 1.0);
    std::vector<Qubit> qubits = {qubit_0, qubit_1, qubit_2};

    Eigen::VectorXcd expected_state;
    expected_state.resize(8);
    expected_state << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);

    std::vector<int> qubit_indices{0, 2};
    circuit.CNOT(qubit_indices);
    REQUIRE(circuit.get_state() == expected_state);
};