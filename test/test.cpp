#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <iostream>
#include <complex>
#include "q_circuit.h"

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
    REQUIRE(fabs(circuit.get_state()(0) - I_H_I_state(0)) < 0.00001);
};

TEST_CASE("Evaluating CNOT gate on non-adjacent qubits") {
    Qubit qubit_0 = Qubit(0.0, 1.0);
    Qubit qubit_1 = Qubit(1.0, 0.0);
    Qubit qubit_2 = Qubit(0.0, 1.0);
    std::vector<Qubit> qubits = {qubit_0, qubit_1, qubit_2};

    Eigen::VectorXcd expected_state;
    expected_state.resize(8);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);

    std::vector<int> qubit_indices{0, 1};   // 101 -> 111
    circuit.CNOT(qubit_indices);
    REQUIRE(circuit.get_state() == expected_state);
};

TEST_CASE("Evaluating SWAP gate on non-adjacent qubits") {
    Qubit qubit_0 = Qubit(1.0, 0.0);
    Qubit qubit_1 = Qubit(0.0, 1.0);
    Qubit qubit_2 = Qubit(0.0, 1.0);
    std::vector<Qubit> qubits = {qubit_0, qubit_1, qubit_2};

    Eigen::VectorXcd expected_state;
    expected_state.resize(8);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);

    std::vector<int> qubit_indices{0, 2};
    circuit.SWAP(qubit_indices);
    REQUIRE(circuit.get_state() == expected_state);
};

TEST_CASE("Testing CCNOT/Toffoli Gate on three qubits") {
    Qubit qubit_0 = Qubit(0.0, 1.0);
    Qubit qubit_1 = Qubit(0.0, 1.0);
    Qubit qubit_2 = Qubit(1.0, 0.0);
    std::vector<Qubit> qubits = {qubit_0, qubit_1, qubit_2};

    Eigen::VectorXcd expected_state;
    expected_state.resize(8);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);
    circuit.CCNOT(0, 1, 2);     // 110 -> 111
    REQUIRE(circuit.get_state() == expected_state);
}

TEST_CASE("Testing CCNOT/Toffoli Gate on four qubits") {
    Qubit qubit_0 = Qubit(0.0, 1.0);
    Qubit qubit_1 = Qubit(1.0, 0.0); 
    Qubit qubit_2 = Qubit(0.0, 1.0);
    Qubit qubit_3 = Qubit(0.0, 1.0); 
    std::vector<Qubit> qubits = {qubit_0, qubit_1, qubit_2, qubit_3};

    Eigen::VectorXcd expected_state;
    expected_state.resize(16);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);
    circuit.CCNOT(0, 2, 3);     // 1011 -> 1010
    REQUIRE(circuit.get_state() == expected_state);
}

TEST_CASE("Quantum Teleportation") {
    Qubit qubit_m = Qubit(1.0/sqrt(3.0), sqrt(2.0/3.0));
    Qubit qubit_1 = Qubit(1.0, 0.0);
    Qubit qubit_2 = Qubit(1.0, 0.0);
    std::vector<Qubit> qubits = {qubit_m, qubit_1, qubit_2};

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);
    circuit.H(1);
    circuit.CNOT(std::vector<int>{1, 2});
    circuit.CNOT(std::vector<int>{0, 1});
    circuit.H(0);
    std::vector<int> results = circuit.measure(std::vector<int>{0, 1});
    std::cout << "Qubit Zero: " << results.at(0) << std::endl;
    std::cout << "Qubit One: " << results.at(1) << std::endl;
    if (results.at(1) == 1) {circuit.X(2);}
    if (results.at(0) == 1) {circuit.Z(2);}
    std::cout << circuit.get_state()<< std::endl;
}

TEST_CASE("Two-bit full adder with CNOT and Toffoli") {
    Qubit qubit_0 = Qubit(0.0, 1.0);
    Qubit qubit_1 = Qubit(0.0, 1.0); 
    Qubit qubit_2 = Qubit(1.0, 0.0);
    Qubit qubit_3 = Qubit(1.0, 0.0); 
    std::vector<Qubit> qubits = {qubit_0, qubit_1, qubit_2, qubit_3};

    Eigen::VectorXcd expected_state;
    expected_state.resize(16);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);
    circuit.CCNOT(0, 1, 3);
    circuit.CNOT(std::vector<int>{0, 1});
    circuit.CCNOT(1, 2, 3);
    circuit.CNOT(std::vector<int>{1, 2});
    circuit.CNOT(std::vector<int>{0, 1});
    REQUIRE(circuit.get_state() == expected_state);
}

TEST_CASE("Default 4_adder from enfield compiler examples") {
    Qubit qubit_0 = Qubit(1.0, 0.0);
    Qubit qubit_1 = Qubit(1.0, 0.0); 
    Qubit qubit_2 = Qubit(1.0, 0.0);
    Qubit qubit_3 = Qubit(1.0, 0.0); 
    std::vector<Qubit> qubits = {qubit_0, qubit_1, qubit_2, qubit_3};

    Eigen::VectorXcd expected_state;
    expected_state.resize(16);
    expected_state << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);
    circuit.CNOT(std::vector<int>{1, 2});
    circuit.H(3);
    circuit.CNOT(std::vector<int>{2, 3});
    circuit.R(M_PI/8.0, 3, 1);
    circuit.CNOT(std::vector<int>{0, 3});
    circuit.R(M_PI/8.0, 3);
    circuit.CNOT(std::vector<int>{2, 3});
    circuit.R(M_PI/8.0, 2);
    circuit.R(M_PI/8.0, 3, 1);
    circuit.CNOT(std::vector<int>{0, 3});
    circuit.CNOT(std::vector<int>{0, 2});
    circuit.R(M_PI/8.0, 3);
    circuit.R(M_PI/8.0, 0);
    circuit.R(M_PI/8.0, 2, 1);
    circuit.H(3);
    REQUIRE(fabs(circuit.get_state()(1) - expected_state(1)) < 0.00001);
}