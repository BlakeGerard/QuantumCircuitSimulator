#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <iostream>
#include <complex>
#include <chrono>
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


TEST_CASE("Evaluating CNOT gate in upward direction") {
    Qubit qubit_0 = Qubit(0.0, 1.0);
    Qubit qubit_1 = Qubit(0.0, 1.0);
    Qubit qubit_2 = Qubit(0.0, 1.0);
    std::vector<Qubit> qubits = {qubit_0, qubit_1, qubit_2};

    Eigen::VectorXcd expected_state;
    expected_state.resize(8);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);

    circuit.CNOT(0, 2);
    REQUIRE(circuit.get_state() == expected_state);
};


TEST_CASE("Evaluating CNOT gate in downward direction") {
    Qubit qubit_0 = Qubit(0.0, 1.0);
    Qubit qubit_1 = Qubit(0.0, 1.0);
    Qubit qubit_2 = Qubit(0.0, 1.0);
    std::vector<Qubit> qubits = {qubit_0, qubit_1, qubit_2};

    Eigen::VectorXcd expected_state;
    expected_state.resize(8);
    expected_state << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);

    circuit.CNOT(2, 0);
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

    circuit.SWAP(0, 2);
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

TEST_CASE("Testing CCNOT/Toffoli Gate on four qubits in upward direction") {
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

TEST_CASE("Testing CCNOT/Toffoli Gate on four qubits in downward direction") {
    Qubit qubit_0 = Qubit(0.0, 1.0);
    Qubit qubit_1 = Qubit(0.0, 1.0); 
    Qubit qubit_2 = Qubit(0.0, 1.0);
    Qubit qubit_3 = Qubit(0.0, 1.0); 
    std::vector<Qubit> qubits = {qubit_0, qubit_1, qubit_2, qubit_3};

    Eigen::VectorXcd expected_state; //1111 ->0111
    expected_state.resize(16);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);
    circuit.CCNOT(3, 2, 0);
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
    circuit.CNOT(1, 2);
    circuit.CNOT(0, 1);
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
    circuit.CNOT(0, 1);
    circuit.CCNOT(1, 2, 3);
    circuit.CNOT(1, 2);
    circuit.CNOT(0, 1);
    REQUIRE(circuit.get_state() == expected_state);
}


TEST_CASE("Default 4_adder from enfield compiler examples") {
    auto start = std::chrono::high_resolution_clock::now(); 
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
    circuit.CNOT(1, 2);
    circuit.H(3);
    circuit.CNOT(2, 3);
    circuit.T(3, 1);
    circuit.CNOT(0, 3);
    circuit.T(3, 0);
    circuit.CNOT(2, 3);
    circuit.T(2, 0);
    circuit.T(3, 1);
    circuit.CNOT(0, 3);
    circuit.CNOT(0, 2);
    circuit.T(3, 0);
    circuit.T(0, 0);
    circuit.T(2, 1);
    circuit.H(3);

    auto end = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    REQUIRE(fabs(circuit.get_state()(1) - expected_state(1)) < 0.00001);
    std::cout << "----4_adder default----" << std::endl;
    std::cout << "Time elapsed: " << duration.count() << " microseconds" << std::endl;
}


TEST_CASE("4_adder optimized by enfield compiler") {
    auto start = std::chrono::high_resolution_clock::now(); 
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
    circuit.U3(M_PI/2.0, 0.0, M_PI, 1);
    circuit.CNOT(3, 2);
    circuit.H(2);
    circuit.H(1);
    circuit.CNOT(1, 2);
    circuit.H(1);
    circuit.H(2);
    circuit.U3(0.0, 0.0, -M_PI/4.0, 1);
    circuit.CNOT(0, 1);
    circuit.U3(0.0, 0.0, M_PI/4.0, 1);
    circuit.H(2);
    circuit.H(1);
    circuit.CNOT(1, 2);
    circuit.H(1);
    circuit.H(2);
    circuit.U3(0.0, 0.0, M_PI/4.0, 2);
    circuit.U3(0.0, 0.0, -M_PI/4.0, 1);
    circuit.CNOT(0, 1);
    circuit.U3(0.0, 0.0, M_PI/4.0, 1);
    circuit.CNOT(0, 2);
    circuit.U3(0.0, 0.0, M_PI/4.0, 0);
    circuit.U3(0.0, 0.0, -M_PI/4.0, 2);
    circuit.U3(M_PI/2.0, 0.0, M_PI, 1);
    
    auto end = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    REQUIRE(fabs(circuit.get_state()(1) - expected_state(1)) < 0.00001);
    std::cout << "----4_adder enfield optimized----" << std::endl;
    std::cout << "Time elapsed: " << duration.count() << " microseconds" << std::endl;
}
