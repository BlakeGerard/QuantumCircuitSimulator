#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <iostream>
#include <complex>
#include "../include/qcs.h"

TEST_CASE("Introduction") {
    std::cout << "Running All QuantumCircuitSimulator Tests." << std::endl;
    std::cout << "We will display the result of one test for brevity." << std::endl;
}

TEST_CASE("Quantum Teleportation") {
    Qubit qubit_m = Qubit(1.0/sqrt(3.0), sqrt(2.0/3.0));
    Qubit qubit_1 = Qubit();
    Qubit qubit_2 = Qubit();
    std::vector<Qubit> qubits = {qubit_m, qubit_1, qubit_2};

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(qubits);
    circuit.H(1);
    circuit.CNOT(1, 2);
    circuit.CNOT(0, 1);
    circuit.H(0);
    std::vector<int> results = circuit.measure({0, 1});
    if (results.at(1) == 1) {circuit.X(2);}
    if (results.at(0) == 1) {circuit.Z(2);}

    std::cout << "Printing quantum teleportation test with |m> = [1.0/sqrt(3.0), sqrt(2.0/3.0)]" << std::endl;
    std::cout << "|" << PSI << "> = " << std::endl << circuit.get_state() << std::endl;
    std::cout << "-----------------" << std::endl;
}

TEST_CASE("Prevent gate operations before qubit initialization") {
    Q_Circuit my_circuit = Q_Circuit();
    REQUIRE_THROWS(my_circuit.H(0));
}

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
    Eigen::VectorXcd expected_state;
    expected_state.resize(8);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits({1, 1, 1});
    circuit.CNOT(0, 2);

    REQUIRE(circuit.get_state() == expected_state);
};


TEST_CASE("Evaluating CNOT gate in downward direction") {
    Eigen::VectorXcd expected_state;
    expected_state.resize(8);
    expected_state << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits({1, 1, 1});
    circuit.CNOT(2, 0);

    REQUIRE(circuit.get_state() == expected_state);
};


TEST_CASE("Evaluating SWAP gate on non-adjacent qubits") {
    Eigen::VectorXcd expected_state;
    expected_state.resize(8);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits({0, 1, 1});
    circuit.SWAP(0, 2);

    REQUIRE(circuit.get_state() == expected_state);
};

TEST_CASE("Testing CCNOT/Toffoli Gate on three qubits") {
    Eigen::VectorXcd expected_state;
    expected_state.resize(8);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits({1, 1, 0});
    circuit.CCNOT(0, 1, 2);     // 110 -> 111

    REQUIRE(circuit.get_state() == expected_state);
}

TEST_CASE("Testing CCNOT/Toffoli Gate on four qubits in upward direction") {
    Eigen::VectorXcd expected_state;
    expected_state.resize(16);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits({1, 0, 1, 1});
    circuit.CCNOT(0, 2, 3);     // 1011 -> 1010

    REQUIRE(circuit.get_state() == expected_state);
}

TEST_CASE("Testing CCNOT/Toffoli Gate on four qubits in downward direction") {
    Eigen::VectorXcd expected_state;
    expected_state.resize(16);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits({1, 1, 1, 1});    //1111 ->0111
    circuit.CCNOT(3, 2, 0);

    REQUIRE(circuit.get_state() == expected_state);
}

TEST_CASE("Two-bit full adder with CNOT and Toffoli") {
    Eigen::VectorXcd expected_state;
    expected_state.resize(16);
    expected_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits({1, 1, 0, 0});
    circuit.CCNOT(0, 1, 3);
    circuit.CNOT(0, 1);
    circuit.CCNOT(1, 2, 3);
    circuit.CNOT(1, 2);
    circuit.CNOT(0, 1);

    REQUIRE(circuit.get_state() == expected_state);
}


TEST_CASE("Default 4_adder from enfield compiler examples") {
    Eigen::VectorXcd expected_state;
    expected_state.resize(16);
    expected_state << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits({0, 0, 0, 0});
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

    REQUIRE(fabs(circuit.get_state()(1) - expected_state(1)) < 0.00001);
}

TEST_CASE("4_adder optimized by enfield compiler") {
    Eigen::VectorXcd expected_state;
    expected_state.resize(16);
    expected_state << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits({0, 0, 0, 0});
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

    REQUIRE(fabs(circuit.get_state()(1) - expected_state(1)) < 0.00001);
}

/*
TEST_CASE("Quantum Fourier Transform for n qubits") {
    int n = 3;
    Q_Circuit circuit = Q_Circuit();
    circuit.add_qubits(n);

    circuit.H(0);
    for (int i = 1; i < n; ++i) {
        for (int j = 1; j <= i; ++j) {
            circuit.R(M_PI / ((double)j * 2.0), (double)i);
        }
        circuit.H(i);
    }
    std::cout << circuit.get_state() << std::endl;
}
*/