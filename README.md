# QuantumCircuitSimulator

A personal exercise in modeling quantum computing circuits. Implemented with Eigen 3.3.7

# Dependencies
This project relies on [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page) for matrix computations.

This project utilizes [Catch2](https://github.com/catchorg/Catch2.git) for testing.

# Available Gates
1. Hadamard (H)
2. Pauli-X (X)
3. Pauli-Y (Y)
4. Pauli-Z (Z)
5. Phase (R)
6. S Gate (S)
7. T Gate (T)
8. U3 Gate (U3), https://arxiv.org/pdf/1707.03429.pdf, p. 5
9. CNOT (CNOT)
10. Controlled-Y (CY)
11. Controlled-Z (CZ)
12. Controlled-Phase (CR)
13. SWAP (SWAP)
14. Toffoli (CCNOT)

# Examples
## Quantum Teleportation with Predefined Qubits
    Q_Circuit circuit = Q_Circuit();

    // For this circuit, we pre-define qubit_m
    // Qubits default to state |0>
    Qubit qubit_m = Qubit(1.0/sqrt(3.0), sqrt(2.0/3.0));
    Qubit qubit_1 = Qubit();
    Qubit qubit_2 = Qubit();    
    std::vector<Qubit> qubits = {qubit_m, qubit_1, qubit_2};

    // Add qubits and apply gates to specified indices
    circuit.add_qubits(qubits);
    circuit.H(1);
    circuit.CNOT(1, 2);
    circuit.CNOT(0, 1);
    circuit.H(0);

    // Measure the results
    std::vector<int> results = circuit.measure({0, 1});
    if (results.at(1) == 1) {circuit.X(2);}
    if (results.at(0) == 1) {circuit.Z(2);}

    // Print the final state of the circuit
    std::cout << circuit.get_state() << std::endl;
