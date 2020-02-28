# QuantumCircuitSimulator

A personal exercise in modeling quantum computing circuits. Implemented with Eigen 3.3.7

# Available Gates
1. Hadamard (H)
2. Pauli-X (X)
3. Pauli-Y (Y)
4. Pauli-Z (Z)
5. Phase (R)
6. S Gate(S)
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
    Qubit qubit_m = Qubit(1.0/sqrt(3.0), sqrt(2.0/3.0));
    Qubit qubit_1 = Qubit(1.0, 0.0);
    Qubit qubit_2 = Qubit(1.0, 0.0);
    std::vector<Qubit> qubits = {qubit_m, qubit_1, qubit_2};

    circuit.add_qubits(qubits);
    circuit.H(1);
    circuit.CNOT(1, 2);
    circuit.CNOT(0, 1);
    circuit.H(0);
    std::vector<int> results = circuit.measure(std::vector<int>{0, 1});
    if (results.at(1) == 1) {circuit.X(2);}
    if (results.at(0) == 1) {circuit.Z(2);}
    std::cout << circuit.get_state()<< std::endl;