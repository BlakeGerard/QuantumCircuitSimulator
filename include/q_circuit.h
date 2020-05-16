/*
        Title: q_circuit.h
        Author: Blake Gerard
        Date: 05/16/2020
        Description: Class to build quantum circuits. 

        The "state" of the circuit is represented as the Kronecker product
        of all qubits contained in each register. Gate applications to particular 
        qubits are tensored with appropriate matrices such that they span the length 
        of the circuit state. A gate can be applied to one or multiple qubits at a time.

        There are three ways to add qubits to the circuit:
            1. Given a circuit "size", generates a state consisting of "size" qubits of the form [1 0]^T
            2. Initialize circuit state as kronecker product of qubits encoding the given bitstring.
            3. Initialize the circuit state as kronecker product of the provided qubit register.
*/

#ifndef Q_CIRCUIT_H
#define Q_CIRCUIT_H

#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <vector>
#include <random>
#include "qubit.h"

#define M_PI 3.14159265358979323846

class Q_Circuit {

    public:
        Q_Circuit();

        // Add qubits to the circuit
        void add_qubits(int size);
        void add_qubits(std::vector<int> bit_string);
        void add_qubits(std::vector<Qubit> qubits);

        // Single-qubit gates
        void H(int qubit_index);
        void H(std::vector<int> qubit_indices);
        void X(int qubit_index);
        void X(std::vector<int> qubit_indices);
        void Y(int qubit_index);
        void Y(std::vector<int> qubit_indices);
        void Z(int qubit_index);
        void Z(std::vector<int> qubit_indices);
        void R(double phi, int qubit_index, int dagger = 0);
        void R(double phi, std::vector<int> qubit_indices, int dagger = 0);
        void S(int qubit_index, int dagger = 0);
        void S(std::vector<int> qubit_indices, int dagger = 0);
        void T(int qubit_index, int dagger = 0);
        void T(std::vector<int> qubit_indices, int dagger = 0);
        void U3(double theta, double phi, double lambda, int qubit_index, int dagger = 0);
        void U3(double theta, double phi, double lambda, std::vector<int> qubit_indices, int dagger = 0);

        // Controlled single-qubit gates
        void CNOT(int control, int target);
        void CY(int control, int target);
        void CZ(int control, int target);
        void CR(double phi, int control, int target);

        // Controlled Controlled single-qubit gates
        void CCNOT(int control1, int control2, int target);

        // SWAP gate
        void SWAP(int qubit1, int qubit2);

        // Measurement and state retrieval
        std::vector<int> measure(std::vector<int> qubit_indices);
        Eigen::VectorXcd get_state();

    private:
        Eigen::VectorXcd state;
        std::mt19937 generator;
        
        void apply_single_qubit_gate(int qubit_index, Eigen::Matrix2cd gate);
        void apply_single_qubit_gate(std::vector<int> qubit_indices, Eigen::Matrix2cd gate);
        void apply_controlled_single_qubit_gate(int control, int target, Eigen::Matrix2cd gate);
        void apply_controlled_two_qubit_gate(int c1, int c2, int target, Eigen::Matrix2cd gate);
        void apply_swap_gate(int qubit1, int qubit2, Eigen::Matrix4d gate);
        int measure_single_qubit(int qubit_index);
        void apply_pre_and_post_identity_matrices(Eigen::MatrixXcd &operation, int control, int target);
};
#endif