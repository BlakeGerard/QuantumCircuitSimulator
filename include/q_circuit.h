#include "qubit.h"
#include <Eigen/Dense>
#include <vector>
#include <random>

class Q_Circuit {
    public:
        Q_Circuit();
        void add_qubits(int size);
        void add_qubits(std::vector<int> bit_string);
        void add_qubits(std::vector<Qubit> qubits);
        void H(int qubit_index);
        void H(std::vector<int> qubit_indices);
        void X(int qubit_index);
        void X(std::vector<int> qubit_indices);
        void Y(int qubit_index);
        void Y(std::vector<int> qubit_indices);
        void Z(int qubit_index);
        void Z(std::vector<int> qubit_indices);
        void CNOT(std::vector<int> qubit_indices);
        void SWAP(std::vector<int> qubit_indices); 
        std::vector<int> measure(std::vector<int> qubit_indices);
        Eigen::VectorXcd get_state();

    private:
        Eigen::VectorXcd state;
        std::mt19937 generator;
        void apply_single_qubit_gate(int qubit_index, Eigen::Matrix2cd gate);
        void apply_single_qubit_gate(std::vector<int> qubit_indices, Eigen::Matrix2cd gate);
        void apply_cnot_gate(std::vector<int> qubit_indices, Eigen::Matrix2cd gate);
        void apply_swap_gate(std::vector<int> qubit_indices, Eigen::Matrix4d gate);
        int measure_single_qubit(int qubit_index);
};