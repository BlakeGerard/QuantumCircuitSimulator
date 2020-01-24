#include </usr/include/eigen3/Eigen/Dense>
#include </usr/include/eigen3/unsupported/Eigen/KroneckerProduct>
#include "qubit.h"
#include <vector>

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
        Eigen::VectorXcd get_state();

    private:
        Eigen::VectorXcd state;
        void apply_single_qubit_gate(int qubit_index, Eigen::Matrix2cd gate);
        void apply_single_qubit_gate(std::vector<int> qubit_indices, Eigen::Matrix2cd gate);
        void apply_two_qubit_gate(std::vector<int> qubit_indices, Eigen::Matrix2cd gate);
};