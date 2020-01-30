#ifndef Q_REGISTER_H
#define Q_REGISTER_H
#include <vector>
#include </usr/include/eigen3/Eigen/Dense>

class Q_Register {
    public:
        Q_Register(int n_qubits);
        void initialize();
        void set_bit_string(std::vector<int> bit_string);
        void set_state(Eigen::VectorXcd new_state);
        Eigen::VectorXcd get_state();

    private:
        int n_qubits;
        Eigen::VectorXcd state;
};

#endif