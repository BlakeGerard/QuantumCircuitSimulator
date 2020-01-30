#include <complex>
#include <iostream>
#include </usr/include/eigen3/unsupported/Eigen/KroneckerProduct>

class State {
    public:
        State() {

        };

        void prepare_state(std::vector<Eigen::Vector2cd> qubits) {
            Eigen::VectorXcd temp_state_vector;
            temp_state_vector.resize(2);
            temp_state_vector << qubits[0];

            for (unsigned int i = 1; i < qubits.size(); ++i) {
                temp_state_vector = Eigen::kroneckerProduct(temp_state_vector, qubits[i]).eval();
            }
            state = temp_state_vector;
        };

        Eigen::VectorXcd get_state() {
            return state;
        };

    private:
        Eigen::VectorXcd state;
};