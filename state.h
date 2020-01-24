#include </usr/include/eigen3/Eigen/Dense>
#include <vector>

class State {
    public:
        State();
        void prepare_state(std::vector<Eigen::Vector2cd> qubits);
        Eigen::VectorXcd get_state();

    private:
        Eigen::VectorXcd state;
};