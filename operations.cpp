#include "operations.h"
#include "q_gate.h"

Eigen::Vector4cd bell_state_entanglement(Qubit *control, Qubit *target) {
     hadamard(control);
     return cnot(control, target);
};

Eigen::VectorXcd quantum_teleportation(Qubit *message, Qubit *local, Qubit *target) {
     Eigen::Vector4cd target_local_entanglement = bell_state_entanglement(local, target);
};