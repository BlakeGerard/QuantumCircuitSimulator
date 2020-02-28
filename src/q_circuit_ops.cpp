#include "q_circuit.h"
#include <Eigen/StdVector>

void Q_Circuit::apply_single_qubit_gate(int qubit_index, Eigen::Matrix2cd gate) {
    Eigen::MatrixXcd operation;

    if (qubit_index == 0) {
        Eigen::MatrixXd post_index_identity;
        post_index_identity.resize(pow(2.0, log2(state.size())-1.0), pow(2.0, log2(state.size())-1.0));
        post_index_identity.setIdentity();
        operation = kroneckerProduct(gate, post_index_identity).eval();

    } else if (qubit_index == log2(state.size()) - 1.0) {
        Eigen::MatrixXd pre_index_identity;
        pre_index_identity.resize(pow(2, qubit_index), pow(2, qubit_index));
        pre_index_identity.setIdentity();
        operation = kroneckerProduct(pre_index_identity, gate).eval();

    } else {
        Eigen::MatrixXd pre_index_identity;
        Eigen::MatrixXd post_index_identity;
        
        pre_index_identity.resize(pow(2.0, qubit_index), pow(2.0, qubit_index));
        pre_index_identity.setIdentity();
        operation = kroneckerProduct(pre_index_identity, gate).eval();

        post_index_identity.resize(pow(2.0, (log2(state.size()) - 1.0) - qubit_index), pow(2.0, (log2(state.size()) - 1.0) - qubit_index));
        post_index_identity.setIdentity();
        operation = kroneckerProduct(operation, post_index_identity).eval();
    }
    state = operation * state;
}; 

void Q_Circuit::apply_single_qubit_gate(std::vector<int> qubit_indices, Eigen::Matrix2cd gate) {
    if (qubit_indices.size() == 1) {
        apply_single_qubit_gate(qubit_indices.at(0), gate);
        return;
    }

    Eigen::MatrixXd identity;
    Eigen::MatrixXcd operation;

    int previous_index = 0;
    int difference = 0;
    for (int i = 0; i < qubit_indices.size(); ++i) {
        difference = qubit_indices.at(i) - previous_index;

        // At first index, need to set operation to either gate or some identity matrix
        if (i == 0) {
        
            // Gate is applied to first qubit
            if (difference == 0) {
                operation = gate;
            } else {
                identity.resize(pow(2, difference), pow(2, difference));
                identity.setIdentity();
                operation = kroneckerProduct(identity, gate).eval(); 
            }
        } else {
            
            // Gate applied to very next qubit
            if (difference == 1) {
                operation = kroneckerProduct(operation, gate).eval();
            } else {
                identity.resize(pow(2, difference - 1), pow(2, difference - 1));
                identity.setIdentity();
                operation = kroneckerProduct(operation, identity).eval();
                operation = kroneckerProduct(operation, gate).eval();
            }

            if (i == qubit_indices.size() - 1 && qubit_indices.at(i) != log2(state.size() - 1)) {
                identity.resize(pow(2, (log2(state.size()) - 1) - qubit_indices.at(i)), pow(2, (log2(state.size()) - 1) - qubit_indices.at(i)));
                identity.setIdentity();
                operation = kroneckerProduct(operation, identity).eval();
            }
        }        
        previous_index = qubit_indices.at(i); 
    }
    state = operation * state;
};

/*
ONLY WORKS WHEN qubit_indices(1) > qubit_indices(0)
*/
void Q_Circuit::apply_controlled_single_qubit_gate(int control, int target, Eigen::Matrix2cd gate) {
    assert(target != control);

    Eigen::Matrix2d zero_proj;
    Eigen::Matrix2d one_proj;
    zero_proj << 1.0, 0.0, 0.0, 0.0;
    one_proj << 0.0, 0.0, 0.0, 1.0;

    Eigen::MatrixXd control_identity;
    Eigen::MatrixXd spanning_identity;
    Eigen::MatrixXcd left_side;
    Eigen::MatrixXcd right_side;
    Eigen::MatrixXcd operation;

    if (control < target) {
        control_identity.resize(pow(2, target - control), pow(2, target - control));
        control_identity.setIdentity();
        left_side = kroneckerProduct(zero_proj, control_identity).eval();

        // Use default CNOT gate or CNOT spanning multiple qubits
        if (target - control != 1) {
            spanning_identity.resize(pow(2, (target - control) - 1), pow(2, (target - control) - 1));
            spanning_identity.setIdentity();
            right_side = kroneckerProduct(one_proj, kroneckerProduct(spanning_identity, gate).eval()).eval();
        } else {
            right_side = kroneckerProduct(one_proj, gate).eval();
        }
        operation = left_side + right_side;
        apply_pre_and_post_identity_matrices(operation, control, target);

    } else {
        control_identity.resize(pow(2, control - target), pow(2, control - target));
        control_identity.setIdentity();
        left_side = kroneckerProduct(control_identity, zero_proj).eval();

        // Use default CNOT gate or CNOT spanning multiple qubits
        if (control - target != 1) {
            spanning_identity.resize(pow(2, (control - target) - 1), pow(2, (control - target) - 1));
            spanning_identity.setIdentity();
            right_side = kroneckerProduct(kroneckerProduct(gate, spanning_identity).eval(), one_proj).eval();
        } else {
            right_side = kroneckerProduct(gate, one_proj).eval();
        }
        operation = left_side + right_side;
        apply_pre_and_post_identity_matrices(operation, control, target);
    }
    state = operation * state;
};

void Q_Circuit::apply_swap_gate(int qubit1, int qubit2, Eigen::Matrix4d gate) {
    assert(qubit2 > qubit1);
    assert(qubit2 != qubit1);

    int difference = qubit2 - qubit1;
    Eigen::MatrixXcd operation;
    operation.resize(pow(2.0, difference + 1), pow(2.0, difference + 1));

    // Qubits are not adjacent, use calculated SWAP gate
    if (difference > 1) {
        Eigen::Matrix2d zero_proj;
        Eigen::Matrix2d one_proj;
        Eigen::Matrix2d two_proj;
        Eigen::Matrix2d three_proj;
        zero_proj << 1.0, 0.0, 0.0, 0.0;
        one_proj << 0.0, 1.0, 0.0, 0.0;
        two_proj << 0.0, 0.0, 1.0, 0.0;
        three_proj << 0.0, 0.0, 0.0, 1.0;

        Eigen::MatrixXd operation_identity;
        operation_identity.resize(pow(2, difference - 1), pow(2, difference - 1));
        operation_identity.setIdentity();

        // Need to resize the shift matrices if there are more than 3 qubits
        Eigen::MatrixXd I0;
        Eigen::MatrixXd I1;
        Eigen::MatrixXd I2;
        Eigen::MatrixXd I3;

        I0 = kroneckerProduct(operation_identity, zero_proj).eval();
        I1 = kroneckerProduct(operation_identity, one_proj).eval();
        I2 = kroneckerProduct(operation_identity, two_proj).eval();
        I3 = kroneckerProduct(operation_identity, three_proj).eval();

        operation.block(0, 0, sqrt(I0.size()), sqrt(I0.size())) = I0;
        operation.block(sqrt(I1.size()), 0, sqrt(I1.size()), sqrt(I1.size())) = I1;
        operation.block(0, sqrt(I2.size()), sqrt(I2.size()), sqrt(I2.size())) = I2;
        operation.block(sqrt(I3.size()), sqrt(I3.size()), sqrt(I3.size()), sqrt(I3.size())) = I3;

    } else {
        operation = gate;
    };

    apply_pre_and_post_identity_matrices(operation, qubit1, qubit2);
    state = operation * state;
};

void Q_Circuit::apply_controlled_two_qubit_gate(int control1, int control2, int target, Eigen::Matrix2cd gate) {
    assert(control1 != control2 != target);
    assert((control1 < control2 < target) || (target < control2 < control1));

    Eigen::Matrix2d one_proj;
    one_proj << 0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix2d identity_2 = Eigen::Matrix2d::Identity();
    Eigen::Matrix2cd gate_modified = gate - identity_2;

    Eigen::MatrixXd spanning_identity;
    Eigen::MatrixXcd gate_configuration;
    Eigen::MatrixXcd operation;

    if (control1 < target) {
        int c2_c1_diff = control2 - control1;
        int t_c2_diff = target - control2;
        int t_c1_diff = target - control1;

        spanning_identity.resize(pow(2, t_c1_diff + 1), pow(2, t_c1_diff + 1));
        spanning_identity.setIdentity();

        operation = spanning_identity;

        if (c2_c1_diff > 1) {
            Eigen::MatrixXd temp_identity;
            temp_identity.resize(pow(2, c2_c1_diff - 1), pow(2, c2_c1_diff - 1));
            temp_identity.setIdentity();
            gate_configuration = kroneckerProduct(one_proj, temp_identity).eval();
        } else {
            gate_configuration = one_proj;
        }

        if (t_c2_diff > 1) {
            Eigen::MatrixXd temp_identity;
            temp_identity.resize(pow(2, t_c2_diff - 1), pow(2, t_c2_diff - 1));
            temp_identity.setIdentity();
            gate_configuration = kroneckerProduct(gate_configuration, one_proj).eval();
            gate_configuration = kroneckerProduct(gate_configuration, temp_identity).eval();
        } else {
            gate_configuration = kroneckerProduct(gate_configuration, one_proj).eval();
        }

        gate_configuration = kroneckerProduct(gate_configuration, gate_modified).eval();
        operation += gate_configuration;

        apply_pre_and_post_identity_matrices(operation, control1, target);
    } else {
        int c1_c2_diff = control1 - control2;
        int c2_t_diff = control2 - target;
        int c1_t_diff = control1 - target;

        spanning_identity.resize(pow(2, c1_t_diff + 1), pow(2, c1_t_diff + 1));
        spanning_identity.setIdentity();

        operation = spanning_identity;
        gate_configuration = gate_modified;

        if (c2_t_diff > 1) {
            Eigen::MatrixXd temp_identity;
            temp_identity.resize(pow(2, c2_t_diff - 1), pow(2, c2_t_diff - 1));
            temp_identity.setIdentity();
            gate_configuration = kroneckerProduct(gate_configuration, temp_identity).eval();
            gate_configuration = kroneckerProduct(gate_configuration, one_proj).eval();
        } else {
            gate_configuration = kroneckerProduct(gate_configuration, one_proj);
        }

        if (c1_c2_diff > 1) {
            Eigen::MatrixXd temp_identity;
            temp_identity.resize(pow(2, c2_t_diff - 1), pow(2, c2_t_diff - 1));
            temp_identity.setIdentity();
            gate_configuration = kroneckerProduct(gate_configuration, temp_identity).eval();
            gate_configuration = kroneckerProduct(gate_configuration, one_proj).eval();
        } else {
            gate_configuration = kroneckerProduct(gate_configuration, one_proj).eval();
        }
        operation += gate_configuration;
        apply_pre_and_post_identity_matrices(operation, control1, target);
    }
    state = operation * state;
}

std::vector<int> Q_Circuit::measure(std::vector<int> qubit_indices) {
    if (qubit_indices.size() == 1) {
        std::vector<int> results;
        int result = measure_single_qubit(qubit_indices.at(0));
        results.push_back(result);
        return results;
    }

    Eigen::Matrix2d zero_proj;
    Eigen::Matrix2d one_proj;
    zero_proj << 1.0, 0.0, 0.0, 0.0;
    one_proj << 0.0, 0.0, 0.0, 1.0;
    std::vector<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d>> projectors;
    projectors.push_back(zero_proj);
    projectors.push_back(one_proj);

    std::vector<double> amplitudes;
    std::vector<Eigen::VectorXcd, Eigen::aligned_allocator<Eigen::VectorXcd>> possible_new_states;
    std::vector<int> results;
    int result;
    int msb = pow(2, qubit_indices.size()) / 2;
    
    for (int i = 0; i < pow(2, qubit_indices.size()); ++i) {
        Eigen::MatrixXd identity;
        Eigen::MatrixXcd projection_matrix;

        int previous_index = 0;
        int difference = 0;
        for (int j = 0; j < qubit_indices.size(); ++j) {
            difference = qubit_indices.at(j) - previous_index;

            int proj_index = (i & (msb >> j) ) >> ((qubit_indices.size() - 1) - j);
            Eigen::Matrix2d current_projector = projectors.at(proj_index);


            // At first index, need to set operation to either gate or some identity matrix
            if (j == 0) {
        
                // Gate is applied to first qubit
                if (difference == 0) {
                    projection_matrix = current_projector;
                } else {
                    identity.resize(pow(2, difference), pow(2, difference));
                    identity.setIdentity();
                    projection_matrix = kroneckerProduct(identity, current_projector).eval(); 
                }

            } else {
            
                // Gate applied to very next qubit
                if (difference == 1) {
                    projection_matrix = kroneckerProduct(projection_matrix, current_projector).eval();
                } else {
                    identity.resize(pow(2, difference - 1), pow(2, difference - 1));
                    identity.setIdentity();
                    projection_matrix = kroneckerProduct(projection_matrix, identity).eval();
                    projection_matrix = kroneckerProduct(projection_matrix, current_projector).eval();
                }

                if (j == qubit_indices.size() - 1 && qubit_indices.at(j) != log2(state.size() - 1)) {
                identity.resize(pow(2, (log2(state.size()) - 1) - qubit_indices.at(j)), pow(2, (log2(state.size()) - 1) - qubit_indices.at(j)));
                identity.setIdentity();
                projection_matrix = kroneckerProduct(projection_matrix, identity).eval();
            }

            }        
            previous_index = qubit_indices.at(j); 
        }
        
        Eigen::VectorXcd state_adjoint = state;
        possible_new_states.push_back(projection_matrix * state);
        std::complex<double> amp = (state_adjoint.adjoint() * projection_matrix) * state;
        amplitudes.push_back(amp.real());
    }

    std::discrete_distribution<int> distribution(amplitudes.begin(), amplitudes.end());
    result = distribution(generator);
    state = possible_new_states.at(result);
    state.normalize();

    for (int k = 0; k < qubit_indices.size(); ++k) {    
        results.push_back((result & (msb >> k) ) >> ((qubit_indices.size() - 1) - k));
    }

    return results;
};

int Q_Circuit::measure_single_qubit(int qubit_index) {
    Eigen::Matrix2d zero_proj;
    Eigen::Matrix2d one_proj;
    zero_proj << 1.0, 0.0, 0.0, 0.0;
    one_proj << 0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix2d gate = zero_proj;
    Eigen::MatrixXcd projection_matrix;

    std::vector<double> amplitudes;
    std::vector<Eigen::VectorXcd, Eigen::aligned_allocator<Eigen::VectorXcd>> possible_new_states;
    std::vector<int> results;
    int result;

    for (int i = 0; i < 2; ++i) {
        if (i == 1) {gate = one_proj;}

        if (qubit_index == 0) {
            Eigen::MatrixXd post_index_identity;
            post_index_identity.resize(pow(2.0, log2(state.size())-1.0), pow(2.0, log2(state.size())-1.0));
            post_index_identity.setIdentity();
            projection_matrix = kroneckerProduct(gate, post_index_identity).eval();

        } else if (qubit_index == log2(state.size()) - 1.0) {
            Eigen::MatrixXd pre_index_identity;
            pre_index_identity.resize(pow(2, qubit_index), pow(2, qubit_index));
            pre_index_identity.setIdentity();
            projection_matrix = kroneckerProduct(pre_index_identity, gate).eval();

        } else {
            Eigen::MatrixXd pre_index_identity;
            Eigen::MatrixXd post_index_identity;
        
            pre_index_identity.resize(pow(2.0, qubit_index), pow(2.0, qubit_index));
            pre_index_identity.setIdentity();
            projection_matrix = kroneckerProduct(pre_index_identity, gate).eval();

            post_index_identity.resize(pow(2.0, (log2(state.size()) - 1.0) - qubit_index), pow(2.0, (log2(state.size()) - 1.0) - qubit_index));
            post_index_identity.setIdentity();
            projection_matrix = kroneckerProduct(projection_matrix, post_index_identity).eval();
        }

        Eigen::VectorXcd state_adjoint = state;
        possible_new_states.push_back(projection_matrix * state);
        std::complex<double> amp = (state_adjoint.adjoint() * projection_matrix) * state;
        amplitudes.push_back(amp.real());
    }

    std::discrete_distribution<int> distribution(amplitudes.begin(), amplitudes.end());
    result = distribution(generator);
    state = possible_new_states.at(result);
    state.normalize();

    return result;
}

void Q_Circuit::apply_pre_and_post_identity_matrices(Eigen::MatrixXcd &operation, int control, int target) {

    // CNOT(2, 0): do things in reverse
    if (control > target) {
        // Pre and post identity matrices both exists
        if (control != log2(state.size()) - 1 && (target != 0)) {

            Eigen::MatrixXd pre_target_index_identity;
            pre_target_index_identity.resize(pow(2, target), pow(2, target));
            pre_target_index_identity.setIdentity();
            operation = kroneckerProduct(pre_target_index_identity, operation).eval();

            Eigen::MatrixXd post_control_index_identity;
            post_control_index_identity.resize(pow(2, (log2(state.size())- 1) - control), pow(2, (log2(state.size())- 1) - control));
            post_control_index_identity.setIdentity();
            operation = kroneckerProduct(operation, post_control_index_identity).eval();
            return;

        // No pre identity, but post exists
        } else if (control != log2(state.size()) - 1 && (target == 0)) {
        
            Eigen::MatrixXd post_control_index_identity;
            post_control_index_identity.resize(pow(2, (log2(state.size())- 1) - control), pow(2, (log2(state.size())- 1) - control));
            post_control_index_identity.setIdentity();
            operation = kroneckerProduct(operation, post_control_index_identity).eval();
            return;

        // No post identity, but pre exists
        } else if (control == log2(state.size()) - 1 && (target != 0)) {
            Eigen::MatrixXd pre_target_index_identity;
            pre_target_index_identity.resize(pow(2, target), pow(2, target));
            pre_target_index_identity.setIdentity();
            operation = kroneckerProduct(pre_target_index_identity, operation).eval();

        // No pre or post identity matrices
        };

    } else {
        // Pre and post identity matrices both exists
        if (control != 0 && target != log2(state.size()) - 1) {

            Eigen::MatrixXd pre_first_index_identity;
            pre_first_index_identity.resize(pow(2, control), pow(2, control));
            pre_first_index_identity.setIdentity();
            operation = kroneckerProduct(pre_first_index_identity, operation).eval();

            Eigen::MatrixXd post_second_index_identity;
            post_second_index_identity.resize(pow(2, (log2(state.size()) - 1) - target), pow(2, (log2(state.size()) - 1) - target));
            post_second_index_identity.setIdentity();
            operation = kroneckerProduct(operation, post_second_index_identity).eval();
            return;
        
        // No pre identity, but post exists
        } else if (control == 0 && target != log2(state.size()) - 1) {
        
            Eigen::MatrixXd post_second_index_identity;
            post_second_index_identity.resize(pow(2, (log2(state.size()) - 1) - target), pow(2, (log2(state.size()) - 1) - target));
            post_second_index_identity.setIdentity();
            operation = kroneckerProduct(operation, post_second_index_identity).eval();
            return;

        // No post identity, but pre exists
        } else if (control != 0 && target == log2(state.size()) - 1) {
        
            Eigen::MatrixXd pre_first_index_identity;
            pre_first_index_identity.resize(pow(2, control), pow(2, control));
            pre_first_index_identity.setIdentity();
            operation = kroneckerProduct(pre_first_index_identity, operation).eval();
            return;

        // No pre or post identity matrices
        };
    }
}