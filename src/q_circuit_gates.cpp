#include <iostream>
#include "q_circuit.h"
#include "gate_set.h"
#include "enum_gates.h"

void Q_Circuit::H(int qubit_index) {
    if (!state.isZero()) {
        apply_single_qubit_gate(qubit_index, gates[Gates::H]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::H(std::vector<int> qubit_indices) {
    if (!state.isZero()) {   
        apply_single_qubit_gate(qubit_indices, gates[Gates::H]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::X(int qubit_index) {
    if (!state.isZero()) {
        apply_single_qubit_gate(qubit_index, gates[Gates::X]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::X(std::vector<int> qubit_indices) {
    if (!state.isZero()) {
        apply_single_qubit_gate(qubit_indices, gates[Gates::X]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::Y(int qubit_index) {
    if (!state.isZero()) {
        apply_single_qubit_gate(qubit_index, gates[Gates::Y]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::Y(std::vector<int> qubit_indices) {
    if (!state.isZero()) {
        apply_single_qubit_gate(qubit_indices, gates[Gates::Y]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::Z(int qubit_index) {
    if (!state.isZero()) {
        apply_single_qubit_gate(qubit_index, gates[Gates::Z]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::Z(std::vector<int> qubit_indices) {
    if (!state.isZero()) {
        apply_single_qubit_gate(qubit_indices, gates[Gates::Z]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }   
}

void Q_Circuit::R(double phi, int qubit_index, int dagger /*=0*/) {
    if (!state.isZero()) {
        Eigen::Matrix2cd gate;
        gate << 1.0, 0.0, 0.0, exp(std::complex<double>(0.0, 1.0) * phi);
        if (dagger) {gate.adjointInPlace();}
        apply_single_qubit_gate(qubit_index, gate);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }   
}

void Q_Circuit::R(double phi, std::vector<int> qubit_indices, int dagger /*=0*/) {
    if (!state.isZero()) {
        Eigen::Matrix2cd gate;
        gate << 1.0, 0.0, 0.0, exp(std::complex<double>(0.0, 1.0) * phi);
        if (dagger) {gate.adjointInPlace();}
        apply_single_qubit_gate(qubit_indices, gate);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::S(int qubit_index, int dagger) {
    if (!state.isZero()) {
        if (dagger) { apply_single_qubit_gate(qubit_index, gates[Gates::S].adjoint()); }
        else        { apply_single_qubit_gate(qubit_index, gates[Gates::S]); }
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::S(std::vector<int> qubit_indices, int dagger) {
    if (!state.isZero()) {
        if (dagger) { apply_single_qubit_gate(qubit_indices, gates[Gates::S].adjoint()); }
        else        { apply_single_qubit_gate(qubit_indices, gates[Gates::S]); }
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::T(int qubit_index, int dagger) {
    if (!state.isZero()) {
        if (dagger) { apply_single_qubit_gate(qubit_index, gates[Gates::T].adjoint()); }
        else        { apply_single_qubit_gate(qubit_index, gates[Gates::T]); }
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::T(std::vector<int> qubit_indices, int dagger) {
    if (!state.isZero()) {
        if (dagger) { apply_single_qubit_gate(qubit_indices, gates[Gates::T].adjoint()); }
        else        { apply_single_qubit_gate(qubit_indices, gates[Gates::T]); }
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

/*
https://arxiv.org/pdf/1707.03429.pdf, p. 5
*/
void Q_Circuit::U3(double theta, double phi, double lambda, int qubit_index, int dagger) {
    if (!state.isZero()) {
        Eigen::Matrix2cd gate;
        std::complex<double> me_i{0.0, 1.0};
        gate << exp((-me_i * (phi+lambda)) / 2.0) * cos(theta/2.0), -exp((-me_i * (phi-lambda)) / 2.0) * sin(theta/2.0),
                exp((me_i * (phi-lambda)) / 2.0) * sin(theta/2.0), exp((me_i * (phi+lambda)) / 2.0) * cos(theta/2.0);
        if (dagger) {gate.adjointInPlace();}
        apply_single_qubit_gate(qubit_index, gate);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

/*
https://arxiv.org/pdf/1707.03429.pdf, p. 5
*/
void Q_Circuit::U3(double theta, double phi, double lambda, std::vector<int> qubit_indices, int dagger) {
    if (!state.isZero()) {
        Eigen::Matrix2cd gate;
        std::complex<double> me_i{0.0, 1.0};
        gate << exp((-me_i * (phi+lambda)) / 2.0) * cos(theta/2.0), -exp((-me_i * (phi-lambda)) / 2.0) * sin(theta/2.0),
                exp((me_i * (phi-lambda)) / 2.0) * sin(theta/2.0), exp((me_i * (phi+lambda)) / 2.0) * cos(theta/2.0);
        if (dagger) {gate.adjointInPlace();}
        apply_single_qubit_gate(qubit_indices, gate);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::CNOT(int control, int target) {
    if (!state.isZero()) {
        apply_controlled_single_qubit_gate(control, target, gates[Gates::X]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::CY(int control, int target) {
    if (!state.isZero()) {
        apply_controlled_single_qubit_gate(control, target, gates[Gates::Y]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::CZ(int control, int target) {
    if (!state.isZero()) {
        apply_controlled_single_qubit_gate(control, target, gates[Gates::Z]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::CR(double phi, int control, int target) {
    if (!state.isZero()) {
        Eigen::Matrix2cd gate;
        gate << 1.0, 0.0, 0.0, exp(std::complex<double>(0.0, 1.0) * phi);
        apply_controlled_single_qubit_gate(control, target, gate);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::SWAP(int qubit1, int qubit2) {
    if (!state.isZero()) {
        apply_swap_gate(qubit1, qubit2, gates[Gates::SWAP]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}

void Q_Circuit::CCNOT(int control1, int control2, int target) {
    if (!state.isZero()) {
        apply_controlled_two_qubit_gate(control1, control2, target, gates[Gates::X]);
    } else {
        throw "Must initialize qubits before calling gate operations.";
    }
}