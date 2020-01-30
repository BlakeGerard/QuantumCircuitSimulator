#include <iostream>
#include <complex>
#include "qubit.h"
#include "q_gate.h"

int main() {
    Qubit* ket_zero = new Qubit(1.0, 0.0);
    Qubit* ket_one = new Qubit(0.0, 1.0);
    Hadamard h = Hadamard();
    Pauli_X px = Pauli_X();
    
    std::cout << ket_zero->get_state() << std::endl;
    px.apply(ket_zero);
    std::cout << "--Pauli_X transformation--" << std::endl;
    std::cout << ket_zero->get_state() << std::endl;
};