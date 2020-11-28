#include <iostream>
#include "include/qcs.h"

int main() {
    Q_Circuit q = Q_Circuit();
    q.add_qubits(3);
    q.H(0);
    std::cout << q.get_state() << std::endl;
}