#ifndef Q_GATE_H
#define Q_GATE_H
#include "qubit.h"
#include </usr/include/eigen3/unsupported/Eigen/KroneckerProduct>

void hadamard(Qubit *qubit);
Eigen::Vector4cd cnot(Qubit *control, Qubit *target);

#endif