/*
        Title: gate_matrices.cpp
        Author: Blake Gerard
        Date: 05/16/2020
        Description: All necessary gate matrices are statically built and stored 
        in the "gate_matrices" vector. Any q_circuit instance will reference this
        vector when a gate is applied (through q_circuit.H, q_circuit.CNOT, etc.).
        We only need to construct and store the supported single-qubit gates, as
        higher-dimensional gates (aside from SWAP) are simply controlled versions
        of single-qubit gates.
*/

#include "gate_matrices.h"

/*
    Static construction of the Hadamard gate as an Eigen matrix.
*/
static Eigen::Matrix2cd H = [] {
    Eigen::Matrix2d gate;        
    gate << 1.0/sqrt(2.0), 1.0/sqrt(2.0), 1.0/sqrt(2.0), -1.0/sqrt(2.0);
    return gate;
}();

/*
    Static construction of the Pauli-X gate as an Eigen matrix.
*/
static Eigen::Matrix2cd X = [] {
    Eigen::Matrix2d gate;
    gate << 0.0, 1.0, 1.0, 0.0;
    return gate;
}();

/*
    Static construction of the Pauli-Y gate as an Eigen matrix.
*/
static Eigen::Matrix2cd Y = [] {
    Eigen::Matrix2cd gate;
    gate << 0.0, std::complex<double>(0.0, -1.0), std::complex<double>(0.0, 1.0), 0.0;
    return gate;
}();

/*
    Static construction of the Pauli-Z gate as an Eigen matrix.
*/
static Eigen::Matrix2cd Z = [] {
    Eigen::Matrix2d gate;
    gate << 1.0, 0.0, 0.0, -1.0;
    return gate;
}();

/*
    Static construction of the S gate as an Eigen matrix.
*/
static Eigen::Matrix2cd S = [] {
    Eigen::Matrix2cd gate;
    gate << 1.0, 0.0, 0.0, exp(std::complex<double>(0.0, 1.0) * (M_PI/2.0));
    return gate;
}();

/*
    Static construction of the T gate as an Eigen matrix.
*/
static Eigen::Matrix2cd T = [] {
    Eigen::Matrix2cd gate;
    gate << 1.0, 0.0, 0.0, exp(std::complex<double>(0.0, 1.0) * (M_PI/4.0));
    return gate;
}();
        
/*
    Static construction of the SWAP gate as an Eigen matrix.
*/
static Eigen::Matrix4cd SWAP = [] {
    Eigen::Matrix4d gate;
    gate << 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,                
            0.0, 1.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 1.0;
    return gate;    
}();

/*
    Vector from which all q_circuit instances pull gates.
*/
static std::vector<Eigen::MatrixXcd, Eigen::aligned_allocator<Eigen::MatrixXcd>> gate_matrices = 
{ H, X, Y, Z, S, T, SWAP };