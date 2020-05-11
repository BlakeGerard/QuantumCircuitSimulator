#include "gate_set.h"

static Eigen::Matrix2cd H = [] {
    Eigen::Matrix2d gate;        
    gate << 1.0/sqrt(2.0), 1.0/sqrt(2.0), 1.0/sqrt(2.0), -1.0/sqrt(2.0);
    return gate;
}();

static Eigen::Matrix2cd X = [] {
    Eigen::Matrix2d gate;
    gate << 0.0, 1.0, 1.0, 0.0;
    return gate;
}();

static Eigen::Matrix2cd Y = [] {
    Eigen::Matrix2cd gate;
    gate << 0.0, std::complex<double>(0.0, -1.0), std::complex<double>(0.0, 1.0), 0.0;
    return gate;
}();

static Eigen::Matrix2cd Z = [] {
    Eigen::Matrix2d gate;
    gate << 1.0, 0.0, 0.0, -1.0;
    return gate;
}();

static Eigen::Matrix2cd S = [] {
    Eigen::Matrix2cd gate;
    gate << 1.0, 0.0, 0.0, exp(std::complex<double>(0.0, 1.0) * (M_PI/2.0));
    return gate;
}();

static Eigen::Matrix2cd T = [] {
    Eigen::Matrix2cd gate;
    gate << 1.0, 0.0, 0.0, exp(std::complex<double>(0.0, 1.0) * (M_PI/4.0));
    return gate;
}();
        
static Eigen::Matrix4cd SWAP = [] {
    Eigen::Matrix4d gate;
    gate << 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,                
            0.0, 1.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 1.0;
    return gate;    
}();

std::vector<Eigen::MatrixXcd, Eigen::aligned_allocator<Eigen::MatrixXcd>> gates = 
{ H, X, Y, Z, S, T, SWAP };